//! Parsing GROMACS simulation output.
//!
//! After an `mdrun`, we use `gmx trjconv` to export the trajectory as a
//! multi-model `.gro` file, then parse each frame here into [`GromacsFrame`].
//! Energy data is extracted from the binary `.edr` file by calling
//! `gmx energy` to convert it to an XVG text file, which is then parsed.

use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{self, BufReader, BufWriter, ErrorKind, Read, Seek, SeekFrom, Write as _},
    path::Path,
    process::{Command, Stdio},
};

use lin_alg::f64::Vec3;

// todo: Consider structs that are "file-like" etc, and impl std traits to read/write etc? (Check examples
// todo in other libs) instead of free-stnading fns.

/// Output thermodynamic / energy properties at a single recorded time step.
///
/// Fields are `Option` because not every simulation records every quantity
/// (e.g. NVT runs have no `pressure`; vacuum runs have no `volume`).
#[derive(Clone, Debug, Default)]
pub struct OutputEnergy {
    /// Simulation time in **ps** — used to match energy records to frames.
    pub time: f64,
    /// Temperature in **K**.
    pub temperature: Option<f32>,
    /// Pressure in **bar**.
    pub pressure: Option<f32>,
    /// Potential energy in **kJ/mol**.
    pub potential_energy: Option<f32>,
    /// Kinetic energy in **kJ/mol**.
    pub kinetic_energy: Option<f32>,
    /// Total energy (potential + kinetic) in **kJ/mol**.
    pub total_energy: Option<f32>,
    /// Simulation box volume in **nm³**.
    pub volume: Option<f32>,
    /// System density in **kg/m³**.
    pub density: Option<f32>,
}

/// Parse a GROMACS `.trr` trajectory file into a sequence of [`GromacsFrame`]
/// values, optionally filtered to a time window (in ps).
///
/// The TRR format is XDR-encoded (big-endian). Each frame contains a header
/// with block sizes that determine both the atom count and the floating-point
/// precision (single vs double) used for that run. Only coordinate blocks are
/// decoded; box, virial, pressure, velocity, and force blocks are skipped.
///
/// Positions are converted from GROMACS native units (nm) to Å on read.
///
/// [Format reference](https://manual.gromacs.org/current/reference-manual/file-formats.html#trr)
pub fn read_trr(
    trr: &Path,
    start_time: Option<f32>,
    end_time: Option<f32>,
) -> io::Result<Vec<GromacsFrame>> {
    const TRR_MAGIC: i32 = 1993;

    let mut r = BufReader::new(File::open(trr)?);
    let mut frames = Vec::new();

    loop {
        // Magic — EOF here is the normal loop termination.
        let magic = {
            let mut b = [0; 4];
            match r.read_exact(&mut b) {
                Err(e) if e.kind() == ErrorKind::UnexpectedEof => break,
                other => other?,
            }
            i32::from_be_bytes(b)
        };
        if magic != TRR_MAGIC {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("TRR bad magic: expected {TRR_MAGIC}, got {magic}"),
            ));
        }

        // Version string: XDR encoding is u32 length + bytes + padding to 4B boundary.
        let ver_len = trr_u32(&mut r)? as usize;
        let ver_pad = (4 - (ver_len % 4)) % 4;
        r.seek(SeekFrom::Current((ver_len + ver_pad) as i64))?;

        // Header — all sizes in bytes.
        let ir_size = trr_i32(&mut r)? as usize;
        let e_size = trr_i32(&mut r)? as usize;
        let box_size = trr_i32(&mut r)? as usize;
        let vir_size = trr_i32(&mut r)? as usize;
        let pres_size = trr_i32(&mut r)? as usize;
        let top_size = trr_i32(&mut r)? as usize;
        let sym_size = trr_i32(&mut r)? as usize;
        let x_size = trr_i32(&mut r)? as usize;
        let v_size = trr_i32(&mut r)? as usize;
        let f_size = trr_i32(&mut r)? as usize;
        let natoms = trr_i32(&mut r)? as usize;
        let _step = trr_i32(&mut r)?;
        let _nre = trr_i32(&mut r)?;

        // Precision: bytes per coordinate component tells us f32 vs f64.
        // Fall back to box_size (9 reals for a 3×3 matrix) if x is absent.
        let double_prec = if x_size > 0 && natoms > 0 {
            x_size == natoms * 3 * 8
        } else if box_size > 0 {
            box_size == 9 * 8
        } else {
            false
        };

        // Time (ps) and lambda — same float width as coordinates.
        let time_ps: f64 = if double_prec {
            trr_f64(&mut r)?
        } else {
            trr_f32(&mut r)? as f64
        };
        let _lambda: f64 = if double_prec {
            trr_f64(&mut r)?
        } else {
            trr_f32(&mut r)? as f64
        };

        // Skip all pre-coordinate data blocks (ir, energy, box, vir, pres, top, sym).
        let pre_skip = ir_size + e_size + box_size + vir_size + pres_size + top_size + sym_size;
        if pre_skip > 0 {
            r.seek(SeekFrom::Current(pre_skip as i64))?;
        }

        // Time-window filter — seek past this frame if it's out of range.
        let in_range = start_time.map_or(true, |t| time_ps >= t as f64)
            && end_time.map_or(true, |t| time_ps <= t as f64);
        if !in_range {
            r.seek(SeekFrom::Current((x_size + v_size + f_size) as i64))?;
            continue;
        }

        // Decode positions (nm → Å).
        let mut atom_posits = Vec::with_capacity(natoms);
        if x_size > 0 {
            for _ in 0..natoms {
                let (x, y, z) = if double_prec {
                    (trr_f64(&mut r)?, trr_f64(&mut r)?, trr_f64(&mut r)?)
                } else {
                    (
                        trr_f32(&mut r)? as f64,
                        trr_f32(&mut r)? as f64,
                        trr_f32(&mut r)? as f64,
                    )
                };

                // These are in nm.
                atom_posits.push(Vec3 { x, y, z });
            }
        }

        // Decode velocities (nm/ps — stored in the same units GROMACS uses internally).
        let mut atom_velocities = Vec::with_capacity(if v_size > 0 { natoms } else { 0 });
        if v_size > 0 {
            for _ in 0..natoms {
                let (x, y, z) = if double_prec {
                    (trr_f64(&mut r)?, trr_f64(&mut r)?, trr_f64(&mut r)?)
                } else {
                    (
                        trr_f32(&mut r)? as f64,
                        trr_f32(&mut r)? as f64,
                        trr_f32(&mut r)? as f64,
                    )
                };
                atom_velocities.push(Vec3 { x, y, z });
            }
        }

        // Skip forces.
        if f_size > 0 {
            r.seek(SeekFrom::Current(f_size as i64))?;
        }

        frames.push(GromacsFrame {
            time: time_ps,
            atom_posits,
            atom_velocities,
            energy: None,
        });
    }

    Ok(frames)
}

/// Write frames to a TRR file, appending if it already exists.
///
/// TRR frames are self-contained, so appending is the canonical way to build
/// up a trajectory incrementally — `gmx traj`, `gmx trjcat`, VMD, and MDAnalysis
/// all read frames sequentially until EOF, with no file-level header to maintain.
///
/// Coordinates are written in GROMACS native units (nm: not Å!).
/// Velocities are written as-is (nm/ps). Both use single-precision (`f32`),
/// which is the GROMACS default. Frames with no velocity data (`atom_velocities`
/// empty or length-mismatched) omit the velocity block (`v_size = 0`).
pub fn write_trr(path: &Path, frames: &[GromacsFrame]) -> io::Result<()> {
    // XDR version string — 12 bytes, already 4B-aligned so no padding needed.
    const VERSION: &[u8] = b"GMX_trn_file";

    let mut w = BufWriter::new(OpenOptions::new().create(true).append(true).open(&*path)?);

    for frame in frames {
        let natoms = frame.atom_posits.len();
        let has_vel = frame.atom_velocities.len() == natoms && natoms > 0;

        let x_size = (natoms * 3 * 4) as i32;
        let v_size = if has_vel { (natoms * 3 * 4) as i32 } else { 0 };

        // Magic
        w.write_all(&1993_i32.to_be_bytes())?;

        // Version string: XDR u32 length + bytes (no padding, len is 4B-aligned)
        w.write_all(&(VERSION.len() as u32).to_be_bytes())?;
        w.write_all(VERSION)?;

        // Header: ir, e, box, vir, pres, top, sym sizes (all 0), then x/v/f,
        // then natoms, step, nre.
        for val in [
            0i32,
            0,
            0,
            0,
            0,
            0,
            0,
            x_size,
            v_size,
            0,
            natoms as i32,
            0,
            0,
        ] {
            w.write_all(&val.to_be_bytes())?;
        }

        // Time (ps) and lambda — single precision, no lambda in our data.
        w.write_all(&(frame.time as f32).to_be_bytes())?;
        w.write_all(&0_f32.to_be_bytes())?;

        // These stay in nm.
        for p in &frame.atom_posits {
            w.write_all(&(p.x as f32).to_be_bytes())?;
            w.write_all(&(p.y as f32).to_be_bytes())?;
            w.write_all(&(p.z as f32).to_be_bytes())?;
        }

        // Velocities: nm/ps, cast to f32.
        if has_vel {
            for v in &frame.atom_velocities {
                w.write_all(&(v.x as f32).to_be_bytes())?;
                w.write_all(&(v.y as f32).to_be_bytes())?;
                w.write_all(&(v.z as f32).to_be_bytes())?;
            }
        }
    }

    Ok(())
}

impl OutputEnergy {
    /// Extract per-frame energy data from a GROMACS binary `.edr` file.
    ///
    /// The EDR format is an XDR-based binary; rather than parsing it directly
    /// this function calls `gmx energy` to convert it to a text XVG file,
    /// then parses that file.
    ///
    /// A broad set of term indices (1–30) is piped to `gmx energy` so that
    /// whatever terms the simulation recorded are captured. GROMACS warns
    /// about out-of-range indices but does not abort.
    ///
    /// Returns `Ok(Vec::new())` when `gmx` is unavailable or the file yields
    /// no usable data, so callers can treat energy as optional.
    pub fn from_edr(path: &Path) -> io::Result<Vec<Self>> {
        let dir = path.parent().unwrap_or(Path::new("."));
        let edr_name = path
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidInput, "invalid EDR path"))?;

        const XVG_OUT: &str = "energy_out.xvg";

        // Feed term indices 1–30 to gmx energy, then "0" to end the
        // selection. Indices beyond the available range produce a warning
        // but are otherwise ignored by GROMACS.
        let selection = b"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 \
                          21 22 23 24 25 26 27 28 29 30 0\n";

        let mut child = match Command::new("gmx")
            .current_dir(dir)
            .args(["energy", "-f", edr_name, "-o", XVG_OUT])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
        {
            Ok(c) => c,
            Err(e) if e.kind() == ErrorKind::NotFound => {
                eprintln!("`gmx` not found; skipping energy parsing");
                return Ok(Vec::new());
            }
            Err(e) => return Err(e),
        };

        if let Some(mut stdin) = child.stdin.take() {
            let _ = stdin.write_all(selection);
        }

        let out = child.wait_with_output()?;
        if !out.status.success() {
            // Non-fatal — GROMACS still writes the XVG for valid terms.
            let msg = String::from_utf8_lossy(&out.stderr);
            eprintln!("gmx energy exited non-zero: {msg}");
        }

        let xvg_path = dir.join(XVG_OUT);
        if !xvg_path.exists() {
            return Ok(Vec::new());
        }

        let text = std::fs::read_to_string(&xvg_path)?;
        parse_xvg(&text)
    }
}

/// A single frame of a GROMACS trajectory. Similar to `DcdFrame`, but includes energy,
/// and velocity, and may use differnet units. (nm vs Å)
#[derive(Clone, Debug)]
pub struct GromacsFrame {
    /// Simulation time in **ps**.
    pub time: f64,
    /// Atom positions in **nm**. Not Å.
    pub atom_posits: Vec<Vec3>,
    /// nm/ps
    pub atom_velocities: Vec<Vec3>,
    /// Energy data for this frame, if the energy recording interval
    /// (`nstenergy`) aligns with the trajectory interval (`nstxout`).
    pub energy: Option<OutputEnergy>,
}

/// The collected output of a `gmx mdrun` run.
#[derive(Clone, Debug, Default)]
pub struct GromacsOutput {
    /// Full text of the GROMACS `.log` file.
    pub log_text: String,
    /// Parsed trajectory frames.
    pub trajectory: Vec<GromacsFrame>,
    /// Number of solute (non-water) atoms per frame. Atoms beyond this index
    /// in each frame's `atom_posits` are water molecules.
    pub solute_atom_count: usize,
}

impl GromacsOutput {
    /// Parse a GROMACS run's output given the raw log text, a multi-frame
    /// `.gro` string produced by `gmx trjconv`, parsed energy records, and
    /// the number of solute atoms.
    ///
    /// Energy records are matched to frames by simulation time; frames whose
    /// time does not align with any energy record will have `energy: None`.
    pub fn new(
        log_text: String,
        // traj_gro: &str,
        mut trajectory: Vec<GromacsFrame>,
        energies: Vec<OutputEnergy>,
        solute_atom_count: usize,
    ) -> io::Result<Self> {
        // let mut trajectory/ = parse_multi_gro(traj_gro)?;
        attach_energies(&mut trajectory, energies);

        Ok(Self {
            log_text,
            trajectory,
            solute_atom_count,
        })
    }
}

/// Parse a multi-model `.gro` file (as produced by `gmx trjconv`) into a
/// sequence of [`GromacsFrame`] values. We prefer to parse trajectories from TRR files, but this
/// can serve as a backup.
///
/// Each model block has the layout:
/// ```text
/// <title — may contain "t= <time_ps>">
/// <natoms>
/// <atom lines…>
/// <box vector line>
/// ```
pub fn parse_multi_gro(text: &str) -> io::Result<Vec<GromacsFrame>> {
    let mut frames = Vec::new();
    let mut lines = text.lines().peekable();

    while lines.peek().is_some() {
        // Title line
        let title = match lines.next() {
            Some(l) => l,
            None => break,
        };

        // Extract time from title if present, e.g. "MD of system, t= 0.00000"
        let time_ps = extract_time(title);

        // Atom count
        let natoms_str = lines.next().ok_or_else(|| {
            io::Error::new(ErrorKind::UnexpectedEof, "Expected atom count after title")
        })?;
        let natoms: usize = natoms_str
            .trim()
            .split_whitespace()
            .next()
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Empty atom count line"))?
            .parse()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse atom count"))?;

        let mut posits = Vec::with_capacity(natoms);

        for _ in 0..natoms {
            let line = lines.next().ok_or_else(|| {
                io::Error::new(ErrorKind::UnexpectedEof, "Unexpected end of atom block")
            })?;

            // GRO fixed-width: x at cols 20-28, y at 28-36, z at 36-44 (nm)
            let x_nm = parse_col(line, 20, 28)?;
            let y_nm = parse_col(line, 28, 36)?;
            let z_nm = parse_col(line, 36, 44)?;

            // These remain in nm.
            posits.push(Vec3 {
                x: x_nm,
                y: y_nm,
                z: z_nm,
            });
        }

        // Box vector line — consume it even if we don't use the value.
        lines.next();

        frames.push(GromacsFrame {
            time: time_ps,
            atom_posits: posits,
            atom_velocities: Vec::new(),
            energy: None,
        });
    }

    Ok(frames)
}

/// Attach energy records to trajectory frames by matching on simulation time.
///
/// The energy recording interval (`nstenergy`) may differ from the trajectory
/// interval (`nstxout`), so some frames legitimately have no energy entry.
/// Times are matched to the nearest millisecond (0.001 ps) to absorb any
/// floating-point rounding in the GRO title and XVG data lines.
fn attach_energies(frames: &mut Vec<GromacsFrame>, energies: Vec<OutputEnergy>) {
    if energies.is_empty() {
        return;
    }
    // Key: time rounded to nearest 0.001 ps expressed as integer (avoids f64 hashing).
    let map: HashMap<i64, OutputEnergy> = energies
        .into_iter()
        .map(|e| ((e.time * 1_000.0).round() as i64, e))
        .collect();

    for frame in frames.iter_mut() {
        let key = (frame.time * 1_000.0).round() as i64;
        frame.energy = map.get(&key).cloned();
    }
}

// ---------------------------------------------------------------------------
// XVG parsing
// ---------------------------------------------------------------------------

/// Parse an XVG file produced by `gmx energy` into a sequence of
/// [`OutputEnergy`] values.
///
/// XVG layout:
/// - Lines starting with `#` are comments — skipped.
/// - Lines starting with `@` are metadata; `@ sN legend "Name"` lines map
///   column indices to energy term names.
/// - All other non-empty lines are data rows: time followed by term values,
///   all space-separated.
fn parse_xvg(text: &str) -> io::Result<Vec<OutputEnergy>> {
    // Column index (0-based after the time column) → GROMACS term name.
    let mut col_names: Vec<String> = Vec::new();
    let mut rows: Vec<Vec<f64>> = Vec::new();

    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed.starts_with('@') {
            if let Some(name) = parse_xvg_legend(trimmed) {
                col_names.push(name);
            }
            continue;
        }
        // Data row: space-separated floats.
        let vals: Vec<f64> = trimmed
            .split_whitespace()
            .filter_map(|t| t.parse().ok())
            .collect();
        if !vals.is_empty() {
            rows.push(vals);
        }
    }

    let energies = rows
        .into_iter()
        .map(|vals| {
            let time = vals.first().copied().unwrap_or(0.0);
            let mut e = OutputEnergy {
                time,
                ..OutputEnergy::default()
            };
            for (col, name) in col_names.iter().enumerate() {
                let val = vals.get(col + 1).copied().map(|v| v as f32);
                match name.as_str() {
                    "Temperature" => e.temperature = val,
                    "Pressure" => e.pressure = val,
                    "Potential" => e.potential_energy = val,
                    n if n.starts_with("Kinetic") => e.kinetic_energy = val,
                    n if n.starts_with("Total") => e.total_energy = val,
                    "Volume" => e.volume = val,
                    "Density" => e.density = val,
                    _ => {}
                }
            }
            e
        })
        .collect();

    Ok(energies)
}

/// Extract the legend name from a `@ sN legend "Name"` XVG directive.
/// Returns `None` for any other `@` line.
fn parse_xvg_legend(line: &str) -> Option<String> {
    // Strip leading '@' and whitespace.
    let rest = line.strip_prefix('@')?.trim();
    // Must start with 's' followed by a column index (digits).
    if !rest.starts_with('s') {
        return None;
    }
    // Strip the sN token to reach the keyword.
    let after_index = rest.trim_start_matches(|c: char| c == 's' || c.is_ascii_digit());
    let rest = after_index.trim();
    // Next token must be "legend".
    let rest = rest.strip_prefix("legend")?.trim();
    // The remainder is the quoted name.
    let name = rest.trim_matches('"');
    if name.is_empty() {
        None
    } else {
        Some(name.to_string())
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Extract simulation time (ps) from a GRO title line such as
/// `"Protein in water, t= 12.00000"`. Returns `0.0` if not found.
fn extract_time(title: &str) -> f64 {
    // Look for "t=" or "t =" anywhere in the title.
    for sep in ["t=", "t ="] {
        if let Some(pos) = title.find(sep) {
            let rest = title[pos + sep.len()..].trim();
            // Take the first whitespace-delimited token.
            if let Some(tok) = rest.split_whitespace().next() {
                if let Ok(v) = tok.parse::<f64>() {
                    return v;
                }
            }
        }
    }
    0.0
}

/// Parse a fixed-width column slice as f64. Returns 0.0 on an empty slice
/// rather than an error, to tolerate short lines at end of file.
fn parse_col(line: &str, start: usize, end: usize) -> io::Result<f64> {
    let len = line.len();
    if start >= len {
        return Ok(0.0);
    }
    let end = end.min(len);
    let s = line[start..end].trim();
    if s.is_empty() {
        return Ok(0.0);
    }
    s.parse::<f64>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, format!("Cannot parse '{s}' as f64")))
}

// ---------------------------------------------------------------------------
// TRR / XDR low-level readers (big-endian IEEE 754)
// ---------------------------------------------------------------------------

fn trr_i32(r: &mut impl Read) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_be_bytes(b))
}

fn trr_u32(r: &mut impl Read) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_be_bytes(b))
}

fn trr_f32(r: &mut impl Read) -> io::Result<f32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(f32::from_be_bytes(b))
}

fn trr_f64(r: &mut impl Read) -> io::Result<f64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_be_bytes(b))
}
