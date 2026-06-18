//! Parsing GROMACS simulation output. See also `trr.rs` for information on that format
//! specifically. Note that GROMACS does not support the DCD format; it uses TRR (full precision),
//! and XTC (Lossily compressed). We currently don't support XTC, but would like to in the future.
//! It is more complicated to read and write.
//!
//! After an `mdrun`, we use `gmx trjconv` to export the trajectory as a
//! multi-model `.gro` file, then parse each frame here into [`GromacsFrame`].
//! Energy data is extracted from the binary `.edr` file by calling
//! `gmx energy` to convert it to an XVG text file, which is then parsed.

use std::{
    collections::HashMap,
    io::{self, ErrorKind, Write as _},
    path::Path,
    process::Stdio,
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

impl OutputEnergy {
    /// Extract per-frame energy data from a GROMACS binary `.edr` file.
    ///
    /// The EDR format is an XDR-based binary; rather than parsing it directly
    /// this function calls `gmx energy` to convert it to a text XVG file,
    /// then parses that file.
    ///
    /// Stable thermodynamic term names are piped to `gmx energy`, avoiding
    /// topology-dependent menu indices.
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

        // Select by stable term names rather than menu indices, which vary with
        // the topology and enabled algorithms.
        let selection =
            b"Temperature\nPressure\nPotential\nKinetic-En.\nTotal-Energy\nVolume\nDensity\n0\n";

        let mut child = match super::gmx_command()
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
    /// Forces in **kJ/(mol·nm)** (GROMACS native units).
    pub atom_forces: Vec<Vec3>,
    /// Energy data for this frame, if the energy recording interval
    /// (`nstenergy`) aligns with the trajectory interval (`nstxout`).
    pub energy: Option<OutputEnergy>,
}

/// The collected output of a `gmx mdrun` run.
#[derive(Clone, Debug, Default)]
pub struct GromacsOutput {
    /// Full text of the GROMACS `.log` file.
    pub log_text: String,
    /// We flag this manually from a rust error code.
    pub setup_failure: bool,
    /// Parsed trajectory frames.
    pub trajectory: Vec<GromacsFrame>,
    /// Number of solute (non-water) atoms per frame. Atoms beyond this index
    /// in each frame's `atom_posits` are water molecules.
    pub solute_atom_count: usize,
    /// Path to the saved input GRO file (mols_in_N.gro), if written.
    pub gro_path: Option<std::path::PathBuf>,
    /// Path to the saved TRR trajectory (traj_N.trr), if written.
    pub trr_path: Option<std::path::PathBuf>,
    /// Path to the saved XTC trajectory (traj_N.xtc), if written.
    pub xtc_path: Option<std::path::PathBuf>,
    /// Final coordinates written by `mdrun`, suitable for continuing in another run.
    pub final_gro_text: Option<String>,
    /// Final topology text, including any counter-ion substitutions made during setup.
    pub final_topology_text: Option<String>,
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
        mut trajectory: Vec<GromacsFrame>,
        energies: Vec<OutputEnergy>,
        solute_atom_count: usize,
    ) -> io::Result<Self> {
        attach_energies(&mut trajectory, energies);

        Ok(Self {
            log_text,
            setup_failure: false,
            trajectory,
            solute_atom_count,
            gro_path: None,
            trr_path: None,
            xtc_path: None,
            final_gro_text: None,
            final_topology_text: None,
        })
    }
}

/// Parse a multi-model `.gro` file (as produced by `gmx trjconv`) into a
/// sequence of [`GromacsFrame`] values. We prefer to parse trajectories from TRR files, but this
/// can serve as a backup, for trajectories stored in `gro` format. This format is larger than TRR and
/// XTC. This reads all frames in the file, unlike our APIs for other traj formats, which read
/// ranges of frames. (From potentially large files)
///
/// Each model block has the layout:
/// ```text
/// <title — may contain "t= <time_ps>">
/// <natoms>
/// <atom lines…>
/// <box vector line>
/// ```
pub fn parse_gro_traj(text: &str) -> io::Result<Vec<GromacsFrame>> {
    let mut frames = Vec::new();
    let mut lines = text.lines().peekable();

    while lines.peek().is_some() {
        let title = match lines.next() {
            Some(l) => l,
            None => break,
        };

        let time_ps = extract_time(title);

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

        // Box vector line - consume it even if we don't use the value.
        lines.next();

        frames.push(GromacsFrame {
            time: time_ps,
            atom_posits: posits,
            atom_velocities: Vec::new(),
            atom_forces: Vec::new(),
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
