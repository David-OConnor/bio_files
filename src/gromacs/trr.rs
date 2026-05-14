//! TRR-specific trajectory reading and writing. See `output.rs` for general
//! trajectory data.

use std::{
    fs::{File, OpenOptions},
    io,
    io::{BufReader, BufWriter, ErrorKind, Read, Seek, SeekFrom, Write},
    path::Path,
    time::Instant,
};

use lin_alg::f64::Vec3;

use crate::{FrameSlice, gromacs::GromacsFrame};

/// This isn't stored directly in the file; transverse frame headers to collect this. Not too
/// slow, as it doesn't load the frame coordinates.
pub struct TrrMetadata {
    pub num_atoms: usize,
    pub num_frames: usize,
    pub start_step: f32,
    pub save_interval_steps: usize,
    pub dt: f32,
    pub end_time: f32,
}

impl TrrMetadata {
    /// Scan all frame headers without decoding coordinate data to collect metadata.
    pub fn read(path: &Path) -> io::Result<Self> {
        let start = Instant::now();
        println!("Starting TRR frame scan to load metadata...");
        const TRR_MAGIC: i32 = 1993;

        let mut r = BufReader::new(File::open(path)?);

        let mut num_atoms = 0;
        let mut num_frames = 0;

        let mut start_step = 0i32;
        let mut first_time = 0.0;
        let mut dt = 0.0;
        let mut save_interval_steps = 0usize;
        let mut end_time = 0.0;

        loop {
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

            // Version string: u32 length + bytes + padding to 4B boundary.
            let ver_len = trr_u32(&mut r)? as usize;
            let ver_pad = (4 - (ver_len % 4)) % 4;
            r.seek(SeekFrom::Current((ver_len + ver_pad) as i64))?;

            // Header fields.
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

            let step = trr_i32(&mut r)?;
            let _nre = trr_i32(&mut r)?;

            let double_prec = if x_size > 0 && natoms > 0 {
                x_size == natoms * 3 * 8
            } else if box_size > 0 {
                box_size == 9 * 8
            } else {
                false
            };

            let time_ps: f64 = if double_prec {
                trr_f64(&mut r)?
            } else {
                trr_f32(&mut r)? as f64
            };
            // lambda — same width, skip.
            if double_prec {
                trr_f64(&mut r)?;
            } else {
                trr_f32(&mut r)?;
            }

            if num_frames == 0 {
                num_atoms = natoms;
                start_step = step;
                first_time = time_ps;
            } else if num_frames == 1 {
                dt = (time_ps - first_time) as f32;
                save_interval_steps = (step - start_step).unsigned_abs() as usize;
            }
            end_time = time_ps as f32;
            num_frames += 1;

            // Seek past all data blocks.
            let skip = ir_size
                + e_size
                + box_size
                + vir_size
                + pres_size
                + top_size
                + sym_size
                + x_size
                + v_size
                + f_size;
            if skip > 0 {
                r.seek(SeekFrom::Current(skip as i64))?;
            }
        }

        let elapsed = start.elapsed().as_millis();
        println!("Scan complete in {elapsed} ms");

        Ok(Self {
            num_atoms,
            num_frames,
            start_step: start_step as f32,
            save_interval_steps,
            dt,
            end_time,
        })
    }
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
pub fn read_trr(trr: &Path, slice: FrameSlice) -> io::Result<Vec<GromacsFrame>> {
    const TRR_MAGIC: i32 = 1993;

    let mut r = BufReader::new(File::open(trr)?);
    let mut frames = Vec::new();
    let mut frame_idx: usize = 0;

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

        // Slice filter — seek past this frame if it's out of range.
        let in_range = match slice {
            FrameSlice::Time { start, end } => {
                start.map_or(true, |t| time_ps >= t) && end.map_or(true, |t| time_ps <= t)
            }
            FrameSlice::Index { start, end } => {
                start.map_or(true, |s| frame_idx >= s) && end.map_or(true, |e| frame_idx <= e)
            }
        };
        frame_idx += 1;
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

        // Decode forces (kJ/(mol·nm) — GROMACS native units).
        let mut atom_forces = Vec::with_capacity(if f_size > 0 { natoms } else { 0 });
        if f_size > 0 {
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
                atom_forces.push(Vec3 { x, y, z });
            }
        }

        frames.push(GromacsFrame {
            time: time_ps,
            atom_posits,
            atom_velocities,
            atom_forces,
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
        let has_forces = frame.atom_forces.len() == natoms && natoms > 0;

        let x_size = (natoms * 3 * 4) as i32;
        let v_size = if has_vel { (natoms * 3 * 4) as i32 } else { 0 };
        let f_size = if has_forces {
            (natoms * 3 * 4) as i32
        } else {
            0
        };

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
            f_size,
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

        // Forces: kJ/(mol·nm), cast to f32.
        if has_forces {
            for f in &frame.atom_forces {
                w.write_all(&(f.x as f32).to_be_bytes())?;
                w.write_all(&(f.y as f32).to_be_bytes())?;
                w.write_all(&(f.z as f32).to_be_bytes())?;
            }
        }
    }

    Ok(())
}

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
