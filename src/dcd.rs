//! For reading and writing the DCD file format, for molecular dynamics simulations.
//! This format is used, for example, by VMD.

use std::{
    fs,
    fs::{File, OpenOptions},
    io,
    io::{BufReader, Read, Seek, SeekFrom, Write},
    path::Path,
    process::Command,
};

use lin_alg::f32::Vec3;

use crate::DensityMap;

pub struct DcdFrame {
    /// fs
    pub time: f64,
    /// Å
    pub atom_posits: Vec<Vec3>,
}

/// Represents a molecular dynamics trajectory, and contains fields specific to DCD files.
/// This is a minimal structure that mainly keeps track of atom positions. It doesn't include velocities,
/// or data like energy, pressure, and temperature.
pub struct DcdTrajectory {
    pub frames: Vec<DcdFrame>,
}

impl DcdTrajectory {
    /// Create or append snapshots to a DCD file. This is a common trajectory/reporter format
    /// used by other software, including OpenMM and VMD.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        if self.frames.is_empty() {
            return Ok(());
        }

        let first = &self.frames[0];

        let mut n_atoms = first.atom_posits.len();

        for s in self.frames.iter().skip(1) {
            let mut n2 = s.atom_posits.len();

            if n2 != n_atoms {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "inconsistent atom counts",
                ));
            }
        }

        let mut f = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true) // todo: QC this.
            .open(path)?;

        let file_len = f.metadata()?.len();

        // write header if new/empty
        if file_len == 0 {
            let nsets = self.frames.len() as i32;
            let istart: i32 = 0;
            let nsavc: i32 = 1;

            let delta: f32 = if self.frames.len() >= 2 {
                (self.frames[1].time - self.frames[0].time) as f32
            } else {
                0.0
            };

            let mut header = Vec::with_capacity(84);
            header.extend_from_slice(b"CORD");

            let mut icntrl = [0i32; 20];
            icntrl[0] = nsets;
            icntrl[1] = istart;
            icntrl[2] = nsavc;
            icntrl[8] = 0;
            icntrl[10] = 0;
            icntrl[11] = 0;
            icntrl[19] = 1;

            for v in icntrl {
                header.extend_from_slice(&v.to_le_bytes());
            }

            header[4 + 36..4 + 40].copy_from_slice(&delta.to_le_bytes());
            write_record(&mut f, &header)?;

            let title = format!("Created by Dynamics  NATOMS={}  NFRAMES={}", n_atoms, nsets);
            let mut line = [0u8; 80];
            let tb = title.as_bytes();
            let n = tb.len().min(80);
            line[..n].copy_from_slice(&tb[..n]);

            let mut title_block = Vec::with_capacity(4 + 80);
            title_block.extend_from_slice(&(1i32).to_le_bytes());
            title_block.extend_from_slice(&line);
            write_record(&mut f, &title_block)?;

            let mut natom_block = Vec::with_capacity(4);
            natom_block.extend_from_slice(&(n_atoms as i32).to_le_bytes());
            write_record(&mut f, &natom_block)?;
        } else {
            // verify header and NATOM; compute current NSET; then append and bump NSET
            f.seek(SeekFrom::Start(0))?;
            let l1 = read_u32_le(&mut f)?;
            let mut hdr = vec![0u8; l1 as usize];

            f.read_exact(&mut hdr)?;
            let l1e = read_u32_le(&mut f)?;
            if l1e != l1 || &hdr[0..4] != b"CORD" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "not a CORD/DCD file",
                ));
            }

            // Current NSET and flags
            let mut icntrl = [0i32; 20];
            for (i, item) in icntrl.iter_mut().enumerate() {
                let off = 4 + i * 4;
                *item = i32::from_le_bytes(hdr[off..off + 4].try_into().unwrap());
            }
            let cur_nset = icntrl[0];

            // Skip the title
            let l2 = read_u32_le(&mut f)?;
            f.seek(SeekFrom::Current(l2 as i64))?;
            let l2e = read_u32_le(&mut f)?;
            if l2e != l2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "corrupt title block",
                ));
            }

            // Read NATOM
            let l3 = read_u32_le(&mut f)?;
            if l3 != 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "unexpected NATOM record length",
                ));
            }

            let mut nb = [0u8; 4];
            f.read_exact(&mut nb)?;

            let natom_existing = i32::from_le_bytes(nb) as usize;
            let l3e = read_u32_le(&mut f)?;
            if l3e != l3 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "corrupt NATOM block",
                ));
            }

            if natom_existing != n_atoms {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "atom count mismatch with existing file",
                ));
            }

            f.seek(SeekFrom::End(0))?;

            let mut xs = vec![0f32; n_atoms];
            let mut ys = vec![0f32; n_atoms];
            let mut zs = vec![0f32; n_atoms];

            for frame in &self.frames {
                let mut i = 0usize;
                let mut push = |v: &Vec<Vec3>| {
                    for p in v {
                        xs[i] = p.x;
                        ys[i] = p.y;
                        zs[i] = p.z;
                        i += 1;
                    }
                };
                push(&frame.atom_posits);

                let xb =
                    unsafe { core::slice::from_raw_parts(xs.as_ptr() as *const u8, xs.len() * 4) };
                let yb =
                    unsafe { core::slice::from_raw_parts(ys.as_ptr() as *const u8, ys.len() * 4) };
                let zb =
                    unsafe { core::slice::from_raw_parts(zs.as_ptr() as *const u8, zs.len() * 4) };

                write_record(&mut f, xb)?;
                write_record(&mut f, yb)?;
                write_record(&mut f, zb)?;
            }

            // Update NSET in header (payload offset = 4-byte marker + 4 for "CORD")
            let new_nset = cur_nset
                .checked_add(self.frames.len() as i32)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "NSET overflow"))?;

            f.seek(SeekFrom::Start(8))?;
            f.write_all(&new_nset.to_le_bytes())?;
            f.flush()?;

            return Ok(());
        }

        let mut xs = vec![0.; n_atoms];
        let mut ys = vec![0.; n_atoms];
        let mut zs = vec![0.; n_atoms];

        for snap in &self.frames {
            let mut i = 0;
            let mut push = |v: &Vec<Vec3>| {
                for p in v {
                    xs[i] = p.x;
                    ys[i] = p.y;
                    zs[i] = p.z;
                    i += 1;
                }
            };
            push(&snap.atom_posits);

            let xb = unsafe { core::slice::from_raw_parts(xs.as_ptr() as *const u8, xs.len() * 4) };
            let yb = unsafe { core::slice::from_raw_parts(ys.as_ptr() as *const u8, ys.len() * 4) };
            let zb = unsafe { core::slice::from_raw_parts(zs.as_ptr() as *const u8, zs.len() * 4) };

            write_record(&mut f, xb)?;
            write_record(&mut f, yb)?;
            write_record(&mut f, zb)?;
        }

        f.flush()
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let f = File::open(path)?;
        let mut r = BufReader::new(f);

        // Header
        let hdr = read_record(&mut r)?;
        if hdr.len() < 84 || &hdr[0..4] != b"CORD" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Not a CORD/DCD file",
            ));
        }
        let mut icntrl = [0i32; 20];
        for (i, item) in icntrl.iter_mut().enumerate() {
            let off = 4 + i * 4;
            *item = i32::from_le_bytes(hdr[off..off + 4].try_into().unwrap());
        }
        let nset_total = icntrl[0] as usize;

        // Delta is at bytes 36..40 after the "CORD"
        let delta = f32::from_le_bytes(hdr[4 + 36..4 + 40].try_into().unwrap()) as f64;

        // Title (ignored)
        let _ = read_record(&mut r)?;

        // NATOM
        let natom_block = read_record(&mut r)?;
        if natom_block.len() != 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Unexpected NATOM block size",
            ));
        }
        let n_atoms = i32::from_le_bytes(natom_block[0..4].try_into().unwrap()) as usize;

        let mut frames = Vec::with_capacity(nset_total);

        for i in 0..nset_total {
            let xb = read_record(&mut r)?;
            let yb = read_record(&mut r)?;
            let zb = read_record(&mut r)?;

            if xb.len() != 4 * n_atoms || yb.len() != 4 * n_atoms || zb.len() != 4 * n_atoms {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Coordinate block size mismatch",
                ));
            }

            let xs = f32s_from_le_bytes(&xb)?;
            let ys = f32s_from_le_bytes(&yb)?;
            let zs = f32s_from_le_bytes(&zb)?;

            let mut atom_posits = Vec::with_capacity(n_atoms);
            for k in 0..n_atoms {
                atom_posits.push(Vec3 {
                    x: xs[k],
                    y: ys[k],
                    z: zs[k],
                });
            }

            frames.push(DcdFrame {
                time: (i as f64) * delta,
                atom_posits,
            });
        }

        Ok(Self { frames })
    }

    /// Converts from a GROMACS XTC file. `[MDTjac](https://www.mdtraj.org/1.9.8.dev0/index.html) must
    /// be installed, and available on the system path. Install with `pip install mdtraj`.
    pub fn load_xtc(path: &Path) -> io::Result<Self> {
        let temp_file = "temp_dcd.dcd";

        let out = Command::new("mdconvert")
            .args(["-o", path.to_str().unwrap(), temp_file])
            .output()?;

        if !out.status.success() {
            let stderr_str = String::from_utf8_lossy(&out.stderr);
            return Err(io::Error::other(format!(
                "Problem parsing DCD file: {}",
                stderr_str
            )));
        }

        let map = Self::load(Path::new(temp_file))?;

        fs::remove_file(Path::new(temp_file))?;

        Ok(map)
    }

    /// Saves this trajectory as a GROMACS XTC file via an intermediate DCD file.
    ///
    /// Requires `mdconvert` from MDTraj to be installed and on PATH:
    /// `pip install mdtraj`
    pub fn save_xtc(&self, out_path: &Path) -> io::Result<()> {
        let temp_file = "temp_dcd.dcd";

        // Write intermediate DCD using our own writer.
        self.save(Path::new(temp_file))?;

        // Convert DCD -> XTC using mdconvert.
        let out = Command::new("mdconvert")
            .args([temp_file, "-o", out_path.to_str().unwrap()])
            .output()?;

        // Always try to remove temp file, even on error.
        let _ = fs::remove_file(Path::new(temp_file));

        if !out.status.success() {
            let stderr_str = String::from_utf8_lossy(&out.stderr);
            return Err(io::Error::other(format!(
                "Problem writing XTC file via mdconvert: {}",
                stderr_str
            )));
        }

        Ok(())
    }
}

/// A wrapper for writing a DCD record: Payload sandwhiched by lenth.
fn write_record<W: Write>(w: &mut W, payload: &[u8]) -> io::Result<()> {
    let len = payload.len() as u32;

    w.write_all(&len.to_le_bytes())?;
    w.write_all(payload)?;
    w.write_all(&len.to_le_bytes())
}

fn read_u32_le<R: Read>(r: &mut R) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_record<R: Read>(r: &mut R) -> io::Result<Vec<u8>> {
    let len = read_u32_le(r)? as usize;
    let mut payload = vec![0u8; len];
    r.read_exact(&mut payload)?;
    let len_end = read_u32_le(r)? as usize;
    if len_end != len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "record length mismatch",
        ));
    }
    Ok(payload)
}

fn f32s_from_le_bytes(b: &[u8]) -> io::Result<Vec<f32>> {
    if !b.len().is_multiple_of(4) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "float block not multiple of 4",
        ));
    }

    let n = b.len() / 4;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let j = 4 * i;
        out.push(f32::from_le_bytes(b[j..j + 4].try_into().unwrap()));
    }
    Ok(out)
}
