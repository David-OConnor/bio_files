//! For reading and writing the [DCD file format](https://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.app.dcdfile.DCDFile.html), for molecular dynamics simulations.
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

#[derive(Clone, Debug)]
pub struct DcdUnitCell {
    pub bounds_low: Vec3,
    pub bounds_high: Vec3,
}

impl DcdUnitCell {
    fn to_dcd_six(&self) -> [f64; 6] {
        let a = (self.bounds_high.x - self.bounds_low.x) as f64;
        let b = (self.bounds_high.y - self.bounds_low.y) as f64;
        let c = (self.bounds_high.z - self.bounds_low.z) as f64;

        // Orthorhombic: angles are 90 degrees.
        // Use the common X-PLOR ordering on disk: [A, gamma, B, beta, alpha, C].
        // For 90/90/90 the permutations don’t change meaning.
        [a, 90.0, b, 90.0, 90.0, c]
    }

    fn from_dcd_six(six: [f64; 6]) -> Self {
        // We only store bounds_low/high, but DCD stores lengths/angles, not an origin.
        // So we reconstruct a box from (0,0,0) to (A,B,C).
        let a = six[0] as f32;
        let b = six[2] as f32;
        let c = six[5] as f32;

        Self {
            bounds_low: Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            bounds_high: Vec3 { x: a, y: b, z: c },
        }
    }
}

#[derive(Clone, Debug)]
pub struct DcdFrame {
    /// fs
    pub time: f64,
    /// Å
    pub atom_posits: Vec<Vec3>,
    /// Also called Periodic box. This is often the bounds of the simulation, with solvents wrapping
    /// around periodic boundary conditions, and long-range forces computed across this.
    pub unit_cell: DcdUnitCell,
}

/// Represents a molecular dynamics trajectory, and contains fields specific to DCD files.
/// This is a minimal structure that mainly keeps track of atom positions. It doesn't include velocities,
/// or data like energy, pressure, and temperature.
#[derive(Clone, Debug)]
pub struct DcdTrajectory {
    pub frames: Vec<DcdFrame>,
}

impl DcdTrajectory {
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

        let has_unitcell = icntrl[19] != 0 && icntrl[10] != 0;

        // Delta is at bytes 36..40 after the "CORD"
        let delta = f32::from_le_bytes(hdr[4 + 36..4 + 40].try_into().unwrap()) as f64;

        skip_title_record(&mut r)?;

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

        let mut unit_cell = DcdUnitCell {
            bounds_low: Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            bounds_high: Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
        };

        for i in 0..nset_total {
            if has_unitcell {
                unit_cell = read_unit_cell_record(&mut r)?;
            }

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

            let istart = icntrl[1] as f64;
            let nsavc = icntrl[2] as f64;

            frames.push(DcdFrame {
                time: (istart + (i as f64) * nsavc) * delta,
                atom_posits,
                unit_cell: unit_cell.clone(),
            });
        }

        Ok(Self { frames })
    }

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
            // .truncate(true) // todo: QC this.
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

            // Writing 1 to index 10 this means "extra block present", and is required when we
            // pass the unit cell.
            icntrl[10] = 1;
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
            if l1 < 84 || l1 > 1024 * 1024 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "unreasonable DCD header size",
                ));
            }
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

            let has_unitcell = icntrl[19] != 0 && icntrl[10] != 0;
            if !has_unitcell {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "existing DCD does not have unit cell blocks enabled (icntrl[10]==0)",
                ));
            }

            // Skip the title
            skip_title_record(&mut f)?;

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
                let mut push = |v: &[Vec3]| {
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

                write_unit_cell_record(&mut f, &frame.unit_cell)?;
                write_record(&mut f, xb)?;
                write_record(&mut f, yb)?;
                write_record(&mut f, zb)?;
            }

            // Update NSET in header (payload offset = 4-byte marker + 4 for "CORD")
            let new_nset = cur_nset
                .checked_add(self.frames.len() as i32)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "NSET overflow"))?;

            // We are at arbitrary place; seek to start and re-read first record marker.
            f.seek(SeekFrom::Start(0))?;
            let l1 = read_u32_le(&mut f)? as u64;
            if l1 < 84 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "unexpected DCD header record length",
                ));
            }

            let mut magic = [0u8; 4];
            f.read_exact(&mut magic)?;
            if &magic != b"CORD" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "not a CORD/DCD file",
                ));
            }

            // file offset of NSET: 4-byte leading length + 4 bytes of "CORD"
            let nset_off = 8u64;
            f.seek(SeekFrom::Start(nset_off))?;
            f.write_all(&new_nset.to_le_bytes())?;

            f.flush()?;

            return Ok(());
        }

        let mut xs = vec![0.; n_atoms];
        let mut ys = vec![0.; n_atoms];
        let mut zs = vec![0.; n_atoms];

        for frame in &self.frames {
            let mut i = 0;
            let mut push = |v: &[Vec3]| {
                for p in v {
                    xs[i] = p.x;
                    ys[i] = p.y;
                    zs[i] = p.z;
                    i += 1;
                }
            };
            push(&frame.atom_posits);

            let xb = unsafe { core::slice::from_raw_parts(xs.as_ptr() as *const u8, xs.len() * 4) };
            let yb = unsafe { core::slice::from_raw_parts(ys.as_ptr() as *const u8, ys.len() * 4) };
            let zb = unsafe { core::slice::from_raw_parts(zs.as_ptr() as *const u8, zs.len() * 4) };

            write_unit_cell_record(&mut f, &frame.unit_cell)?;
            write_record(&mut f, xb)?;
            write_record(&mut f, yb)?;
            write_record(&mut f, zb)?;
        }

        f.flush()
    }

    /// Converts from a GROMACS XTC file. [MDTraj](https://www.mdtraj.org/1.9.8.dev0/index.html) must
    /// be installed, and available on the system path. Install with `pip install mdtraj`.
    pub fn load_xtc(path: &Path) -> io::Result<Self> {
        let temp_file = "temp_dcd.dcd";

        let out = Command::new("mdconvert")
            .args(["-o", temp_file, path.to_str().unwrap()])
            .output()?;

        if !out.status.success() {
            let stderr_str = String::from_utf8_lossy(&out.stderr);
            return Err(io::Error::other(format!(
                "Problem parsing XTC file: {}",
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

fn write_unit_cell_record<W: Write>(w: &mut W, unit_cell: &DcdUnitCell) -> io::Result<()> {
    let six = unit_cell.to_dcd_six();

    let mut payload = [0u8; 48];
    for (i, v) in six.iter().enumerate() {
        let b = v.to_le_bytes();
        payload[i * 8..i * 8 + 8].copy_from_slice(&b);
    }

    write_record(w, &payload)
}

fn read_unit_cell_record<R: Read>(r: &mut R) -> io::Result<DcdUnitCell> {
    let b = read_record(r)?;
    if b.len() != 48 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unexpected unit cell record size (expected 48 bytes)",
        ));
    }

    let mut six = [0f64; 6];
    for i in 0..6 {
        let j = i * 8;
        six[i] = f64::from_le_bytes(b[j..j + 8].try_into().unwrap());
    }

    Ok(DcdUnitCell::from_dcd_six(six))
}

fn skip_title_record<R: Read>(r: &mut R) -> io::Result<()> {
    let b = read_record(r)?;
    if b.len() < 4 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "title record too short",
        ));
    }
    let ntitle = i32::from_le_bytes(b[0..4].try_into().unwrap());
    if ntitle < 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid NTITLE"));
    }
    // Expected size is 4 + 80*ntitle. Some writers may pad; you can allow >=.
    let expected = 4usize + (ntitle as usize) * 80;
    if b.len() < expected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "truncated title record",
        ));
    }
    Ok(())
}
