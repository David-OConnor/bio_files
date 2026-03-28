//! For reading and writing the [DCD file format](https://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.app.dcdfile.DCDFile.html), for molecular dynamics simulations.
//! This format is used, for example, by VMD.

use std::{
    env, fs,
    fs::{File, OpenOptions},
    io,
    io::{BufReader, Read, Seek, SeekFrom, Write},
    path::Path,
    process::Command,
};

use lin_alg::f32::Vec3;

// ---------------------------------------------------------------------------
// Embedded Python helpers for XTC I/O via MDTraj
// ---------------------------------------------------------------------------
//
// Wire format shared between Rust and Python (all little-endian):
//   Header : n_frames (i32), n_atoms (i32)
//   Per frame:
//     time_ps (f32)
//     box_a, box_b, box_c (f32, Å — orthorhombic diagonal lengths)
//     xyz: n_atoms × 3 × f32 (Å, interleaved x0,y0,z0, x1,y1,z1, …)
//
// Positions and box lengths are in Å in Rust; the Python scripts convert
// to/from MDTraj's native nm.

/// Python script executed by `write_xtc`.  Reads the wire-format binary
/// written by Rust and appends/creates an XTC file via MDTraj.
const WRITE_XTC_PY: &str = r#"
import sys, struct
import numpy as np
try:
    from mdtraj.formats import XTCTrajectoryFile
except ImportError:
    sys.exit("mdtraj not found. Install with: pip install mdtraj")

mode, xtc_path, bin_path = sys.argv[1], sys.argv[2], sys.argv[3]
ANG_TO_NM = 0.1

with open(bin_path, 'rb') as bf:
    n_frames, n_atoms = struct.unpack('<ii', bf.read(8))
    xyz  = np.empty((n_frames, n_atoms, 3), dtype=np.float32)
    time_arr = np.empty(n_frames, dtype=np.float32)
    box_arr  = np.zeros((n_frames, 3, 3), dtype=np.float32)
    for i in range(n_frames):
        t, a, b, c = struct.unpack('<ffff', bf.read(16))
        coords = np.frombuffer(bf.read(n_atoms * 12), dtype='<f4').reshape(n_atoms, 3)
        time_arr[i]    = t
        box_arr[i,0,0] = a * ANG_TO_NM
        box_arr[i,1,1] = b * ANG_TO_NM
        box_arr[i,2,2] = c * ANG_TO_NM
        xyz[i]         = coords * ANG_TO_NM

step_arr = np.arange(n_frames, dtype=np.int32)
with XTCTrajectoryFile(xtc_path, mode) as f:
    f.write(xyz, time=time_arr, step=step_arr, box=box_arr)
"#;

/// Python script executed by `read_xtc`.  Reads an XTC file via MDTraj,
/// filters to the requested time window, and writes matching frames in the
/// wire-format binary that Rust then parses.
const READ_XTC_PY: &str = r#"
import sys, struct
import numpy as np
try:
    from mdtraj.formats import XTCTrajectoryFile
except ImportError:
    sys.exit("mdtraj not found. Install with: pip install mdtraj")

xtc_path, bin_path = sys.argv[1], sys.argv[2]
start = None if sys.argv[3] == 'None' else float(sys.argv[3])
end   = None if sys.argv[4] == 'None' else float(sys.argv[4])
NM_TO_ANG = 10.0

frames = []
with XTCTrajectoryFile(xtc_path, 'r') as f:
    while True:
        xyz, time_arr, _step, box_arr = f.read(n_frames=1)
        if len(xyz) == 0:
            break
        t = float(time_arr[0])
        if start is not None and t < start:
            continue
        if end is not None and t > end:
            break  # frames are chronological — nothing later can match
        a = float(box_arr[0,0,0]) * NM_TO_ANG if box_arr is not None else 0.0
        b = float(box_arr[0,1,1]) * NM_TO_ANG if box_arr is not None else 0.0
        c = float(box_arr[0,2,2]) * NM_TO_ANG if box_arr is not None else 0.0
        coords = (xyz[0] * NM_TO_ANG).astype(np.float32)
        frames.append((t, a, b, c, coords))

n_frames = len(frames)
n_atoms  = frames[0][4].shape[0] if n_frames > 0 else 0
with open(bin_path, 'wb') as bf:
    bf.write(struct.pack('<ii', n_frames, n_atoms))
    for t, a, b, c, coords in frames:
        bf.write(struct.pack('<ffff', t, a, b, c))
        bf.write(coords.flatten().astype('<f4').tobytes())
"#;

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
    /// Load all frames from a DCD file. Equivalent to `read_dcd(path, None, None)`.
    pub fn load(path: &Path) -> io::Result<Self> {
        Ok(Self {
            frames: read_dcd(path, None, None)?,
        })
    }

    /// Create or append snapshots to a DCD file. Equivalent to `write_dcd(path, &self.frames)`.
    /// This is a common trajectory/reporter format used by other software, including OpenMM and VMD.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        write_dcd(path, &self.frames)
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

/// Read frames from a DCD file, optionally filtered to a time window (in ps).
///
/// Frames whose time falls outside `[start_time, end_time]` are skipped by
/// seeking past their coordinate blocks, so only matching frames are decoded.
/// Pass `None` for either bound to leave that end open.
///
/// Returns `Ok(Vec::new())` for an empty file.
pub fn read_dcd(
    path: &Path,
    start_time: Option<f32>,
    end_time: Option<f32>,
) -> io::Result<Vec<DcdFrame>> {
    let f = File::open(path)?;
    let mut r = BufReader::new(f);

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
    let istart = icntrl[1] as f64;
    let nsavc = icntrl[2] as f64;
    let has_unitcell = icntrl[19] != 0 && icntrl[10] != 0;
    let delta = f32::from_le_bytes(hdr[4 + 36..4 + 40].try_into().unwrap()) as f64;

    skip_title_record(&mut r)?;

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
            x: 0.,
            y: 0.,
            z: 0.,
        },
        bounds_high: Vec3 {
            x: 0.,
            y: 0.,
            z: 0.,
        },
    };

    for i in 0..nset_total {
        let time = (istart + (i as f64) * nsavc) * delta;
        let in_range = start_time.map_or(true, |t| time >= t as f64)
            && end_time.map_or(true, |t| time <= t as f64);

        if !in_range {
            // Seek past unit cell record if present.
            if has_unitcell {
                let len = read_u32_le(&mut r)?;
                r.seek(SeekFrom::Current(len as i64 + 4))?;
            }
            // Seek past X, Y, Z records.
            for _ in 0..3 {
                let len = read_u32_le(&mut r)?;
                r.seek(SeekFrom::Current(len as i64 + 4))?;
            }
            continue;
        }

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

        frames.push(DcdFrame {
            time,
            atom_posits,
            unit_cell: unit_cell.clone(),
        });
    }

    Ok(frames)
}

/// Write frames to a DCD file, creating it if it does not exist or appending
/// to an existing file.
///
/// On first creation the full CORD header, title, and NATOM blocks are written.
/// On subsequent calls only the coordinate frames are appended and the NSET
/// counter in the header is updated.
///
/// DCD stores times implicitly via `istart`, `nsavc`, and `delta` in the header
/// (unlike TRR which embeds a timestamp per frame).  For newly created files
/// `istart = 0`, `nsavc = 1`, and `delta` is derived from the first two frames'
/// time difference so that frame times round-trip correctly.
pub fn write_dcd(path: &Path, frames: &[DcdFrame]) -> io::Result<()> {
    if frames.is_empty() {
        return Ok(());
    }

    let n_atoms = frames[0].atom_posits.len();
    for s in frames.iter().skip(1) {
        if s.atom_posits.len() != n_atoms {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "inconsistent atom counts across frames",
            ));
        }
    }

    let mut f = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(path)?;

    let file_len = f.metadata()?.len();

    if file_len == 0 {
        // New file: write CORD header, title, and NATOM.
        let nsets = frames.len() as i32;
        let delta: f32 = if frames.len() >= 2 {
            (frames[1].time - frames[0].time) as f32
        } else {
            0.0
        };

        let mut header = Vec::with_capacity(84);
        header.extend_from_slice(b"CORD");

        let mut icntrl = [0i32; 20];
        icntrl[0] = nsets;
        icntrl[1] = 0; // istart
        icntrl[2] = 1; // nsavc
        // icntrl[10] = 1: extra block present (unit cell); icntrl[19] = 1: has unit cell.
        icntrl[10] = 1;
        icntrl[19] = 1;
        for v in icntrl {
            header.extend_from_slice(&v.to_le_bytes());
        }
        header[4 + 36..4 + 40].copy_from_slice(&delta.to_le_bytes());
        write_record(&mut f, &header)?;

        let title = format!(
            "Created by bio_files  NATOMS={}  NFRAMES={}",
            n_atoms, nsets
        );
        let mut line = [0u8; 80];
        let tb = title.as_bytes();
        line[..tb.len().min(80)].copy_from_slice(&tb[..tb.len().min(80)]);
        let mut title_block = Vec::with_capacity(4 + 80);
        title_block.extend_from_slice(&1i32.to_le_bytes());
        title_block.extend_from_slice(&line);
        write_record(&mut f, &title_block)?;

        let mut natom_block = Vec::with_capacity(4);
        natom_block.extend_from_slice(&(n_atoms as i32).to_le_bytes());
        write_record(&mut f, &natom_block)?;
    } else {
        // Existing file: verify header and NATOM, seek to end, then update NSET.
        f.seek(SeekFrom::Start(0))?;
        let l1 = read_u32_le(&mut f)?;
        if !(84..=1024 * 1024).contains(&l1) {
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

        let mut icntrl = [0i32; 20];
        for (i, item) in icntrl.iter_mut().enumerate() {
            let off = 4 + i * 4;
            *item = i32::from_le_bytes(hdr[off..off + 4].try_into().unwrap());
        }
        let cur_nset = icntrl[0];

        if icntrl[19] == 0 || icntrl[10] == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "existing DCD does not have unit cell blocks enabled",
            ));
        }

        skip_title_record(&mut f)?;

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
                "atom count mismatch with existing DCD file",
            ));
        }

        f.seek(SeekFrom::End(0))?;

        let mut xs = vec![0f32; n_atoms];
        let mut ys = vec![0f32; n_atoms];
        let mut zs = vec![0f32; n_atoms];
        for frame in frames {
            for (i, p) in frame.atom_posits.iter().enumerate() {
                xs[i] = p.x;
                ys[i] = p.y;
                zs[i] = p.z;
            }
            let xb = unsafe { core::slice::from_raw_parts(xs.as_ptr() as *const u8, xs.len() * 4) };
            let yb = unsafe { core::slice::from_raw_parts(ys.as_ptr() as *const u8, ys.len() * 4) };
            let zb = unsafe { core::slice::from_raw_parts(zs.as_ptr() as *const u8, zs.len() * 4) };
            write_unit_cell_record(&mut f, &frame.unit_cell)?;
            write_record(&mut f, xb)?;
            write_record(&mut f, yb)?;
            write_record(&mut f, zb)?;
        }

        let new_nset = cur_nset
            .checked_add(frames.len() as i32)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "NSET overflow"))?;

        // NSET lives at byte offset 8: 4-byte leading length marker + 4 bytes "CORD".
        f.seek(SeekFrom::Start(8))?;
        f.write_all(&new_nset.to_le_bytes())?;
        f.flush()?;

        return Ok(());
    }

    // New file: write all frames.
    let mut xs = vec![0f32; n_atoms];
    let mut ys = vec![0f32; n_atoms];
    let mut zs = vec![0f32; n_atoms];
    for frame in frames {
        for (i, p) in frame.atom_posits.iter().enumerate() {
            xs[i] = p.x;
            ys[i] = p.y;
            zs[i] = p.z;
        }
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

/// Read frames from an XTC file, optionally filtered to a time window (ps).
///
/// Requires Python 3 and `mdtraj` (`pip install mdtraj`) to be available on
/// the system PATH.  The XTC is read frame-by-frame inside the Python helper
/// so only matching frames are decoded.
///
/// Returns `Ok(Vec::new())` when no frames fall within the requested window.
pub fn read_xtc(
    path: &Path,
    start_time: Option<f32>,
    end_time: Option<f32>,
) -> io::Result<Vec<DcdFrame>> {
    let pid = std::process::id();
    let tmp = env::temp_dir();
    let script_path = tmp.join(format!("bio_read_xtc_{pid}.py"));
    let bin_path = tmp.join(format!("bio_xtc_out_{pid}.bin"));

    fs::write(&script_path, READ_XTC_PY)?;

    let start_str = start_time
        .map(|t| t.to_string())
        .unwrap_or_else(|| "None".to_string());
    let end_str = end_time
        .map(|t| t.to_string())
        .unwrap_or_else(|| "None".to_string());

    let out = run_python(&[
        script_path.to_str().unwrap(),
        path.to_str().unwrap(),
        bin_path.to_str().unwrap(),
        &start_str,
        &end_str,
    ]);

    let _ = fs::remove_file(&script_path);

    let out = out?;
    if !out.status.success() {
        let _ = fs::remove_file(&bin_path);
        return Err(io::Error::other(format!(
            "read_xtc failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }

    let data = fs::read(&bin_path)?;
    let _ = fs::remove_file(&bin_path);

    parse_xtc_bin(&data)
}

/// Write frames to an XTC file, creating it if it does not exist or appending
/// to an existing one.
///
/// Requires Python 3 and `mdtraj` (`pip install mdtraj`) to be available on
/// the system PATH.  Each call spawns Python once; no persistent process is
/// needed.  MDTraj opens XTC in `'a'` (append) mode when the file already
/// exists, which is safe because the XTC format is a plain sequence of
/// self-contained frames.
pub fn write_xtc(path: &Path, frames: &[DcdFrame]) -> io::Result<()> {
    if frames.is_empty() {
        return Ok(());
    }

    let n_atoms = frames[0].atom_posits.len();
    for f in frames.iter().skip(1) {
        if f.atom_posits.len() != n_atoms {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "inconsistent atom counts across frames",
            ));
        }
    }

    let pid = std::process::id();
    let tmp = env::temp_dir();
    let script_path = tmp.join(format!("bio_write_xtc_{pid}.py"));
    let bin_path = tmp.join(format!("bio_xtc_in_{pid}.bin"));

    // Serialize frames into the wire-format binary.
    {
        let mut bf = File::create(&bin_path)?;
        bf.write_all(&(frames.len() as i32).to_le_bytes())?;
        bf.write_all(&(n_atoms as i32).to_le_bytes())?;
        for frame in frames {
            let uc = &frame.unit_cell;
            let a = uc.bounds_high.x - uc.bounds_low.x;
            let b = uc.bounds_high.y - uc.bounds_low.y;
            let c = uc.bounds_high.z - uc.bounds_low.z;
            bf.write_all(&(frame.time as f32).to_le_bytes())?;
            bf.write_all(&a.to_le_bytes())?;
            bf.write_all(&b.to_le_bytes())?;
            bf.write_all(&c.to_le_bytes())?;
            for p in &frame.atom_posits {
                bf.write_all(&p.x.to_le_bytes())?;
                bf.write_all(&p.y.to_le_bytes())?;
                bf.write_all(&p.z.to_le_bytes())?;
            }
        }
    }

    fs::write(&script_path, WRITE_XTC_PY)?;

    let mode = if path.exists() { "a" } else { "w" };
    let out = run_python(&[
        script_path.to_str().unwrap(),
        mode,
        path.to_str().unwrap(),
        bin_path.to_str().unwrap(),
    ]);

    let _ = fs::remove_file(&script_path);
    let _ = fs::remove_file(&bin_path);

    let out = out?;
    if !out.status.success() {
        return Err(io::Error::other(format!(
            "write_xtc failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }

    Ok(())
}

/// Try `python3` first; fall back to `python` (needed on some Windows installs).
fn run_python(args: &[&str]) -> io::Result<std::process::Output> {
    match Command::new("python3").args(args).output() {
        Err(e) if e.kind() == io::ErrorKind::NotFound => Command::new("python").args(args).output(),
        other => other,
    }
}

/// Parse the wire-format binary produced by the read-XTC Python helper.
fn parse_xtc_bin(data: &[u8]) -> io::Result<Vec<DcdFrame>> {
    if data.len() < 8 {
        return Ok(Vec::new());
    }

    let n_frames = i32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    let n_atoms = i32::from_le_bytes(data[4..8].try_into().unwrap()) as usize;

    let frame_bytes = 16 + n_atoms * 12; // time+box (16) + xyz (n_atoms*3*4)
    let expected = 8 + n_frames * frame_bytes;
    if data.len() < expected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "XTC binary too short: expected {} bytes, got {}",
                expected,
                data.len()
            ),
        ));
    }

    let mut frames = Vec::with_capacity(n_frames);
    let mut off = 8usize;

    for _ in 0..n_frames {
        let time = f32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as f64;
        off += 4;
        let a = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
        off += 4;
        let b = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
        off += 4;
        let c = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
        off += 4;

        let mut atom_posits = Vec::with_capacity(n_atoms);
        for _ in 0..n_atoms {
            let x = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
            off += 4;
            let y = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
            off += 4;
            let z = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
            off += 4;
            atom_posits.push(Vec3 { x, y, z });
        }

        frames.push(DcdFrame {
            time,
            atom_posits,
            unit_cell: DcdUnitCell {
                bounds_low: Vec3 {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                },
                bounds_high: Vec3 { x: a, y: b, z: c },
            },
        });
    }

    Ok(frames)
}

/// A wrapper for writing a DCD record: Payload sandwiched by length.
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
    for (i, v) in six.iter_mut().enumerate() {
        let j = i * 8;
        *v = f64::from_le_bytes(b[j..j + 8].try_into().unwrap());
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
