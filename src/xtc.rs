use std::{env, fs, fs::File, io, io::Write, path::Path, process::Command};

use lin_alg::f32::Vec3;

use crate::dcd::{DcdFrame, DcdUnitCell};

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
