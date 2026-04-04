use std::{
    fs,
    fs::OpenOptions,
    io::{self, Write},
    path::Path,
    process::Command,
};

use crate::{
    FrameSlice,
    dcd::{DcdFrame, DcdMetadata, DcdTrajectory, read_dcd},
};

// ---------------------------------------------------------------------------
// Internal helper
// ---------------------------------------------------------------------------

/// Run `mdconvert input -o output`, returning an error if the process fails or
/// is not found.  mdtraj installs `mdconvert` as a console-script alongside
/// the Python interpreter, so it is always on PATH when mdtraj is installed —
/// regardless of which Python alias (python / python3 / py) the OS exposes.
fn run_mdconvert(input: &Path, output: &Path) -> io::Result<()> {
    let out = Command::new("mdconvert")
        .arg(input)
        .arg("-o")
        .arg(output)
        .arg("--force")
        .output()
        .map_err(|e| {
            if e.kind() == io::ErrorKind::NotFound {
                io::Error::other("mdconvert not found. Install mdtraj: pip install mdtraj")
            } else {
                e
            }
        })?;

    if !out.status.success() {
        return Err(io::Error::other(format!(
            "mdconvert failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Lightweight metadata about an XTC file.
pub struct XtcMetadata {
    pub num_atoms: usize,
    pub num_frames: usize,
    /// Time of the first frame, in ps.
    pub start_time: f32,
    /// Time of the last frame, in ps.
    pub end_time: f32,
    /// Time step between consecutive saved frames, in ps.  0 if < 2 frames.
    pub dt: f32,
}

impl XtcMetadata {
    /// Read metadata from an XTC file by converting it to a temporary DCD and
    /// reading that file's header.  Requires `mdconvert` (installed with mdtraj).
    pub fn read(path: &Path) -> io::Result<Self> {
        let pid = std::process::id();
        let tmp = std::env::temp_dir();
        let dcd_path = tmp.join(format!("bio_xtc_meta_{pid}.dcd"));

        run_mdconvert(path, &dcd_path)?;

        let md = DcdMetadata::read(&dcd_path);
        let _ = fs::remove_file(&dcd_path);
        let md = md?;

        // DcdMetadata: start_step = istart (step number), dt = delta (ps/step),
        // save_interval_steps = nsavc.  Actual start/end time in ps:
        //   t_start = istart * dt
        //   t_end   = (istart + (n_frames-1) * nsavc) * dt  ← already in md.end_time
        let start_time = md.start_step * md.dt;
        let dt = if md.num_frames > 1 {
            md.save_interval_steps as f32 * md.dt
        } else {
            0.0
        };

        Ok(Self {
            num_atoms: md.num_atoms,
            num_frames: md.num_frames,
            start_time,
            end_time: md.end_time,
            dt,
        })
    }
}

// ---------------------------------------------------------------------------
// Read / write
// ---------------------------------------------------------------------------

/// Read frames from an XTC file, filtered by a [`FrameSlice`].
///
/// Converts the XTC to a temporary DCD via `mdconvert`, reads the DCD with the
/// requested slice, then removes the temporary file.
pub fn read_xtc(path: &Path, slice: FrameSlice) -> io::Result<Vec<DcdFrame>> {
    let pid = std::process::id();
    let tmp = std::env::temp_dir();
    let dcd_path = tmp.join(format!("bio_xtc_read_{pid}.dcd"));

    run_mdconvert(path, &dcd_path)?;

    let frames = read_dcd(&dcd_path, slice);
    let _ = fs::remove_file(&dcd_path);
    frames
}

/// Write (append) frames to an XTC file.
///
/// Converts the frames to a temporary DCD (pure Rust, no subprocess), then
/// calls `mdconvert` to produce a temporary XTC chunk, and binary-appends
/// that chunk to `path`.  The XTC format is a plain sequence of self-contained
/// frame records, so appending is always safe.
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
    let tmp = std::env::temp_dir();
    let dcd_path = tmp.join(format!("bio_xtc_write_{pid}.dcd"));
    let xtc_tmp = tmp.join(format!("bio_xtc_write_{pid}.xtc"));

    // Write frames to a temporary DCD (pure-Rust, no subprocess).
    let traj = DcdTrajectory {
        frames: frames.to_vec(),
    };
    traj.save(&dcd_path)?;

    // Convert DCD → XTC chunk.
    let result = run_mdconvert(&dcd_path, &xtc_tmp);
    let _ = fs::remove_file(&dcd_path);
    result?;

    // Binary-append the new chunk to the output XTC.
    let chunk = fs::read(&xtc_tmp);
    let _ = fs::remove_file(&xtc_tmp);

    let mut out = OpenOptions::new().create(true).append(true).open(path)?;
    out.write_all(&chunk?)?;

    Ok(())
}
