//! Shared MD output folder naming: `md_out/{prefix}_{YYYY-MM-DD_HH-MM-SS}`.

use std::{
    fs,
    io::{self, ErrorKind},
    path::{Path, PathBuf},
    thread,
    time::Duration,
};

pub const MD_OUT_DIR: &str = "md_out";

/// Create a new run folder, e.g. `md_out/gromacs_2026-11-12_22-47-22`. If that name
/// exists, waits for the next second.
pub fn create_run_dir(prefix: &str) -> io::Result<PathBuf> {
    create_run_dir_in(Path::new(MD_OUT_DIR), prefix)
}

pub fn create_run_dir_in(root: &Path, prefix: &str) -> io::Result<PathBuf> {
    fs::create_dir_all(root)?;

    loop {
        let stamp = chrono::Local::now().format("%Y-%m-%d_%H-%M-%S");
        let dir = root.join(format!("{prefix}_{stamp}"));

        match fs::create_dir(&dir) {
            Ok(()) => return Ok(dir),
            Err(e) if e.kind() == ErrorKind::AlreadyExists => {
                thread::sleep(Duration::from_millis(100))
            }
            Err(e) => return Err(e),
        }
    }
}
