//! This example demonstrates syntax for saving and loading MD trajectories in
//! [DCD](https://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.app.dcdfile.DCDFile.html) and XTC formats.

use std::path::Path;

use bio_files::{
    FrameSlice,
    dcd::DcdTrajectory,
    xtc::{read_xtc, write_xtc},
};

fn main() {
    let path_dcd = Path::new("traj.dcd");
    let path_xtc = Path::new("traj.xtc");

    let traj = DcdTrajectory::load(path_dcd).unwrap();

    // XTC files load as the same frame type. Use `FrameSlice` to load a subset.
    let all = FrameSlice::Time {
        start: None,
        end: None,
    };
    let _frames_xtc = read_xtc(path_xtc, all).unwrap();

    for frame in &traj.frames {
        println!(
            "Time: {:.1} ps unit cell: {:?}",
            frame.time, frame.unit_cell
        );
        println!("Positions:");
        for _posit in &frame.atom_posits {
            // ...
        }
    }

    // To save:
    traj.save(path_dcd).unwrap();
    write_xtc(path_xtc, &traj.frames).unwrap();

    // If using the `dynamics` lib:
    // let snapshots = dynamics::load_snapshots_from_file(path_dcd).unwrap();
    // Or
    // let snaps = dynamics::Snapshot::from_dcd(&traj);

    // To convert a Snapshot:
    // snap.to_dcd(cell, false), // Unit cell, if we should write water atoms.

    // You will likely use a similar adapter for your application or library's native snapshot/frame
    // format.
}
