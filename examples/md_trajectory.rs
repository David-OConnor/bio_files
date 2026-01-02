//! This example demonstrates syntax for saving and loading MD trajectories in
//! [DCD](https://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.app.dcdfile.DCDFile.html) and XTC formats.

use std::path::Path;

use bio_files::dcd::DcdTrajectory;

fn main() {
    let path_dcd = Path::new("traj.dcd");
    let path_xtc = Path::new("traj.xtc");

    let traj = DcdTrajectory::load(path_dcd).unwrap();
    // Or, if you have MDTraj installed, load XTC files:
    let traj = DcdTrajectory::load_xtc(path_xtc).unwrap();

    for frame in &traj.frames {
        println!("Time: {} unit cell: {:?}", frame.time, frame.unit_cell);
        println!("Positions:");
        for posit in &frame.atom_posits {
            // ...
        }
    }

    // To save:
    traj.save(path_dcd).unwrap();
    traj.save_xtc(path_xtc).unwrap();
}
