//! [Ab initio Molecular Dynamics](https://www.faccts.de/docs/orca/6.1/manual/contents/moleculardynamics/moldyn.html)

use std::{
    io,
    path::{Path, PathBuf},
};

use crate::{Xyz, load_xyz_trajectory, orca::make_inp_block};

/// [Thermostat](https://www.faccts.de/docs/orca/6.1/manual/contents/moleculardynamics/moldyn.html#thermostat)
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Thermostat {
    Berendensen,
    Csvr,
    Nhc,
    None,
}

impl Thermostat {
    pub fn keyword(self) -> String {
        match self {
            Self::Berendensen => "Berendensen",
            Self::Csvr => "CSVR",
            Self::Nhc => "NHC",
            Self::None => "None",
        }
        .to_string()
    }
}

/// [Ab initio Molecular Dynamics Command List](https://www.faccts.de/docs/orca/6.1/manual/contents/moleculardynamics/moldyn.html#command-list)
#[derive(Clone, Debug)]
pub struct Dynamics {
    /// [Timestep](https://www.faccts.de/docs/orca/6.1/manual/contents/moleculardynamics/moldyn.html#timestep) fs
    pub timestep: f32,
    /// [Initvel](https://www.faccts.de/docs/orca/6.1/manual/contents/moleculardynamics/moldyn.html#initvel) Kelvin
    pub init_vel: f32,
    /// [Thermostat](https://www.faccts.de/docs/orca/6.1/manual/contents/moleculardynamics/moldyn.html#thermostat)
    pub thermostat: Thermostat,
    /// Kelvin
    pub thermostat_temp: f32,
    /// fs. 10-100 is reasonable.
    pub thermostat_timecon: f32,
    /// [Dump](https://www.faccts.de/docs/orca/6.1/manual/contents/moleculardynamics/moldyn.html#dump)
    pub traj_out_dir: PathBuf,
    pub steps: u32,
    // todo: Epxand this
}

impl Dynamics {
    /// Note: We can ommit units from the strings in favor of default units, if we wish.
    pub fn make_inp(&self) -> String {
        let mut contents = vec![
            ("Timestep", format!("{:.1}_fs", self.timestep)),
            ("Initvel", format!("{:.1}_K", self.init_vel)),
            (
                "Thermostat",
                format!(
                    "{} {:.1}_K Timecon {:.1}_fs",
                    self.thermostat.keyword(),
                    self.thermostat_temp,
                    self.thermostat_timecon
                ),
            ),
            (
                "Dump",
                format!(
                    "Position Stride 1 Filename \"{}\"",
                    self.traj_out_dir.to_str().unwrap()
                ),
            ),
            ("Run", self.steps.to_string()),
        ];

        make_inp_block("md", &contents, &[])
    }
}

#[derive(Clone, Debug)]
pub struct DynamicsOutput {
    pub text: String,
    pub trajectory: Vec<Xyz>,
}

impl DynamicsOutput {
    pub fn new(traj_path: &Path, text: String) -> io::Result<Self> {
        let trajectory = load_xyz_trajectory(traj_path)?;
        Ok(Self { text, trajectory })
    }
}
