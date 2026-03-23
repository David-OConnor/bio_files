//! GROMACS MDP (Molecular Dynamics Parameters) file types and generation.
//!
//! See the [GROMACS manual — MDP options](https://manual.gromacs.org/current/user-guide/mdp-options.html).
//!
//! Fields left as `None` are omitted from the generated file, deferring to GROMACS defaults —
//! the same philosophy used by the ORCA module's `Option` fields.

use std::fmt;

/// Integration algorithm.
///
/// [GROMACS manual: run-control](https://manual.gromacs.org/current/user-guide/mdp-options.html#run-control)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Integrator {
    /// Leap-frog MD integrator (default).
    #[default]
    Md,
    /// Velocity Verlet.
    MdVv,
    /// Steep-descent energy minimization.
    Steep,
    /// Conjugate gradient energy minimization.
    Cg,
}

impl Integrator {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Md => "md",
            Self::MdVv => "md-vv",
            Self::Steep => "steep",
            Self::Cg => "cg",
        }
    }
}

impl fmt::Display for Integrator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.keyword())
    }
}

/// Temperature-coupling algorithm.
///
/// [GROMACS manual: temperature-coupling](https://manual.gromacs.org/current/user-guide/mdp-options.html#temperature-coupling)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Thermostat {
    No,
    Berendsen,
    /// Velocity-rescaling thermostat; produces the canonical ensemble.
    #[default]
    VRescale,
    /// Nosé–Hoover; canonical ensemble.
    NoseHoover,
}

impl Thermostat {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::No => "no",
            Self::Berendsen => "berendsen",
            Self::VRescale => "v-rescale",
            Self::NoseHoover => "nose-hoover",
        }
    }
}

impl fmt::Display for Thermostat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.keyword())
    }
}

/// Pressure-coupling algorithm.
///
/// [GROMACS manual: pressure-coupling](https://manual.gromacs.org/current/user-guide/mdp-options.html#pressure-coupling)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Barostat {
    #[default]
    No,
    Berendsen,
    /// Parrinello–Rahman; NPT ensemble.
    ParrinelloRahman,
}

impl Barostat {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::No => "no",
            Self::Berendsen => "Berendsen",
            Self::ParrinelloRahman => "Parrinello-Rahman",
        }
    }
}

impl fmt::Display for Barostat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.keyword())
    }
}

/// Coulomb interaction treatment.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum CoulombType {
    /// Particle Mesh Ewald — default for periodic systems.
    #[default]
    Pme,
    CutOff,
    ReactionField,
}

impl CoulombType {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Pme => "PME",
            Self::CutOff => "Cut-off",
            Self::ReactionField => "Reaction-Field",
        }
    }
}

impl fmt::Display for CoulombType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.keyword())
    }
}

/// Van der Waals interaction treatment.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum VdwType {
    #[default]
    CutOff,
    Pme,
}

impl VdwType {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::CutOff => "Cut-off",
            Self::Pme => "PME",
        }
    }
}

impl fmt::Display for VdwType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.keyword())
    }
}

/// Periodic boundary conditions.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Pbc {
    #[default]
    Xyz,
    No,
    Xy,
}

impl Pbc {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Xyz => "xyz",
            Self::No => "no",
            Self::Xy => "xy",
        }
    }
}

impl fmt::Display for Pbc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.keyword())
    }
}

/// Bond constraint algorithm.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum ConstraintAlgorithm {
    /// LINCS — recommended for MD.
    #[default]
    Lincs,
    /// SHAKE — slower but supports angle constraints.
    Shake,
}

impl ConstraintAlgorithm {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Lincs => "LINCS",
            Self::Shake => "SHAKE",
        }
    }
}

/// Which bonds to constrain.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Constraints {
    None,
    /// Constrain bonds involving H atoms; allows longer timesteps.
    #[default]
    HBonds,
    AllBonds,
}

impl Constraints {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::HBonds => "h-bonds",
            Self::AllBonds => "all-bonds",
        }
    }
}

/// Complete set of GROMACS MDP parameters for classical MD.
///
/// All parameters correspond directly to
/// [GROMACS MDP options](https://manual.gromacs.org/current/user-guide/mdp-options.html).
/// Units follow GROMACS conventions (ps, nm, kJ/mol, K, bar).
#[derive(Clone, Debug)]
pub struct MdpParams {
    // --- Run control ---
    pub integrator: Integrator,
    /// Total number of integration steps.
    pub nsteps: u64,
    /// Integration timestep in **ps** (e.g. `0.002` = 2 fs).
    pub dt: f32,

    // --- Output control ---
    /// Write full-precision coordinates to `.trr` every N steps (0 = never).
    pub nstxout: u32,
    /// Write velocities to `.trr` every N steps (0 = never).
    pub nstvout: u32,
    /// Write compressed coordinates to `.xtc` every N steps (0 = never).
    pub nstxout_compressed: u32,
    /// Write energies to `.edr` every N steps.
    pub nstenergy: u32,
    /// Write to `.log` every N steps.
    pub nstlog: u32,

    // --- Non-bonded interactions ---
    pub coulombtype: CoulombType,
    /// Coulomb cutoff in **nm**.
    pub rcoulomb: f32,
    pub vdwtype: VdwType,
    /// VdW cutoff in **nm**.
    pub rvdw: f32,
    /// PME Fourier grid spacing in **nm**. `None` → GROMACS default (0.12).
    pub fourierspacing: Option<f32>,

    // --- Temperature coupling ---
    pub thermostat: Thermostat,
    /// Time constant (ps) for each coupling group. Must match length of `ref_t`.
    pub tau_t: Vec<f32>,
    /// Reference temperature (K) for each coupling group.
    pub ref_t: Vec<f32>,

    // --- Pressure coupling ---
    pub barostat: Barostat,
    /// Pressure coupling time constant (ps).
    pub tau_p: f32,
    /// Reference pressure (bar).
    pub ref_p: f32,
    /// Isothermal compressibility (bar⁻¹), e.g. `4.5e-5` for water.
    pub compressibility: f32,

    // --- Periodic boundary ---
    pub pbc: Pbc,

    // --- Velocity generation ---
    /// Generate initial velocities from a Maxwell–Boltzmann distribution.
    pub gen_vel: bool,
    /// Temperature (K) for velocity generation.
    pub gen_temp: f32,
    /// Random seed for velocity generation; `-1` lets GROMACS choose.
    pub gen_seed: i32,

    // --- Constraints ---
    pub constraints: Constraints,
    pub constraint_algorithm: ConstraintAlgorithm,
}

impl Default for MdpParams {
    fn default() -> Self {
        Self {
            integrator: Integrator::Md,
            nsteps: 50_000,
            dt: 0.002,
            nstxout: 0,
            nstvout: 0,
            nstxout_compressed: 500,
            nstenergy: 500,
            nstlog: 500,
            coulombtype: CoulombType::Pme,
            rcoulomb: 1.0,
            vdwtype: VdwType::CutOff,
            rvdw: 1.0,
            fourierspacing: None,
            thermostat: Thermostat::VRescale,
            tau_t: vec![0.1],
            ref_t: vec![310.],
            barostat: Barostat::No,
            tau_p: 2.0,
            ref_p: 1.0,
            compressibility: 4.5e-5,
            pbc: Pbc::Xyz,
            gen_vel: true,
            gen_temp: 310.,
            gen_seed: -1,
            constraints: Constraints::HBonds,
            constraint_algorithm: ConstraintAlgorithm::Lincs,
        }
    }
}

impl MdpParams {
    /// Render as a GROMACS `.mdp` file string.
    pub fn make_mdp(&self) -> String {
        let mut s = String::from("; GROMACS MDP file — generated by Molchanica\n\n");

        // Run control
        s.push_str("; Run control\n");
        s.push_str(&format!(
            "integrator              = {}\n",
            self.integrator.keyword()
        ));
        s.push_str(&format!("nsteps                  = {}\n", self.nsteps));
        s.push_str(&format!("dt                      = {}\n", self.dt));

        // Output
        s.push_str("\n; Output control\n");
        s.push_str(&format!("nstxout                 = {}\n", self.nstxout));
        s.push_str(&format!("nstvout                 = {}\n", self.nstvout));
        s.push_str(&format!(
            "nstxout-compressed      = {}\n",
            self.nstxout_compressed
        ));
        s.push_str(&format!("nstenergy               = {}\n", self.nstenergy));
        s.push_str(&format!("nstlog                  = {}\n", self.nstlog));

        // Non-bonded
        s.push_str("\n; Non-bonded interactions\n");
        s.push_str(&format!(
            "coulombtype             = {}\n",
            self.coulombtype.keyword()
        ));
        s.push_str(&format!("rcoulomb                = {}\n", self.rcoulomb));
        s.push_str(&format!(
            "vdw-type                = {}\n",
            self.vdwtype.keyword()
        ));
        s.push_str(&format!("rvdw                    = {}\n", self.rvdw));
        if let Some(fs) = self.fourierspacing {
            s.push_str(&format!("fourierspacing          = {fs}\n"));
        }

        // Temperature coupling
        s.push_str("\n; Temperature coupling\n");
        s.push_str(&format!(
            "tcoupl                  = {}\n",
            self.thermostat.keyword()
        ));
        if self.thermostat != Thermostat::No {
            let n = self.ref_t.len().max(1);
            let tc_grps = (0..n).map(|_| "System").collect::<Vec<_>>().join(" ");
            let tau_t = self
                .tau_t
                .iter()
                .map(|v| format!("{v}"))
                .collect::<Vec<_>>()
                .join(" ");
            let ref_t = self
                .ref_t
                .iter()
                .map(|v| format!("{v}"))
                .collect::<Vec<_>>()
                .join(" ");
            s.push_str(&format!("tc-grps                 = {tc_grps}\n"));
            s.push_str(&format!("tau-t                   = {tau_t}\n"));
            s.push_str(&format!("ref-t                   = {ref_t}\n"));
        }

        // Pressure coupling
        s.push_str("\n; Pressure coupling\n");
        s.push_str(&format!(
            "pcoupl                  = {}\n",
            self.barostat.keyword()
        ));
        if self.barostat != Barostat::No {
            s.push_str(&format!("tau-p                   = {}\n", self.tau_p));
            s.push_str(&format!("ref-p                   = {}\n", self.ref_p));
            s.push_str(&format!(
                "compressibility         = {}\n",
                self.compressibility
            ));
        }

        // PBC
        s.push_str("\n; Periodic boundary conditions\n");
        s.push_str(&format!(
            "pbc                     = {}\n",
            self.pbc.keyword()
        ));

        // Velocity generation
        s.push_str("\n; Velocity generation\n");
        s.push_str(&format!(
            "gen-vel                 = {}\n",
            if self.gen_vel { "yes" } else { "no" }
        ));
        if self.gen_vel {
            s.push_str(&format!("gen-temp                = {}\n", self.gen_temp));
            s.push_str(&format!("gen-seed                = {}\n", self.gen_seed));
        }

        // Constraints
        s.push_str("\n; Constraints\n");
        s.push_str(&format!(
            "constraints             = {}\n",
            self.constraints.keyword()
        ));
        if self.constraints != Constraints::None {
            s.push_str(&format!(
                "constraint-algorithm    = {}\n",
                self.constraint_algorithm.keyword()
            ));
        }

        s
    }
}
