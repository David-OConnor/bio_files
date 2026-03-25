//! GROMACS MDP (Molecular Dynamics Parameters) file types and generation. This is responsible
//! for general MD configuration.
//!
//! See the [GROMACS manual — MDP options](https://manual.gromacs.org/current/user-guide/mdp-options.html).
//! [Example MDP file](https://manual.gromacs.org/2026.1/reference-manual/file-formats.html#mdp)
//!
//! Fields left as `None` are omitted from the generated file, deferring to GROMACS defaults.

use std::fmt;

/// Integration algorithm.
///
/// [GROMACS manual: run-control](https://manual.gromacs.org/current/user-guide/mdp-options.html#run-control)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Integrator {
    /// Leap-frog MD integrator (default in gromacsαθθs).
    #[default]
    Md,
    /// Velocity Verlet. Slower than `md`, but perhaps more rigorous.
    MdVv,
    MdVvAvek,
    /// Stochastic dynamics (Langevin). `tcoupl` must be `no`; friction is
    /// specified via `tau-t` (the inverse friction constant, ps) per group.
    /// Also known as Langevin, with built-in stochastic thermostat.
    Sd,
    Bd,
    /// Steep-descent energy minimization. (Good default for energy minimization)
    Steep,
    /// Conjugate gradient energy minimization.
    Cg,
    LBfgs,
    Nm,
    Tpi,
    Tpic,
    Mimic,
}

impl Integrator {
    pub fn keyword(self) -> &'static str {
        use Integrator::*;
        match self {
            Md => "md",
            MdVv => "md-vv",
            MdVvAvek => "md-vv-avek",
            Bd => "bd",
            Sd => "sd",
            Steep => "steep",
            Cg => "cg",
            LBfgs => "l-bfgs",
            Nm => "nm",
            Tpi => "tpi",
            Tpic => "tpic",
            Mimic => "mimic",
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
    /// Nosé–Hoover; canonical ensemble.
    NoseHoover,
    Andersen,
    AndersenMassive,
    /// Velocity-rescaling thermostat; produces the canonical ensemble.
    #[default]
    VRescale,
}

impl Thermostat {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::No => "no",
            Self::NoseHoover => "nose-hoover",
            Self::Andersen => "andersen",
            Self::AndersenMassive => "andersen-massive",
            Self::VRescale => "v-rescale",
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
    No,
    /// Not recommended, by GROMACS.
    Berendsen,
    #[default]
    /// Exponential relaxation pressure coupling with time constant tau-p, including a stochastic term to enforce correct volume fluctuations.
    CRescale,
    /// Parrinello–Rahman; NPT ensemble.
    ParrinelloRahman,
    /// Martyna-Tuckerman-Tobias-Klein implementation, only useable with integrator=md-vv or integrator=md-vv-avek, very similar to Parrinello-Rahman.
    Mtkk,
}

impl Barostat {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::No => "no",
            Self::CRescale => "C-rescale",
            Self::Berendsen => "Berendsen",
            Self::ParrinelloRahman => "Parrinello-Rahman",
            Self::Mtkk => "MTTK",
        }
    }
}

impl fmt::Display for Barostat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.keyword())
    }
}

/// Pressure-coupling algorithm.
///
/// [GROMACS manual: pressure-coupling](https://manual.gromacs.org/current/user-guide/mdp-options.html#pressure-coupling)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum PressureCouplingType {
    #[default]
    ///  todo: What is the gromacs default?
    Isotropic,
    Semiisotropic,
    Anisotropic,
    SurfaceTension,
}

impl PressureCouplingType {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Isotropic => "isotropic",
            Self::Semiisotropic => "semiisotropic",
            Self::Anisotropic => "anisotropic",
            Self::SurfaceTension => "surface-tension",
        }
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
    /// All directions
    #[default]
    Xyz,
    No,
    /// X and Y directions only
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

/// Complete set of GROMACS Molecular Dynamics Parameters (MDP) for classical MD.
///
/// All parameters correspond directly to
/// [GROMACS MDP options](https://manual.gromacs.org/current/user-guide/mdp-options.html).
/// Similar in concept to `dyanmics::MdConfig`.
///
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
    pub pcoupl: Barostat,
    pub pcoupltype: PressureCouplingType,
    /// Pressure coupling time constant (ps).
    pub tau_p: f32,
    /// Reference pressure (bar).
    /// todo: All pcoupltype values other than isotropic require multiple ref_p values. Wrapped enum on
    /// todo: PressureCouplingType
    pub ref_p: f32,
    /// Isothermal compressibility (bar⁻¹), e.g. `4.5e-5` for water.
    /// todo: Same note on ref_p about pressure coupling type applies here.
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
    /// Valid values are 3-12.
    pub pme_order: u8,
}

/// We use the GROMACS defaults here.
impl Default for MdpParams {
    fn default() -> Self {
        Self {
            integrator: Integrator::Md,
            nsteps: 0,
            dt: 0.001,
            nstxout: 0,
            nstvout: 0,
            nstxout_compressed: 0,
            nstenergy: 1_000,
            nstlog: 1_000,
            coulombtype: CoulombType::default(),
            rcoulomb: 1.0,
            vdwtype: VdwType::default(),
            rvdw: 1.0,
            fourierspacing: None,
            thermostat: Thermostat::default(),
            tau_t: vec![1.0], // todo: 0.1?
            ref_t: vec![300.],
            pcoupl: Barostat::default(),
            pcoupltype: PressureCouplingType::default(),
            tau_p: 5.0,
            ref_p: 1.0,
            compressibility: 4.5e-5,
            pbc: Pbc::default(),
            gen_vel: true,
            gen_temp: 300.,
            gen_seed: -1,
            constraints: Constraints::default(),
            constraint_algorithm: ConstraintAlgorithm::default(),
            pme_order: 4,
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
        // For the `sd` integrator, tcoupl must be `no`, but tc-grps / tau-t / ref-t
        // are still required: tau-t is the inverse friction constant (1/γ, ps).
        let needs_tc_params =
            self.thermostat != Thermostat::No || self.integrator == Integrator::Sd;
        if needs_tc_params {
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
            self.pcoupl.keyword()
        ));

        if self.pcoupl != Barostat::No {
            // todo: QC pcoupltype.
            s.push_str(&format!(
                "pcoupltype              = {}\n",
                self.pcoupltype.keyword()
            ));
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

        // pme-order is only relevant when PME is used for Coulomb or VdW.
        if self.coulombtype == CoulombType::Pme || self.vdwtype == VdwType::Pme {
            s.push_str(&format!("pme-order               = {}\n", self.pme_order));
        }

        s
    }
}
