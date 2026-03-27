//! GROMACS MDP (Molecular Dynamics Parameters) file types and generation. This is responsible
//! for general MD configuration.
//!
//! See the [GROMACS manual — MDP options](https://manual.gromacs.org/current/user-guide/mdp-options.html).
//! [Example MDP file](https://manual.gromacs.org/2026.1/reference-manual/file-formats.html#mdp)
//!
//! Fields left as `None` are omitted from the generated file, deferring to GROMACS defaults.

use std::{
    collections::HashMap,
    fmt, io,
    io::{Error, ErrorKind},
    path::Path,
};

use crate::gromacs::save_txt_to_file;

/// Abramowitz & Stegun 7.1.26 — max error 1.5×10⁻⁷.
fn erfc_approx(x: f32) -> f32 {
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let p = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    p * (-x * x).exp()
}

/// Invert erfc by bisection. Only called at parse time so performance is irrelevant.
fn erfc_inv_approx(y: f32) -> f32 {
    let (mut lo, mut hi) = (0.0_f32, 6.0_f32);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if erfc_approx(mid) > y {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Helper function to append a key and value to an .mdp file.
fn append_inp<T: ToString>(inp: &mut String, k: &str, v: T) {
    inp.push_str(&format!("{k:<25}= {}\n", v.to_string()));
}

/// Float wrapper that formats with scientific notation when `|v| < 0.001`.
/// This matches GROMACS conventions: `compressibility = 4.5e-5`, `ewald-rtol = 1e-5`, etc.
struct Gf(f32);

impl fmt::Display for Gf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let v = self.0;
        if v != 0.0 && v.abs() < 0.001 {
            write!(f, "{:e}", v)
        } else {
            write!(f, "{v}")
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PmeConfig {
    /// PME Fourier grid spacing in **nm**.
    pub fourierspacing: f32,
    /// Valid values are 3-12.
    pub order: u8,
    /// Ewald splitting parameter α in **nm⁻¹**.  Controls the balance between
    /// real-space and reciprocal-space work.  When writing an MDP file, this is
    /// converted to `ewald-rtol = erfc(α × rcoulomb)`.
    ///
    /// Default (3.12 nm⁻¹) gives ewald-rtol ≈ 1×10⁻⁵ at the GROMACS default
    /// rcoulomb of 1.0 nm, matching the GROMACS default behaviour.
    pub alpha: f32,
    pub rtol_lj: f32,
    pub epsilon_surface: Option<f32>,
}

impl Default for PmeConfig {
    fn default() -> Self {
        Self {
            fourierspacing: 0.12,
            order: 4,
            alpha: 3.12, // erfc(3.12 × 1.0 nm) ≈ 1e-5 — matches GROMACS default ewald-rtol
            rtol_lj: 1e-3,
            epsilon_surface: None,
        }
    }
}

impl PmeConfig {
    /// Render PME parameters to MDP key-value lines.
    /// `rcoulomb` (nm) is required to derive `ewald-rtol = erfc(alpha × rcoulomb)`.
    pub fn make_inp(&self, rcoulomb: f32) -> String {
        let mut res = String::new();

        append_inp(&mut res, "fourierspacing", Gf(self.fourierspacing));
        append_inp(&mut res, "pme-order", self.order);
        append_inp(
            &mut res,
            "ewald-rtol",
            Gf(erfc_approx(self.alpha * rcoulomb)),
        );
        append_inp(&mut res, "ewald-rtol-lj", Gf(self.rtol_lj));

        if let Some(v) = &self.epsilon_surface {
            append_inp(&mut res, "epsilon-surface", Gf(*v));
        }

        res
    }
}

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

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum FreeEnergyCalculations {
    #[default]
    No,
    Yes,
    Expanded,
}

impl FreeEnergyCalculations {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::No => "no",
            Self::Yes => "yes",
            Self::Expanded => "expanded",
        }
    }
}

/// [https://manual.gromacs.org/current/user-guide/mdp-options.html#energy-minimization](MDP Guide: Energy min)
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "encode", derive(bincode::Encode, bincode::Decode))]
pub struct EnergyMinimization {
    /// (10.0) [kJ mol-1 nm-1] the minimization is converged when the maximum force is smaller
    /// than this value
    pub emtol: f32,
    /// (0.01) [nm] initial step-size
    pub emstep: f32,
    /// (1000) [steps] interval of performing 1 steepest descent step while doing conjugate
    /// gradient energy minimization.
    pub nstcgsteep: u16,
    /// (10) Number of correction steps to use for L-BFGS minimization. A higher number
    /// is (at least theoretically) more accurate, but slower.
    pub nbfgscorr: u8,
}

// GROMACS defaults.
impl Default for EnergyMinimization {
    fn default() -> Self {
        Self {
            emtol: 10.,
            emstep: 0.01,
            nstcgsteep: 1_000,
            nbfgscorr: 0,
        }
    }
}

impl EnergyMinimization {
    pub fn make_inp(&self) -> String {
        let mut res = String::new();

        append_inp(&mut res, "emtol", Gf(self.emtol));
        append_inp(&mut res, "emstep", Gf(self.emstep));
        append_inp(&mut res, "nstcgsteep", self.nstcgsteep);
        append_inp(&mut res, "nbfgscorr", self.nbfgscorr);

        res
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
#[derive(Clone, Debug, PartialEq)]
pub enum Barostat {
    No,
    /// Not recommended, by GROMACS.
    Berendsen(BarostatCfg),
    /// Exponential relaxation pressure coupling with time constant tau-p, including a stochastic term to enforce correct volume fluctuations.
    CRescale(BarostatCfg),
    /// Parrinello–Rahman; NPT ensemble.
    ParrinelloRahman(BarostatCfg),
    /// Martyna-Tuckerman-Tobias-Klein implementation, only useable with integrator=md-vv or integrator=md-vv-avek, very similar to Parrinello-Rahman.
    Mtkk(BarostatCfg),
}

impl Default for Barostat {
    fn default() -> Self {
        Self::CRescale(BarostatCfg::default())
    }
}

impl Barostat {
    pub fn make_inp(&self) -> String {
        let mut res = String::new();

        match self {
            Self::No => {
                append_inp(&mut res, "pcoupl", "no");
            }
            Self::Berendsen(cfg) => {
                append_inp(&mut res, "pcoupl", "Berendsen");
                res.push_str(&cfg.make_inp());
            }
            Self::CRescale(cfg) => {
                append_inp(&mut res, "pcoupl", "C-rescale");
                res.push_str(&cfg.make_inp());
            }
            Self::ParrinelloRahman(cfg) => {
                append_inp(&mut res, "pcoupl", "Parrinello-Rahman");
                res.push_str(&cfg.make_inp());
            }
            Self::Mtkk(cfg) => {
                append_inp(&mut res, "pcoupl", "MTTK");
                res.push_str(&cfg.make_inp());
            }
        }

        res
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BarostatCfg {
    pub pcoupltype: PressureCouplingType,
    /// Pressure coupling time constant (ps).
    pub tau_p: f32,
}

impl BarostatCfg {
    pub fn make_inp(&self) -> String {
        let mut res = String::new();

        res.push_str(&self.pcoupltype.make_inp());
        append_inp(&mut res, "tau-p", self.tau_p);

        res
    }
}

impl Default for BarostatCfg {
    fn default() -> Self {
        Self {
            pcoupltype: Default::default(),
            tau_p: 5.,
        }
    }
}

/// Pressure-coupling algorithm.
///
/// [GROMACS manual: pressure-coupling](https://manual.gromacs.org/current/user-guide/mdp-options.html#pressure-coupling)
#[derive(Clone, Debug, PartialEq)]
pub enum PressureCouplingType {
    Isotropic {
        /// Isothermal compressibility (bar⁻¹), e.g. `4.5e-5` for water.
        compressibility: f32,
        /// Reference pressure (bar).
        ref_p: f32,
    },
    Semiisotropic {
        compressibility_xy: f32,
        compressibility_z: f32,
        ref_p_xy: f32,
        ref_p_z: f32,
    },
    Anisotropic {
        compressibility: [f32; 6],
        ref_p: [f32; 6],
    },
    SurfaceTension {
        compressibility_xy: f32,
        compressibility_z: f32,
        ref_p_sfc_tension: f32,
        ref_p_z: f32,
    },
}

impl Default for PressureCouplingType {
    fn default() -> Self {
        Self::Isotropic {
            ref_p: 1.0,
            compressibility: 4.5e-5,
        }
    }
}

impl PressureCouplingType {
    pub fn make_inp(&self) -> String {
        let mut res = String::new();

        match self {
            Self::Isotropic {
                compressibility,
                ref_p,
            } => {
                append_inp(&mut res, "pcoupltype", "isotropic");
                append_inp(&mut res, "ref-p", ref_p);
                append_inp(&mut res, "compressibility", Gf(*compressibility));
            }
            Self::Semiisotropic {
                compressibility_xy,
                compressibility_z,
                ref_p_xy,
                ref_p_z,
            } => {
                append_inp(&mut res, "pcoupltype", "semiisotropic");
                append_inp(
                    &mut res,
                    "ref-p",
                    format!("{} {}", Gf(*ref_p_xy), Gf(*ref_p_z)),
                );
                append_inp(
                    &mut res,
                    "compressibility",
                    format!("{} {}", Gf(*compressibility_xy), Gf(*compressibility_z)),
                );
            }
            Self::Anisotropic {
                compressibility,
                ref_p,
            } => {
                append_inp(&mut res, "pcoupltype", "anisotropic");
                append_inp(
                    &mut res,
                    "ref-p",
                    ref_p
                        .iter()
                        .map(|v| Gf(*v).to_string())
                        .collect::<Vec<_>>()
                        .join(" "),
                );
                append_inp(
                    &mut res,
                    "compressibility",
                    compressibility
                        .iter()
                        .map(|v| Gf(*v).to_string())
                        .collect::<Vec<_>>()
                        .join(" "),
                );
            }
            Self::SurfaceTension {
                compressibility_xy,
                compressibility_z,
                ref_p_sfc_tension,
                ref_p_z,
            } => {
                append_inp(&mut res, "pcoupltype", "surface-tension");
                append_inp(
                    &mut res,
                    "ref-p",
                    format!("{} {}", Gf(*ref_p_sfc_tension), Gf(*ref_p_z)),
                );
                append_inp(
                    &mut res,
                    "compressibility",
                    format!("{} {}", Gf(*compressibility_xy), Gf(*compressibility_z)),
                );
            }
        }

        res
    }
}

/// Coulomb interaction treatment.
#[derive(Clone, Debug, PartialEq)]
pub enum CoulombType {
    /// Particle Mesh Ewald — default for periodic systems.
    Pme(PmeConfig),
    CutOff,
    ReactionField,
}

impl Default for CoulombType {
    fn default() -> Self {
        Self::Pme(PmeConfig::default())
    }
}

impl CoulombType {
    pub fn make_inp(&self, rcoulomb: f32) -> String {
        let mut res = String::new();

        match self {
            Self::Pme(pme) => {
                append_inp(&mut res, "coulombtype", "PME");
                res.push_str(&pme.make_inp(rcoulomb));
            }
            Self::CutOff => append_inp(&mut res, "coulombtype", "Cut-off"),
            Self::ReactionField => append_inp(&mut res, "coulombtype", "Reaction-Field"),
        }

        res
    }

    pub fn keyword(&self) -> &'static str {
        match self {
            Self::Pme(_) => "PME",
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
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ConstraintAlgorithm {
    /// LINCS — recommended for MD.
    Lincs { order: u8, iter: u8 },
    /// SHAKE — slower but supports angle constraints.
    Shake { tol: f32 },
}

impl Default for ConstraintAlgorithm {
    fn default() -> Self {
        Self::Lincs { order: 4, iter: 1 }
        //Note: Shake tol default is 0.0001
    }
}

impl ConstraintAlgorithm {
    pub fn make_inp(self) -> String {
        let mut res = String::new();

        match self {
            Self::Lincs { order, iter } => {
                append_inp(&mut res, "constraint-algorithm", "LINCS");
                append_inp(&mut res, "lincs-order", order);
                append_inp(&mut res, "lincs-iter", iter);
            }
            Self::Shake { tol } => {
                append_inp(&mut res, "constraint-algorithm", "SHAKE");
                append_inp(&mut res, "shake-tol", Gf(tol));
            }
        }

        res
    }
}

/// Which bonds to constrain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraints {
    None,
    /// Constrain bonds involving H atoms; allows longer timesteps.
    HBonds(ConstraintAlgorithm),
    AllBonds(ConstraintAlgorithm),
    HAngles(ConstraintAlgorithm),
    AllAngles(ConstraintAlgorithm),
}

impl Default for Constraints {
    fn default() -> Self {
        Self::HAngles(ConstraintAlgorithm::default())
    }
}

impl Constraints {
    pub fn make_inp(&self) -> String {
        let mut res = String::new();

        match self {
            Self::None => {
                append_inp(&mut res, "constraints", "none");
            }
            Self::HBonds(ca) => {
                append_inp(&mut res, "constraints", "h-bonds");
                res.push_str(&ca.make_inp());
            }
            Self::AllBonds(ca) => {
                append_inp(&mut res, "constraints", "all-bonds");
                res.push_str(&ca.make_inp());
            }
            Self::HAngles(ca) => {
                append_inp(&mut res, "constraints", "h-angles");
                res.push_str(&ca.make_inp());
            }
            Self::AllAngles(ca) => {
                append_inp(&mut res, "constraints", "all-angles");
                res.push_str(&ca.make_inp());
            }
        }

        res
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
    /// None means use GROMACS defaults; don't write to the file.
    pub energy_minimization: Option<EnergyMinimization>,

    // --- Temperature coupling ---
    pub thermostat: Thermostat,
    /// Time constant (ps) for each coupling group. Must match length of `ref_t`.
    pub tau_t: Vec<f32>,
    /// Reference temperature (K) for each coupling group.
    pub ref_t: Vec<f32>,
    pub pcoupl: Barostat,

    // --- Periodic boundary ---
    pub pbc: Pbc,

    // --- Velocity generation ---
    /// Generate initial velocities from a Maxwell–Boltzmann distribution.
    pub gen_vel: bool,
    /// Temperature (K) for velocity generation.
    pub gen_temp: f32,
    /// Random seed for velocity generation; `-1` lets GROMACS choose.
    pub gen_seed: Option<i32>,
    // --- Constraints ---
    pub constraints: Constraints,
    // ---
    pub free_energy_calculations: FreeEnergyCalculations,
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
            thermostat: Thermostat::default(),
            tau_t: vec![1.0], // todo: 0.1?
            ref_t: vec![300.],
            pcoupl: Barostat::default(),
            pbc: Pbc::default(),
            gen_vel: true,
            gen_temp: 300.,
            gen_seed: None,
            constraints: Constraints::default(),
            energy_minimization: Default::default(),
            free_energy_calculations: Default::default(),
        }
    }
}

impl MdpParams {
    /// Render as a GROMACS `.mdp` file string.
    pub fn to_mdp_str(&self) -> String {
        let mut s = String::from("; GROMACS MDP - generated by Bio Files\n\n");

        // Run control
        s.push_str("; Run control\n");
        append_inp(&mut s, "integrator", self.integrator.keyword());
        append_inp(&mut s, "nsteps", self.nsteps);
        append_inp(&mut s, "dt", Gf(self.dt));

        // Output
        s.push_str("\n; Output control\n");
        append_inp(&mut s, "nstxout", self.nstxout);
        append_inp(&mut s, "nstvout", self.nstvout);
        append_inp(&mut s, "nstxout-compressed", self.nstxout_compressed);
        append_inp(&mut s, "nstenergy", self.nstenergy);
        append_inp(&mut s, "nstlog", self.nstlog);

        // Non-bonded
        s.push_str("\n; Non-bonded interactions\n");
        s.push_str(&self.coulombtype.make_inp(self.rcoulomb));

        append_inp(&mut s, "rcoulomb", Gf(self.rcoulomb));
        append_inp(&mut s, "vdw-type", self.vdwtype.keyword());
        append_inp(&mut s, "rvdw", Gf(self.rvdw));

        // Temperature coupling
        s.push_str("\n; Temperature coupling\n");
        append_inp(&mut s, "tcoupl", self.thermostat.keyword());
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
                .map(|v| Gf(*v).to_string())
                .collect::<Vec<_>>()
                .join(" ");

            let ref_t = self
                .ref_t
                .iter()
                .map(|v| Gf(*v).to_string())
                .collect::<Vec<_>>()
                .join(" ");
            append_inp(&mut s, "tc-grps", tc_grps);
            append_inp(&mut s, "tau-t", tau_t);
            append_inp(&mut s, "ref-t", ref_t);
        }

        // Pressure coupling
        s.push_str("\n; Pressure coupling\n");
        s.push_str(&self.pcoupl.make_inp());

        // PBC
        s.push_str("\n; Periodic boundary conditions\n");
        append_inp(&mut s, "pbc", self.pbc.keyword());

        // Velocity generation
        s.push_str("\n; Velocity generation\n");
        append_inp(&mut s, "gen-vel", if self.gen_vel { "yes" } else { "no" });
        if self.gen_vel {
            append_inp(&mut s, "gen-temp", Gf(self.gen_temp));

            if let Some(v) = &self.gen_seed {
                append_inp(&mut s, "gen-seed", v);
            }
        }

        // Constraints
        s.push_str("\n; Constraints\n");
        s.push_str(&self.constraints.make_inp());

        if let Some(em) = &self.energy_minimization {
            s.push_str("\n; Energy minimization\n");

            s.push_str(&em.make_inp());
        }

        append_inp(
            &mut s,
            "free-energy",
            self.free_energy_calculations.keyword(),
        );

        s
    }

    pub fn from_mdp_str(data: &str) -> io::Result<Self> {
        fn inv(k: &str, e: impl std::fmt::Display) -> io::Error {
            Error::new(ErrorKind::InvalidData, format!("{k}: {e}"))
        }

        fn get_s<'a>(map: &'a HashMap<String, String>, k: &str) -> Option<&'a str> {
            map.get(k).map(String::as_str)
        }

        fn f32_or(map: &HashMap<String, String>, k: &str, d: f32) -> io::Result<f32> {
            match map.get(k) {
                None => Ok(d),
                Some(v) => v.parse().map_err(|e| inv(k, e)),
            }
        }

        fn u8_or(map: &HashMap<String, String>, k: &str, d: u8) -> io::Result<u8> {
            match map.get(k) {
                None => Ok(d),
                Some(v) => v.parse().map_err(|e| inv(k, e)),
            }
        }

        fn u16_or(map: &HashMap<String, String>, k: &str, d: u16) -> io::Result<u16> {
            match map.get(k) {
                None => Ok(d),
                Some(v) => v.parse().map_err(|e| inv(k, e)),
            }
        }

        fn u32_or(map: &HashMap<String, String>, k: &str, d: u32) -> io::Result<u32> {
            match map.get(k) {
                None => Ok(d),
                Some(v) => v.parse().map_err(|e| inv(k, e)),
            }
        }

        fn u64_or(map: &HashMap<String, String>, k: &str, d: u64) -> io::Result<u64> {
            match map.get(k) {
                None => Ok(d),
                Some(v) => v.parse().map_err(|e| inv(k, e)),
            }
        }

        fn i32_or(map: &HashMap<String, String>, k: &str, d: i32) -> io::Result<i32> {
            match map.get(k) {
                None => Ok(d),
                Some(v) => v.parse().map_err(|e| inv(k, e)),
            }
        }

        fn f32_vec(map: &HashMap<String, String>, k: &str, d: Vec<f32>) -> io::Result<Vec<f32>> {
            match map.get(k) {
                None => Ok(d),
                Some(v) => v
                    .split_whitespace()
                    .map(|s| s.parse::<f32>().map_err(|e| inv(k, e)))
                    .collect(),
            }
        }

        // Parse all key = value pairs; strip inline comments and lowercase keys.
        let mut map = HashMap::<String, String>::new();
        for line in data.lines() {
            let line = line.split(';').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            if let Some((k, v)) = line.split_once('=') {
                map.insert(k.trim().to_ascii_lowercase(), v.trim().to_owned());
            }
        }

        let integrator_s = get_s(&map, "integrator")
            .unwrap_or("md")
            .to_ascii_lowercase();

        let integrator = match integrator_s.as_str() {
            "md" => Integrator::Md,
            "md-vv" => Integrator::MdVv,
            "md-vv-avek" => Integrator::MdVvAvek,
            "bd" => Integrator::Bd,
            "sd" => Integrator::Sd,
            "steep" => Integrator::Steep,
            "cg" => Integrator::Cg,
            "l-bfgs" => Integrator::LBfgs,
            "nm" => Integrator::Nm,
            "tpi" => Integrator::Tpi,
            "tpic" => Integrator::Tpic,
            "mimic" => Integrator::Mimic,
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown integrator: {o}"),
                ));
            }
        };

        let nsteps = u64_or(&map, "nsteps", 0)?;
        let dt = f32_or(&map, "dt", 0.001)?;
        let nstxout = u32_or(&map, "nstxout", 0)?;
        let nstvout = u32_or(&map, "nstvout", 0)?;
        let nstxout_compressed = u32_or(&map, "nstxout-compressed", 0)?;
        let nstenergy = u32_or(&map, "nstenergy", 1_000)?;
        let nstlog = u32_or(&map, "nstlog", 1_000)?;

        // Parse rcoulomb early — needed to convert ewald-rtol → alpha.
        let rcoulomb = f32_or(&map, "rcoulomb", 1.0)?;

        let coulomb_s = get_s(&map, "coulombtype")
            .unwrap_or("pme")
            .to_ascii_lowercase();

        let coulombtype = match coulomb_s.as_str() {
            "pme" => {
                let d = PmeConfig::default();
                let fourierspacing = f32_or(&map, "fourierspacing", d.fourierspacing)?;
                let order = u8_or(&map, "pme-order", d.order)?;
                // Parse ewald-rtol (GROMACS default 1e-5), then invert to alpha (nm⁻¹).
                let rtol = f32_or(&map, "ewald-rtol", 1e-5_f32)?;
                let rtol_lj = f32_or(&map, "ewald-rtol-lj", d.rtol_lj)?;
                let epsilon_surface = map
                    .get("epsilon-surface")
                    .map(|v| v.parse::<f32>().map_err(|e| inv("epsilon-surface", e)))
                    .transpose()?;
                CoulombType::Pme(PmeConfig {
                    fourierspacing,
                    order,
                    alpha: erfc_inv_approx(rtol) / rcoulomb,
                    rtol_lj,
                    epsilon_surface,
                })
            }
            "cut-off" => CoulombType::CutOff,
            "reaction-field" => CoulombType::ReactionField,
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown coulombtype: {o}"),
                ));
            }
        };

        let vdw_s = get_s(&map, "vdw-type")
            .unwrap_or("cut-off")
            .to_ascii_lowercase();

        let vdwtype = match vdw_s.as_str() {
            "cut-off" => VdwType::CutOff,
            "pme" => VdwType::Pme,
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown vdw-type: {o}"),
                ));
            }
        };
        let rvdw = f32_or(&map, "rvdw", 1.0)?;

        let tcoupl_s = get_s(&map, "tcoupl")
            .unwrap_or("v-rescale")
            .to_ascii_lowercase();

        let thermostat = match tcoupl_s.as_str() {
            "no" => Thermostat::No,
            "nose-hoover" => Thermostat::NoseHoover,
            "andersen" => Thermostat::Andersen,
            "andersen-massive" => Thermostat::AndersenMassive,
            "v-rescale" => Thermostat::VRescale,
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown tcoupl: {o}"),
                ));
            }
        };
        let tau_t = f32_vec(&map, "tau-t", vec![1.0])?;
        let ref_t = f32_vec(&map, "ref-t", vec![300.0])?;

        let pcoupl_s = get_s(&map, "pcoupl").unwrap_or("no").to_ascii_lowercase();

        let pcoupltype_s = get_s(&map, "pcoupltype")
            .unwrap_or("isotropic")
            .to_ascii_lowercase();

        let pcoupltype = match pcoupltype_s.as_str() {
            "isotropic" => PressureCouplingType::Isotropic {
                compressibility: f32_or(&map, "compressibility", 4.5e-5)?,
                ref_p: f32_or(&map, "ref-p", 1.0)?,
            },
            "semiisotropic" => {
                let comp = f32_vec(&map, "compressibility", vec![4.5e-5, 4.5e-5])?;
                let rp = f32_vec(&map, "ref-p", vec![1.0, 1.0])?;
                PressureCouplingType::Semiisotropic {
                    compressibility_xy: comp.first().copied().unwrap_or(4.5e-5),
                    compressibility_z: comp.get(1).copied().unwrap_or(4.5e-5),
                    ref_p_xy: rp.first().copied().unwrap_or(1.0),
                    ref_p_z: rp.get(1).copied().unwrap_or(1.0),
                }
            }
            "anisotropic" => {
                let comp = f32_vec(&map, "compressibility", vec![4.5e-5; 6])?;
                let rp = f32_vec(&map, "ref-p", vec![1.0; 6])?;
                let mut compressibility = [4.5e-5_f32; 6];
                let mut ref_p = [1.0_f32; 6];
                for (i, v) in comp.iter().enumerate().take(6) {
                    compressibility[i] = *v;
                }
                for (i, v) in rp.iter().enumerate().take(6) {
                    ref_p[i] = *v;
                }
                PressureCouplingType::Anisotropic {
                    compressibility,
                    ref_p,
                }
            }
            "surface-tension" => {
                let comp = f32_vec(&map, "compressibility", vec![4.5e-5, 4.5e-5])?;
                let rp = f32_vec(&map, "ref-p", vec![0.0, 1.0])?;
                PressureCouplingType::SurfaceTension {
                    compressibility_xy: comp.first().copied().unwrap_or(4.5e-5),
                    compressibility_z: comp.get(1).copied().unwrap_or(4.5e-5),
                    ref_p_sfc_tension: rp.first().copied().unwrap_or(0.0),
                    ref_p_z: rp.get(1).copied().unwrap_or(1.0),
                }
            }
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown pcoupltype: {o}"),
                ));
            }
        };

        let tau_p = f32_or(&map, "tau-p", 5.0)?;
        let barostat_cfg = BarostatCfg { pcoupltype, tau_p };

        let pcoupl = match pcoupl_s.as_str() {
            "no" => Barostat::No,
            "berendsen" => Barostat::Berendsen(barostat_cfg),
            "c-rescale" => Barostat::CRescale(barostat_cfg),
            "parrinello-rahman" => Barostat::ParrinelloRahman(barostat_cfg),
            "mttk" => Barostat::Mtkk(barostat_cfg),
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown pcoupl: {o}"),
                ));
            }
        };

        let pbc_s = get_s(&map, "pbc").unwrap_or("xyz").to_ascii_lowercase();
        let pbc = match pbc_s.as_str() {
            "xyz" => Pbc::Xyz,
            "no" => Pbc::No,
            "xy" => Pbc::Xy,
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown pbc: {o}"),
                ));
            }
        };

        let gen_vel_s = get_s(&map, "gen-vel").unwrap_or("no").to_ascii_lowercase();
        let gen_vel = matches!(gen_vel_s.as_str(), "yes" | "true" | "1");
        let gen_temp = f32_or(&map, "gen-temp", 300.0)?;

        let mut gen_seed = Some(i32_or(&map, "gen-seed", -1)?);
        if gen_seed == Some(-1) {
            gen_seed = None;
        }

        // Handle both the correct spelling and the typo in ConstraintAlgorithm::make_inp.
        let ca_s = map
            .get("constraint-algorithm")
            .or_else(|| map.get("constraint-algorithm"))
            .map(String::as_str)
            .unwrap_or("lincs")
            .to_ascii_lowercase();

        let constraint_algorithm = match ca_s.as_str() {
            "lincs" => ConstraintAlgorithm::Lincs {
                order: u8_or(&map, "lincs-order", 4)?,
                iter: u8_or(&map, "lincs-iter", 1)?,
            },
            "shake" => ConstraintAlgorithm::Shake {
                tol: f32_or(&map, "shake-tol", 0.0001)?,
            },
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown constraint-algorithm: {o}"),
                ));
            }
        };

        let constr_s = get_s(&map, "constraints")
            .unwrap_or("h-bonds")
            .to_ascii_lowercase();

        let constraints = match constr_s.as_str() {
            "none" => Constraints::None,
            "h-bonds" => Constraints::HBonds(constraint_algorithm),
            "all-bonds" => Constraints::AllBonds(constraint_algorithm),
            "h-angles" => Constraints::HAngles(constraint_algorithm),
            "all-angles" => Constraints::AllAngles(constraint_algorithm),
            o => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("unknown constraints: {o}"),
                ));
            }
        };

        let em_keys = ["emtol", "emstep", "nstcgsteep", "nbfgscorr"];
        let has_em = em_keys.iter().any(|k| map.contains_key(*k));
        let energy_minimization = if has_em {
            let d = EnergyMinimization::default();
            Some(EnergyMinimization {
                emtol: f32_or(&map, "emtol", d.emtol)?,
                emstep: f32_or(&map, "emstep", d.emstep)?,
                nstcgsteep: u16_or(&map, "nstcgsteep", d.nstcgsteep)?,
                nbfgscorr: u8_or(&map, "nbfgscorr", d.nbfgscorr)?,
            })
        } else {
            None
        };

        let free_energy_calculations = match get_s(&map, "free-energy").unwrap_or("no") {
            "yes" => FreeEnergyCalculations::Yes,
            "expanded" => FreeEnergyCalculations::Expanded,
            _ => FreeEnergyCalculations::No,
        };

        Ok(Self {
            integrator,
            nsteps,
            dt,
            nstxout,
            nstvout,
            nstxout_compressed,
            nstenergy,
            nstlog,
            coulombtype,
            rcoulomb,
            vdwtype,
            rvdw,
            thermostat,
            tau_t,
            ref_t,
            pcoupl,
            pbc,
            gen_vel,
            gen_temp,
            gen_seed,
            constraints,
            energy_minimization,
            free_energy_calculations,
        })
    }

    /// Save to a .mdp file
    pub fn save(&self, path: &Path) -> io::Result<()> {
        save_txt_to_file(path, &self.to_mdp_str())
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        Self::from_mdp_str(&std::fs::read_to_string(path)?)
    }
}
