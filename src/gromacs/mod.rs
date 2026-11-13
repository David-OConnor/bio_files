//! Interop with [GROMACS](https://www.gromacs.org/) molecular dynamics software.
//! [GROMACS User guide](https://manual.gromacs.org/current/user-guide/index.html)
//! [GROMACS reference Manual](https://manual.gromacs.org/2026.1/reference-manual/index.html)
//! [GROMACS manual on Zenodo](https://zenodo.org/records/18886967)
//!
//! ## Overview
//!
//! This module constructs GROMACS input files (`.mdp`, `.gro`, `.top`), runs
//! the simulation via the `gmx` command-line tool, and parses the resulting
//! trajectory into a portable format.
//!
//! ## Workflow
//!
//! 1. Build a [`GromacsInput`] from molecule atoms/bonds, force-field params,
//!    and an [`MdpParams`] config.
//! 2. Call [`GromacsInput::save`] to write inputs to a directory, or
//!    [`GromacsInput::run`] to execute the full pipeline and collect output.
//!
//! The `run` workflow:
//! - Writes inputs and runs GROMACS in a new `md_out/gromacs_{timestamp}/` directory.
//! - Keeps the native files and console logs, and parses the trajectory.
//!
//! Units throughout follow the rest of the codebase (Å for positions).
//! Conversions to GROMACS units (nm) happen inside this module.

pub mod ff_params;
pub mod forcefield;
pub mod gro;
pub mod mdp;
pub mod output;
pub mod rtp;
pub mod solvate;
pub mod top_parse;
pub mod topology;
pub mod trr;

use std::{
    fs::{self, File, OpenOptions},
    io::{self, ErrorKind, Write},
    path::{Path, PathBuf},
    process::{Command, Output, Stdio},
    sync::OnceLock,
};

pub use ff_params::GromacsFfParams;
pub use mdp::{MdpParams, OutputControl};
pub use output::{GromacsFrame, GromacsOutput, OutputEnergy};
pub use rtp::Rtp;
use solvate::Solvent;
pub use top_parse::GromacsTopology;
pub use topology::MoleculeTopology;
use trr::read_trr;

use crate::{
    AtomGeneric, BondGeneric, FrameSlice, md_out::create_run_dir, md_params::ForceFieldParams,
};

pub(in crate::gromacs) const GRO_NAME: &str = "conf.gro";
pub(in crate::gromacs) const TOP_NAME: &str = "topo.top";
pub(in crate::gromacs) const MDP_NAME: &str = "md.mdp";

pub(in crate::gromacs) const CUSTOM_SOLVENT_GRO: &str = "solvent.gro";

// This must match the 4-site water template in the GROMACS data dir.
pub(in crate::gromacs) const GMX_BUILTIN_TIP4P_WATER: &str = "tip4p.gro";

pub(in crate::gromacs) const SOLVATED_NAME: &str = "mols_in_solvated.gro";
pub(in crate::gromacs) const IONIZED_NAME: &str = "ionized.gro";

pub(in crate::gromacs) const TRR_NAME: &str = "traj.trr";
pub(in crate::gromacs) const XTC_NAME: &str = "traj.xtc";
pub(in crate::gromacs) const GRO_OUT_NAME: &str = "confout.gro";
pub(in crate::gromacs) const ENERGY_OUT_NAME: &str = "energy.edr";

const LOG_NAME: &str = "md.log";

/// One molecule entry passed to a [`GromacsInput`].
///
/// Mirrors the per-molecule data used by the `dynamics` crate, adapted for
/// GROMACS input generation.
#[derive(Clone, Debug)]
pub struct MoleculeInput {
    /// Short identifier used as the molecule name in the topology (e.g. `"MOL"`, `"PEP"`).
    pub name: String,
    pub atoms: Vec<AtomGeneric>,
    pub bonds: Vec<BondGeneric>,
    // todo: We don't want to have a copy of the whole Gaff2 etc for each input.
    /// Molecule-specific force-field parameters (e.g. GAFF2 for a ligand).
    pub ff_params: Option<ForceFieldParams>,
    /// Number of copies to list in `[ molecules ]`.
    pub count: usize,
    /// Optional atom positions for each explicitly packed copy.
    ///
    /// When omitted, each copy uses the positions stored in `atoms`. When set,
    /// the outer length must equal `count`, and each inner vector must contain
    /// one position for every atom in this molecule topology.
    pub copy_atom_posits: Option<Vec<Vec<lin_alg::f64::Vec3>>>,
}

/// GROMACS simulation input. This contains everything required to launch a simulation.
///
/// Analogous to `orca::OrcaInput`.
#[derive(Clone, Debug)]
pub struct GromacsInput {
    /// MD simulation parameters (timestep, thermostat, etc.).
    pub mdp: MdpParams,
    /// All molecules that participate in the simulation, in order.
    pub molecules: Vec<MoleculeInput>,
    /// Simulation box dimensions in **nm** (`x, y, z`).
    /// `None` → box line written as `0 0 0` (GROMACS will error on `grompp` —
    /// provide a box when running production MD).
    pub box_nm: Option<(f64, f64, f64)>,
    /// Optional lower corner of the input coordinate frame in **Å**.
    ///
    /// When set, solute coordinates are translated from `[origin_a, ...]` into
    /// GROMACS' `[0, box_nm]` frame instead of being recentered by centroid.
    /// Use this for deliberately offset structures such as boundary layers.
    pub coordinate_origin_a: Option<lin_alg::f64::Vec3>,
    /// System-wide (global) force-field parameters used as a fall-back for
    /// any per-molecule terms that are absent from `MoleculeInput::ff_params`.
    pub ff_global: Option<ForceFieldParams>,
    /// Non-bonded conventions, water and ions. Supplied molecule parameters
    /// must belong to this family (selection does not re-type molecules).
    pub topology_config: forcefield::TopologyConfig,
    /// When set, `run()` will call `gmx solvate` to fill the box with water
    /// before preprocessing, and include the model's topology in the `.top`.
    pub solvent: Option<Solvent>,
    /// Optional complete starting structure in GRO format.
    ///
    /// Use this for callers that have already built a solvated or otherwise
    /// pre-packed system and only want this module to own the GROMACS file/run
    /// pipeline. When set, `run()` skips `gmx solvate` even if `solvent` is set;
    /// the solvent marker can still be used to include solvent topology such as
    /// OPC water.
    pub initial_gro: Option<String>,
    /// Extra `[ molecules ]` entries to append to the generated topology.
    ///
    /// This pairs with `initial_gro` for prebuilt systems containing molecule
    /// types whose topology is supplied by `solvent`, for example `SOL` from
    /// Amber OPC water.
    pub extra_molecule_counts: Vec<(String, usize)>,
    /// Optional complete topology text for a continued run.
    ///
    /// Use this with `initial_gro` when resuming from a previous `mdrun`, so
    /// counter-ion substitutions and molecule counts remain identical.
    pub topology_override: Option<String>,
    /// Additional arguments appended to the production `gmx mdrun` command.
    ///
    /// For example, a caller driving a shrinking box can use `["-ntmpi", "1"]`
    /// to disable domain decomposition as the cell becomes smaller.
    pub mdrun_extra_args: Vec<String>,
    /// Skip automatic counter-ion insertion. Continued runs should enable this
    /// after reusing the topology emitted by the first run.
    pub skip_counterion_insertion: bool,
    pub minimize_energy: bool,
}

impl Default for GromacsInput {
    fn default() -> Self {
        Self {
            mdp: MdpParams::default(),
            molecules: Vec::new(),
            box_nm: None,
            coordinate_origin_a: None,
            ff_global: None,
            topology_config: Default::default(),
            solvent: None,
            initial_gro: None,
            extra_molecule_counts: Vec::new(),
            topology_override: None,
            mdrun_extra_args: Vec::new(),
            skip_counterion_insertion: false,
            minimize_energy: true,
        }
    }
}

impl GromacsInput {
    /// Generate the `.mdp` file contents.
    pub fn make_mdp(&self) -> String {
        self.mdp.to_mdp_str()
    }

    /// Generate the GROMACS `.top` topology file.
    pub fn make_top(&self) -> io::Result<String> {
        if let Some(top) = &self.topology_override {
            return Ok(top.clone());
        }
        let mol_tops: Vec<MoleculeTopology<'_>> = self
            .molecules
            .iter()
            .map(|m| MoleculeTopology {
                name: &m.name,
                atoms: &m.atoms,
                bonds: &m.bonds,
                ff_mol: m.ff_params.as_ref(),
                count: m.count,
            })
            .collect();

        let solvent_mol_tops: Vec<MoleculeTopology<'_>> = match &self.solvent {
            Some(Solvent::Custom(template)) => template
                .topology_molecules
                .iter()
                .map(|m| MoleculeTopology {
                    name: &m.name,
                    atoms: &m.atoms,
                    bonds: &m.bonds,
                    ff_mol: m.ff_params.as_ref(),
                    count: m.count,
                })
                .collect(),
            _ => Vec::new(),
        };

        let mut top = topology::make_top_with_config(
            &mol_tops,
            &solvent_mol_tops,
            self.ff_global.as_ref(),
            self.solvent.as_ref(),
            &self.topology_config,
        )?;

        for (name, count) in &self.extra_molecule_counts {
            top.push_str(&format!("{:<14}  {}\n", name, count));
        }

        Ok(top)
    }

    /// Total number of solute (non-water) atoms across all molecule copies.
    pub fn solute_atom_count(&self) -> usize {
        self.molecules.iter().map(|m| m.atoms.len() * m.count).sum()
    }

    /// Write all input files (`conf.gro`, `topol.top`, `md.mdp`) to `dir`.
    /// When a water model is set, also writes the solvent box
    /// so that `gmx solvate -cs <model>.gro` works without requiring the file
    /// to be pre-installed in the GROMACS data directory.
    pub fn save(&self, dir: &Path) -> io::Result<()> {
        self.save_with_mdp(dir, &self.make_mdp(), None)
    }

    /// Save a production MDP verbatim, including settings outside MdpParams.
    pub fn save_with_mdp(&self, dir: &Path, mdp: &str, index: Option<&Path>) -> io::Result<()> {
        fs::create_dir_all(dir)?;

        save_txt_to_file(dir.join(MDP_NAME), mdp)?;
        if let Some(index) = index {
            let dest = dir.join("index.ndx");
            if index.canonicalize().ok() != dest.canonicalize().ok() {
                fs::copy(index, dest)?;
            }
        }
        if let Some(gro_text) = &self.initial_gro {
            save_txt_to_file(dir.join(GRO_NAME), gro_text)?;
        } else {
            save_txt_to_file(
                dir.join(GRO_NAME),
                &gro::make_gro_with_origin(
                    &self.molecules,
                    &self.box_nm,
                    self.coordinate_origin_a,
                )?,
            )?;
        }
        let topology = match &self.topology_override {
            Some(topology) => topology.clone(),
            None => self.make_top()?,
        };
        save_txt_to_file(dir.join(TOP_NAME), &topology)?;

        if self.initial_gro.is_none()
            && let Some(Solvent::Custom(template)) = &self.solvent
        {
            save_txt_to_file(dir.join(CUSTOM_SOLVENT_GRO), &template.gro_text)?;
        }

        Ok(())
    }

    /// Run GROMACS and keep its files in a new `md_out/gromacs_{timestamp}/` directory.
    pub fn run(&self) -> io::Result<GromacsOutput> {
        self.run_with_mdp(&self.make_mdp(), None)
    }

    /// Run with a production MDP and optional saved index file.
    pub fn run_with_mdp(&self, mdp: &str, index: Option<&Path>) -> io::Result<GromacsOutput> {
        let dir = create_run_dir("gromacs")?;

        self.run_in(&dir, mdp, index)
            .map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", dir.display())))
    }

    fn run_in(&self, dir: &Path, mdp: &str, index: Option<&Path>) -> io::Result<GromacsOutput> {
        self.save_with_mdp(dir, mdp, index)?;

        let solute_atom_count = self.solute_atom_count();
        let run_command = |args: &[&str]| run_gmx_logged(dir, args, None);

        // Solvation: fill the box with water before preprocessing.
        // We don't specify `box` (aka cell/sim box); it's present in the solute coordinate file (`-cp`).
        // -cs specifies the solute template. `tip4p.gro` is included with GROMACS; we use that.
        // [Docs for GMX solvate](https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html#gmx-solvate)
        let structure_gro = if self.initial_gro.is_none() {
            if let Some(ref wm) = self.solvent
                && let Some(solvent_gro) = wm.prepare_gro(
                    dir,
                    self.box_nm,
                    if self.topology_config.water().om == 0. {
                        "spc216.gro"
                    } else {
                        GMX_BUILTIN_TIP4P_WATER
                    },
                )?
            {
                run_command(&[
                    "solvate",
                    "-cp",
                    GRO_NAME,
                    "-cs",
                    solvent_gro,
                    "-o",
                    SOLVATED_NAME,
                    "-p",
                    TOP_NAME,
                ])?;
                SOLVATED_NAME
            } else {
                GRO_NAME
            }
        } else {
            GRO_NAME
        };

        // Add counter-ions if the solute carries a net charge.
        // Only meaningful when solvent is present (ions replace water molecules).
        // `gmx genion -neutral` adds exactly the number of Na+/Cl- needed.
        let net_q: f32 = self
            .molecules
            .iter()
            .map(|m| {
                m.atoms
                    .iter()
                    .map(|a| a.partial_charge.unwrap_or(0.0))
                    .sum::<f32>()
                    * m.count as f32
            })
            .sum();

        let structure_gro =
            if !self.skip_counterion_insertion && self.solvent.is_some() && net_q.abs() >= 0.5 {
                // grompp needs a valid MDP to build ions.tpr; reuse the EM MDP.
                save_txt_to_file(dir.join("ions.mdp"), em_mdp_str())?;
                run_command(&[
                    "grompp",
                    "-f",
                    "ions.mdp",
                    "-c",
                    structure_gro,
                    "-p",
                    TOP_NAME,
                    "-o",
                    "ions.tpr",
                    "-po",
                    "ions_out.mdp",
                    "-maxwarn",
                    "5",
                ])?;

                // genion asks interactively which group to replace; we select SOL.
                run_gmx_logged(
                    dir,
                    &[
                        "genion",
                        "-s",
                        "ions.tpr",
                        "-o",
                        IONIZED_NAME,
                        "-p",
                        TOP_NAME,
                        "-pname",
                        "NA",
                        "-nname",
                        "CL",
                        "-neutral",
                    ],
                    Some(b"SOL\n"),
                )?;
                IONIZED_NAME
            } else {
                structure_gro
            };

        // Energy minimization — removes bad contacts
        let md_input_gro: String = if self.minimize_energy {
            // A preset configuration file.
            save_txt_to_file(dir.join("em.mdp"), em_mdp_str())?;

            run_command(&[
                "grompp",
                "-f",
                "em.mdp",
                "-c",
                structure_gro,
                "-p",
                TOP_NAME,
                "-o",
                "em.tpr",
                "-po",
                "em_out.mdp",
                "-maxwarn",
                "5",
            ])?;

            run_command(&["mdrun", "-deffnm", "em"])?;

            "em.gro".to_string()
        } else {
            structure_gro.to_string()
        };

        // grompp: Preprocessor. Creates a [binary] TPR file from the generated input files.
        // [Data on TPR](https://manual.gromacs.org/2026.1/reference-manual/file-formats.html#tpr)
        // This format contains everything GROMACS needs to run MD.
        let mut grompp_args = vec![
            "grompp",
            "-f",
            MDP_NAME,
            "-c",
            &md_input_gro,
            "-p",
            TOP_NAME,
            "-o",
            "topol.tpr",
            "-maxwarn",
            "5",
        ];
        if index.is_some() {
            grompp_args.extend(["-n", "index.ndx"]);
        }
        if self.initial_gro.is_some() {
            grompp_args.extend(["-r", &md_input_gro]);
        }
        run_command(&grompp_args)?;

        // mdrun: Run MD.
        let mut mdrun_args = vec![
            "mdrun",
            "-s",
            "topol.tpr",
            "-o",
            TRR_NAME,
            "-x",
            XTC_NAME,
            "-c",
            GRO_OUT_NAME,
            "-e",
            ENERGY_OUT_NAME,
            "-g",
            LOG_NAME,
        ];
        mdrun_args.extend(self.mdrun_extra_args.iter().map(String::as_str));
        run_command(&mdrun_args)?;

        let trr = dir.join(TRR_NAME);
        let xtc = dir.join(XTC_NAME);
        let gro = dir.join(&md_input_gro);

        let log_text = read_text(dir.join(LOG_NAME)).unwrap_or_default();
        let trr_frames = if trr.exists() {
            read_trr(
                &trr,
                FrameSlice::Time {
                    start: None,
                    end: None,
                },
            )?
        } else {
            Vec::new()
        };

        // Energy data is optional: if gmx energy fails or is unavailable, run()
        // still returns a valid trajectory - frames just have `energy: None`.
        let energies =
            OutputEnergy::from_edr_impl(&dir.join(ENERGY_OUT_NAME), true).unwrap_or_default();

        let mut result = GromacsOutput::new(log_text, trr_frames, energies, solute_atom_count)?;
        result.trr_path = trr.exists().then_some(trr);
        result.xtc_path = xtc.exists().then_some(xtc);
        result.gro_path = gro.exists().then_some(gro);
        result.final_gro_text = read_text(dir.join(GRO_OUT_NAME)).ok();
        result.final_topology_text = read_text(dir.join(TOP_NAME)).ok();

        Ok(result)
    }
}

/// Minimal MDP for a steep-descent energy minimization pass.
///
/// Run before production MD to remove bad contacts — especially important
/// when the solvent was placed on a regular lattice by `opc_water_box_gro`.
/// `emtol = 1000` kJ/mol/nm is intentionally loose: we only need to eliminate
/// clashes, not find the true minimum.  `constraints = none` lets GROMACS move
/// hydrogen atoms freely; SETTLE still applies to water via `[ settles ]`.
fn em_mdp_str() -> &'static str {
    "; Energy minimization — generated by Bio Files\n\
     integrator               = steep\n\
     nsteps                   = 5000\n\
     emtol                    = 1000.0\n\
     emstep                   = 0.01\n\
     \n\
     cutoff-scheme            = Verlet\n\
     coulombtype              = PME\n\
     fourierspacing           = 0.16\n\
     rcoulomb                 = 1.0\n\
     vdw-type                 = Cut-off\n\
     rvdw                     = 1.0\n\
     \n\
     pbc                      = xyz\n\
     constraints              = none\n"
}

fn save_txt_to_file(path: impl AsRef<Path>, text: &str) -> io::Result<()> {
    let mut f = File::create(path)?;
    write!(f, "{text}")
}

fn read_text(path: impl AsRef<Path>) -> io::Result<String> {
    fs::read_to_string(path)
}

static GMX_RESOLVER: OnceLock<fn() -> Option<PathBuf>> = OnceLock::new();

/// Sets how to locate `gmx`, for applications that install it off `PATH`. Called before each
/// run; `None` falls back to `PATH`. Only the first call has an effect.
pub fn set_gmx_resolver(resolver: fn() -> Option<PathBuf>) {
    let _ = GMX_RESOLVER.set(resolver);
}

fn gmx_command() -> Command {
    let executable = GMX_RESOLVER
        .get()
        .and_then(|resolve| resolve())
        .unwrap_or_else(|| PathBuf::from("gmx"));

    let mut cmd = Command::new(&executable);
    cmd.env("GMX_MAXBACKUP", "-1");
    if std::env::var_os("GMXLIB").is_none()
        && let Some(prefix) = executable.parent().and_then(Path::parent)
    {
        let top = prefix.join("share/gromacs/top");
        if top.join("residuetypes.dat").is_file() {
            cmd.env("GMXLIB", top);
        }
    }

    cmd
}

/// Run a `gmx` subcommand, returning an error on failure.
pub fn run_gmx(dir: &Path, args: &[&str]) -> io::Result<()> {
    check_gmx_status(args, &gmx_output(dir, args, None)?)
}

/// Like [`run_gmx`] but supplies `stdin_data` to the process.
pub fn run_gmx_stdin(dir: &Path, args: &[&str], stdin_data: &[u8]) -> io::Result<()> {
    check_gmx_status(args, &gmx_output(dir, args, Some(stdin_data))?)
}

fn run_gmx_logged(dir: &Path, args: &[&str], stdin_data: Option<&[u8]>) -> io::Result<()> {
    let out = gmx_output(dir, args, stdin_data)?;
    save_gmx_logs(dir, &out)?;

    check_gmx_status(args, &out)
}

fn gmx_output(dir: &Path, args: &[&str], stdin_data: Option<&[u8]>) -> io::Result<Output> {
    let mut cmd = gmx_command();
    cmd.current_dir(dir).args(args);

    let out = if let Some(stdin_data) = stdin_data {
        cmd.stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        cmd.spawn().and_then(|mut child| {
            let stdin_result = child.stdin.take().unwrap().write_all(stdin_data);
            let out = child.wait_with_output();

            if let Err(e) = stdin_result
                && e.kind() != ErrorKind::BrokenPipe
            {
                return Err(e);
            }

            out
        })
    } else {
        cmd.output()
    };

    match out {
        Err(e) if e.kind() == ErrorKind::NotFound => Err(io::Error::new(
            ErrorKind::NotFound,
            "`gmx` executable not found on the system PATH",
        )),
        result => result,
    }
}

fn save_gmx_logs(dir: &Path, out: &Output) -> io::Result<()> {
    for (name, bytes) in [("stdout.log", &out.stdout), ("stderr.log", &out.stderr)] {
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(dir.join(name))?
            .write_all(bytes)?;
    }

    Ok(())
}

fn check_gmx_status(args: &[&str], out: &Output) -> io::Result<()> {
    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        return Err(io::Error::other(format!(
            "gmx {} failed: {}",
            args.first().unwrap_or(&"?"),
            stderr,
        )));
    }

    Ok(())
}
