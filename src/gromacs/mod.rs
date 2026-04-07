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
//! 2. Call [`GromacsInput::save_input_files`] to write inputs to a directory, or
//!    [`GromacsInput::run`] to execute the full pipeline and collect output.
//!
//! The `run` workflow:
//! - Creates a temporary directory `gromacs_out/`
//! - Writes `conf.gro`, `topol.top`, `md.mdp`
//! - Runs `gmx grompp` to produce `topol.tpr`
//! - Runs `gmx mdrun` to produce `traj.trr`, `confout.gro`, `energy.edr`, `md.log`
//! - Converts `traj.trr` to a multi-frame `traj.gro` via `gmx trjconv`
//! - Parses the trajectory and log, then removes the temporary directory
//!
//! Units throughout follow the rest of the codebase (Å for positions).
//! Conversions to GROMACS units (nm) happen inside this module.

pub mod gro;
pub mod mdp;
pub mod output;
pub mod solvate;
pub mod topology;

use std::{
    fmt::Write as _,
    fs,
    fs::File,
    io::{self, ErrorKind, Write},
    path::Path,
    process::{Command, Stdio},
};

pub use mdp::{MdpParams, OutputControl};
pub use output::{GromacsFrame, GromacsOutput, OutputEnergy};
use solvate::Solvent;
pub use topology::MoleculeTopology;

use crate::{
    AtomGeneric, BondGeneric, FrameSlice,
    gromacs::{gro::make_gro, output::read_trr},
    md_params::ForceFieldParams,
};

// Used for creating intermediate files
const TEMP_DIR: &str = "gromacs_out";
// Permanent output directory — numbered files are written here.
const MD_OUT_DIR: &str = "md_out";

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
}

/// GROMACS simulation input.
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
    /// System-wide (global) force-field parameters used as a fall-back for
    /// any per-molecule terms that are absent from `MoleculeInput::ff_params`.
    pub ff_global: Option<ForceFieldParams>,
    /// When set, `run()` will call `gmx solvate` to fill the box with water
    /// before preprocessing, and include the model's topology in the `.top`.
    pub solvent: Option<Solvent>,
    pub minimize_energy: bool,
}

impl Default for GromacsInput {
    fn default() -> Self {
        Self {
            mdp: MdpParams::default(),
            molecules: Vec::new(),
            box_nm: None,
            ff_global: None,
            solvent: None,
            minimize_energy: true,
        }
    }
}

// ---------------------------------------------------------------------------
// File generation
// ---------------------------------------------------------------------------

impl GromacsInput {
    /// Generate the `.mdp` file contents.
    pub fn make_mdp(&self) -> String {
        self.mdp.to_mdp_str()
    }

    /// Generate the GROMACS `.top` topology file.
    pub fn make_top(&self) -> io::Result<String> {
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

        topology::make_top(&mol_tops, self.ff_global.as_ref(), self.solvent.as_ref())
    }

    /// Total number of solute (non-water) atoms across all molecule copies.
    pub fn solute_atom_count(&self) -> usize {
        self.molecules.iter().map(|m| m.atoms.len() * m.count).sum()
    }

    /// Write all input files (`conf.gro`, `topol.top`, `md.mdp`) to `dir`.
    /// When a water model is set, also writes the solvent box (e.g. `opc.gro`)
    /// so that `gmx solvate -cs <model>.gro` works without requiring the file
    /// to be pre-installed in the GROMACS data directory.
    pub fn save(&self, dir: &Path) -> io::Result<()> {
        fs::create_dir_all(dir)?;

        save_txt_to_file(dir.join(MDP_NAME), &self.make_mdp())?;
        save_txt_to_file(dir.join(GRO_NAME), &make_gro(&self.molecules, &self.box_nm))?;
        save_txt_to_file(dir.join(TOP_NAME), &self.make_top()?)?;

        if let Some(Solvent::Custom(template)) = &self.solvent {
            save_txt_to_file(dir.join(CUSTOM_SOLVENT_GRO), &template.to_gro())?;
        }

        Ok(())
    }

    /// Run the full GROMACS pipeline and return parsed output.
    /// Requires `gmx` (GROMACS) to be available on the system `PATH`, and
    /// `.gro`, `.top`, and `.mdp` files already created.
    /// [Docs: grompp (preprocessor)](https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html)
    /// [Docs: mdrun](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html)
    ///
    /// Steps:
    /// 1. Write inputs to a temporary directory.
    /// 2. `gmx grompp` → `topol.tpr`
    /// 3. `gmx mdrun` → trajectory + log
    /// 4. `gmx trjconv` → multi-frame `.gro`
    /// 5. Parse and return [`GromacsOutput`].
    /// 6. Remove the temporary directory.
    pub fn run(&self) -> io::Result<GromacsOutput> {
        let dir = Path::new(TEMP_DIR);
        self.save(dir)?;

        let solute_atom_count = self.solute_atom_count();

        // Find the lowest run index N for which no output files exist yet.
        let out_dir = Path::new(MD_OUT_DIR);
        fs::create_dir_all(out_dir)?;
        let run_n = (1_usize..)
            .find(|&n| {
                !out_dir.join(format!("traj_{n}.trr")).exists()
                    && !out_dir.join(format!("traj_{n}.xtc")).exists()
            })
            .unwrap_or(1);

        // todo: Support custom solvents.

        // Solvation: fill the box with water before preprocessing.
        // We don't specify `box` (aka cell/sim box); it's present in the solute coordinate file (`-cp`).
        // -cs specifies the solute template. `opc.gro` is included with GROMACS; we use that.
        // [Docs for GMX solvate](https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html#gmx-solvate)
        let structure_gro = if let Some(ref wm) = self.solvent {
            run_gmx(
                dir,
                &[
                    "solvate",
                    "-cp",
                    GRO_NAME,
                    "-cs",
                    wm.gro_filename(),
                    "-o",
                    SOLVATED_NAME,
                    "-p",
                    TOP_NAME,
                ],
            )?;
            SOLVATED_NAME
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

        let structure_gro = if self.solvent.is_some() && net_q.abs() >= 0.5 {
            // grompp needs a valid MDP to build ions.tpr; reuse the EM MDP.
            save_txt_to_file(dir.join("ions.mdp"), em_mdp_str())?;
            run_gmx(
                dir,
                &[
                    "grompp",
                    "-f",
                    "ions.mdp",
                    "-c",
                    structure_gro,
                    "-p",
                    TOP_NAME,
                    "-o",
                    "ions.tpr",
                    "-maxwarn",
                    "5",
                ],
            )?;

            // genion asks interactively which group to replace; we select SOL.
            run_gmx_stdin(
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
                b"SOL\n",
            )?;
            IONIZED_NAME
        } else {
            structure_gro
        };

        // Energy minimization — removes bad contacts
        let md_input_gro: String = if self.minimize_energy {
            // A preset configuration file.
            save_txt_to_file(dir.join("em.mdp"), em_mdp_str())?;

            run_gmx(
                dir,
                &[
                    "grompp",
                    "-f",
                    "em.mdp",
                    "-c",
                    structure_gro,
                    "-p",
                    TOP_NAME,
                    "-o",
                    "em.tpr",
                    "-maxwarn",
                    "5",
                ],
            )?;

            run_gmx(
                dir,
                &[
                    "mdrun",
                    "-s",
                    "em.tpr",
                    "-c",
                    "em.gro",
                    "-e",
                    ENERGY_OUT_NAME,
                    "-g",
                    "em.log",
                ],
            )?;

            "em.gro".to_string()
        } else {
            structure_gro.to_string()
        };

        // grompp: Preprocessor. Creates a [binary] TPR file from the generated input files.
        // [Data on TPR](https://manual.gromacs.org/2026.1/reference-manual/file-formats.html#tpr)
        // This format contains everything GROMACS needs to run MD.
        run_gmx(
            dir,
            &[
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
            ],
        )?;

        // mdrun: Run MD.
        run_gmx(
            dir,
            &[
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
            ],
        )?;

        // Copy output files to numbered paths in the permanent output directory.
        let trr_src = dir.join(TRR_NAME);
        let trr_dest = out_dir.join(format!("traj_{run_n}.trr"));
        if trr_src.exists() {
            fs::copy(&trr_src, &trr_dest)?;
        }

        let xtc_src = dir.join(XTC_NAME);
        let xtc_dest = if xtc_src.exists() {
            let dest = out_dir.join(format!("traj_{run_n}.xtc"));
            fs::copy(&xtc_src, &dest)?;
            Some(dest)
        } else {
            None
        };

        let gro_src = dir.join(&md_input_gro);
        let gro_dest = out_dir.join(format!("mols_in_{run_n}.gro"));
        if gro_src.exists() {
            fs::copy(&gro_src, &gro_dest)?;
        }

        let log_text = read_text(dir.join("md.log")).unwrap_or_default();
        let trr_frames = read_trr(
            &trr_src,
            FrameSlice::Time {
                start: None,
                end: None,
            },
        )?;

        // Energy data is optional: if gmx energy fails or is unavailable, run()
        // still returns a valid trajectory — frames just have `energy: None`.
        let energies = OutputEnergy::from_edr(&dir.join(ENERGY_OUT_NAME)).unwrap_or_default();

        let mut result = GromacsOutput::new(log_text, trr_frames, energies, solute_atom_count)?;
        result.trr_path = if trr_dest.exists() {
            Some(trr_dest)
        } else {
            None
        };
        result.xtc_path = xtc_dest;
        result.gro_path = if gro_dest.exists() {
            Some(gro_dest)
        } else {
            None
        };

        Ok(result)
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

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

/// Run a `gmx` sub-command, returning an error if `gmx` is not found or the
/// process exits non-zero. This is primarily called by `run`, but can be called directly
/// by applications.
pub fn run_gmx(dir: &Path, args: &[&str]) -> io::Result<()> {
    let mut cmd = Command::new("gmx");
    cmd.current_dir(dir).args(args);

    let out = match cmd.output() {
        Ok(o) => o,
        Err(e) if e.kind() == ErrorKind::NotFound => {
            return Err(io::Error::new(
                ErrorKind::NotFound,
                "`gmx` executable not found on the system PATH",
            ));
        }
        Err(e) => return Err(e),
    };

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

/// Like [`run_gmx`] but supplies `stdin_data` to the process (for interactive
/// group-selection prompts such as `gmx trjconv`). This is primarily called by `run`, but can be called directly
/// by applications.
pub fn run_gmx_stdin(dir: &Path, args: &[&str], stdin_data: &[u8]) -> io::Result<()> {
    let mut cmd = Command::new("gmx");
    cmd.current_dir(dir)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let mut child = match cmd.spawn() {
        Ok(c) => c,
        Err(e) if e.kind() == ErrorKind::NotFound => {
            return Err(io::Error::new(
                ErrorKind::NotFound,
                "`gmx` executable not found on the system PATH",
            ));
        }
        Err(e) => return Err(e),
    };

    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(stdin_data)?;
    }

    let out = child.wait_with_output()?;

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
