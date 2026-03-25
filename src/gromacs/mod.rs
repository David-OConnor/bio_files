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
//! - Creates a temporary directory `gromacs_temp/`
//! - Writes `conf.gro`, `topol.top`, `md.mdp`
//! - Runs `gmx grompp` to produce `topol.tpr`
//! - Runs `gmx mdrun` to produce `traj.trr`, `confout.gro`, `energy.edr`, `md.log`
//! - Converts `traj.trr` to a multi-frame `traj.gro` via `gmx trjconv`
//! - Parses the trajectory and log, then removes the temporary directory
//!
//! Units throughout follow the rest of the codebase (Å for positions).
//! Conversions to GROMACS units (nm) happen inside this module.

pub mod mdp;
pub mod output;
pub mod topology;

use std::{
    fmt::Write as _,
    fs,
    fs::File,
    io::{self, ErrorKind, Write},
    path::Path,
    process::{Command, Stdio},
};

pub use mdp::MdpParams;
pub use output::{GromacsFrame, GromacsOutput};
pub use topology::MoleculeTopology;

/// Water model to use for explicit solvation via `gmx solvate`.
///
/// Selecting a model causes the pipeline to:
/// 1. Include the model's atom-type and molecule-type sections in the topology.
/// 2. Run `gmx solvate` to fill the simulation box with water before `grompp`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WaterModel {
    /// OPC 4-point rigid water (recommended with Amber ff19SB).
    /// Requires GROMACS 2022+ which ships `opc.gro` in its data directory.
    Opc,
}

impl WaterModel {
    /// Pre-equilibrated water-box filename for `gmx solvate -cs`.
    /// GROMACS resolves this from its data directory automatically.
    fn gro_filename(self) -> &'static str {
        match self {
            Self::Opc => "opc.gro",
        }
    }
}

use crate::{AtomGeneric, BondGeneric, md_params::ForceFieldParams};

// Used for creating intermediate files
const TEMP_DIR: &str = "gromacs_temp";

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
    pub water_model: Option<WaterModel>,
}

impl Default for GromacsInput {
    fn default() -> Self {
        Self {
            mdp: MdpParams::default(),
            molecules: Vec::new(),
            box_nm: None,
            ff_global: None,
            water_model: None,
        }
    }
}

// ---------------------------------------------------------------------------
// File generation
// ---------------------------------------------------------------------------

impl GromacsInput {
    /// Generate the `.mdp` file contents.
    pub fn make_mdp(&self) -> String {
        self.mdp.make_mdp()
    }

    /// Generate a GROMACS `.gro` structure file from all molecule atoms.
    /// [Example format](https://manual.gromacs.org/2026.1/reference-manual/file-formats.html#gro)
    ///
    /// Positions are converted from Å (internal) to nm (GROMACS). This has similar data to that stored
    /// in various molecule formats (SDF, Mol2, mmCIF etc), but also includes velocities,
    pub fn make_gro(&self) -> String {
        let total_atoms: usize = self.molecules.iter().map(|m| m.atoms.len() * m.count).sum();

        let mut s = String::from("Generated by Molchanica\n");
        let _ = writeln!(s, "{}", total_atoms);

        let mut atom_serial = 1;
        let mut res_serial = 1;

        for mol in &self.molecules {
            for _copy in 0..mol.count {
                for atom in &mol.atoms {
                    // atom name: prefer type_in_res, then type_in_res_general, then element+index
                    let atom_name = atom
                        .type_in_res
                        .as_ref()
                        .map(|t| t.to_string())
                        .or_else(|| atom.type_in_res_general.clone())
                        .unwrap_or_else(|| format!("{}{}", atom.element.to_letter(), atom_serial));

                    let x_nm = atom.posit.x / 10.0;
                    let y_nm = atom.posit.y / 10.0;
                    let z_nm = atom.posit.z / 10.0;

                    // GRO fixed-width format:
                    // resid(5) resname(5) atom(5) serial(5) x(8.3) y(8.3) z(8.3)
                    let _ = writeln!(
                        s,
                        "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                        res_serial,
                        &mol.name[..mol.name.len().min(5)],
                        &atom_name[..atom_name.len().min(5)],
                        atom_serial % 100_000,
                        x_nm,
                        y_nm,
                        z_nm,
                    );
                    atom_serial += 1;
                }
                res_serial += 1;
            }
        }

        // Box vector line (nm)
        let (bx, by, bz) = self.box_nm.unwrap_or((0.0, 0.0, 0.0));
        let _ = writeln!(s, "{:>10.5}{:>10.5}{:>10.5}", bx, by, bz);

        s
    }

    /// Generate the GROMACS `.top` topology file.
    pub fn make_top(&self) -> String {
        let mol_tops: Vec<topology::MoleculeTopology<'_>> = self
            .molecules
            .iter()
            .map(|m| topology::MoleculeTopology {
                name: &m.name,
                atoms: &m.atoms,
                bonds: &m.bonds,
                ff_mol: m.ff_params.as_ref(),
                count: m.count,
            })
            .collect();

        topology::make_top(&mol_tops, self.ff_global.as_ref(), self.water_model)
    }

    /// Total number of solute (non-water) atoms across all molecule copies.
    pub fn solute_atom_count(&self) -> usize {
        self.molecules.iter().map(|m| m.atoms.len() * m.count).sum()
    }

    /// Write all input files (`conf.gro`, `topol.top`, `md.mdp`) to `dir`.
    /// These must be in place prior to running MD with `grompp` and `mdrun`.
    pub fn save_input_files(&self, dir: &Path) -> io::Result<()> {
        fs::create_dir_all(dir)?;

        write_text(dir.join("conf.gro"), &self.make_gro())?;
        write_text(dir.join("topol.top"), &self.make_top())?;
        write_text(dir.join("md.mdp"), &self.make_mdp())?;

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
        self.save_input_files(dir)?;

        let solute_atom_count = self.solute_atom_count();

        // Optional solvation: fill the box with water before preprocessing.
        let structure_gro = if let Some(wm) = self.water_model {
            run_gmx(
                dir,
                &[
                    "solvate",
                    "-cp",
                    "conf.gro",
                    "-cs",
                    wm.gro_filename(),
                    "-o",
                    "solvated.gro",
                    "-p",
                    "topol.top",
                ],
            )?;
            "solvated.gro"
        } else {
            "conf.gro"
        };

        // grompp: Preprocessor. Creates a [binary] TPR file from the generated input files.
        // [Data on TPR](https://manual.gromacs.org/2026.1/reference-manual/file-formats.html#tpr)
        // This format contains everything GROMACS needs to run MD.
        run_gmx(
            dir,
            &[
                "grompp",
                "-f",
                "md.mdp",
                "-c",
                structure_gro,
                "-p",
                "topol.top",
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
                "traj.trr",
                "-x",
                "traj.xtc",
                "-c",
                "confout.gro",
                "-e",
                "energy.edr",
                "-g",
                "md.log",
            ],
        )?;

        // --- gmx trjconv: export trajectory as multi-frame GRO ---
        // `-pbc mol` re-wraps molecules into the box; `-b 0` starts from t=0.
        // The echo provides the "System" group selection for the interactive prompt.
        run_gmx_stdin(
            dir,
            &[
                "trjconv",
                "-f",
                "traj.trr",
                "-s",
                "topol.tpr",
                "-o",
                "traj.gro",
                "-pbc",
                "mol",
            ],
            b"0\n", // select group 0 = "System"
        )?;

        let log_text = read_text(dir.join("md.log")).unwrap_or_default();
        let traj_gro = read_text(dir.join("traj.gro"))?;

        let result = GromacsOutput::new(log_text, &traj_gro, solute_atom_count)?;

        // Remove the folder containing temporary (e.g. input and output) files.

        // todo: Temp leaving in
        // fs::remove_dir_all(dir)?;

        Ok(result)
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn write_text(path: impl AsRef<Path>, text: &str) -> io::Result<()> {
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
