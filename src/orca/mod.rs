#![allow(non_camel_case_types)]

//! For interop with the [ORCA quantum chemistry software package](https://www.faccts.de/orca/). For example,
//! This module constructs inputs for ORCA, and can either save them to disk, or run them directly
//!  if ORCA is available in the system PATH. It can display and parse the output.
//!
//! The workflow when running commands directly:
//!   - Create a temporary directly inside the working directory
//!   - Save a `.inp` script to that dir
//!   - Launch that script. Note that it will add temporary and output files to that dir
//!   - Collect the output via `stdout` and `stderr`
//!   - Delete the temporary directory
//!
//! This currently supports a limited subset of ORCA's functionality.
//!
//! We make use of `Option` on input fields, allowing them to be `None` instead of setting
//! a default value. Fields marked as None  mean *use ORCA's defaults*; this is more maintainable
//! and simpler for the user to understand than maintaining a duplicate set of default values.
//! It also has implications for writing vs omiting items in the input string.
//!
//! Most structs and enums, and many fields, link to the relevant section of the [ORCA Manual](https://www.faccts.de/docs/orca/6.1/manual),
//! in lieu of direct documentation.
//!
//! While the `bio-files` library in general works in Python via bindings, we have not enabled the Orca
//! module in Python, because FACCTS provides its own [high-quality ORCA Python Interface library](https://www.faccts.de/docs/opi/1.0/docs/)

pub mod basis_sets;
pub mod charges;
pub mod dynamics;
pub mod geom;
pub mod method;
mod plots;
pub mod scf;
pub mod solvation;

use std::{
    fmt::Display,
    fs,
    fs::File,
    io::{self, ErrorKind, Write},
    path::Path,
    process::Command,
};

use basis_sets::BasisSet;
use lin_alg::f64::Vec3;
use method::{Method, MethodSection};
use scf::Scf;
use solvation::{Solvator, SolvatorImplicit};

use crate::{
    AtomGeneric,
    orca::{
        charges::{AtomChargeData, ChargesOutput},
        dynamics::{Dynamics, DynamicsOutput},
        geom::Geom,
        plots::Plots,
    },
};

/// A helper. The &str and String use reflects how we use this in practie,
/// e.g. with &str literals vs format!().
fn make_inp_block(block_name: &str, contents: &[(&str, String)], keywords: &[&str]) -> String {
    let mut r = String::new();

    r.push('%');
    r.push_str(block_name);
    for kw in keywords {
        r.push(' ');
        r.push_str(kw);
    }
    r.push('\n');

    for (k, v) in contents {
        r.push_str(&format!("    {} {}\n", k, v));
    }

    r.push_str("end");
    r
}

/// [Counterpoise Corrections](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/counterpoise.html?q=dftminis&n=0#counterpoise-corrections)
/// See Table 2.52
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum GcpOption {
    HfMinis,
    HfSv,
    Hf631Gd,
    HfSvp,
    HfTz,
    DftMinis,
    DftSv,
    Dft631Gd,
    DftLanl,
    DftVsP_,
    DftSvp,
    DftTz,
    File,
}

impl GcpOption {
    pub fn keyword(self) -> String {
        match self {
            Self::HfMinis => "hf/minis",
            Self::HfSv => "hf/sv",
            Self::Hf631Gd => "hf/631gd",
            Self::HfSvp => "hf/svp",
            Self::HfTz => "hf/tz",
            Self::DftMinis => "dft/minis",
            Self::DftSv => "dft/sv",
            Self::Dft631Gd => "dft/631gd",
            Self::DftLanl => "dft/lanl",
            Self::DftVsP_ => "dft/sv(p)",
            Self::DftSvp => "dft/svp",
            Self::DftTz => "dft/tz",
            Self::File => "file",
        }
        .to_string()
    }
}

/// [Geometry optimization thresholds](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html?q=tightopt&n=0#geometry-optimization-thresholds)
/// Just the general ones, for their use as keywords. Only include this if intending to perform
/// geometry optimization, for example, instead of calculating single-point energies.
#[derive(Clone, Copy, Default, PartialEq, Debug)]
pub enum GeomOptThresh {
    Loose,
    /// Default in Orca, e.g. if the keyword is omitted.
    #[default]
    Opt,
    Tight,
    VeryTight,
}

impl GeomOptThresh {
    pub fn keyword(self) -> String {
        match self {
            Self::Loose => "LooseOpt".to_string(),
            Self::Opt => "Opt".to_string(),
            Self::Tight => "TightOpt".to_string(),
            Self::VeryTight => "VeryTightOpt".to_string(),
        }
    }
}

impl Display for GeomOptThresh {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let v = match self {
            Self::Loose => "Loose",
            Self::Opt => "Opt",
            Self::Tight => "Tight",
            Self::VeryTight => "Very tight",
        };

        write!(f, "{}", v)
    }
}

/// This seems to be implicit, but of critical importance: This controls the primary mode of operation.
/// the nature of the computation's results depend on this.
#[derive(Clone, Debug, Default)]
pub enum Task {
    /// The default ORCA mode of operation; no associated keyword.
    #[default]
    SinglePoint,
    GeometryOptimization((GeomOptThresh, Option<Geom>)),
    /// Minimal Basis Iterative Stockholder. Can be used for force field parameterization.
    /// [Charge tutorial](https://www.faccts.de/docs/orca/5.0/tutorials/prop/charges.html)
    /// [MBIS Charges](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html?q=mbis&n=0#mbis-charges)
    MbisCharges,
    MolDynamics(Dynamics),
    // todo: Others A/R
}

impl Display for Task {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let v = match self {
            Self::SinglePoint => "Single point energy",
            Self::GeometryOptimization(_) => "Optimize geometry",
            Self::MbisCharges => "MBIS charges",
            Self::MolDynamics(_) => "Mol dynamics (Ab-initio)",
        };

        write!(f, "{v}")
    }
}

/// Misc other keywords not including method and basis set.
/// Note that we ommit some that are part of other fields we have, as with solvents.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Keyword {
    // Mbis,
    Freq,
    NumericalGradient,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/disp.html
    D4Dispersion,
    /// AKA GOAT. https://www.faccts.de/docs/orca/6.0/tutorials/prop/goat.html
    ConformerSearch,
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/counterpoise.html
    Gcp(GcpOption),
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/symmetry.html
    UseSymmetry,
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html
    AnFreq,
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html
    NumFreq,
}

impl Keyword {
    /// Prefixed with an !, starts the .inp file.
    pub fn keyword(self) -> String {
        match self {
            // Self::Mbis => "MBIS".to_string(),
            Self::Freq => "FREQ".to_string(),
            Self::NumericalGradient => "NUMGRAD".to_string(),
            Self::D4Dispersion => "D4".to_string(),
            Self::ConformerSearch => "GOAT".to_string(),
            Self::Gcp(option) => format!("GCP({}", option.keyword()),
            Self::UseSymmetry => "UseSymmetry".to_string(),
            Self::AnFreq => "AnFreq".to_string(),
            Self::NumFreq => "NumFreq".to_string(),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum LocalizationMethod {
    #[default]
    PipekMezey,
    FosterBoys,
}

impl LocalizationMethod {
    pub fn keyword(self) -> String {
        match self {
            Self::PipekMezey => "PM",
            Self::FosterBoys => "FB",
        }
        .to_owned()
    }
}

#[derive(Clone, Debug)]
pub struct BondLocalization {
    pub method: LocalizationMethod,
}

/// [ORCA and Symmetry](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/symmetry.html)
#[derive(Clone, Debug, Default)]
pub struct Symmetry {
    pub sym_thresh: Option<f32>,
    pub prefer_c2v: Option<bool>,
    pub point_group: Option<String>,
    // todo: Fill out the remaining items using the table A/R.
}

impl Symmetry {
    pub fn make_inp(&self) -> String {
        let mut contents = vec![("UseSymmetry", "true".to_owned())];

        if let Some(v) = self.sym_thresh {
            contents.push(("SymThresh", format!("{v:.6}")));
        }

        if let Some(v) = self.prefer_c2v {
            contents.push(("PreferC2v", format!("{v:?}")));
        }

        if let Some(v) = &self.point_group {
            contents.push(("PreferC2v", v.clone()));
        }

        make_inp_block("sym", &contents, &[])
    }
}

/// [General Structure of the Input File](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html)
/// Any fields marked as `Optional here`
#[derive(Debug, Clone, Default)]
pub struct OrcaInput {
    pub task: Task,
    // For now, We keep `Method` separate from the optional `method_section` part;
    // the former is required, and writes to the first line. We could, alternatively,
    // have it write to the `method` key in the `%method` block we write to  in `method_section`.
    pub method: Method,
    pub method_section: Option<MethodSection>, // Sets the %method section.
    pub basis_set: BasisSet,
    // /// If None, calculate single point energies as the default mode.
    // pub opt_mode: Option<GeomOptThresh>,
    pub keywords: Vec<Keyword>,
    pub atoms: Vec<AtomGeneric>,
    /// todo: Ref [this list of input blocks from the docs](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html);
    pub solvator: Option<Solvator>,
    pub solvator_implicit: Option<SolvatorImplicit>,
    pub bond_localization: Option<BondLocalization>,
    pub scf: Option<Scf>,
    pub symmetry: Option<Symmetry>, // todo: Combine into task?
    // pub dynamics: Option<Dynamics>,
    pub plots: Option<Plots>,
    // todo: A/R: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/stabilityanalysis.html
    // pub shark: Option<Shark>,
}

impl OrcaInput {
    // todo: Keywords?
    pub fn new(method: Method, basis_set: BasisSet, atoms: &[AtomGeneric]) -> Self {
        Self {
            method,
            basis_set,
            atoms: atoms.to_vec(),
            ..Self::default()
        }
    }

    /// Create an .inp string for input into ORCA.
    pub fn make_inp(&self) -> String {
        let mut result = String::new();

        // --- Initial line ---
        result.push_str(&format!("!{}", self.method.keyword()));

        if !self.method.is_composite() {
            result.push_str(&format!(" {}", self.basis_set.keyword()));
        }

        match &self.task {
            Task::SinglePoint => {} // No action or keyword.
            Task::GeometryOptimization((thresh, geom)) => {
                // todo: Integrate our geom block here?
                result.push_str(&format!(" {}", thresh.keyword()));

                if let Some(v) = geom {
                    result.push('\n');
                    result.push_str(&v.make_inp());
                }
            }
            Task::MbisCharges => {
                result.push_str(" MBIS");
            }
            Task::MolDynamics(md) => {
                result.push_str(&format!(" {}", md.make_inp()));
            }
        }

        for kw in &self.keywords {
            if kw == &Keyword::D4Dispersion && self.method.is_composite() {
                continue;
            }

            result.push_str(&format!(" {}", kw.keyword()));
        }

        // --- Blocks ---
        if let Some(v) = &self.method_section {
            result.push('\n');
            result.push_str(&v.make_inp());
        }

        // todo: Generalization over these blocks A/R with a macro.
        if let Some(v) = &self.solvator {
            result.push('\n');
            result.push_str(&v.make_inp());
        }

        if let Some(v) = &self.solvator_implicit {
            result.push('\n');
            result.push_str(&v.make_inp());
        }

        if let Some(v) = &self.bond_localization {
            result.push('\n');
            result.push_str(&make_inp_block(
                "loc",
                &[("locmet", v.method.keyword())],
                &[],
            ));
        }

        if let Some(v) = &self.scf {
            result.push('\n');
            result.push_str(&v.make_inp());
        }

        if let Some(v) = &self.symmetry {
            result.push('\n');
            result.push_str(&v.make_inp());
        }

        if let Some(v) = &self.plots {
            result.push('\n');
            result.push_str(&v.make_inp());
        }

        result.push_str("\n\n* xyz 0 1\n");

        // --- Atoms ---
        for atom in &self.atoms {
            result.push_str(&format!(
                "{:<2} {:>12.5} {:>12.5} {:>12.5}\n",
                atom.element.to_letter(),
                atom.posit.x,
                atom.posit.y,
                atom.posit.z
            ));
        }

        result.push('*');

        result
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut file = File::create(path)?;
        let text = self.make_inp();

        write!(file, "{text}")
    }

    /// Run this command in Orca, and collect the output. Requires `orca` to be available
    /// on the system PATH environment variable.
    /// todo: Outputs a string for now; adjust this as required into a custom output struct
    pub fn run(&self) -> io::Result<OrcaOutput> {
        let dir = Path::new("orca_temp");
        fs::create_dir_all(dir)?;

        let file_name = "temp_orca_input.inp";
        let path = dir.join(Path::new(file_name));
        self.save(&path)?;

        let cmd_out = match Command::new("orca")
            .current_dir(dir)
            .args([file_name])
            .output()
        {
            Ok(out) => out,
            Err(e) if e.kind() == ErrorKind::NotFound => {
                // Orca binary not found on PATH
                fs::remove_dir_all(dir)?;

                return Err(io::Error::new(
                    ErrorKind::NotFound,
                    "`orca` executable not found in the system PATH",
                ));
            }
            Err(e) => return Err(e),
        };

        if !cmd_out.status.success() {
            let stderr_str = String::from_utf8_lossy(&cmd_out.stderr);
            fs::remove_dir_all(dir)?;

            return Err(io::Error::other(format!(
                "Problem reading out temporary ORCA file: {}",
                stderr_str
            )));
        }

        // Convert stdout bytes to String
        let result_text = String::from_utf8(cmd_out.stdout)
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;

        if !result_text.contains("****ORCA TERMINATED NORMALLY****") {
            return Err(io::Error::other(format!(
                "ORCA terminated abnormally: {result_text}"
            )));
        }

        let result = match &self.task {
            Task::SinglePoint => {
                // todo
                OrcaOutput::Text(result_text)
            }
            Task::MolDynamics(md) => {
                let out = dir.join(&md.traj_out_dir);
                let out = DynamicsOutput::new(&out, result_text)?;
                OrcaOutput::Dynamics(out)
            }
            Task::MbisCharges => {
                let out = ChargesOutput::new(result_text)?;
                OrcaOutput::Charges(out)
            }
            Task::GeometryOptimization(_) => {
                let out = GeometryOutput::new(result_text)?;
                OrcaOutput::Geometry(out)
            }
        };

        // Remove the entire temporary directory.
        fs::remove_dir_all(dir)?;

        Ok(result)
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum TerminationStatus {
    Normal,
    Error, // todo a/r
}

#[derive(Debug, Clone)]
pub enum OrcaOutput {
    Text(String),
    Dynamics(DynamicsOutput),
    Charges(ChargesOutput),
    /// E.g. from geometry optimization.
    Geometry(GeometryOutput),
    // termination_status: TerminationStatus,
}

// impl OrcaOutput {
//     /// Create output by parsing Orca's stdout text.
//     pub fn new(data: &str) -> Self {
//         let mut termination_status: TerminationStatus = TerminationStatus::Error;
//         if data.contains("****ORCA TERMINATED NORMALLY****") {
//             termination_status = TerminationStatus::Error
//         }
//
//         Self { termination_status }
//     }
// }

#[derive(Debug, Clone)]
pub struct GeometryOutput {
    pub text: String,
    pub posits: Vec<Vec3>,
}

impl GeometryOutput {
    pub fn new(text: String) -> io::Result<Self> {
        let mut posits = Vec::new();

        // 1. Find the start of the "FINAL ENERGY EVALUATION" section.
        // We look for this first to ensure we aren't grabbing initial or intermediate steps.
        let final_eval_marker = "*** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***";
        let section_start = text.find(final_eval_marker).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                "Final stationary point not reached yet",
            )
        })?;

        // 2. Only search within the text AFTER that marker
        let remaining_text = &text[section_start..];
        let coord_header = "CARTESIAN COORDINATES (ANGSTROEM)";

        let header_pos = remaining_text.find(coord_header).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                "Could not find coordinates in final section",
            )
        })?;

        // 3. Iterate through lines starting from the coordinate header
        let mut lines = remaining_text[header_pos..].lines();

        // Skip the header line and the "-----------------" separator line
        lines.next();
        lines.next();

        for line in lines {
            let trimmed = line.trim();

            // ORCA usually ends these blocks with a line of dashes or an empty line
            if trimmed.is_empty() || trimmed.starts_with('-') {
                break;
            }

            let parts: Vec<&str> = trimmed.split_whitespace().collect();

            // Format: Symbol  X  Y  Z
            if parts.len() >= 4 {
                let x = parts[1]
                    .parse::<f64>()
                    .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;
                let y = parts[2]
                    .parse::<f64>()
                    .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;
                let z = parts[3]
                    .parse::<f64>()
                    .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;

                posits.push(Vec3 { x, y, z });
            }
        }

        if posits.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "Coordinate block was empty or malformed",
            ));
        }

        Ok(Self { text, posits })
    }
}
