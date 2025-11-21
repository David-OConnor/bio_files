#![allow(non_camel_case_types)]

//! For interop with ORCA files. For example, `inp` and `out` files. Can construct inputs for ORCA,
//! run these if ORCA is available in the system PATH, and display and parse the output.
//!
//! This currently supports a limited subset of ORCA's functionality.
//!
//! We make heavy use of `Option` on input fields, allowing them to be None instead of setting
//! a default value. Fields marked as None  mean *use ORCA's defaults*; this is more maintainable
//! and simpler for the user to understand than maintaining a duplicate set of default values.
//! It also has implications for writing vs ommiting items in the input string.

pub mod basis_sets;
pub mod method;
pub mod scf;
pub mod solvation;

use std::{
    fs,
    fs::File,
    io::{self, ErrorKind, Write},
    path::Path,
    process::Command,
};

use basis_sets::BasisSet;
use method::{Method, MethodSection};
use scf::Scf;
use solvation::{Solvator, SolvatorImplicit};

use crate::AtomGeneric;

/// A helper. The &str and String use reflects how we use this in practie,
/// e.g. with &str literals vs format!().
fn make_inp_block(block_name: &str, contents: &[(&str, String)]) -> String {
    let mut r = format!("%{block_name}\n");

    for (k, v) in contents {
        r.push_str(&format!("    {} {}\n", k, v));
    }

    r.push_str("end");
    r
}

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/counterpoise.html
/// Table 2.52
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

/// Misc other keywords not including method and basis set.
/// Note that we ommit some that are part of other fields we have, as with solvents.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Keyword {
    /// Minimal Basis Iterative Stockholder. Can be used for force field parameterization.
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/charges.html
    Mbis,
    OptimizeGeometry,
    // OptimizeHydrogens,
    TightOptimization,
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
}

impl Keyword {
    /// Prefixed with an !, starts the .inp file.
    pub fn keyword(self) -> String {
        match self {
            Self::Mbis => "MBIS".to_string(),
            Self::OptimizeGeometry => "OPT".to_string(),
            // Self::OptimizeHydrogens => "OPT",
            Self::TightOptimization => "TIGHTOPT".to_string(),
            Self::Freq => "FREQ".to_string(),
            Self::NumericalGradient => "NUMGRAD".to_string(),
            Self::D4Dispersion => "D4".to_string(),
            Self::ConformerSearch => "GOAT".to_string(),
            Self::Gcp(option) => format!("GCP({}", option.keyword()),
            Self::UseSymmetry => "UseSymmetry".to_string(),
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

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/symmetry.html
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

        make_inp_block("sym", &contents)
    }
}

#[derive(Debug, Clone, Default)]
pub struct OrcaInput {
    /// For now, We keep `Method` separate from the optional `method_section` part;
    /// the former is required, and writes to the first line. We could, alternatively,
    /// have it write to the `method` key in the `%method` block we write to  in `method_section`.
    pub method: Method, // Sets the first line
    pub method_section: Option<MethodSection>, // Sets the %method section.
    pub basis_set: BasisSet,
    pub keywords: Vec<Keyword>,
    pub atoms: Vec<AtomGeneric>,
    /// todo: Ref [this list of input blocks from the docs](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html);
    pub solvator: Option<Solvator>,
    pub solvator_implicit: Option<SolvatorImplicit>,
    pub bond_localization: Option<BondLocalization>,
    pub scf: Option<Scf>,
    pub symmetry: Option<Symmetry>,
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
        result.push_str(&format!(
            "!{} {}",
            self.method.keyword(),
            self.basis_set.keyword()
        ));

        for kw in &self.keywords {
            result.push_str(&format!(" {}", kw.keyword()));
        }

        // --- Blocks ---
        if let Some(m) = &self.method_section {
            result.push('\n');
            result.push_str(&m.make_inp());
        }

        // todo: Generalization over these blocks A/R with a macro.
        if let Some(solvator) = &self.solvator {
            result.push('\n');
            result.push_str(&solvator.make_inp());
        }

        if let Some(solvator) = &self.solvator_implicit {
            result.push('\n');
            result.push_str(&solvator.make_inp());
        }

        if let Some(loc) = &self.bond_localization {
            result.push('\n');
            result.push_str(&make_inp_block("loc", &[("locmet", loc.method.keyword())]));
        }

        if let Some(scf) = &self.scf {
            result.push('\n');
            result.push_str(&scf.make_inp());
        }

        if let Some(sym) = &self.symmetry {
            result.push('\n');
            result.push_str(&sym.make_inp());
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
    pub fn run(&self) -> io::Result<String> {
        // pub fn run(&self) -> io::Result<OrcaOutput> {
        let dir = Path::new("orca_temp");
        fs::create_dir_all(dir)?;

        let file_name = "temp_orca_input.inp";
        let path = Path::new(file_name);
        self.save(path)?;

        let out = match Command::new("orca")
            .current_dir(dir)
            .args([file_name])
            .output()
        {
            Ok(out) => out,
            Err(e) if e.kind() == ErrorKind::NotFound => {
                // Orca binary not found on PATH
                fs::remove_file(path)?;
                return Err(io::Error::new(
                    ErrorKind::NotFound,
                    "`orca` executable not found in the system PATH",
                ));
            }
            Err(e) => return Err(e),
        };

        if !out.status.success() {
            let stderr_str = String::from_utf8_lossy(&out.stderr);
            fs::remove_file(path)?;
            return Err(io::Error::other(format!(
                "Problem reading out temporary orca file: {}",
                stderr_str
            )));
        }

        // Remove the entire temporary directory.
        fs::remove_dir(dir)?;
        // fs::remove_file(path)?;

        // Convert stdout bytes to String
        let result =
            String::from_utf8(out.stdout).map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;

        Ok(result)
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum TerminationStatus {
    Normal,
    Error, // todo a/r
}

#[derive(Debug, Clone)]
pub struct OrcaOutput {
    termination_status: TerminationStatus,
    // atoms: Vec<AtomGeneric>,
}

impl OrcaOutput {
    /// Create output by parsing Orca's stdout text.
    pub fn new(data: &str) -> Self {
        let mut termination_status: TerminationStatus = TerminationStatus::Error;
        if data.contains("****ORCA TERMINATED NORMALLY****") {
            termination_status = TerminationStatus::Error
        }

        Self { termination_status }
    }
}
