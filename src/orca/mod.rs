#![allow(non_camel_case_types)]

//! For interop with ORCA files. For example, `inp` and `out` files. Can construct inputs for ORCA,
//! run these if ORCA is available in the system PATH, and display and parse the output.
//!
//! This currently supports a limited subset of ORCA's functionality.
//!
//! Fields marked as None generally mean *use ORCA's defaults*.

pub mod basis_sets;
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
use scf::Scf;
use solvation::{Solvator, SolvatorImplicit};

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

use crate::AtomGeneric;
/// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum Method {
    #[default]
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#hartree-fock-hf
    HartreeFock,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#density-functional-theory-dft
    /// todo: Multiple DFT functionals.
    Dft,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#mp2-perturbation-theory
    Mp2Perturbation,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#spin-component-scaled-mp2-scs-mp2
    SpinComponentScaledMp2,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#orbital-optimized-mp2-oo-mp2
    OrbitalOptimzedMp2,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#regularized-mp2
    RegularlizedMp2,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#double-hybrid-dft-dh-dft
    DoubleHybridDft,
    TripleHybridDft, // todo: QC this one.
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#coupled-cluster-cc
    CoupledCluster,
    Xtb,
    // todo: Support other semiemperical methods; XTB2 is just one.
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#semiempirical-methods-sqm
    SemiEmpericalMethods,

    // --- Gradient corrected functions (GGA) -----
    // [Docs, Table 3.2](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#gradient-corrected-functionals-meta-ggas)
    BP86,
    BLYP,
    OLYP,
    GLYP,
    XLYP,
    PW91,
    MPWPW,
    MPWLYP,
    PBE,
    RPBE,
    REVPBE,
    RPW86PBE,
    PWP,
    B97_3C,
    B97M_V,
    V97M_D3BJ,
    B97M_D4,
    SCANFUNC,
    RSCAN,
    R2SCAN,
    TPSS,
    REVTPSS,
    R2SCAN_3C,
    // --- End GGAs
    // --- Global hybrid functions. [Docs Table 3.4](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#global-hybrid-functionals)
    B1LYP,
    B3LYP,
    B3LYP_G,
    O3LYP,
    X3LYP,
    B1P86,
    B3PW91,
    PW1PW,
    MPW1PW,
    MPW1LYP,
    PBE0,
    REVPBE0,
    REVPBE38,
    BHANDHLYP,
    M06,
    M062X,
    PW6B95,
    TPSSH,
    TPSS0,
    R2SCANH,
    R2SCAN0,
    R2SCAN50,
    PBEH_3C,
    B3LYP_3C,
    // --- End HFX
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/fod.html
    FractionalOccupationDensity,
    None, // todo: I believe this is right?
}

impl Method {
    /// Prefixed with an !, starts the .inp file.
    pub fn keyword(self) -> String {
        match self {
            Self::HartreeFock => "HF",
            Self::Dft => "B3LYP",
            Self::Mp2Perturbation => "RI-MP2", // Or "DLPNO-MP2"
            Self::SpinComponentScaledMp2 => "RI-SCS-MP2",
            Self::OrbitalOptimzedMp2 => "OO-RI-MP2",
            Self::RegularlizedMp2 => "RI-MP2",
            Self::DoubleHybridDft => "B2PLYP",
            Self::TripleHybridDft => "B3LYP-D3(BJ)", // todo: QC this one
            Self::CoupledCluster => "DLPNO-CCSD(T)",
            Self::Xtb => "XTB", // todo?
            Self::SemiEmpericalMethods => "XTB2",
            // GGA
            Self::BP86 => "BP86",
            Self::BLYP => "BLYP",
            Self::OLYP => "OLYP",
            Self::GLYP => "GLYP",
            Self::XLYP => "XLYP",
            Self::PW91 => "PW91",
            Self::MPWPW => "MPWPW",
            Self::MPWLYP => "MPWLYP",
            Self::PBE => "PBE",
            Self::RPBE => "RPBE",
            Self::REVPBE => "REVPBE",
            Self::RPW86PBE => "RPW86PBE",
            Self::PWP => "PWP",
            Self::B97_3C => "B97-3C",
            Self::B97M_V => "B97M-V",
            Self::V97M_D3BJ => "V97M-D3BJ",
            Self::B97M_D4 => "B97M-D4",
            Self::SCANFUNC => "SCANFUNC",
            Self::RSCAN => "RSCAN",
            Self::R2SCAN => "R2SCAN",
            Self::TPSS => "TPSS",
            Self::REVTPSS => "REVTPSS",
            Self::R2SCAN_3C => "R2SCAN-3C",
            // end GGA
            // Start Global Hybrid functionals
            Self::B1LYP => "B1LYP",
            Self::B3LYP => "B3LYP",
            Self::B3LYP_G => "B3LYP-G",
            Self::O3LYP => "O3LYP",
            Self::X3LYP => "X3LYP",
            Self::B1P86 => "B1P86",
            Self::B3PW91 => "B3PW91",
            Self::PW1PW => "PW1PW",
            Self::MPW1PW => "MPW1PW",
            Self::MPW1LYP => "MPW1LYP",
            Self::PBE0 => "PBE0",
            Self::REVPBE0 => "REVPBE0",
            Self::REVPBE38 => "REVPBE38",
            Self::BHANDHLYP => "BHANDHLYP",
            Self::M06 => "M06",
            Self::M062X => "M062X",
            Self::PW6B95 => "PW6B95",
            Self::TPSSH => "TPSSH",
            Self::TPSS0 => "TPSS0",
            Self::R2SCANH => "R2SCANH",
            Self::R2SCAN0 => "R2SCAN0",
            Self::R2SCAN50 => "R2SCAN50",
            Self::PBEH_3C => "PBEH-3C",
            Self::B3LYP_3C => "B3LYP-3C",
            // End Global hybrid functionals
            Self::FractionalOccupationDensity => "FOD",
            Self::None => "",
        }
        .to_string()
    }
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
    pub method: Method,
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
        //
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
