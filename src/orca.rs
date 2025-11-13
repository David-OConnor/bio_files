//! For interop with ORCA files. For example, `inp` and `out` files.
//!
//! todo: Consider an XYZ format module, to assist with ORA input?

use std::{
    fs,
    fs::File,
    io::{self, ErrorKind, Write},
    path::Path,
    process::Command,
};

use crate::AtomGeneric;

// todo: Rename A/R
/// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum SinglePointEnergyType {
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
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#coupled-cluster-cc
    CoupledCluster,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#semiempirical-methods-sqm
    SemiEmpericalMethods,
}

impl SinglePointEnergyType {
    /// Prefixed with an !, starts the .inp file.
    pub fn command(self) -> String {
        match self {
            Self::HartreeFock => "HF",
            Self::Dft => "B3LYP",
            Self::Mp2Perturbation => "RI-MP2", // Or "DLPNO-MP2"
            Self::SpinComponentScaledMp2 => "RI-SCS-MP2",
            Self::OrbitalOptimzedMp2 => "OO-RI-MP2",
            Self::RegularlizedMp2 => "RI-MP2",
            Self::DoubleHybridDft => "B2PLYP",
            Self::CoupledCluster => "DLPNO-CCSD(T)",
            Self::SemiEmpericalMethods => "XTB2",
        }
        .to_string()
    }
}

#[derive(Debug, Clone, Default)]
pub struct OrcaInput {
    pub energy_type: SinglePointEnergyType,
    pub atoms: Vec<AtomGeneric>,
    /// Optimize geometry
    pub opt: bool,
    /// Calculate vibrational frequencies
    pub freq: bool,
    /// Request a numerical gradient
    pub numerical_gradient: bool,
    /// Optimize positions of hydrogen atoms only
    pub optimize_hydrogens: bool,
}

impl OrcaInput {
    pub fn new(energy_type: SinglePointEnergyType, atoms: &[AtomGeneric]) -> Self {
        Self {
            energy_type,
            atoms: atoms.to_vec(),
            opt: false,
            freq: false,
            numerical_gradient: false,
            optimize_hydrogens: false,
        }
    }

    /// Create an .inp string for input into ORCA.
    pub fn make_inp(&self) -> String {
        let mut result = String::new();

        result.push_str(&format!("!{} DEF2-SVP", self.energy_type.command()));

        if self.opt {
            result.push_str(" OPT")
        }
        if self.freq {
            result.push_str(" FREQ")
        }
        if self.numerical_gradient {
            result.push_str(" NUMGRAD")
        }

        if self.optimize_hydrogens {
            result.push_str("\n\n%geom\n  optimizehydrogens true\nend\n\n")
        }

        result.push_str("\n* xyz 0 1\n");

        for atom in &self.atoms {
            result.push_str(&format!(
                "{:<2} {:>12.5} {:>12.5} {:>12.5}\n",
                atom.element.to_letter(),
                atom.posit.x,
                atom.posit.y,
                atom.posit.z
            ));
        }

        result.push_str("*");

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
        let file_name = "temp_orca_input.inp";
        let path = Path::new(file_name);
        self.save(&path)?;

        let out = match Command::new("orca").args([file_name]).output() {
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

        fs::remove_file(path)?;

        // Convert stdout bytes to String
        let result =
            String::from_utf8(out.stdout).map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;

        Ok(result)
    }
}

#[derive(Debug, Clone)]
pub struct OrcaOutput {
    // atoms: Vec<AtomGeneric>,
}
