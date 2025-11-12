//! For interop with ORCA files. For example, `inp` and `out` files.
//!
//! todo: Consider an XYZ format module, to assist with ORA input?

use crate::AtomGeneric;

// todo: Rename A/R
/// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum SinglePointEnergyType {
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

#[derive(Debug, Clone)]
pub struct OrcaInput {
    entry_type: SinglePointEnergyType,
    atoms: Vec<AtomGeneric>,
    opt: bool,
}

impl OrcaInput {
    pub fn new(entry_type: SinglePointEnergyType, atoms: &[AtomGeneric]) -> Self {
        Self {
            entry_type,
            atoms: atoms.to_vec(),
            opt: false,
        }
    }

    pub fn make_file(&self) -> String {
        let mut result = String::new();

        result.push_str(&format!("!{} DEF2-SVP", self.entry_type.command()));

        if self.opt {
            result.push_str(" OPT")
        }

        result.push_str("\n* xyz 0 1\n");

        // todo: Look at other files to see how you do this. E.g. write! on file objects, < > in formatting etc.
        for atom in &self.atoms {
            result.push_str(&format!(
                "{}    {:>5.4}    {:>5.4}    {:>10.4}  \n",
                atom.element.to_letter(),
                atom.posit.x,
                atom.posit.y,
                atom.posit.z
            ));
        }

        result.push_str("*");

        result
    }
}

#[derive(Debug, Clone)]
pub struct OrcaOutput {
    // atoms: Vec<AtomGeneric>,
}
