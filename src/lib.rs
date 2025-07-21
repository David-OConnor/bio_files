#![allow(confusable_idents)]
#![allow(mixed_script_confusables)]

//! The `generic` label in names in this module are to differentiate from ones used in more specific
//! applications.

pub mod ab1;
pub mod map;
pub mod mol2;
pub mod sdf;

pub mod amber_params;
mod cif_sf;
pub mod dat;
pub mod frcmod;
mod mmcif;
mod mmcif_aux;
mod mtz;

use std::{
    fmt,
    fmt::{Display, Formatter},
    io,
    io::ErrorKind,
    str::FromStr,
};

pub use ab1::*;
use lin_alg::f64::Vec3;
pub use map::*;
pub use mmcif::*;
pub use mol2::*;
use na_seq::{AminoAcid, AtomTypeInRes, Element};
pub use sdf::*;

#[derive(Clone, Debug, Default)]
pub struct AtomGeneric {
    pub serial_number: u32,
    pub posit: Vec3,
    pub element: Element,
    /// e.g. "CG1", "CA", "O", "C", "HA", "CD", "C9" etc.
    pub type_in_res: Option<AtomTypeInRes>,
    /// E.g. "c6", "ca", "n3", "ha", "h0" etc, as seen in Mol2 files from AMBER.
    /// e.g.: "ha": hydrogen attached to an aromatic carbon.
    /// "ho": hydrogen on a hydroxyl oxygen
    /// "n3": sp³ nitrogen with three substituents
    /// "c6": sp² carbon in a pure six-membered aromatic ring (new in GAFF2; lets GAFF distinguish
    /// a benzene carbon from other aromatic caca carbons)
    /// For proteins, this appears to be teh same as for `name`.
    /// Internal term: "Type 3".
    pub force_field_type: Option<String>,
    pub occupancy: Option<f32>,
    pub partial_charge: Option<f32>,
    pub hetero: bool,
}

#[derive(Clone, Debug)]
pub struct BondGeneric {
    pub bond_type: String, // todo: Enum
    pub atom_0_sn: u32,
    pub atom_1_sn: u32,
}

#[derive(Debug, Clone)]
pub enum ResidueType {
    AminoAcid(AminoAcid),
    Water,
    Other(String),
}

impl Default for ResidueType {
    fn default() -> Self {
        Self::Other(String::new())
    }
}

impl ResidueType {
    /// Parses from the "name" field in common text-based formats lik CIF, PDB, and PDBQT.
    pub fn from_str(name: &str) -> Self {
        if name.to_uppercase() == "HOH" {
            ResidueType::Water
        } else {
            match AminoAcid::from_str(name) {
                Ok(aa) => ResidueType::AminoAcid(aa),
                Err(_) => ResidueType::Other(name.to_owned()),
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct ResidueGeneric {
    /// We use serial number of display, search etc, and array index to select. Residue serial number is not
    /// unique in the molecule; only in the chain.
    pub serial_number: u32,
    pub res_type: ResidueType,
    /// Serial number
    pub atom_sns: Vec<u32>,
}

#[derive(Debug, Clone)]
pub struct ChainGeneric {
    pub id: String,
    // todo: Do we want both residues and atoms stored here? It's an overconstraint.
    /// Serial number
    pub residue_sns: Vec<u32>,
    /// Serial number
    pub atom_sns: Vec<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    Coil,
}

#[derive(Clone, Debug)]
/// See note elsewhere regarding serial numbers vs indices: In your downstream applications, you may
/// wish to convert sns to indices, for faster operations.
pub struct BackboneSS {
    /// Atom serial numbers.
    pub start_sn: u32,
    pub end_sn: u32,
    pub sec_struct: SecondaryStructure,
}

#[derive(Clone, Copy, PartialEq, Debug)]
/// The method used to find a given molecular structure. This data is present in mmCIF files
/// as the `_exptl.method` field.
pub enum ExperimentalMethod {
    XRayDiffraction,
    ElectronDiffraction,
    NeutronDiffraction,
    /// i.e. Cryo-EM
    ElectronMicroscopy,
    SolutionNmr,
}

impl ExperimentalMethod {
    /// E.g. for displaying in the space-constrained UI.
    pub fn to_str_short(&self) -> String {
        match self {
            Self::XRayDiffraction => "X-ray",
            Self::NeutronDiffraction => "ND",
            Self::ElectronDiffraction => "ED",
            Self::ElectronMicroscopy => "EM",
            Self::SolutionNmr => "NMR",
        }
        .to_owned()
    }
}

impl Display for ExperimentalMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let val = match self {
            Self::XRayDiffraction => "X-Ray diffraction",
            Self::NeutronDiffraction => "Neutron diffraction",
            Self::ElectronDiffraction => "Electron diffraction",
            Self::ElectronMicroscopy => "Electron microscopy",
            Self::SolutionNmr => "Solution NMR",
        };
        write!(f, "{val}")
    }
}

impl FromStr for ExperimentalMethod {
    type Err = io::Error;

    /// Parse an mmCIF‐style method string into an ExperimentalMethod.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let normalized = s.to_lowercase();
        let s = normalized.trim();
        let method = match s {
            "x-ray diffraction" => ExperimentalMethod::XRayDiffraction,
            "neutron diffraction" => ExperimentalMethod::NeutronDiffraction,
            "electron diffraction" => ExperimentalMethod::ElectronDiffraction,
            "electron microscopy" => ExperimentalMethod::ElectronMicroscopy,
            "solution nmr" => ExperimentalMethod::SolutionNmr,
            other => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Error parsing experimental method: {other}"),
                ));
            }
        };
        Ok(method)
    }
}
