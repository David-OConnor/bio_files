#![allow(confusable_idents)]
#![allow(mixed_script_confusables)]

//! The `generic` label in names in this module are to differentiate from ones used in more specific
//! applications.

pub mod ab1;
pub mod map;
pub mod mol2;
pub mod pdbqt;
pub mod sdf;

pub mod amber_params;
pub mod dat;
pub mod frcmod;
mod mmcif;
mod mmcif_aux;

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
pub use pdbqt::Pdbqt;
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
    /// "n3": sp³ nitrogen with three substitutes
    /// "c6": sp² carbon in a pure six-membered aromatic ring (new in GAFF2; lets GAFF distinguish
    /// a benzene carbon from other aromatic caca carbons)
    /// For proteins, this appears to be the same as for `name`.
    pub force_field_type: Option<String>,
    /// An atom-centered electric charge, used in molecular dynamics simulations.
    /// These are sometimes loaded from Mol2 or SDF files, and sometimes added after.
    pub partial_charge: Option<f32>,
    /// Indicates, in proteins, that the atom isn't part of an amino acid. E.g., water or
    /// ligands.
    pub hetero: bool,
    pub occupancy: Option<f32>,
}

impl Display for AtomGeneric {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let ff_type = match &self.force_field_type {
            Some(f) => f,
            None => "None",
        };

        let q = match &self.partial_charge {
            Some(q_) => format!("{q_:.3}"),
            None => "None".to_string(),
        };

        write!(
            f,
            "Atom {}: {}, {}. {:?}, ff: {ff_type}, q: {q}",
            self.serial_number,
            self.element.to_letter(),
            self.posit,
            self.type_in_res,
        )?;

        if self.hetero {
            write!(f, ", Het")?;
        }

        Ok(())
    }
}

/// These are the Mol2 standard types, unless otherwise noted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
    Amide,
    Dummy,
    Unknown,
    NotConnected,
    /// mmCIF, rare
    Quadruple,
    /// mmCIF. Distinct from aromatic; doesn't need to be a classic ring.
    Delocalized,
    /// mmCif; mostly for macromolecular components
    PolymericLink,
}

impl Display for BondType {
    /// Return the exact MOL2 bond-type token as an owned `String`.
    /// (Use `&'static str` if you never need it allocated.)
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Single => "1",
            Self::Double => "2",
            Self::Triple => "3",
            Self::Aromatic => "ar",
            Self::Amide => "am",
            Self::Dummy => "du",
            Self::Unknown => "un",
            Self::NotConnected => "nc",
            Self::Quadruple => "quad",
            Self::Delocalized => "delo",
            Self::PolymericLink => "poly",
        };

        write!(f, "{name}")
    }
}

impl BondType {
    /// SDF format uses a truncated set, and does things like mark every other
    /// aromatic bond as double.
    pub fn to_str_sdf(&self) -> String {
        match self {
            Self::Single | Self::Double | Self::Triple => *self,
            _ => Self::Single,
        }
        .to_string()
    }
}

impl FromStr for BondType {
    type Err = io::Error;

    /// Can ingest from mol2, SDF, and mmCIF formats.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "1" | "sing" => Ok(BondType::Single),
            "2" | "doub" => Ok(BondType::Double),
            "3" | "trip" => Ok(BondType::Triple),
            "4" | "ar" | "arom" => Ok(BondType::Aromatic),
            "am" => Ok(BondType::Amide),
            "du" => Ok(BondType::Dummy),
            "un" => Ok(BondType::Unknown),
            "nc" => Ok(BondType::NotConnected),
            "quad" => Ok(BondType::Quadruple),
            "delo" => Ok(BondType::Delocalized),
            "poly" => Ok(BondType::PolymericLink),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid BondType: {s}"),
            )),
        }
    }
}

#[derive(Clone, Debug)]
pub struct BondGeneric {
    pub bond_type: BondType,
    pub atom_0_sn: u32,
    pub atom_1_sn: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ResidueType {
    AminoAcid(AminoAcid),
    Water,
    Other(String),
}

impl Display for ResidueType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let name = match &self {
            ResidueType::Other(n) => n.clone(),
            ResidueType::Water => "Water".to_string(),
            ResidueType::AminoAcid(aa) => aa.to_string(),
        };

        write!(f, "{name}")
    }
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
    pub end: ResidueEnd,
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ResidueEnd {
    Internal,
    NTerminus,
    CTerminus,
    /// Not part of a protein/polypeptide.
    Hetero,
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
