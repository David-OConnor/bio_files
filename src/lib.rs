//! The `generic` label in names in this module are to differentiate from more specific
//! ones used in e.g. Daedalus.

pub mod ab1;
pub mod map;
pub mod mol2;
pub mod sdf;

pub mod amber_params;
mod cif_sf;
pub mod dat;
pub mod frcmod;
mod mtz;

use std::str::FromStr;

pub use ab1::*;
use lin_alg::f64::Vec3;
pub use map::*;
pub use mol2::*;
use na_seq::{AminoAcid, Element};

#[derive(Clone, Debug, Default)]
pub struct AtomGeneric {
    pub serial_number: usize,
    pub posit: Vec3,
    pub element: Element,
    /// e.g. "CG1", "CA", "O", "C" etc.
    pub name: Option<String>,
    // pub role: Option<AtomRole>,
    // pub residue: Option<usize>,
    pub occupancy: Option<f32>,
    pub partial_charge: Option<f32>,
    /// E.g. "c6", "ca", "n3", "ha", "h0" etc, as seen in Mol2 files from AMBER.
    /// e.g.: "ha": hydrogen attached to an aromatic carbon.
    /// "ho": hydrogen on a hydroxyl oxygen
    /// "n3": sp³ nitrogen with three substituents
    /// "c6": 	sp² carbon in a pure six-membered aromatic ring (new in GAFF2; lets GAFF distinguish
    /// a benzene carbon from other aromatic caca carbons)
    pub force_field_atom_type: Option<String>,
}

#[derive(Clone, Debug)]
pub struct BondGeneric {
    pub bond_type: String, // todo: Enum
    pub atom_0: usize,
    pub atom_1: usize,
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
    pub serial_number: isize, // pdbtbx uses isize. Negative allowed?
    pub res_type: ResidueType,
    pub atoms: Vec<usize>, // Atom index
}

#[derive(Debug, Clone)]
pub struct Chain {
    pub id: String,
    // todo: Do we want both residues and atoms stored here? It's an overconstraint.
    pub residues: Vec<usize>,
    /// Indexes
    pub atoms: Vec<usize>,
    // todo: Perhaps vis would make more sense in a separate UI-related place.
    pub visible: bool,
}
