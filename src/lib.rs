#![allow(confusable_idents)]
#![allow(mixed_script_confusables)]

//! The `generic` label in names in this module are to differentiate from ones used in more specific
//! applications.

pub mod gro;
pub mod mmcif;
pub mod mol2;
pub mod pdbqt;
pub mod sdf;

pub mod ab1;
pub mod map;

pub mod dat;
pub mod frcmod;
pub mod md_params;
pub mod orca;

pub mod amber_typedef;
mod bond_inference;
pub mod cif_sf;
pub mod dcd;
mod mmcif_aux;
pub mod mol_templates;
pub mod prmtop;
pub mod xyz;

use std::{
    fmt,
    fmt::{Display, Formatter},
    io,
    io::ErrorKind,
    str::FromStr,
};

pub use ab1::*;
pub use bond_inference::create_bonds;
pub use gro::*;
use lin_alg::f64::Vec3;
pub use map::*;
pub use mmcif::*;
pub use mol2::*;
use na_seq::{AminoAcid, AtomTypeInRes, Element};
pub use pdbqt::Pdbqt;
pub use sdf::*;
pub use xyz::*;

// todo: SHould this be in na_seq?
/// Common lipid types, as defined in Amber params. For phospholipics, chains and head groups
/// are separate entries, and can be combined.
/// The repr assumes ingested from lipid21 in a deterministic (e.g. alphabetcial) order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(usize)]
pub enum LipidStandard {
    /// Arachidonoyl chain (20:4 ω-6).
    Ar,
    /// Cholesterol
    Chl,
    /// Docosahexaenoyl chain (22:6 ω-3).
    Dha,
    /// Linoleoyl chain (18:2 ω-6). (Amber tag: LAL)
    Lal,
    /// Myristoyl chain (14:0).
    My,
    /// Oleoyl chain (18:1 ω-9).
    Ol,
    /// Palmitoyl chain (16:0).
    Pa,
    /// Phosphatidylcholine headgroup (PC).
    Pc,
    /// Phosphatidylethanolamine headgroup (PE).
    Pe,
    /// Phosphatidylglycerol headgroup (PG / PGR). Note: ChemForma currently uses "PG".
    Pgr,
    /// Phosphatidylglycerol sulfate / related PG variant (Amber tag: PGS).
    Pgs,
    /// Phosphate head (Amber tag: "PH-"; protonation/charge variant used in Amber lipids).
    Ph,
    /// Phosphatidylserine headgroup (PS).
    Ps,
    /// Stearoyl chain (18:0). (Amber tag: SA)
    Sa,
    /// Sphingomyelin (SPM).
    Spm,
    /// Stearoyl chain (18:0). (Amber tag: ST)
    St,
    /// Phosphatidylinositol headgroup (PI). *Not in lib21.dat.*
    Pi,
    /// Cardiolipin (diphosphatidylglycerol). *Not in lib21.dat.*
    Cardiolipin,
}

impl Display for LipidStandard {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Ar => "Ar",
            Self::Chl => "Chl",
            Self::Dha => "Dha",
            Self::Lal => "Lal",
            Self::My => "My",
            Self::Ol => "Ol",
            Self::Pa => "Pa",
            Self::Pc => "Pc",
            Self::Pe => "Pe",
            Self::Pgr => "Pgr",
            Self::Pgs => "Pgs",
            Self::Ph => "Ph", // corresponds to "PH-"
            Self::Ps => "Ps",
            Self::Sa => "Sa",
            Self::Spm => "Spm",
            Self::St => "St",
            Self::Pi => "Pi",                   // not in lib21.dat
            Self::Cardiolipin => "Cardiolipin", // not in lib21.dat
        };

        write!(f, "{name}")
    }
}

impl FromStr for LipidStandard {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_uppercase().as_ref() {
            // lib21.dat tags
            "AR" => Self::Ar,
            "CHL" => Self::Chl,
            "DHA" => Self::Dha,
            "LAL" => Self::Lal,
            "MY" => Self::My,
            "OL" => Self::Ol,
            "PA" => Self::Pa,
            "PC" => Self::Pc,
            "PE" => Self::Pe,
            "PGR" => Self::Pgr,
            "PGS" => Self::Pgs,
            "PH-" => Self::Ph,
            "PS" => Self::Ps,
            "SA" => Self::Sa,
            "SPM" => Self::Spm,
            "ST" => Self::St,
            // Common aliases / non-lib21 extras
            "PI" => Self::Pi, // not in lib21.dat
            "CARDIOLIPIN" | "CL" | "CDL" | "CDL2" => Self::Cardiolipin, // tood; Not sure.
            _ => {
                return Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("Unknown lipid standard: '{s}'"),
                ));
            }
        })
    }
}

// // todo: Move this to NA/Seq too likely (even more likely than LipidStandard)
// #[derive(Clone, PartialEq, Debug)]
// pub enum AtomTypeInLipid {
//     // todo: Remove this wrapping, and use a plain string if it makmes sense
//     Val(String),
//     // todo: This may be the wrong column
//     // Ca,
//     // Cb,
//     // Cc,
//     // Cd,
//     // Hb,
//     // Hl,
//     // He,
//     // Ho,
//     // Hx,
//     // Pa,
//     // Oc,
//     // Oh,
//     // Op,
//     // Os,
//     // Ot,
//     // H(String),
//     // /// E.g. ligands and water molecules.
//     // Hetero(String),
// }
//
// impl Display for AtomTypeInLipid {
//     fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
//         match self {
//             Self::Val(s) => write!(f, "{s}"),
//         }
//     }
// }
//
// impl FromStr for AtomTypeInLipid {
//     type Err = io::Error;
//
//     fn from_str(s: &str) -> Result<Self, Self::Err> {
//         Ok(Self::Val(s.to_owned()))
//     }
// }

/// This represents an atom, and can be used for various purposes. It is used in various format-specific
/// molecules in this library. You may wish to augment the data here with a custom application-specific
/// format.
#[derive(Clone, Debug, Default)]
pub struct AtomGeneric {
    /// A unique identifier for this atom, within its molecule. This may originate from data in
    /// mmCIF files, Mol2, SDF files, etc.
    pub serial_number: u32,
    pub posit: Vec3,
    pub element: Element,
    /// This identifier will be unique within a given residue. For example, within an
    /// amino acid on a protein. Different residues will have different sets of these.
    /// e.g. "CG1", "CA", "O", "C", "HA", "CD", "C9" etc.
    /// todo: This setup might be protein/aa specific.
    pub type_in_res: Option<AtomTypeInRes>,
    /// There are too many variants of this (with different numbers) for lipids, nucleic
    /// acids etc to use an enum effectively.
    pub type_in_res_general: Option<String>,
    /// Used by Amber and other force fields to apply the correct molecular dynamics parameters for
    /// this atom.
    /// E.g. "c6", "ca", "n3", "ha", "h0" etc, as seen in Mol2 files from AMBER.
    /// e.g.: "ha": hydrogen attached to an aromatic carbon.
    /// "ho": hydrogen on a hydroxyl oxygen
    /// "n3": sp³ nitrogen with three substitutes
    /// "c6": sp² carbon in a pure six-membered aromatic ring (new in GAFF2; lets GAFF distinguish
    /// a benzene carbon from other aromatic caca carbons)
    /// For proteins, this appears to be the same as for `name`.
    pub force_field_type: Option<String>,
    /// An atom-centered electric charge, used in molecular dynamics simulations. In elementary charge units.
    /// These are sometimes loaded from Amber-provided Mol2 or SDF files, and sometimes added after.
    /// We get partial charge for ligands from (e.g. Amber-provided) Mol files, so we load it from the atom, vice
    /// the loaded FF params. Convert to appropriate units prior to running dynamics.
    pub partial_charge: Option<f32>,
    /// Indicates, in proteins, that the atom isn't part of an amino acid. E.g., water or
    /// ligands.
    pub hetero: bool,
    pub occupancy: Option<f32>,
    /// Used by mmCIF files to store alternate conformations. If this isn't None, there may
    /// be, for example, an "A" and "B" variant of this atom at slightly different positions.
    pub alt_conformation_id: Option<String>,
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
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Single => "Single",
            Self::Double => "Double",
            Self::Triple => "Triple",
            Self::Aromatic => "Aromatic",
            Self::Amide => "Amide",
            Self::Dummy => "Dummy",
            Self::Unknown => "Unknown",
            Self::NotConnected => "Not connected",
            Self::Quadruple => "Quadruple",
            Self::Delocalized => "Delocalized",
            Self::PolymericLink => "Polymeric link",
        };
        write!(f, "{name}")
    }
}

impl BondType {
    pub fn order(self) -> f32 {
        match self {
            Self::Aromatic => 1.5,
            Self::Double => 2.,
            Self::Triple => 3.,
            Self::Quadruple => 4.,
            _ => 1.,
        }
    }

    /// A shorthand, visual string.
    pub fn to_visual_str(&self) -> String {
        match self {
            Self::Single => "-",
            Self::Double => "=",
            Self::Triple => "≡",
            Self::Aromatic => "=–",
            Self::Amide => "-am-",
            Self::Dummy => "-",
            Self::Unknown => "-un-",
            Self::NotConnected => "-nc-",
            Self::Quadruple => "-#-",
            Self::Delocalized => "-delo-",
            Self::PolymericLink => "-poly-",
        }
        .to_string()
    }

    /// Return the exact MOL2 bond-type token as an owned `String`.
    /// (Use `&'static str` if you never need it allocated.)
    pub fn to_mol2_str(&self) -> String {
        match self {
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
        }
        .to_string()
    }

    /// SDF format uses a truncated set, and does things like mark every other
    /// aromatic bond as double.
    pub fn to_str_sdf(&self) -> String {
        match self {
            Self::Single | Self::Double | Self::Triple => *self,
            Self::Aromatic => return "4".to_string(),
            _ => Self::Single,
        }
        .to_mol2_str()
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
    /// You may wish to augment these serial numbers with atom indices in downstream
    /// applications, for lookup speed.
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

/// Can be used for amino acid, and nucleotide sequences.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
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
