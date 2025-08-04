//! For opening Mol2 files. These are common molecular descriptions for ligands.
//! [This unofficial resource](https://chemicbook.com/2021/02/20/mol2-file-format-explained-for-beginners-part-2.html)
//! descripts the format.

use std::{
    fmt,
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
    str::FromStr,
};

use lin_alg::f64::Vec3;
use na_seq::{AtomTypeInRes, Element};

use crate::{AtomGeneric, BondGeneric};

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum MolType {
    Small,
    Bipolymer,
    Protein,
    NucleicAcid,
    Saccharide,
}

impl MolType {
    pub fn to_str(self) -> String {
        match self {
            Self::Small => "SMALL",
            Self::Bipolymer => "BIPOLYMER",
            Self::Protein => "PROTEIN",
            Self::NucleicAcid => "NUCLEIC_ACID",
            Self::Saccharide => "SACCHARIDE",
        }
        .to_owned()
    }
}

impl FromStr for MolType {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "SMALL" => Ok(MolType::Small),
            "BIPOLYMER" => Ok(MolType::Bipolymer),
            "PROTEIN" => Ok(MolType::Protein),
            "NUCLEIC_ACID" => Ok(MolType::NucleicAcid),
            "SACCHARIDE" => Ok(MolType::Saccharide),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid MolType: {s}"),
            )),
        }
    }
}

#[derive(Clone, PartialEq, Debug)]
pub enum ChargeType {
    None,
    DelRe,
    Gasteiger,
    GastHuck,
    Huckel,
    Pullman,
    Gauss80,
    Ampac,
    Mulliken,
    Dict,
    MmFf94,
    User,
    Amber,
    Other(String),
}

impl fmt::Display for ChargeType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ChargeType::None => write!(f, "NO_CHARGES"),
            ChargeType::DelRe => write!(f, "DEL_RE"),
            ChargeType::Gasteiger => write!(f, "GASTEIGER"),
            ChargeType::GastHuck => write!(f, "GAST_HUCK"),
            ChargeType::Huckel => write!(f, "HUCKEL"),
            ChargeType::Pullman => write!(f, "PULLMAN"),
            ChargeType::Gauss80 => write!(f, "GAUSS80_CHARGES"),
            ChargeType::Ampac => write!(f, "AMPAC_CHARGES"),
            ChargeType::Mulliken => write!(f, "MULLIKEN_CHARGES"),
            ChargeType::Dict => write!(f, "DICT_CHARGES"),
            ChargeType::MmFf94 => write!(f, "MMFF94_CHARGES"),
            ChargeType::User => write!(f, "USER_CHARGES"),
            ChargeType::Amber => write!(f, "ABCG2"),
            ChargeType::Other(v) => write!(f, "{}", v),
        }
    }
}

impl FromStr for ChargeType {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "NO_CHARGES" => Ok(ChargeType::None),
            "DEL_RE" => Ok(ChargeType::DelRe),
            "GASTEIGER" => Ok(ChargeType::Gasteiger),
            "GAST_HUCK" => Ok(ChargeType::GastHuck),
            "HUCKEL" => Ok(ChargeType::Huckel),
            "PULLMAN" => Ok(ChargeType::Pullman),
            "GAUSS80_CHARGES" => Ok(ChargeType::Gauss80),
            "AMPAC_CHARGES" => Ok(ChargeType::Ampac),
            "MULLIKEN_CHARGES" => Ok(ChargeType::Mulliken),
            "DICT_CHARGES" => Ok(ChargeType::Dict),
            "MMFF94_CHARGES" => Ok(ChargeType::MmFf94),
            "USER_CHARGES" => Ok(ChargeType::User),
            "ABCG2" => Ok(ChargeType::Amber),
            "AMBER FF14SB" => Ok(ChargeType::Amber),
            _ => Ok(ChargeType::Other(s.to_owned())),
        }
    }
}

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
}

impl BondType {
    /// Return the exact MOL2 bond-type token as an owned `String`.
    /// (Use `&'static str` if you never need it allocated.)
    pub fn to_str(self) -> String {
        match self {
            Self::Single => "1",
            Self::Double => "2",
            Self::Triple => "3",
            Self::Aromatic => "ar",
            Self::Amide => "am",
            Self::Dummy => "du",
            Self::Unknown => "un",
            Self::NotConnected => "nc",
        }
        .to_owned()
    }
}

impl FromStr for BondType {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "1" => Ok(BondType::Single),
            "2" => Ok(BondType::Double),
            "3" => Ok(BondType::Triple),
            "ar" => Ok(BondType::Aromatic),
            "am" => Ok(BondType::Amide),
            "du" => Ok(BondType::Dummy),
            "un" => Ok(BondType::Unknown),
            "nc" => Ok(BondType::NotConnected),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid BondType: {s}"),
            )),
        }
    }
}

#[derive(Debug)]
pub struct Mol2 {
    pub ident: String,
    pub mol_type: MolType,
    pub charge_type: ChargeType,
    pub comment: Option<String>,
    // pub metadata: HashMap<String, String>,
    pub atoms: Vec<AtomGeneric>,
    pub bonds: Vec<BondGeneric>,
}

impl Mol2 {
    /// From a string of a Mol2 text file.
    pub fn new(text: &str) -> io::Result<Self> {
        // todo: For these `new` methods in general that take a &str param: Should we use
        // todo R: Reed + Seek instead, and pass a Cursor or File object? Probably doesn't matter.
        // todo Either way, we should keep it consistent between the files.
        let lines: Vec<&str> = text.lines().collect();

        // Example Mol2 header:
        // "
        // @<TRIPOS>MOLECULE
        // 5287969
        // 48 51
        // SMALL
        // USER_CHARGES
        // ****
        // Charges calculated by ChargeFW2 0.1, method: SQE+qp
        // @<TRIPOS>ATOM
        // "

        if lines.len() < 5 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines to parse a MOL2 header",
            ));
        }

        let mut atoms = Vec::new();
        let mut bonds = Vec::new();

        let mut in_atom_section = false;
        let mut in_bond_section = false;

        for line in &lines {
            if line.trim().is_empty() {
                continue;
            }

            let upper = line.to_uppercase();
            if upper.contains("<TRIPOS>ATOM") {
                in_atom_section = true;
                in_bond_section = false;
                continue;
            }

            if upper.contains("<TRIPOS>BOND") {
                in_atom_section = false;
                in_bond_section = true;
                continue;
            }

            if upper.contains("@<TRIPOS>SUBSTRUCTURE") {
                // todo: As required. Example:
                //    1 SER     2 RESIDUE           4 A     SER     1 ROOT
                //      2 VAL    13 RESIDUE           4 A     VAL     2
                //      3 PRO    29 RESIDUE           4 A     PRO     2
                in_atom_section = false;
                in_bond_section = false;
                continue;
            }

            if upper.contains("@<TRIPOS>SET") {
                // todo: As required. Example:
                // ANCHOR          STATIC     ATOMS    <user>   **** Anchor Atom Set
                // 63 127 1110 128 129 610 130 131 132 133 134 740 135 53 741 54 55 617 1482 612 57 58 1485 60 1487 1488 742 1489 743 59 614 1075 611 1486 1076 1481 1077 613 1078 1079 615 616 744 1081 56 618 61 745 1080 1483 738 1074 739 1103 746 1104 1484 1105 1106 1107 1108 1109 1102 1082
                // RIGID           STATIC     BONDS    <user>   **** Rigid Bond Set
                // 56 280 58 59 281 60 61 62 63 64 65 671 672 673 674 332 675 676 484 485 677 678 284 333 24 480 481 26 282 482 334 483 28 486 283 30 31 285 335 487 286 337 338 336 493 494 495 496 27 497 25 499 500 331 279 498 29
                in_atom_section = false;
                in_bond_section = false;
                continue;
            }

            // atom_id atom_name x y z atom_type [subst_id[subst_name [charge [status_bit]]]]
            // Where:
            //
            // atom_id (integer): The ID number of the atom at the time the file was created. This is provided for reference only and is not used when the .mol2 file is read into any mol2 parser software
            // atom_name (string): The name of the atom
            // x (real): The x coordinate of the atom
            // y (real): The y coordinate of the atom
            // z (real): The z coordinate of the atom
            // atom_type (string): The SYBYL atom type for the atom
            // subst_id (integer): The ID number of the substructure containing the atom
            // subst_name (string): The name of the substructure containing the atom
            // charge (real): The charge associated with the atom
            // status_bit (string): The internal SYBYL status bits associated with the atom. These should never be set by the user. Valid status bits are DSPMOD, TYPECOL, CAP, BACKBONE, DICT, ESSENTIAL, WATER, and DIRECT

            if in_atom_section {
                let cols: Vec<&str> = line.split_whitespace().collect();

                if cols.len() < 5 {
                    return Err(io::Error::new(
                        ErrorKind::InvalidData,
                        "Not enough columns on Mol2 atom row",
                    ));
                }

                let serial_number = cols[0].parse::<u32>().map_err(|_| {
                    io::Error::new(ErrorKind::InvalidData, "Could not parse serial number")
                })?;

                // Col 1: e.g. "H", "HG22" etc. Col 5: "C.3", "N.p13" etc.
                let mut atom_name = cols[1].to_owned();
                if let Some((before_dot, _after_dot)) = atom_name.split_once('.') {
                    atom_name = before_dot.to_string();
                }

                let element = match Element::from_letter(&atom_name) {
                    Ok(l) => l,
                    Err(e) => {
                        if atom_name.len() > 1 {
                            // It might be something like "c3", "c1" etc."
                            Element::from_letter(&atom_name[0..1])?
                        } else {
                            return Err(e);
                        }
                    }
                };

                let x = cols[2].parse::<f64>().map_err(|_| {
                    io::Error::new(ErrorKind::InvalidData, "Could not parse X coordinate")
                })?;
                let y = cols[3].parse::<f64>().map_err(|_| {
                    io::Error::new(ErrorKind::InvalidData, "Could not parse Y coordinate")
                })?;
                let z = cols[4].parse::<f64>().map_err(|_| {
                    io::Error::new(ErrorKind::InvalidData, "Could not parse Z coordinate")
                })?;

                let charge = if cols.len() >= 9 {
                    cols[8].parse::<f32>().unwrap_or_default()
                } else {
                    0.
                };

                // todo: ALso, parse the charge type at line 4.
                let partial_charge = if charge.abs() < 0.000001 {
                    None
                } else {
                    Some(charge)
                };

                let type_in_res = if !atom_name.is_empty() {
                    Some(AtomTypeInRes::Hetero(atom_name.to_string()))
                } else {
                    None
                };

                atoms.push(AtomGeneric {
                    serial_number,
                    type_in_res,
                    posit: Vec3 { x, y, z }, // or however you store coordinates
                    element,
                    occupancy: None,
                    partial_charge,
                    force_field_type: Some(cols[5].to_string()),
                    hetero: true,
                });
            }

            if in_bond_section {
                let cols: Vec<&str> = line.split_whitespace().collect();

                if cols.len() < 4 {
                    return Err(io::Error::new(
                        ErrorKind::InvalidData,
                        "Not enough columns when parsing Mol2 bonds.",
                    ));
                }

                let atom_0_sn = cols[1].parse::<u32>().map_err(|_| {
                    io::Error::new(
                        ErrorKind::InvalidData,
                        format!("Could not parse atom 0 in bond: {}", cols[1]),
                    )
                })?;

                let atom_1_sn = cols[2].parse::<u32>().map_err(|_| {
                    io::Error::new(
                        ErrorKind::InvalidData,
                        format!("Could not parse atom 1 in bond: {}", cols[2]),
                    )
                })?;

                // 1 = single
                // 2 = double
                // 3 = triple
                // am = amide
                // ar = aromatic
                // du = dummy
                // un = unknown (cannot be determined from the parameter tables)
                // nc = not connected
                bonds.push(BondGeneric {
                    bond_type: cols[3].to_owned(),
                    atom_0_sn,
                    atom_1_sn,
                });
            }
        }

        // Note: This may not be the identifier we think of.
        let ident = lines[1].to_owned();
        let mol_type = MolType::from_str(lines[3])?;
        let charge_type = ChargeType::from_str(lines[4])?;

        // todo: Multi-line comments are supported by Mol2.
        let comment = if lines[5].contains("****") {
            None
        } else {
            Some(lines[5].to_owned())
        };

        Ok(Self {
            ident,
            mol_type,
            charge_type,
            atoms,
            bonds,
            comment,
        })
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        //todo: Fix this so it outputs mol2 instead of sdf.
        let mut file = File::create(path)?;

        // There is a subtlety here. Add that to your parser as well. There are two values
        // todo in the files we have; this top ident is not the DB id.
        writeln!(file, "@<TRIPOS>MOLECULE")?;
        writeln!(file, "{}", self.ident)?;
        writeln!(file, "{} {}", self.atoms.len(), self.bonds.len())?;
        writeln!(file, "{}", self.mol_type.to_str())?;
        writeln!(file, "{}", self.charge_type)?;

        // //  todo: Multi-line comments are supported by Mol2
        // let comment = match &self.comment {
        //     Some(c) => &c,
        //     None => "****",
        // };

        // **** Means a non-optional field is empty.
        // writeln!(file, "{comment}")?;
        // Optional line (comments, molecule weight, etc.)

        writeln!(file, "")?;
        writeln!(file, "")?;

        writeln!(file, "@<TRIPOS>ATOM")?;
        for (i, atom) in self.atoms.iter().enumerate() {
            let type_in_res = match &atom.type_in_res {
                Some(n) => n.to_string(),
                None => atom.element.to_letter(),
            };

            let ff_type = match &atom.force_field_type {
                Some(f) => f.to_owned(),
                // Not ideal, but will do as a placeholder.
                None => atom.element.to_letter().to_lowercase(),
            };

            // todo: A/R
            // let res_name = String::new();
            // for res in &self.

            writeln!(
                file,
                "{:>7} {:<8} {:>10.4} {:>10.4} {:>10.4} {:<6} {:>5} {:<8} {:>9.6}",
                atom.serial_number,
                type_in_res,
                atom.posit.x,
                atom.posit.y,
                atom.posit.z,
                ff_type,
                "1",        // Assumes 1 residue.
                self.ident, // todo: This should really be the residue information.
                atom.partial_charge.unwrap_or_default()
            )?;
        }

        writeln!(file, "@<TRIPOS>BOND")?;

        for (i, bond) in self.bonds.iter().enumerate() {
            let bond_type = if bond.bond_type.is_empty() {
                "1"
            } else {
                &bond.bond_type
            };

            writeln!(
                file,
                "{:>6}{:>6}{:>6} {:<3}",
                i + 1,
                bond.atom_0_sn,
                bond.atom_1_sn,
                bond_type,
            )?;
        }

        Ok(())
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        let data_str: String = String::from_utf8(buffer)
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid UTF8"))?;

        Self::new(&data_str)
    }
}
