//! A [`winnow`](https://docs.rs/winnow) re-implementation of [`crate::mol2`]. It opens Mol2 files
//! (the [Tripos Mol2 spec](https://zhanggroup.org/DockRMSD/mol2.pdf) variant) and is otherwise
//! API-compatible with `mol2.rs`: the [`Mol2`] struct, [`MolType`], [`ChargeType`], and the
//! read/write/`From` surface all mirror the originals.
//!
//! Only the *parsing* direction uses winnow; winnow has no serialization support, so the
//! `write_to`/`save` path is identical to `mol2.rs`. The shared winnow primitives ([`run`],
//! [`ws_token`], ...) and the pharmacophore helpers are reused from [`crate::sdf_winnow`],
//! mirroring how `mol2.rs` pulls those helpers from `sdf.rs`.

use std::{
    collections::HashMap,
    fmt, fs,
    fs::File,
    io,
    io::{BufWriter, ErrorKind, Write},
    path::Path,
    str::FromStr,
};

use bio_apis::amber_geostd;
use lin_alg::f64::Vec3;
use na_seq::AtomTypeInRes;
use winnow::combinator::opt;

use crate::{
    AtomGeneric, BondGeneric, BondType, PharmacophoreFeatureGeneric, el_from_atom_name,
    sdf_winnow::{
        Sdf, format_pharmacophore_features, parse_pharmacophore_features, run, ws_f32, ws_f64,
        ws_token, ws_u32,
    },
};

// For our custom format addition.
const PHARMACOPHORE_TAG: &str = "@<BIO_FILES>PHARMACOPHORE";

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

/// This implements the [Tripos Mol2 spec](https://zhanggroup.org/DockRMSD/mol2.pdf).
/// It's a format used for small organic molecules, and may include force field names and partial
/// charges, making it suitable for molecular dynamics applications. This struct will likely
/// be used as an intermediate format, and converted to something application-specific.
#[derive(Clone, Debug)]
pub struct Mol2 {
    pub ident: String,
    pub metadata: HashMap<String, String>,
    pub atoms: Vec<AtomGeneric>,
    pub bonds: Vec<BondGeneric>,
    pub mol_type: MolType,
    pub charge_type: ChargeType,
    /// Note: I have not observed these in Mol2 files in the wild, but
    /// we are using them in Molchanica; setting up in a simimlar way to
    /// how they're stored in SDF.
    pub pharmacophore_features: Vec<PharmacophoreFeatureGeneric>,
    pub comment: Option<String>,
}

/// A parsed Mol2 `@<TRIPOS>ATOM` row:
/// `atom_id atom_name x y z atom_type [subst_id subst_name charge ...]`.
/// The charge is column 9 (index 8) and only present when the optional trailing
/// `subst_id`/`subst_name`/`charge` triple is.
fn atom_row(line: &str) -> io::Result<AtomGeneric> {
    let (serial_number, raw_name, x, y, z, ff_type, trailing): (
        u32,
        &str,
        f64,
        f64,
        f64,
        &str,
        Option<(&str, &str, f32)>,
    ) = run(
        (
            ws_u32,
            ws_token,
            ws_f64,
            ws_f64,
            ws_f64,
            ws_token,
            opt((ws_token, ws_token, ws_f32)),
        ),
        line,
        "Mol2 atom row",
    )?;

    // Col 1: e.g. "H", "HG22" etc. Col 5: "C.3", "N.p13" etc. Drop the SYBYL
    // subtype suffix after the dot.
    let atom_name = match raw_name.split_once('.') {
        Some((before_dot, _)) => before_dot,
        None => raw_name,
    };

    let element = el_from_atom_name(atom_name);

    let charge = trailing.map(|(_, _, q)| q).unwrap_or(0.);
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

    Ok(AtomGeneric {
        serial_number,
        type_in_res,
        posit: Vec3 { x, y, z },
        element,
        partial_charge,
        force_field_type: Some(ff_type.to_string()),
        hetero: true,
        ..Default::default()
    })
}

/// A parsed Mol2 `@<TRIPOS>BOND` row: `bond_id atom_0 atom_1 bond_type`.
fn bond_row(line: &str) -> io::Result<BondGeneric> {
    let (_bond_id, atom_0_sn, atom_1_sn, bond_type): (u32, u32, u32, &str) =
        run((ws_u32, ws_u32, ws_u32, ws_token), line, "Mol2 bond row")?;

    Ok(BondGeneric {
        bond_type: BondType::from_str(bond_type)?,
        atom_0_sn,
        atom_1_sn,
    })
}

impl Mol2 {
    /// From a string of a Mol2 text file.
    pub fn new(text: &str) -> io::Result<Self> {
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
        let mut pharmacophore_rows = Vec::new();

        let mut in_atom_section = false;
        let mut in_bond_section = false;
        let mut in_bond_pharmacophore_section = false;
        let mut metadata = HashMap::new();
        let mut pending_metadata_key: Option<String> = None;
        let mut pending_metadata_lines: Vec<String> = Vec::new();

        for line in &lines {
            if line.trim().is_empty() {
                continue;
            }

            let upper = line.to_uppercase();
            if upper.contains("<TRIPOS>ATOM") {
                in_atom_section = true;
                in_bond_section = false;
                in_bond_pharmacophore_section = false;
                if let Some(key) = pending_metadata_key.take() {
                    metadata.insert(key, pending_metadata_lines.join("\n"));
                    pending_metadata_lines.clear();
                }
                continue;
            }

            if upper.contains("<TRIPOS>BOND") {
                in_atom_section = false;
                in_bond_section = true;
                in_bond_pharmacophore_section = false;
                if let Some(key) = pending_metadata_key.take() {
                    metadata.insert(key, pending_metadata_lines.join("\n"));
                    pending_metadata_lines.clear();
                }
                continue;
            }

            if upper.contains("@<TRIPOS>SUBSTRUCTURE") {
                // todo: As required. Example:
                //    1 SER     2 RESIDUE           4 A     SER     1 ROOT
                //      2 VAL    13 RESIDUE           4 A     VAL     2
                //      3 PRO    29 RESIDUE           4 A     PRO     2
                in_atom_section = false;
                in_bond_section = false;
                in_bond_pharmacophore_section = false;
                if let Some(key) = pending_metadata_key.take() {
                    metadata.insert(key, pending_metadata_lines.join("\n"));
                    pending_metadata_lines.clear();
                }
                continue;
            }

            if upper.contains("@<TRIPOS>SET") {
                // todo: As required. Example:
                // ANCHOR          STATIC     ATOMS    <user>   **** Anchor Atom Set
                // RIGID           STATIC     BONDS    <user>   **** Rigid Bond Set
                in_atom_section = false;
                in_bond_section = false;
                in_bond_pharmacophore_section = false;
                continue;
            }

            if upper.contains(PHARMACOPHORE_TAG) {
                in_atom_section = false;
                in_bond_section = false;
                in_bond_pharmacophore_section = true;
                if let Some(key) = pending_metadata_key.take() {
                    metadata.insert(key, pending_metadata_lines.join("\n"));
                    pending_metadata_lines.clear();
                }
                continue;
            }

            // Our custom metadata parsing. Match all @ lines that don't match one above.
            if line.starts_with('@') && !upper.contains("<TRIPOS>") {
                in_atom_section = false;
                in_bond_section = false;
                in_bond_pharmacophore_section = false;
                if let Some(key) = pending_metadata_key.take() {
                    metadata.insert(key, pending_metadata_lines.join("\n"));
                    pending_metadata_lines.clear();
                }
                pending_metadata_key = Some(line.trim_start_matches('@').trim().to_owned());
                continue;
            }

            if pending_metadata_key.is_some() {
                pending_metadata_lines.push(line.to_string());
                continue;
            }

            // atom_id atom_name x y z atom_type [subst_id[subst_name [charge [status_bit]]]]
            // See `atom_row` for the field breakdown.
            if in_atom_section {
                atoms.push(atom_row(line)?);
            }

            if in_bond_section {
                bonds.push(bond_row(line)?);
            }

            if in_bond_pharmacophore_section {
                pharmacophore_rows.push(line.to_owned());
            }
        }

        // Flush any trailing custom metadata section.
        if let Some(key) = pending_metadata_key.take() {
            metadata.insert(key, pending_metadata_lines.join("\n"));
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

        let pharmacophore_features = if pharmacophore_rows.is_empty() {
            Vec::new()
        } else {
            parse_pharmacophore_features(&pharmacophore_rows)?
        };

        Ok(Self {
            ident,
            metadata,
            mol_type,
            charge_type,
            atoms,
            bonds,
            comment,
            pharmacophore_features,
        })
    }

    pub fn write_to(&self, w: &mut impl Write) -> io::Result<()> {
        // There is a subtlety here. Add that to your parser as well. There are two values
        // todo in the ws we have; this top ident is not the DB id.
        writeln!(w, "@<TRIPOS>MOLECULE")?;
        writeln!(w, "{}", self.ident)?;
        writeln!(w, "{} {}", self.atoms.len(), self.bonds.len())?;
        writeln!(w, "{}", self.mol_type.to_str())?;
        writeln!(w, "{}", self.charge_type)?;

        // **** Means a non-optional field is empty.
        writeln!(w)?;
        writeln!(w)?;

        writeln!(w, "@<TRIPOS>ATOM")?;
        for atom in &self.atoms {
            let type_in_res = match &atom.type_in_res {
                Some(n) => n.to_string(),
                None => atom.element.to_letter(),
            };

            let ff_type = match &atom.force_field_type {
                Some(f) => f.to_owned(),
                // Not ideal, but will do as a placeholder.
                None => atom.element.to_letter().to_lowercase(),
            };

            writeln!(
                w,
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

        writeln!(w, "@<TRIPOS>BOND")?;

        for (i, bond) in self.bonds.iter().enumerate() {
            writeln!(
                w,
                "{:>6}{:>6}{:>6} {:<3}",
                i + 1,
                bond.atom_0_sn,
                bond.atom_1_sn,
                bond.bond_type.to_mol2_str(),
            )?;
        }

        if !self.pharmacophore_features.is_empty() {
            writeln!(w, "{PHARMACOPHORE_TAG}")?;
            let v = format_pharmacophore_features(&self.pharmacophore_features);
            write!(w, "{v}")?;
        }

        // Unofficial way of writing metadata
        for (k, v) in &self.metadata {
            writeln!(w, "\n@{k}")?;
            writeln!(w, "{v}")?;
        }

        Ok(())
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);
        self.write_to(&mut w)
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let data_str = fs::read_to_string(path)?;
        Self::new(&data_str)
    }

    /// Download  rom our Amber Geostd DB using a PDBe/Amber ID.
    pub fn load_amber_geostd(ident: &str) -> io::Result<Self> {
        let data_str = amber_geostd::load_mol2(ident)
            .map_err(|e| io::Error::other(format!("Error loading: {e:?}")))?;
        Self::new(&data_str)
    }
}

impl From<Sdf> for Mol2 {
    fn from(m: Sdf) -> Self {
        Self {
            ident: m.ident.clone(),
            metadata: m.metadata.clone(),
            atoms: m.atoms.clone(),
            bonds: m.bonds.clone(),
            mol_type: MolType::Small,
            charge_type: ChargeType::None,
            pharmacophore_features: m.pharmacophore_features,
            comment: None,
        }
    }
}
