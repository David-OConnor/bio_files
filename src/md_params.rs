#![allow(clippy::excessive_precision)]

//! Contains parameters used in Amber Forcefields. For details on these formats,
//! see the [Amber Reference Manual](https://ambermd.org/doc12/Amber25.pdf), section
//! 15: Reading and modifying Amber parameter files.
//!
//! Called by both the `dat`, and `frcmod` modules. These formats share line formats, but
//! arrange them in different ways.
//!
//! For ligands, `atom_type` is a "Type 3". For proteins/AAs, we are currently treating it
//! as a type 1, but we're unclear on this.

use crate::{LipidStandard, ResidueEnd};
use na_seq::Element::Hydrogen;
use na_seq::{AminoAcidGeneral, AtomTypeInRes, Nucleotide};
use std::collections::HashSet;
use std::{
    collections::HashMap,
    fs::File,
    io::{self, ErrorKind, Read},
    path::Path,
    str::FromStr,
};

/// Data for a MASS entry: e.g. "CT 12.01100" with optional comment.
#[derive(Debug, Clone)]
pub struct MassParams {
    pub atom_type: String,
    /// Atomic mass units (Daltons)
    pub mass: f32,
    // /// ATPOL: Atomic polarizability (Å^3).
    // /// Intended for Slater–Kirkwood or future polarizable models, and unused by Amber (?)
    // pub polarizability: f32,
    pub comment: Option<String>,
}

impl MassParams {
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.split_whitespace().collect();

        // Allow Skipping ATPOL which we don't currently use, and is sometimes missing.
        if cols.len() < 2 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Not enough cols (Mass) when parsing line {line}"),
            ));
        }

        let atom_type = cols[0].to_string();
        let mass = parse_float(cols[1])?;

        // Note: This skips comments where this is a missing col[2].
        let mut comment = None;
        if cols.len() >= 4 {
            comment = Some(cols[3..].join(" "));
        }

        Ok(Self {
            atom_type,
            mass,
            comment,
        })
    }
}

/// Amber RM 2025, 15.1.6
/// Data for a BOND entry: e.g. "CT-CT  310.0    1.526" with optional comment
/// Length between 2 covalently bonded atoms.
#[derive(Debug, Clone)]
pub struct BondStretchingParams {
    pub atom_types: (String, String),
    /// Force constant. (Similar to a spring constant). kcal/mol/Å²
    pub k_b: f32,
    /// Equilibrium bond length. Å
    pub r_0: f32,
    pub comment: Option<String>,
}

impl BondStretchingParams {
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Not enough cols (Bond) when parsing line {line}"),
            ));
        }

        let (atom_types, col1_i) = get_atom_types(&cols);
        let atom_types = (atom_types[0].to_owned(), atom_types[1].to_owned());

        let k = parse_float(cols[col1_i])?;
        let r_0 = parse_float(cols[col1_i + 1])?;

        // We ignore the remaining cols for now: Source, # of ref geometries used to fit,
        // and RMS deviation of the fit.

        let mut comment = None;
        if cols.len() >= col1_i + 2 {
            comment = Some(cols[col1_i + 2..].join(" "));
        }

        Ok(Self {
            atom_types,
            k_b: k,
            r_0,
            comment,
        })
    }
}

/// Amber RM 2025, 15.1.6
/// Data for an ANGLE entry: e.g. "CT-CT-CT  63.0    109.5" with optional comment
/// Angle between 3 linear covalently-bonded atoms (2 bonds)
#[derive(Debug, Clone)]
pub struct AngleBendingParams {
    pub atom_types: (String, String, String),
    /// Force constant. kcal/mol/rad²
    pub k: f32,
    /// In radians.
    pub theta_0: f32,
    pub comment: Option<String>,
}

impl AngleBendingParams {
    /// Parse a single valence-angle record from a GAFF/Amber `.dat` or `.frcmod` file.
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Not enough cols (Angle) when parsing line {line}"),
            ));
        }

        let (atom_types, col1_i) = get_atom_types(&cols);
        let atom_types = (
            atom_types[0].to_owned(),
            atom_types[1].to_owned(),
            atom_types[2].to_owned(),
        );

        let k = parse_float(cols[col1_i])?;
        let angle = parse_float(cols[col1_i + 1])?.to_radians();

        // We ignore the remaining cols for now: Source, # of ref geometries used to fit,
        // and RMS deviation of the fit.

        let mut comment = None;
        if cols.len() >= col1_i + 2 {
            comment = Some(cols[col1_i + 2..].join(" "));
        }

        Ok(Self {
            atom_types,
            k,
            theta_0: angle,
            comment,
        })
    }
}

/// Also known as Torsion angle.
///
/// Angle between 4 linear covalently-bonded atoms ("proper"), or 3 atoms in a hub-and-spoke
/// configuration, with atom 3 as the hub ("improper"). In either case, this is the angle between the planes of
/// atoms 1-2-3, and 2-3-4. (Rotation around the 2-3 bond)
#[derive(Debug, Clone, Default)]
pub struct DihedralParams {
    /// "ca", "n", "cd", "sh" etc.
    pub atom_types: (String, String, String, String),
    /// Scaling factor used for barrier height.
    /// "Splits the torsion term into individual contributions for
    /// each pair of atoms involved in the torsion."
    /// Always 1 for improper dihedrals. (Not present in the Amber files for improper)
    pub divider: u8,
    /// Also known as V_n. kcal/mol.
    pub barrier_height: f32,
    /// Phase, in radians. Often 0 or τ/2. Maximum energy
    /// is encountered at this value, and other values implied by periodicity.
    /// For example, if this is 0, and periodicity is 3, there is no torsion
    /// force applied for dihedral angles 0, τ/3, and 2τ/3.
    pub phase: f32,
    /// An integer, relative to a full rotation; there is a minimum once every
    /// this/τ radians.
    ///
    /// "If the torsion definition has a "negative" periodicity (-2 in the case above), it tells
    /// programs reading the parameter file that additional terms are present for that
    /// particular connectivity.
    pub periodicity: u8,
    pub comment: Option<String>,
}

impl DihedralParams {
    /// For both FRCMOD, and Dat. For both proper, and improper. Returns `true` if improper.
    pub fn from_line(line: &str) -> io::Result<(Self, bool)> {
        let cols: Vec<_> = line.split_whitespace().collect();

        if cols.len() < 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Not enough cols (Dihedral) when parsing line {line}"),
            ));
        }

        let (atom_types, mut col1_i) = get_atom_types(&cols);
        let atom_types = (
            atom_types[0].to_owned(),
            atom_types[1].to_owned(),
            atom_types[2].to_owned(),
            atom_types[3].to_owned(),
        );

        let mut improper = true;
        let mut integer_divisor = 1; // Default, for dihedral.
        // Determine if an improper or not, prescense of decimal in col 1. This means it's improper,
        // as we're skipping the integer.

        if !cols[col1_i].contains(".") {
            integer_divisor = parse_float(cols[col1_i])? as u8;
            col1_i += 1;
            improper = false;
        }

        let barrier_height_vn = parse_float(cols[col1_i])?;
        let phase = parse_float(cols[col1_i + 1])?.to_radians();

        // A negative periodicity in Amber params indicates that there are additional terms
        // are present. We ignore those for now.
        let periodicity = parse_float(cols[col1_i + 2])?.abs() as u8;

        // We ignore the remaining cols for now: Source, # of ref geometries used to fit,
        // and RMS deviation of the fit.

        let mut comment = None;
        if cols.len() >= col1_i + 3 {
            comment = Some(cols[col1_i + 3..].join(" "));
        }

        Ok((
            Self {
                atom_types,
                divider: integer_divisor,
                barrier_height: barrier_height_vn,
                phase,
                periodicity,
                comment,
            },
            improper,
        ))
    }
}

#[derive(Debug, Clone)]
/// Represents Lennard Jones parameters. This approximate Pauli Exclusion (i.e. exchange interactions)
/// with Van Der Waals ones. Note: Amber stores Rmin / 2 in Å. (This is called R star). We convert to σ, which can
/// be used in more general LJ formulas. The relation: R_min = 2^(1/6) σ. σ = 2 R_star / 2^(1/6)
/// Amber RM, section 15.1.7
pub struct LjParams {
    pub atom_type: String,
    /// σ. derived from Van der Waals radius, Å. Note that Amber parameter files use R_min,
    /// vice σ. The value in this field is σ, which we compute when parsing.
    pub sigma: f32,
    /// Energy, kcal/mol. (Represents depth of the potential well).
    /// σ(i, j) = 0.5 * (σ_i + σ_j)
    /// ε(i, j) = sqrt(ε_i * ε_j)
    pub eps: f32,
}

impl LjParams {
    /// Parse a single van-der-Waals (Lennard-Jones) parameter line in a dat file, e.g. `parm19.dat`
    /// from Amber.
    pub fn from_line(line: &str) -> io::Result<Self> {
        // todo: QC this factor of 2!
        // 1.122 is 2^(1/6)
        // todo: We're getting conflicting information on if we should
        // todo use a factor of 2, or 4 as the prefix here.
        const SIGMA_FACTOR: f32 = 2. / 1.122_462_048_309_373;

        let cols: Vec<_> = line.split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Not enough cols (Lennard Jones) when parsing line {line}"),
            ));
        }

        let atom_type = cols[0].to_string();
        let r_star = parse_float(cols[1])?;
        let eps = parse_float(cols[2])?;

        let sigma = r_star * SIGMA_FACTOR;

        Ok(Self {
            atom_type,
            sigma,
            eps,
        })
    }
}

#[derive(Clone, Debug)]
pub struct ChargeParamsProtein {
    /// For proteins. The residue-specific ID. We use this value to map forcefield type
    /// to atoms loaded from mmCIF etc; these will have this `type_in_res`, but not
    /// an Amber ff type. We apply the charge here to the atom based on its `type_in_res` and AA type,
    /// and apply its FF type.
    ///
    /// Once the FF type is applied here, we can map other params, e.g. Vdw, and bonded terms to it.
    pub type_in_res: AtomTypeInRes,
    /// "XC", "H1" etc.
    pub ff_type: String,
    /// Partial charge. Units of elementary charge.
    pub charge: f32,
}

/// See notes on `ChargeParams`; equivalent here. For lipids, nucleic acids etc.
#[derive(Clone, Debug)]
pub struct ChargeParams {
    pub type_in_res: String,
    pub ff_type: String,
    pub charge: f32,
}

/// Top-level lib, dat, or frcmod data. We store the name-tuples in fields, vice as HashMaps here,
/// for parsing flexibility.
///
/// Note that we don't include partial charges here, as they come from Mol2 files; this struct
/// is for data parsed from DAT, FRCMOD etc files.
#[derive(Debug, Default)]
pub struct ForceFieldParamsVec {
    /// Length between 2 covalently bonded atoms.
    pub bond: Vec<BondStretchingParams>,
    /// Angle between 3 linear covalently-bonded atoms (2 bonds)
    pub angle: Vec<AngleBendingParams>,
    /// Angle between 4 linear covalently-bonded atoms (3 bonds). This is
    /// the angle between the planes of atoms 1-2-3, and 2-3-4. (Rotation around the 2-3 bond)
    pub dihedral: Vec<DihedralParams>,
    /// Angle between 4 covalently-bonded atoms (3 bonds), in a hub-and-spoke
    /// arrangement. The third atom is the hub. This is the angle between the planes of
    /// atoms 1-2-3, and 2-3-4. Note that these are generally only included for planar configurations,
    /// and always hold a planar dihedral shape. (e.g. τ/2 with symmetry 2)
    pub improper: Vec<DihedralParams>,
    pub mass: Vec<MassParams>,
    pub lennard_jones: Vec<LjParams>,
    pub remarks: Vec<String>,
}

/// This variant of forcefield parameters offers the fastest lookups. Unlike the Vec and Hashmap
/// based parameter structs, the indices are provincial
/// to specific sets of atoms, bonds, sets of 3 atoms etc. For a description of fields, see `ForceFieldParams`, or the individual
/// param-type structs here.
///
/// Note: The single-atom fields of `mass` and `partial_charges` are omitted: They're part of our
/// `AtomDynamics` struct.`
///
/// We can use this for fast lookups in rust applications and libraries, or as a precursor to writing
/// parameter files for other MD engines like GROMACS.
#[derive(Clone, Debug, Default)]
pub struct ForceFieldParamsIndexed {
    pub mass: HashMap<usize, MassParams>,
    pub bond_stretching: HashMap<(usize, usize), BondStretchingParams>,
    /// Any bond to Hydrogen if configured as constrained. (Distance^2 in Å, 1 / mass in Daltons)
    pub bond_rigid_constraints: HashMap<(usize, usize), (f32, f32)>,
    pub angle: HashMap<(usize, usize, usize), AngleBendingParams>,
    pub dihedral: HashMap<(usize, usize, usize, usize), Vec<DihedralParams>>,
    pub improper: HashMap<(usize, usize, usize, usize), Vec<DihedralParams>>,
    /// We use this to determine which 1-2 exclusions to apply for non-bonded forces. We use this
    /// instead of `bond_stretching`, because `bond_stretching` omits bonds to Hydrogen, which we need
    /// to account when applying exclusions.
    pub bonds_topology: HashSet<(usize, usize)>,
    pub lennard_jones: HashMap<usize, LjParams>,
}

/// Force field parameters, e.g. from Amber. Similar to `ForceFieldParams` but
/// with Hashmap-based keys (of atom-name tuples) for fast look-ups. See that struct
/// for a description of each field.
///
/// For descriptions of each field and the units used, reference the structs in bio_files, of which
/// this uses internally.
#[derive(Clone, Debug, Default)]
pub struct ForceFieldParams {
    /// Length between 2 covalently bonded atoms.
    pub bond: HashMap<(String, String), BondStretchingParams>,
    /// Angle between 3 linear covalently-bonded atoms (2 bonds)
    pub angle: HashMap<(String, String, String), AngleBendingParams>,
    /// Angle between 4 linear covalently-bonded atoms (3 bonds). This is
    /// the angle between the planes of atoms 1-2-3, and 2-3-4. (Rotation around the 2-3 bond)
    /// This is a Vec, as there can be multiple terms for proper dihedrals. (Negative
    /// periodicity is a flag meaning there are follow-on terms)
    pub dihedral: HashMap<(String, String, String, String), Vec<DihedralParams>>,
    /// Angle between 4 covalently-bonded atoms (3 bonds), in a hub-and-spoke
    /// arrangement. The third atom is the hub. This is the angle between the planes of
    /// atoms 1-2-3, and 2-3-4. Note that these are generally only included for planar configurations,
    /// and always hold a planar dihedral shape. (e.g. τ/2 with symmetry 2)
    /// It's possible, but unlikely there can be more than one improper term
    pub improper: HashMap<(String, String, String, String), Vec<DihedralParams>>,
    pub mass: HashMap<String, MassParams>,
    pub lennard_jones: HashMap<String, LjParams>,
}

impl ForceFieldParams {
    /// Restructures params so the `atom_type` fields are arranged as HashMap keys, for faster
    /// lookup.
    pub fn new(params: &ForceFieldParamsVec) -> Self {
        let mut result = Self::default();

        for val in &params.mass {
            result.mass.insert(val.atom_type.clone(), val.clone());
        }

        for val in &params.bond {
            result.bond.insert(val.atom_types.clone(), val.clone());
        }

        for val in &params.angle {
            result.angle.insert(val.atom_types.clone(), val.clone());
        }

        // Insert, or append, as required. There can be multiple proper dihedral terms.
        for val in &params.dihedral {
            result
                .dihedral
                .entry(val.atom_types.clone())
                .and_modify(|v| v.push(val.clone()))
                .or_insert_with(|| vec![val.clone()]);
        }

        for val in &params.improper {
            result
                .improper
                .entry(val.atom_types.clone())
                .and_modify(|v| v.push(val.clone()))
                .or_insert_with(|| vec![val.clone()]);
        }

        for val in &params.lennard_jones {
            result
                .lennard_jones
                .insert(val.atom_type.clone(), val.clone());
        }

        result
    }

    /// Merge `other` into a copy of `self`. Entries already present in `self` win;
    /// entries only in `other` are added. This lets you build a comprehensive global
    /// fallback from multiple force-field tables (e.g. GAFF2 + ff19SB).
    pub fn merge_with(&self, other: &Self) -> Self {
        let mut result = other.clone();
        result.mass.extend(self.mass.clone());
        result.lennard_jones.extend(self.lennard_jones.clone());
        result.bond.extend(self.bond.clone());
        result.angle.extend(self.angle.clone());
        // For dihedrals/impropers, self's entry replaces other's entirely (same key).
        result.dihedral.extend(self.dihedral.clone());
        result.improper.extend(self.improper.clone());
        result
    }

    /// A convenience wrapper.
    pub fn from_frcmod(text: &str) -> io::Result<Self> {
        Ok(Self::new(&ForceFieldParamsVec::from_frcmod(text)?))
    }

    /// A convenience wrapper.
    pub fn from_dat(text: &str) -> io::Result<Self> {
        Ok(Self::new(&ForceFieldParamsVec::from_dat(text)?))
    }

    /// A convenience wrapper.
    pub fn load_frcmod(path: &Path) -> io::Result<Self> {
        Ok(Self::new(&ForceFieldParamsVec::load_frcmod(path)?))
    }

    /// A convenience wrapper.
    pub fn load_dat(path: &Path) -> io::Result<Self> {
        Ok(Self::new(&ForceFieldParamsVec::load_dat(path)?))
    }

    /// For the `get_` methods below. Expand possible wildcard forms of an atom type, keeping priority order:
    /// 1. Exact atom name
    /// 2. Pattern with same first letter and '*'
    /// 3. Global wildcard "X"
    fn wildcard_variants(atom: &str) -> Vec<String> {
        let mut out = Vec::new();
        out.push(atom.to_string()); // exact

        if !atom.is_empty() {
            let first = atom.chars().next().unwrap();
            // Only add meaningful ones like C*, N*, O*, etc.
            if first.is_ascii_alphabetic() {
                out.push(format!("{}*", first));
            }
        }
        out.push("X".to_string());
        out
    }

    /// A utility function that handles proper and improper dihedral data,
    /// tries both atom orders, and falls back to wildcard (“X”) matches on
    /// the outer atoms when an exact hit is not found.
    pub fn get_bond(
        &self,
        atom_types: &(String, String),
        wildcard_allowed: bool,
    ) -> Option<&BondStretchingParams> {
        let (a_variants, b_variants) = if wildcard_allowed {
            (
                Self::wildcard_variants(&atom_types.0),
                Self::wildcard_variants(&atom_types.1),
            )
        } else {
            (
                vec![atom_types.0.to_string()],
                vec![atom_types.1.to_string()],
            )
        };

        // Priority: exact before partial before X
        for a in &a_variants {
            for b in &b_variants {
                // try both orders
                for &(k0, k1) in &[(a, b), (b, a)] {
                    let key = (k0.clone(), k1.clone());
                    if let Some(hit) = self.bond.get(&key) {
                        return Some(hit);
                    }
                }
            }
        }
        None
    }

    // todo: YOu may need to augment all these helps with support for "C*", "N*" etc.

    /// A utility function that handles proper and improper dihedral data,
    /// tries both atom orders, and falls back to wildcard (“X”) matches on
    /// the outer atoms when an exact hit is not found.
    pub fn get_valence_angle(
        &self,
        atom_types: &(String, String, String),
        wildcard_allowed: bool,
    ) -> Option<&AngleBendingParams> {
        let (a_variants, b_variants, c_variants) = if wildcard_allowed {
            (
                Self::wildcard_variants(&atom_types.0),
                Self::wildcard_variants(&atom_types.1),
                Self::wildcard_variants(&atom_types.2),
            )
        } else {
            (
                vec![atom_types.0.to_string()],
                vec![atom_types.1.to_string()],
                vec![atom_types.2.to_string()],
            )
        };

        // Try combinations in both directions (a-b-c and c-b-a)
        for a in &a_variants {
            for b in &b_variants {
                for c in &c_variants {
                    for &(k0, k1, k2) in &[(a, b, c), (c, b, a)] {
                        let key = (k0.clone(), k1.clone(), k2.clone());
                        if let Some(hit) = self.angle.get(&key) {
                            return Some(hit);
                        }
                    }
                }
            }
        }
        None
    }

    /// A utility function that handles proper and improper dihedral data,
    /// tries both atom orders, and falls back to wildcard (“X”) matches on
    /// the outer atoms when an exact hit is not found.
    ///
    /// We return multiple, as there can be multiple dihedrals for a single atom type set;
    /// we add them during computations.
    pub fn get_dihedral(
        &self,
        atom_types: &(String, String, String, String),
        proper: bool,
        wildcard_allowed: bool,
    ) -> Option<&Vec<DihedralParams>> {
        let (a_variants, b_variants, c_variants, d_variants) = if wildcard_allowed {
            (
                Self::wildcard_variants(&atom_types.0),
                Self::wildcard_variants(&atom_types.1),
                Self::wildcard_variants(&atom_types.2),
                Self::wildcard_variants(&atom_types.3),
            )
        } else {
            (
                vec![atom_types.0.to_string()],
                vec![atom_types.1.to_string()],
                vec![atom_types.2.to_string()],
                vec![atom_types.3.to_string()],
            )
        };

        for a in &a_variants {
            for b in &b_variants {
                for c in &c_variants {
                    for d in &d_variants {
                        for &(k0, k1, k2, k3) in &[(a, b, c, d), (d, c, b, a)] {
                            let key = (k0.clone(), k1.clone(), k2.clone(), k3.clone());
                            let hit = if proper {
                                self.dihedral.get(&key)
                            } else {
                                self.improper.get(&key)
                            };
                            if let Some(h) = hit {
                                return Some(h);
                            }
                        }
                    }
                }
            }
        }
        None
    }
}

/// Helper to deal with spaces in the FF-type col, while still allowing col separation
/// by whitespace.
/// Note: it appears the whitespace is due to the spacing being padded to 2 chars each.
pub(crate) fn get_atom_types(cols: &[&str]) -> (Vec<String>, usize) {
    let mut atom_types = cols[0].to_string();
    let mut col_1_i = 1;

    for col in &cols[1..] {
        if col.parse::<f32>().is_ok() || col_1_i >= 4 {
            break; // This prevents adding negative integers and comments into this.
        }
        if col.starts_with("-") {
            atom_types += col;
            col_1_i += 1;
        }
    }

    let atom_types: Vec<_> = atom_types.split("-").map(|v| v.to_owned()).collect();

    (atom_types, col_1_i)
}
/// Helper to prevent repetition
fn parse_float(v: &str) -> io::Result<f32> {
    v.parse()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, format!("Invalid float: {v}")))
}

/// Parse a lib (partial charge, FF type map), then convert keys to amino acids, and use our peptide-specific
/// type-in-res enum instead of strings.
pub fn parse_lib_peptide(
    text: &str,
) -> io::Result<HashMap<AminoAcidGeneral, Vec<ChargeParamsProtein>>> {
    let parsed = parse_lib(text)?;

    let mut result = HashMap::new();

    for (tag, v) in parsed {
        // This currently fails on alternate variants like ASSH for ASP that's protonated.
        // other examples are LYS/LYN. todo: Impl if you need.
        let Ok(aa) = AminoAcidGeneral::from_str(&tag) else {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Unable to parse amino acid from lib",
            ));
        };

        let mut v_prot = Vec::new();
        for param in v {
            v_prot.push(ChargeParamsProtein {
                type_in_res: AtomTypeInRes::from_str(&param.type_in_res)?,
                ff_type: param.ff_type,
                charge: param.charge,
            });
        }

        result.insert(aa, v_prot);
    }

    Ok(result)
}

/// Parse a lib (partial charge, FF type map), then convert keys to lipid standards.
///
/// Load charge data from Amber's `amino19.lib`, `aminoct12.lib`, `aminont12.lib`, and similar.
/// This provides partial charges for all amino acids, as well as a mapping between atom type in residue,
/// e.g. "C1", "NA" etc, to amber force field type, e.g. "XC".
/// See [Amber RM](https://ambermd.org/doc12/Amber25.pdf), section 13.2: Residue naming conventions,
/// for info on the protenation variants, and their 3-letter identifiers.
pub fn parse_lib_lipid(text: &str) -> io::Result<HashMap<LipidStandard, Vec<ChargeParams>>> {
    let parsed = parse_lib(text)?;

    let mut result = HashMap::new();
    for (tag, v) in parsed {
        let Ok(s) = LipidStandard::from_str(&tag) else {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Unable to parse lipid from lib: {tag}"),
            ));
        };
        result.insert(s, v);
    }

    Ok(result)
}

/// For DNA or RNA.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct NucleotideTemplate {
    pub nt: Nucleotide,
    pub end: ResidueEnd,
}

impl FromStr for NucleotideTemplate {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // OL24.lib uses a "D" prefix we can skip.
        // RNA.lib has no such prefix.
        let s = s.to_owned().replace("D", "");

        let first = match s.chars().next() {
            Some(c) => c,
            None => {
                return Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("Nucleotide in template too short: {s}"),
                ));
            }
        };

        let nt = match first {
            'A' => Nucleotide::A,
            'C' => Nucleotide::C,
            'T' | 'U' => Nucleotide::T,
            'G' => Nucleotide::G,
            _ => {
                return Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("Unrecognized nucleotide in template: {s}"),
                ));
            }
        };
        let end = if s.ends_with('5') {
            ResidueEnd::NTerminus
        } else if s.ends_with('3') {
            ResidueEnd::CTerminus
        } else if s.ends_with('N') {
            ResidueEnd::Hetero // "neutral" / not part of a chain?
        } else {
            ResidueEnd::Internal
        };

        Ok(Self { nt, end })
    }
}

/// Parse a lib (partial charge, FF type map), then convert keys to nucleotides.
/// Note: For RNA, assume T = U here? // todo
pub fn parse_lib_nucleic_acid(
    text: &str,
) -> io::Result<HashMap<NucleotideTemplate, Vec<ChargeParams>>> {
    let parsed = parse_lib(text)?;

    let mut result = HashMap::new();
    for (tag, v) in parsed {
        // OHE is a special cap for RNA on a 5' terminal phosphate.
        // Skip it for now.
        if &tag == "OHE" {
            continue;
        }
        let Ok(s) = NucleotideTemplate::from_str(&tag) else {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Unable to parse nucleotide from lib: {tag}"),
            ));
        };
        result.insert(s, v);
    }

    Ok(result)
}

// todo: This is DRY with the parse_amino_charges fn above. Fix it. Too much repetition for too little diff.
/// A general function for parsing partial charge, and maps of atom type to FF type/name.
/// We post-process these with molecule-type specific keys.
pub fn parse_lib(text: &str) -> io::Result<HashMap<String, Vec<ChargeParams>>> {
    enum Mode {
        Scan,                    // not inside an atoms table
        InAtoms { res: String }, // currently reading atom lines for this residue
    }

    let mut state = Mode::Scan;
    let mut result: HashMap<String, Vec<ChargeParams>> = HashMap::new();

    let lines: Vec<&str> = text.lines().collect();

    for line in lines {
        let ltrim = line.trim_start();

        // Section headers
        if let Some(rest) = ltrim.strip_prefix("!entry.") {
            state = Mode::Scan;

            if let Some((tag, tail)) = rest.split_once('.') {
                // We only care about "<RES>.unit.atoms table"
                if tail.starts_with("unit.atoms table") {
                    state = Mode::InAtoms {
                        res: tag.to_owned(),
                    };

                    result.entry(tag.to_owned()).or_default(); // make sure map key exists
                }
            }
            continue;
        }

        // If inside atoms table, parse data line
        if let Mode::InAtoms { ref res } = state {
            // tables end when we hit an empty line or a comment
            if ltrim.is_empty() || ltrim.starts_with('!') {
                state = Mode::Scan;
                continue;
            }

            let mut tokens = Vec::<&str>::new();
            let mut in_quote = false;
            let mut start = 0usize;
            let bytes = ltrim.as_bytes();
            for (i, &b) in bytes.iter().enumerate() {
                match b {
                    b'"' => in_quote = !in_quote,
                    b' ' | b'\t' if !in_quote => {
                        if start < i {
                            tokens.push(&ltrim[start..i]);
                        }
                        start = i + 1;
                    }
                    _ => {}
                }
            }
            if start < ltrim.len() {
                tokens.push(&ltrim[start..]);
            }

            let type_in_res = tokens[0].trim_matches('"').to_string();
            let ff_type = tokens[1].trim_matches('"').to_string();
            let charge = parse_float(tokens.last().unwrap())?;

            result.get_mut(res).unwrap().push(ChargeParams {
                type_in_res,
                ff_type,
                charge,
            });
        }
    }

    Ok(result)
}

pub fn load_amino_charges(
    path: &Path,
) -> io::Result<HashMap<AminoAcidGeneral, Vec<ChargeParamsProtein>>> {
    let mut file = File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    let data_str: String = String::from_utf8(buffer)
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid UTF8"))?;

    parse_lib_peptide(&data_str)
}

// // todo: C+P from dynamics, for use in GROMACS functionlity. Currently doesn't
// // todo have the ff_q_maps. Do we need those for GROMACS? Should we just move this data struct here
// // todo from dynamics?
// #[derive(Default, Clone, Debug)]
// /// A set of general parameters that aren't molecule-specific. E.g. from GAFF2, OL3, RNA, or amino19.
// /// These are used as a baseline, and in some cases, overridden by molecule-specific parameters.
// pub struct FfParamSet2 {
//     pub peptide: Option<ForceFieldParams>,
//     pub small_mol: Option<ForceFieldParams>,
//     pub dna: Option<ForceFieldParams>,
//     pub rna: Option<ForceFieldParams>,
//     pub lipids: Option<ForceFieldParams>,
//     pub carbohydrates: Option<ForceFieldParams>,
//     // /// In addition to charge, this also contains the mapping of res type to FF type; required to map
//     // /// other parameters to protein atoms. E.g. from `amino19.lib`, and its N and C-terminus variants.
//     // pub peptide_ff_q_map: Option<ProtFfChargeMapSet>,
//     // pub lipid_ff_q_map: Option<LipidFfChargeMap>,
//     // // todo: QC these types; lipid as place holder. See how they parse.
//     // pub dna_ff_q_map: Option<NucleicAcidFfChargeMap>,
//     // pub rna_ff_q_map: Option<NucleicAcidFfChargeMap>,
// }

/// Associate loaded Force field data (e.g. from Amber) into the atom indices used in a specific
/// dynamics sim. This handles combining general and molecule-specific parameter sets, and converting
/// between atom name, and the specific indices of the atoms we're using.
///
/// This code is straightforward if params are available; much of the logic here is related to handling
/// missing parameters.
impl ForceFieldParamsIndexed {
    pub fn new(
        params: &ForceFieldParams,
        atoms: &[AtomDynamics],
        adjacency_list: &[Vec<usize>],
        h_constrained: bool,
    ) -> io::Result<Self> {
        let mut result = Self::default();

        for (i, atom) in atoms.iter().enumerate() {
            let ff_type = &atom.force_field_type;

            // Mass
            if let Some(mass) = params.mass.get(ff_type) {
                result.mass.insert(i, mass.clone());
            } else if ff_type.starts_with("C") {
                match params.mass.get("C") {
                    Some(m) => {
                        result.mass.insert(i, m.clone());
                        println!("Using C fallback mass for {ff_type}");
                    }
                    None => {
                        return Err(ParamError::new(&format!(
                            "\nMD failure: Missing mass params for {ff_type}"
                        )));
                    }
                }
            } else if ff_type.starts_with("N") {
                match params.mass.get("N") {
                    Some(m) => {
                        result.mass.insert(i, m.clone());
                        println!("Using N fallback mass for {ff_type}");
                    }
                    None => {
                        return Err(ParamError::new(&format!(
                            "\nMD failure: Missing mass params for {ff_type}"
                        )));
                    }
                }
            } else if ff_type.starts_with("O") {
                match params.mass.get("O") {
                    Some(m) => {
                        result.mass.insert(i, m.clone());
                        println!("Using O fallback mass for {ff_type}");
                    }
                    None => {
                        return Err(ParamError::new(&format!(
                            "\nMD failure: Missing mass params for {ff_type}"
                        )));
                    }
                }
            } else {
                result.mass.insert(
                    i,
                    MassParams {
                        atom_type: "".to_string(),
                        mass: atom.element.atomic_weight(),
                        comment: None,
                    },
                );

                eprintln!("\nMissing mass params on {atom}; using element default.");
            }

            // Lennard-Jones / van der Waals
            if let Some(vdw) = params.lennard_jones.get(ff_type) {
                result.lennard_jones.insert(i, vdw.clone());
                // If the key is missing for the given FF type in our loaded data, check for certain
                // special cases.
            } else {
                // The mass values for all 4 of these are present in frcmod.ff19sb.
                if ff_type == "2C" || ff_type == "3C" || ff_type == "C8" {
                    result
                        .lennard_jones
                        .insert(i, params.lennard_jones.get("CT").unwrap().clone());
                } else if ff_type == "CO" {
                    result
                        .lennard_jones
                        .insert(i, params.lennard_jones.get("C").unwrap().clone());
                } else if ff_type == "OXT" {
                    result
                        .lennard_jones
                        .insert(i, params.lennard_jones.get("O2").unwrap().clone());
                } else if ff_type.starts_with("N") {
                    result
                        .lennard_jones
                        .insert(i, params.lennard_jones.get("N").unwrap().clone());
                    println!("Using N fallback VdW for {atom}");
                } else if ff_type.starts_with("O") {
                    result
                        .lennard_jones
                        .insert(i, params.lennard_jones.get("O").unwrap().clone());
                    println!("Using O fallback LJ for {atom}");
                } else {
                    println!("\nMissing LJ params for {atom}; setting to 0.");
                    // 0. no interaction.
                    // todo: If this is "CG" etc, fall back to other carbon params instead.
                    result.lennard_jones.insert(
                        i,
                        LjParams {
                            atom_type: "".to_string(),
                            sigma: 0.,
                            eps: 0.,
                        },
                    );
                }

                // return Err(ParamError::new(&format!(
                //     "MD failure: Missing Van der Waals params for {ff_type}"
                // )));
            }
        }

        // Map from serial number fo index, for bonds and the atoms they point ot.
        // let mut index_map = HashMap::new();
        // for (i, atom) in atoms.iter().enumerate() {
        //     index_map.insert(atom.serial_number, i);
        // }

        // Bond lengths.
        // for bond in bonds {
        for (i0, neighbors) in adjacency_list.iter().enumerate() {
            for &i1 in neighbors {
                if i0 >= i1 {
                    continue; // Only add each bond once.
                }

                let type_0 = &atoms[i0].force_field_type;
                let type_1 = &atoms[i1].force_field_type;

                let data = params.get_bond(&(type_0.clone(), type_1.clone()), true);

                let Some(data) = data else {
                    // todo: Consider removing this, and return an error.
                    eprintln!(
                        "\nMissing bond parameters for {type_0}-{type_1} on {} - {}. Using a safe default.",
                        atoms[i0], atoms[i1]
                    );
                    result.bond_stretching.insert(
                        (i0.min(i1), i0.max(i1)),
                        BondStretchingParams {
                            atom_types: (String::new(), String::new()),
                            k_b: 300.,
                            r_0: (atoms[i0].posit - atoms[i1].posit).magnitude(),
                            comment: None,
                        },
                    );
                    continue;
                };
                let data = data.clone();

                // If using fixed hydrogens, don't add these to our bond stretching params;
                // add to a separate hydrogen rigid param variable.
                if h_constrained && (atoms[i0].element == Hydrogen || atoms[i1].element == Hydrogen)
                {
                    // Set up inverse mass using params directly, so we don't have to have
                    // mass loaded into the atom directly yet.
                    let ff_type_0 = &atoms[i0].force_field_type;
                    let ff_type_1 = &atoms[i1].force_field_type;

                    // Mass
                    let (Some(mass_0), Some(mass_1)) =
                        (params.mass.get(ff_type_0), params.mass.get(ff_type_1))
                    else {
                        return Err(ParamError::new(&format!(
                            "MD failure: Missing mass params for {ff_type_0} or {ff_type_1}"
                        )));
                    };

                    // `bonds_topology` exists separately from `bond_params` specifically so we can
                    // account for bonds to H in exclusions.
                    // We will populate inverse mass in a second loop.
                    let inv_mass = 1. / mass_0.mass + 1. / mass_1.mass;

                    result
                        .bond_rigid_constraints
                        .insert((i0, i1), (data.r_0.powi(2), inv_mass));
                    result.bonds_topology.insert((i0, i1));
                    continue;
                }

                result.bond_stretching.insert((i0, i1), data);
                result.bonds_topology.insert((i0, i1));
            }
        }

        // Valence angles: Every connection between 3 atoms bonded linearly.
        for (ctr, neighbors) in adjacency_list.iter().enumerate() {
            if neighbors.len() < 2 {
                continue;
            }
            for (&n0, &n1) in neighbors.iter().tuple_combinations() {
                let type_n0 = &atoms[n0].force_field_type;
                let type_ctr = &atoms[ctr].force_field_type;
                let type_n1 = &atoms[n1].force_field_type;

                let data = params
                    .get_valence_angle(&(type_n0.clone(), type_ctr.clone(), type_n1.clone()), true);

                let Some(data) = data else {
                    // This comes up with the Hydrogen bound to NB in His. I don't know what to make of it.
                    // I'm not sure exactly why I can't find this CR-NB-H angle, but
                    // try subbing the NA variant:
                    // "CR-NA-H     50.0      120.00    AA his,    changed based on NMA nmodes"
                    if (type_n0 == "H" || type_n1 == "H")
                        && type_ctr == "NB"
                        && (type_n0 == "CR"
                            || type_n1 == "CR"
                            || type_n0 == "CV"
                            || type_n1 == "CV")
                    {
                        // todo: Get to the bottom of this.
                        println!(
                            "His HB: Skipping valence angle. For now, inserting a dummy one with no force."
                        );
                        result.angle.insert(
                            (n0, ctr, n1),
                            AngleBendingParams {
                                atom_types: (String::new(), String::new(), String::new()),
                                k: 0.,
                                theta_0: 0.,
                                comment: None,
                            },
                        );
                        continue;

                        // if let Some(v) = params.get_valence_angle(&(
                        //     type_n0.clone(),
                        //     "NA".to_string(),
                        //     type_n1.clone(),
                        // )) {
                        //     println!("Substituting NA for NB in valence angle params");
                        //     result.angle.insert((n0, ctr, n1), v.clone());
                        //     continue;
                        // }
                    }

                    return Err(ParamError::new(&format!(
                        "\nMD failure: Missing valence angle params for {type_n0}-{type_ctr}-{type_n1}. (sns: {} - {} - {})",
                        atoms[n0].serial_number, atoms[ctr].serial_number, atoms[n1].serial_number,
                    )));
                    // }
                };
                let data = data.clone();

                result.angle.insert((n0, ctr, n1), data);
            }
        }

        // Proper and improper dihedral angles.
        let mut seen = HashSet::<(usize, usize, usize, usize)>::new();

        // Proper dihedrals: Atoms 1-2-3-4 bonded linearly
        for (i1, nbr_j) in adjacency_list.iter().enumerate() {
            for &i2 in nbr_j {
                if i1 >= i2 {
                    continue;
                } // handle each j-k bond once

                for &i0 in adjacency_list[i1].iter().filter(|&&x| x != i2) {
                    for &i3 in adjacency_list[i2].iter().filter(|&&x| x != i1) {
                        if i0 == i3 {
                            continue;
                        }

                        // Canonicalise so (i1, i2) is always (min, max)
                        let idx_key = if i1 < i2 {
                            (i0, i1, i2, i3)
                        } else {
                            (i3, i2, i1, i0)
                        };
                        if !seen.insert(idx_key) {
                            continue;
                        }

                        let type_0 = &atoms[i0].force_field_type;
                        let type_1 = &atoms[i1].force_field_type;
                        let type_2 = &atoms[i2].force_field_type;
                        let type_3 = &atoms[i3].force_field_type;

                        let key = (
                            type_0.clone(),
                            type_1.clone(),
                            type_2.clone(),
                            type_3.clone(),
                        );

                        if let Some(dihes) = params.get_dihedral(&key, true, true) {
                            let mut dihes = dihes.clone();

                            for d in &mut dihes {
                                // Divide here; then don't do it during the dynamics run. Optimization.
                                d.barrier_height /= d.divider as f32;
                                d.divider = 1;
                            }
                            result.dihedral.insert(idx_key, dihes);
                        // }
                        // else if allow_missing_dihedral_params {
                        //     // Default of no constraint
                        //     result.dihedral.insert(
                        //         idx_key,
                        //         vec![DihedralParams {
                        //             atom_types: key,
                        //             divider: 1,
                        //             barrier_height: 0.,
                        //             phase: 0.,
                        //             periodicity: 1,
                        //             comment: None,
                        //         }],
                        //     );
                        } else {
                            return Err(ParamError::new(&format!(
                                "\nMD failure: Missing dihedral params for {type_0}-{type_1}-{type_2}-{type_3}. (atom0 sn: {})",
                                atoms[i0].serial_number
                            )));
                        }
                    }
                }
            }
        }

        // Improper dihedrals 2-1-3-4. Atom 3 is the hub, with the other 3 atoms bonded to it.
        // The order of the others in the angle calculation affects the sign of the result.
        // Generally only for planar configs.
        //
        // Note: The sattelites are expected to be in alphabetical order, re their FF types.
        // So, for the hub of "ca" with sattelites of "ca", "ca", and "os", the correct combination
        // to look for in the params is "ca-ca-ca-os"
        for (ctr, satellites) in adjacency_list.iter().enumerate() {
            if satellites.len() < 3 {
                continue;
            }

            // Unique unordered triples of neighbours
            for a in 0..satellites.len() - 2 {
                for b in a + 1..satellites.len() - 1 {
                    for d in b + 1..satellites.len() {
                        let (sat0, sat1, sat2) = (satellites[a], satellites[b], satellites[d]);

                        let idx_key = (sat0, sat1, ctr, sat2); // Order is fixed; no swap
                        if !seen.insert(idx_key) {
                            continue;
                        }

                        let type_0 = &atoms[sat0].force_field_type;
                        let type_1 = &atoms[sat1].force_field_type;
                        let type_ctr = &atoms[ctr].force_field_type;
                        let type_2 = &atoms[sat2].force_field_type;

                        // Sort satellites alphabetically; required to ensure we don't miss combinations.
                        let mut sat_types = [type_0.clone(), type_1.clone(), type_2.clone()];
                        sat_types.sort();

                        let key = (
                            sat_types[0].clone(),
                            sat_types[1].clone(),
                            type_ctr.clone(),
                            sat_types[2].clone(),
                        );

                        // In the case of improper, unlike all other param types, we are allowed to
                        // have missing values. Impropers are only, by Amber convention, for planar
                        // hub and spoke setups, so non-planar ones will be omitted. These may occur,
                        // for example, at ring intersections.
                        if let Some(dihes) = params.get_dihedral(&key, false, true) {
                            let mut dihes = dihes.clone();
                            for d in &mut dihes {
                                // Generally, there is no divisor for impropers, but set it up here
                                // to be more general.
                                d.barrier_height /= d.divider as f32;
                                d.divider = 1;
                            }
                            result.improper.insert(idx_key, dihes);
                        }
                    }
                }
            }
        }

        Ok(result)
    }
}
