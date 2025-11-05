//! Contains parameters used in Amber Forcefields. For details on these formats,
//! see the [Amber Reference Manual](https://ambermd.org/doc12/Amber25.pdf), section
//! 15: Reading and modifying Amber parameter files.
//!
//! Called by both the `dat`, and `frcmod` modules. These formats share line formats, but
//! arrange them in different ways.
//!
//! For ligands, `atom_type` is a "Type 3". For proteins/AAs, we are currently treating it
//! as a type 1, but we're unclear on this.

use std::{
    collections::HashMap,
    fs::File,
    io::{self, ErrorKind, Read},
    path::Path,
    str::FromStr,
};

use lin_alg::f64::Vec3;
use na_seq::{AminoAcidGeneral, AtomTypeInRes, Element};

use crate::{AtomGeneric, BondGeneric, BondType, LipidStandard};


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
    /// Parse a single van-der-Waals (Lennard-Jones) parameter line.
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
pub struct ChargeParams {
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

/// See notes on `ChargeParams`; equivalent here.
#[derive(Clone, Debug)]
pub struct ChargeParamsLipid {
    // pub type_in_res: AtomTypeInLipid,
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

        if atom.len() > 0 {
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
    pub fn get_bond(&self, atom_types: &(String, String)) -> Option<&BondStretchingParams> {
        let a_variants = Self::wildcard_variants(&atom_types.0);
        let b_variants = Self::wildcard_variants(&atom_types.1);

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
    ) -> Option<&AngleBendingParams> {
        let a_variants = Self::wildcard_variants(&atom_types.0);
        let b_variants = Self::wildcard_variants(&atom_types.1);
        let c_variants = Self::wildcard_variants(&atom_types.2);

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
    ) -> Option<&Vec<DihedralParams>> {
        let a_variants = Self::wildcard_variants(&atom_types.0);
        let b_variants = Self::wildcard_variants(&atom_types.1);
        let c_variants = Self::wildcard_variants(&atom_types.2);
        let d_variants = Self::wildcard_variants(&atom_types.3);

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

/// Load charge data from Amber's `amino19.lib`, `aminoct12.lib`, `aminont12.lib`, and similar.
/// This provides partial charges for all amino acids, as well as a mapping between atom type in residue,
/// e.g. "C1", "NA" etc, to amber force field type, e.g. "XC".
/// See [Amber RM](https://ambermd.org/doc12/Amber25.pdf), section 13.2: Residue naming conventions,
/// for info on the protenation variants, and their 3-letter identifiers.
pub fn parse_amino_charges(text: &str) -> io::Result<HashMap<AminoAcidGeneral, Vec<ChargeParams>>> {
    enum Mode {
        Scan,                              // not inside an atoms table
        InAtoms { res: AminoAcidGeneral }, // currently reading atom lines for this residue
    }

    let mut state = Mode::Scan;
    let mut result: HashMap<AminoAcidGeneral, Vec<ChargeParams>> = HashMap::new();

    let lines: Vec<&str> = text.lines().collect();

    for line in lines {
        let ltrim = line.trim_start();

        // Section headers
        if let Some(rest) = ltrim.strip_prefix("!entry.") {
            state = Mode::Scan;

            if let Some((tag, tail)) = rest.split_once('.') {
                // We only care about "<RES>.unit.atoms table"
                if tail.starts_with("unit.atoms table") {
                    // This currently fails on alternate variants like ASSH for ASP that's protonated.
                    // other examples are LYS/LYN. todo: Impl if you need.
                    let Ok(aa) = AminoAcidGeneral::from_str(tag) else {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Unable to parse AA from lib",
                        ));
                    };

                    state = Mode::InAtoms { res: aa };

                    result.entry(aa).or_default(); // make sure map key exists
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
                type_in_res: AtomTypeInRes::from_str(&type_in_res)?,
                ff_type,
                charge,
            });
        }
    }

    Ok(result)
}

// todo: This is DRY with the parse_amino_charges fn above. Fix it. Too much repetition for too little diff.
pub fn parse_lipid_charges(
    text: &str,
) -> io::Result<HashMap<LipidStandard, Vec<ChargeParamsLipid>>> {
    enum Mode {
        Scan,                           // not inside an atoms table
        InAtoms { res: LipidStandard }, // currently reading atom lines for this residue
    }

    let mut state = Mode::Scan;
    let mut result: HashMap<LipidStandard, Vec<ChargeParamsLipid>> = HashMap::new();

    let lines: Vec<&str> = text.lines().collect();

    for line in lines {
        let ltrim = line.trim_start();

        // Section headers
        if let Some(rest) = ltrim.strip_prefix("!entry.") {
            state = Mode::Scan;

            if let Some((tag, tail)) = rest.split_once('.') {
                // We only care about "<RES>.unit.atoms table"
                if tail.starts_with("unit.atoms table") {
                    let Ok(aa) = LipidStandard::from_str(tag) else {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            format!("Unable to parse lipid from lib: {tag}"),
                        ));
                    };

                    state = Mode::InAtoms { res: aa };

                    result.entry(aa).or_default(); // make sure map key exists
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

            result.get_mut(res).unwrap().push(ChargeParamsLipid {
                // type_in_res: AtomTypeInLipid::from_str(&type_in_res)?,
                type_in_res,
                ff_type,
                charge,
            });
        }
    }

    Ok(result)
}

pub fn load_amino_charges(path: &Path) -> io::Result<HashMap<AminoAcidGeneral, Vec<ChargeParams>>> {
    let mut file = File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    let data_str: String = String::from_utf8(buffer)
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid UTF8"))?;

    parse_amino_charges(&data_str)
}

// todo: Consider an amino acid variant of this too!
/// Creates a set of lipid atom and bonds for all items in a `.lib` file, e.g. `lipid21.lib` from Amber.
/// The hashmap key is the lipid name, e.g. "AR", "CHL" etc.
pub fn load_lipid_templates(
    text: &str,
) -> io::Result<HashMap<String, (Vec<AtomGeneric>, Vec<BondGeneric>)>> {
    #[derive(Clone)]
    struct AtomRow {
        name: String,    // per-residue atom name (eg. "C116")
        ff_type: String, // Amber/GAFF atom type (eg. "cD")
        z: u32,          // atomic number
        q: f64,          // partial charge
    }

    #[derive(Default)]
    struct Work {
        atoms: Vec<AtomRow>,
        positions: Vec<Vec3>,
        bonds: Vec<(u32, u32, u32)>, // (a1, a2, flag)
    }

    enum Sec {
        Atoms,
        Positions,
        Connectivity,
    }

    let mut out: HashMap<String, (Vec<AtomGeneric>, Vec<BondGeneric>)> = HashMap::new();
    let mut cur_key: Option<String> = None;
    let mut cur_sec: Option<Sec> = None;
    let mut cur: Work = Work::default();

    let element_from_z = |z: u32| -> Element {
        match z {
            1 => Element::Hydrogen,
            6 => Element::Carbon,
            7 => Element::Nitrogen,
            8 => Element::Oxygen,
            9 => Element::Fluorine,
            15 => Element::Phosphorus,
            16 => Element::Sulfur,
            17 => Element::Chlorine,
            35 => Element::Bromine,
            53 => Element::Iodine,
            _ => Element::Tellurium, // fallback; adjust if you have a dedicated Unknown
        }
    };

    let bond_from_flag = |f: u32| -> BondType {
        match f {
            1 => BondType::Single,
            2 => BondType::Double,
            3 => BondType::Triple,
            // Add/adjust if you encode aromatic or order via other bits
            _ => BondType::Single,
        }
    };

    let mut finalize = |key: Option<String>, work: &mut Work| {
        if let Some(k) = key {
            if !work.atoms.is_empty() {
                let n = work.atoms.len();

                // Ensure positions length
                if work.positions.len() < n {
                    work.positions.extend(
                        std::iter::repeat_with(|| Vec3::new(0.0, 0.0, 0.0))
                            .take(n - work.positions.len()),
                    );
                }

                // Build atoms
                let atoms: Vec<AtomGeneric> = work
                    .atoms
                    .iter()
                    .enumerate()
                    .map(|(i, ar)| AtomGeneric {
                        serial_number: (i as u32) + 1,
                        posit: work.positions[i],
                        element: element_from_z(ar.z),
                        type_in_res: None,
                        type_in_res_lipid: Some(ar.name.clone()),
                        force_field_type: Some(ar.ff_type.clone()),
                        partial_charge: Some(ar.q as f32),
                        hetero: false,
                        occupancy: None,
                        alt_conformation_id: None,
                    })
                    .collect();

                // Build bonds
                let bonds: Vec<BondGeneric> = work
                    .bonds
                    .iter()
                    .map(|&(a1, a2, fl)| BondGeneric {
                        bond_type: bond_from_flag(fl),
                        atom_0_sn: a1,
                        atom_1_sn: a2,
                    })
                    .collect();

                out.insert(k, (atoms, bonds));
            }
            *work = Work::default();
        }
    };

    for line in text.lines() {
        let l = line.trim_start();

        if let Some(rest) = l.strip_prefix("!entry.") {
            // New entry section header
            // Forms we care about:
            // !entry.<NAME>.unit.atoms table ...
            // !entry.<NAME>.unit.positions table ...
            // !entry.<NAME>.unit.connectivity table ...
            // Extract key = between "!entry." and ".unit."
            if let Some(dot_unit_idx) = rest.find(".unit.") {
                let key = &rest[..dot_unit_idx];
                let after_unit = &rest[dot_unit_idx + ".unit.".len()..];

                // If switching to a new NAME (key) and we already had a current key, finalize the previous lipid
                let key_changed = match &cur_key {
                    Some(k) => k != key,
                    None => true,
                };
                if key_changed {
                    finalize(cur_key.take(), &mut cur);
                    cur_key = Some(key.to_string());
                }

                if after_unit.starts_with("atoms table") {
                    cur_sec = Some(Sec::Atoms);
                    continue;
                } else if after_unit.starts_with("positions table") {
                    cur_sec = Some(Sec::Positions);
                    continue;
                } else if after_unit.starts_with("connectivity table") {
                    cur_sec = Some(Sec::Connectivity);
                    continue;
                } else {
                    // Section we don't parse; leave state None so we skip rows
                    cur_sec = None;
                    continue;
                }
            } else {
                // Some other entry; ignore
                cur_sec = None;
                continue;
            }
        }

        match cur_sec {
            Some(Sec::Atoms) => {
                // Expect:  "NAME" "TYPE" 0 1 131073 1 6 0.039100
                // Grab the two quoted fields first.
                let bytes = l.as_bytes();
                let mut qpos = Vec::with_capacity(4);
                for (i, &b) in bytes.iter().enumerate() {
                    if b == b'"' {
                        qpos.push(i);
                    }
                    if qpos.len() == 4 {
                        break;
                    }
                }
                if qpos.len() == 4 {
                    let name = &l[qpos[0] + 1..qpos[1]];
                    let ff_type = &l[qpos[2] + 1..qpos[3]];
                    let tail = &l[qpos[3] + 1..];
                    // Remaining numeric tokens; we need elmnt (5th) and chg (6th) after the two quoted tokens
                    // Layout: typex resx flags seq elmnt chg
                    let nums: Vec<&str> = tail.split_whitespace().collect();
                    if nums.len() >= 6 {
                        let elmnt = nums[4].parse::<u32>().unwrap_or(0);
                        let chg = nums[5].parse::<f64>().unwrap_or(0.0);
                        cur.atoms.push(AtomRow {
                            name: name.to_string(),
                            ff_type: ff_type.to_string(),
                            z: elmnt,
                            q: chg,
                        });
                    }
                }
            }
            Some(Sec::Positions) => {
                // Expect three floats per line
                let mut it = l.split_whitespace();
                if let (Some(x), Some(y), Some(z)) = (it.next(), it.next(), it.next())
                    && let (Ok(xv), Ok(yv), Ok(zv)) =
                        (x.parse::<f64>(), y.parse::<f64>(), z.parse::<f64>())
                {
                    cur.positions.push(Vec3::new(xv, yv, zv));
                }
            }
            Some(Sec::Connectivity) => {
                // Expect three ints: atom1x atom2x flags
                let mut it = l.split_whitespace();
                if let (Some(a1), Some(a2), Some(flg)) = (it.next(), it.next(), it.next())
                    && let (Ok(a1v), Ok(a2v), Ok(fv)) =
                        (a1.parse::<u32>(), a2.parse::<u32>(), flg.parse::<u32>())
                {
                    cur.bonds.push((a1v, a2v, fv));
                }
            }
            None => {}
        }
    }

    // Finalize last lipid
    finalize(cur_key.take(), &mut cur);

    Ok(out)
}
