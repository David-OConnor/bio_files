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
    io::{self, ErrorKind},
    str::FromStr,
};

use na_seq::{AminoAcid, AtomTypeInRes, AminoAcidGeneral};

/// Data for a MASS entry: e.g. "CT 12.01100" with optional comment
#[derive(Debug, Clone)]
pub struct MassParams {
    pub atom_type: String,
    /// AMU
    pub mass: f32,
    // /// ATPOL: Atomic polarizability (Å^3).
    // /// Intended for Slater–Kirkwood or future polarizable models, and unused by Amber (?)
    // pub polarizability: f32,
    pub comment: Option<String>,
}

impl MassParams {
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        // Allow Skipping ATPOL which we don't currently use, and is sometimes missing.
        if cols.len() < 2 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Mass))",
            ));
        }

        let atom_type = cols[0].to_string();
        let mass = parse_float(cols[1])?;

        // Skipping polarizability; in for historical reasons?
        // // Note: This skips comments where this is a missing col[2].
        // let mut polarizability = 0.0;
        // if cols.len() >= 3 {
        //     polarizability = parse_float(cols[2])?;
        // }

        // Note: This skips comments where this is a missing col[2].
        let mut comment = None;
        if cols.len() >= 4 {
            comment = Some(cols[3..].join(" "));
        }

        Ok(Self {
            atom_type,
            mass,
            // polarizability,
            comment,
        })
    }
}

/// Amber RM 2025, 15.1.6
/// Data for a BOND entry: e.g. "CT-CT  310.0    1.526" with optional comment
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
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Bond).",
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
#[derive(Debug, Clone)]
pub struct AngleBendingParams {
    pub atom_types: (String, String, String),
    /// Force constant. kcal/mol/rad²
    pub k: f32,
    /// In degrees.
    pub theta_0: f32,
    pub comment: Option<String>,
}

impl AngleBendingParams {
    /// Parse a single valence-angle record from a GAFF/Amber `.dat` or `.frcmod` file.
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Angle).",
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

/// Also known as Torsion. Data for both proper, and improper dihedral data.
#[derive(Debug, Clone, Default)]
pub struct DihedralParams {
    /// "ca", "n", "cd", "sh" etc.
    pub atom_types: (String, String, String, String),
    /// Scaling factor used for barrier height.
    /// "Splits the torsion term into individual contributions for
    /// each pair of atoms involved in the torsion."
    /// Always 1 for improper dihedrdals. (Not present in the Amber files for improper)
    pub divider: u8,
    /// Also known as V_n., kcal/mol/rad²
    pub barrier_height: f32,
    /// Equilibrium angle, or phase, in radians. Often 0 or τ/2. Maximum (?) energy
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
    pub periodicity: i8,
    pub comment: Option<String>,
}

impl DihedralParams {
    /// For both FRCMOD, and Dat. For both proper, and improper. Returns `true` if improper.
    pub fn from_line(line: &str) -> io::Result<(Self, bool)> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols. (Dihedral)",
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
        let periodicity = parse_float(cols[col1_i + 2])? as i8;

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
/// Amber RM, section 15.1.7
pub struct VdwParams {
    pub atom_type: String,
    /// σ. derived from Van der Waals radius, Å. Note that Amber parameter files use R_min,
    /// vice σ. This value is σ, which we compute when parsing.
    ///
    /// R_min (i, j) = 0.5(R_min_i + R_min_j)
    /// σ_min (i, j) = 0.5(σ_min_i + σ_min_j)
    pub sigma: f32,
    /// Energy, kcal/mol. (Represents depth of the potential well).
    /// ε(i, j) = sqrt(ε_i * ε_j)
    pub eps: f32,
}

impl VdwParams {
    /// Parse a single van-der-Waals (Lennard-Jones) parameter line.
    pub fn from_line(line: &str) -> io::Result<Self> {
        const SIGMA_FACTOR: f32 = 1.122_462_048_309_373; // 2^(1/6)

        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Vdw).",
            ));
        }

        let atom_type = cols[0].to_string();
        let r_min = parse_float(cols[1])?;
        let eps = parse_float(cols[2])?;

        let sigma = 2.0 * r_min / SIGMA_FACTOR;

        Ok(Self {
            atom_type,
            sigma,
            eps,
        })
    }
}

#[derive(Debug)]
pub struct ChargeParams {
    /// The residue-specific ID. We use this value to map forcefield type
    /// to atoms loaded from mmCIF etc; these will have this `type_in_res`, but not
    /// an Amber ff type. We apply the charge here to the atom based on its `type_in_res` and AA type,
    /// and apply its FF type.
    ///
    /// Once the FF type is applied here, we can map other params, e.g. Vdw, and bonded terms to it.
    pub type_in_res: AtomTypeInRes,
    /// "XC", "H1" etc.
    pub ff_type: String,
    pub charge: f32, // partial charge (q_i)
}

// impl ChargeParams {
//     /// Parse a single van-der-Waals (Lennard-Jones) parameter line.
//     /// Note: This comes from a .lib file; not .dat or .frcmod. This differs
//     /// from the other parsings inthis file.
//     pub fn from_line(line: &str) -> io::Result<Self> {
//         let cols: Vec<_> = line.trim().split_whitespace().collect();
//
//         if cols.len() < 3 {
//             return Err(io::Error::new(
//                 ErrorKind::InvalidData,
//                 "Not enough cols (Charge).",
//             ));
//         }
//
//         let atom_type = cols[0].to_string();
//         let r_min = parse_float(cols[1])?;
//         let eps = parse_float(cols[2])?;
//
//         Ok(Self {
//             atom_type,
//             atom_type,
//             charge
//         })
//     }
// }

/// Top-level dat or frcmod data. We store the name-tuples in fields, vice as HashMaps here,
/// for parsing flexibility.
///
/// Note that we don't include partial charges here, as they come from Mol2 files; this struct
/// is for data parsed from DAT, FRCMOD etc files.
#[derive(Debug, Default)]
pub struct ForceFieldParams {
    pub mass: Vec<MassParams>,
    pub bond: Vec<BondStretchingParams>,
    pub angle: Vec<AngleBendingParams>,
    pub dihedral: Vec<DihedralParams>,
    pub improper: Vec<DihedralParams>,
    // pub partial_charge: Vec<PartialChargeData>,
    pub van_der_waals: Vec<VdwParams>,
    pub remarks: Vec<String>,
}

/// Force field parameters, e.g. from Amber. Similar to that in `bio_files`, but
/// with Hashmap-based keys (of atom-name tuples) for fast look-ups.
///
/// For descriptions of each field and the units used, reference the structs in bio_files, of which
/// this uses internally.
#[derive(Clone, Debug, Default)]
pub struct ForceFieldParamsKeyed {
    pub mass: HashMap<String, MassParams>,
    pub bond: HashMap<(String, String), BondStretchingParams>,
    pub angle: HashMap<(String, String, String), AngleBendingParams>,
    pub dihedral: HashMap<(String, String, String, String), DihedralParams>,
    pub dihedral_improper: HashMap<(String, String, String, String), DihedralParams>,
    pub van_der_waals: HashMap<String, VdwParams>,
    // todo: Partial charges here A/R. Removed here for now, since we get them from
    // Mol2 etc files, so their state is associated that way.
    // pub partial_charges: HashMap<String, f32>,
}

impl ForceFieldParamsKeyed {
    pub fn new(params: &ForceFieldParams) -> Self {
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

        for val in &params.dihedral {
            result.dihedral.insert(val.atom_types.clone(), val.clone());
        }

        for val in &params.improper {
            result
                .dihedral_improper
                .insert(val.atom_types.clone(), val.clone());
        }

        for val in &params.van_der_waals {
            result
                .van_der_waals
                .insert(val.atom_type.clone(), val.clone());
        }

        result
    }

    /// A utility function that handles proper and improper dihedral data,
    /// tries both atom orders, and falls back to wildcard (“X”) matches on
    /// the outer atoms when an exact hit is not found.
    pub fn get_dihedral(
        &self,
        atom_types: &(String, String, String, String),
        proper: bool, // todo: Experimenting.
    ) -> Option<&DihedralParams> {
        let (a, b, c, d) = (
            atom_types.0.clone(),
            atom_types.1.clone(),
            atom_types.2.clone(),
            atom_types.3.clone(),
        );

        // Candidate keys ordered from most-specific to most-generic.
        let mut keys: Vec<(String, String, String, String)> = Vec::with_capacity(8);

        // Exact
        keys.push((a.clone(), b.clone(), c.clone(), d.clone()));
        keys.push((d.clone(), c.clone(), b.clone(), a.clone()));

        // One X
        keys.push(("X".into(), b.clone(), c.clone(), d.clone()));
        keys.push(("X".into(), c.clone(), b.clone(), a.clone()));
        keys.push((a.clone(), b.clone(), c.clone(), "X".into()));
        keys.push((d.clone(), c.clone(), b.clone(), "X".into()));

        // Xs on both ends
        keys.push(("X".into(), b.clone(), c.clone(), "X".into()));
        keys.push(("X".into(), c.clone(), b.clone(), "X".into()));

        // Proper
        if proper {
            for k in &keys {
                if let Some(data) = self.dihedral.get(k) {
                    return Some(data);
                }
            }
        } else {
            // Improper
            for k in &keys {
                if let Some(data) = self.dihedral_improper.get(k) {
                    return Some(data);
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
        Scan,                       // not inside an atoms table
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
                            "Unable to parse AA from lib"),
                        );
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
