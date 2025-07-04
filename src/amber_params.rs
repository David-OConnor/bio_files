//! Contains parameters used in Amber Forcefields.
//! Called by both the `dat`, and `frcmod` modules. These formats share line formats, but
//! arrange them in different ways.
//!
//! Note that this is a bit fragile, but there may not be a better way given the choice to A: In .dat
//! files not divide into explicit sections, and B: Using a whitespace-separated col format, with white
//! space in some of the cols.

use std::{
    collections::HashMap,
    io::{self, ErrorKind},
};

/// Data for a MASS entry: e.g. "CT 12.01100" with optional comment
#[derive(Debug, Clone)]
pub struct MassData {
    pub ff_type: String,
    /// Daltons?
    pub mass: f32,
    pub comment: Option<String>,
}

impl MassData {
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Mass))",
            ));
        }

        let ff_type = cols[0].to_string();
        let mass = cols[1]
            .parse()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid mass"))?;

        // todo: Skipping the col after mass for now.

        let comment = cols[3..].join(" ");

        Ok(Self {
            ff_type,
            mass,
            comment: Some(comment),
        })
    }
}

/// Data for a BOND entry: e.g. "CT-CT  310.0    1.526" with optional comment
#[derive(Debug, Clone)]
pub struct BondData {
    pub ff_types: (String, String),
    /// Force constant. (Similar to a spring constant). kcal/mol/Å²
    pub k: f32,
    /// Equilibrium bond length. Å
    pub r_0: f32,
    pub comment: Option<String>,
}

impl BondData {
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Bond).",
            ));
        }

        let (ff_types, col1_i) = get_ff_types(&cols);
        let ff_types = (ff_types[0].to_owned(), ff_types[1].to_owned());

        let k = parse_float(cols[col1_i])?;
        let r_0 = parse_float(cols[col1_i + 1])?;

        // We ignore the remaining cols for now: Source, # of ref geometries used to fit,
        // and RMS deviation of the fit.

        let mut comment = None;
        if cols.len() >= col1_i + 2 {
            comment = Some(cols[col1_i + 2..].join(" "));
        }

        Ok(Self {
            ff_types,
            k,
            r_0,
            comment,
        })
    }
}

/// Data for an ANGLE entry: e.g. "CT-CT-CT  63.0    109.5" with optional comment
#[derive(Debug, Clone)]
pub struct AngleData {
    pub ff_types: (String, String, String),
    /// Force constant. kcal/mol/rad²
    pub k: f32,
    /// In degrees.
    pub angle: f32,
    pub comment: Option<String>,
}

impl AngleData {
    /// Parse a single valence-angle record from a GAFF/Amber `.dat` or `.frcmod` file.
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Angle).",
            ));
        }

        let (ff_types, col1_i) = get_ff_types(&cols);
        let ff_types = (
            ff_types[0].to_owned(),
            ff_types[1].to_owned(),
            ff_types[2].to_owned(),
        );

        let k = parse_float(cols[col1_i])?;
        let angle = parse_float(cols[col1_i + 1])?;

        // We ignore the remaining cols for now: Source, # of ref geometries used to fit,
        // and RMS deviation of the fit.

        let mut comment = None;
        if cols.len() >= col1_i + 2 {
            comment = Some(cols[col1_i + 2..].join(" "));
        }

        Ok(Self {
            ff_types,
            k,
            angle,
            comment,
        })
    }
}

/// Data for both proper, and improper dihedral data
#[derive(Debug, Clone, Default)]
pub struct DihedralData {
    /// "ca", "n", "cd", "sh" etc.
    pub ff_types: (String, String, String, String),
    /// Aka idivf. 	Scaling factor for barrier height (divide Vn by this).
    /// Always 1 for improper dihedrdals.
    pub integer_divisor: u8,
    /// AKA PK-V_n., kcal/mol/rad²
    pub barrier_height_vn: f32,
    /// Equilibrium angle, or phase. Degrees. Often 0 or τ/2.
    pub phase: f32,
    /// An integer, but uses decimals in the file format.
    pub periodicity: i8,
    pub comment: Option<String>,
}

impl DihedralData {
    /// For both FRCMOD, and Dat. For both proper, and improper. Returns `true` if improper.
    pub fn from_line(line: &str) -> io::Result<(Self, bool)> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols. (Dihedral)",
            ));
        }

        let (ff_types, mut col1_i) = get_ff_types(&cols);
        let ff_types = (
            ff_types[0].to_owned(),
            ff_types[1].to_owned(),
            ff_types[2].to_owned(),
            ff_types[3].to_owned(),
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
        let phase = parse_float(cols[col1_i + 1])?;
        let periodicity = parse_float(cols[col1_i + 2])? as i8;

        // We ignore the remaining cols for now: Source, # of ref geometries used to fit,
        // and RMS deviation of the fit.

        let mut comment = None;
        if cols.len() >= col1_i + 3 {
            comment = Some(cols[col1_i + 3..].join(" "));
        }

        Ok((
            Self {
                ff_types,
                integer_divisor,
                barrier_height_vn,
                phase,
                periodicity,
                comment,
            },
            improper,
        ))
    }
}

#[derive(Debug, Clone)]
pub struct VdwData {
    pub ff_type: String,
    pub r_star: f32,
    pub eps: f32,
}

impl VdwData {
    /// Parse a single van-der-Waals (Lennard-Jones) parameter line.
    pub fn from_line(line: &str) -> io::Result<Self> {
        let cols: Vec<_> = line.trim().split_whitespace().collect();

        if cols.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough cols (Vdw).",
            ));
        }

        let ff_type = cols[0].to_string();

        let r_star = parse_float(cols[1])?;
        let eps = parse_float(cols[2])?;

        Ok(Self {
            ff_type,
            r_star,
            eps,
        })
    }
}

#[derive(Debug, Clone)]
pub struct PartialChargeData {
    pub ff_type: String,
    pub charge: f32,
    pub comment: Option<String>,
}

/// Top-level dat or frcmod data. We store the name-tuples in fields, vice as HashMaps here,
/// for parsing flexibility.
#[derive(Debug, Default)]
pub struct ForceFieldParams {
    pub mass: Vec<MassData>,
    pub bond: Vec<BondData>,
    pub angle: Vec<AngleData>,
    pub dihedral: Vec<DihedralData>,
    pub improper: Vec<DihedralData>,
    pub partial_charge: Vec<PartialChargeData>,
    pub van_der_waals: Vec<VdwData>,
    pub remarks: Vec<String>,
}

/// Force field parameters, e.g. from Amber. Similar to that in `bio_files`, but
/// with Hashmap-based keys (of atom-name tuples) for fast look-ups.
///
/// For descriptions of each field and the units used, reference the structs in bio_files, of which
/// this uses internally.
#[derive(Clone, Debug, Default)]
pub struct ForceFieldParamsKeyed {
    pub mass: HashMap<String, MassData>,
    pub bond: HashMap<(String, String), BondData>,
    pub angle: HashMap<(String, String, String), AngleData>,
    pub dihedral: HashMap<(String, String, String, String), DihedralData>,
    pub dihedral_improper: HashMap<(String, String, String, String), DihedralData>,
    pub van_der_waals: HashMap<String, VdwData>,
    // todo: Partial charges here A/R. Note that we also assign them to the atom directly.
    pub partial_charges: HashMap<String, f32>,
}

impl ForceFieldParamsKeyed {
    pub fn new(params: &ForceFieldParams) -> Self {
        let mut result = Self::default();

        for val in &params.mass {
            result.mass.insert(val.ff_type.clone(), val.clone());
        }

        for val in &params.bond {
            result.bond.insert(val.ff_types.clone(), val.clone());
        }

        for val in &params.angle {
            result.angle.insert(val.ff_types.clone(), val.clone());
        }

        for val in &params.dihedral {
            result.dihedral.insert(val.ff_types.clone(), val.clone());
        }

        for val in &params.improper {
            result
                .dihedral_improper
                .insert(val.ff_types.clone(), val.clone());
        }

        for val in &params.van_der_waals {
            result
                .van_der_waals
                .insert(val.ff_type.clone(), val.clone());
        }

        result
    }
}

/// Helper to deal with spaces in the FF-type col, while still allowing col separation
/// by whitespace.
pub(crate) fn get_ff_types(cols: &[&str]) -> (Vec<String>, usize) {
    let mut ff_types = cols[0].to_string();
    let mut col_1_i = 1;

    for col in &cols[1..] {
        if col.parse::<f32>().is_ok() || col_1_i >= 4 {
            break; // This prevents adding negative integers and comments into this.
        }
        if col.starts_with("-") {
            ff_types += col;
            col_1_i += 1;
        }
    }

    let ff_types: Vec<_> = ff_types.split("-").map(|v| v.to_owned()).collect();

    (ff_types, col_1_i)
}
/// Helper to prevent repetition
fn parse_float(v: &str) -> io::Result<f32> {
    v.parse()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, format!("Invalid float: {v}")))
}
