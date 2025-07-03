///! For operating on frcmod files, which describe Amber force fields for small molecules.
use std::{
    fs::File,
    io::{self, ErrorKind, Read, Write},
    path::Path,
};

/// Data for a MASS entry: e.g. "CT 12.01100" with optional comment
#[derive(Debug, Clone)]
pub struct MassData {
    pub ff_type: String,
    /// Daltons?
    pub mass: f32,
    pub comment: Option<String>,
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

/// Data for a DIHEDRAL (proper) entry
#[derive(Debug, Clone, Default)]
pub struct DihedralData {
    /// "ca", "n", "cd", "sh" etc.
    pub ff_types: (String, String, String, String),
    /// Aka idivf. 	Scaling factor for barrier height (divide Vn by this)
    pub scaling_factor: u8,
    /// kcal/mol/rad²
    pub barrier_height_vn: f32,
    /// Equilibrium angle, or phase. Degrees. Often 0 or τ/2.
    pub gamma: f32,
    /// An integer, but uses decimals in the file format.
    pub periodicity: i8,
    pub notes: Option<String>,
    pub penalty_score: f32,
}

/// Used to envorce planarity or chirality. For when one atom is out-of-plane.
/// todo: Remove this struct in liu of `DihedralData`?
#[derive(Debug, Clone)]
pub struct ImproperDihedralData {
    pub ff_types: (String, String, String, String),
    /// kcal/mol/rad²
    pub barrier_height_vn: f32,
    /// Equilibrium angle, or phase. Degrees. Often 0 or τ/2.
    pub gamma: f32,
    pub periodicity: i8,
    pub comment: Option<String>,
}

#[derive(Debug, Clone)]
pub struct VdwData {
    pub ff_type: String,
    pub sigma: f32,
    pub eps: f32,
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
    pub improper: Vec<ImproperDihedralData>,
    pub partial_charge: Vec<PartialChargeData>,
    pub van_der_waals: Vec<VdwData>,
    pub remarks: Vec<String>,
}

#[derive(Debug, PartialEq)]
enum Section {
    Remark,
    Mass,
    Bond,
    Angle,
    Dihedral,
    Improper,
    Nonbond,
    None,
}

impl ForceFieldParams {
    /// From a string of a FRCMOD text file.
    pub fn from_frcmod(text: &str) -> io::Result<Self> {
        let mut result = Self::default();

        let lines: Vec<&str> = text.lines().collect();

        let mut section = Section::Remark;

        for raw in lines {
            let line = raw.trim();
            if line.is_empty() {
                continue;
            }
            match line {
                "MASS" => {
                    section = Section::Mass;
                    continue;
                }
                "BOND" => {
                    section = Section::Bond;
                    continue;
                }
                "ANGLE" => {
                    section = Section::Angle;
                    continue;
                }
                "DIHE" | "DIHEDRAL" => {
                    section = Section::Dihedral;
                    continue;
                }
                "IMPROPER" => {
                    section = Section::Improper;
                    continue;
                }
                "NONBON" => {
                    section = Section::Nonbond;
                    continue;
                }
                _ => {}
            }

            match section {
                Section::Remark => {
                    result.remarks.push(line.to_owned());
                }
                Section::Mass => {
                    let mut parts = line
                        .splitn(3, char::is_whitespace)
                        .filter(|s| !s.is_empty());
                    let atom = parts.next().ok_or_else(|| {
                        io::Error::new(ErrorKind::InvalidData, "MASS missing atom type")
                    })?;
                    let mass_s = parts.next().ok_or_else(|| {
                        io::Error::new(ErrorKind::InvalidData, "MASS missing mass")
                    })?;
                    let mass = mass_s
                        .parse::<f32>()
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid mass"))?;
                    let comment = parts
                        .next()
                        .map(|s| s.trim_start_matches('!').trim().to_string());
                    result.mass.push(MassData {
                        ff_type: atom.to_string(),
                        mass,
                        comment,
                    });
                }
                Section::Bond => {
                    let mut parts = line.split_whitespace();
                    let pair = parts.next().ok_or_else(|| {
                        io::Error::new(ErrorKind::InvalidData, "BOND missing pair")
                    })?;
                    let mut atoms = pair.split('-');
                    let a1 = atoms.next().unwrap().to_string();
                    let a2 = atoms.next().unwrap().to_string();

                    let k = parts
                        .next()
                        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "BOND missing k"))?
                        .parse::<f32>()
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid BOND k"))?;

                    let length = parts
                        .next()
                        .ok_or_else(|| {
                            io::Error::new(ErrorKind::InvalidData, "BOND missing length")
                        })?
                        .parse::<f32>()
                        .map_err(|_| {
                            io::Error::new(ErrorKind::InvalidData, "Invalid BOND length")
                        })?;

                    let comment = parts.next().map(|s| s.to_string());
                    result.bond.push(BondData {
                        ff_types: (a1, a2),
                        k,
                        r_0: length,
                        comment,
                    });
                }
                Section::Angle => {
                    let mut parts = line.split_whitespace();
                    let triple = parts.next().ok_or_else(|| {
                        io::Error::new(ErrorKind::InvalidData, "ANGLE missing triple")
                    })?;
                    let mut atoms = triple.split('-');

                    let a1 = atoms.next().unwrap().to_string();
                    let a2 = atoms.next().unwrap().to_string();
                    let a3 = atoms.next().unwrap().to_string();

                    let k = parts
                        .next()
                        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "ANGLE missing k"))?
                        .parse::<f32>()
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid ANGLE k"))?;

                    let angle = parts
                        .next()
                        .ok_or_else(|| {
                            io::Error::new(ErrorKind::InvalidData, "ANGLE missing angle")
                        })?
                        .parse::<f32>()
                        .map_err(|_| {
                            io::Error::new(ErrorKind::InvalidData, "Invalid ANGLE angle")
                        })?.to_radians();

                    let comment = parts.next().map(|s| s.to_string());

                    result.angle.push(AngleData {
                        ff_types: (a1, a2, a3),
                        k,
                        angle,
                        comment,
                    });
                }
                Section::Dihedral => {
                    let cols: Vec<&str> = line.split_whitespace().collect();

                    if cols.len() < 5 {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            format!("Dihedral line is too short: {:?}", line),
                        ));
                    }

                    let ff_types: Vec<&str> = cols[0].split('-').collect();
                    let scaling_factor = cols[1].parse().unwrap_or(1);
                    let barrier_height_vn = cols[2].parse().unwrap_or(0.0);
                    let gamma = cols[3].parse::<f32>().unwrap_or(0.0).to_radians();

                    /// Integer, but often represented as a float, e.g. "1.000" in the files Amber
                    /// generates.
                    let periodicity: f32 = cols[4].parse().unwrap_or(0.0);

                    // let mut notes_raw = cols[5].trim().to_string();

                    let mut notes = String::new();
                    for col in &cols[5..cols.len() - 1] {
                        notes += &format!("{col} ");
                    }
                    //
                    // let mut penalty = 0.0;
                    // if let Some(idx) = notes_raw.to_lowercase().find("penalty score=") {
                    //     let note_part = notes_raw[..idx].trim().to_string();
                    //     let score_part = notes_raw[idx..]
                    //         .trim_start_matches(|c: char| !c.is_digit(10) && c != '-' && c != '.')
                    //         .trim();
                    //     penalty = score_part.parse::<f32>().unwrap_or(0.0);
                    //     notes_raw = note_part;
                    // }
                    // let notes = if notes_raw.is_empty() {
                    //     None
                    // } else {
                    //     Some(notes_raw)
                    // };

                    result.dihedral.push(DihedralData {
                        ff_types: (
                            ff_types.get(0).unwrap_or(&"").to_string(),
                            ff_types.get(1).unwrap_or(&"").to_string(),
                            ff_types.get(2).unwrap_or(&"").to_string(),
                            ff_types.get(3).unwrap_or(&"").to_string(),
                        ),
                        scaling_factor,
                        barrier_height_vn,
                        gamma,
                        periodicity: periodicity as i8,
                        notes: if notes.len() > 0 { Some(notes) } else { None },
                        /// FOr now, we don't parse penalty.
                        penalty_score: 0.,
                    });
                }
                Section::Improper => {
                    let tokens: Vec<&str> = raw
                        .splitn(5, char::is_whitespace)
                        .filter(|s| !s.is_empty())
                        .collect();
                    if tokens.len() < 5 {
                        continue;
                    }
                    let ff_types: Vec<&str> = tokens[0].split('-').collect();
                    let barrier_height_vn = tokens[1].parse::<f32>().unwrap_or(0.0);
                    let gamma = tokens[2].parse::<f32>().unwrap_or(0.0).to_radians();

                    let per = tokens[3].parse::<f32>().unwrap_or(0.0) as i8;
                    let comment = tokens.get(4).map(|s| s.to_string());

                    result.improper.push(ImproperDihedralData {
                        ff_types: (
                            ff_types.get(0).unwrap_or(&"").to_string(),
                            ff_types.get(1).unwrap_or(&"").to_string(),
                            ff_types.get(2).unwrap_or(&"").to_string(),
                            ff_types.get(3).unwrap_or(&"").to_string(),
                        ),
                        barrier_height_vn,
                        gamma,
                        periodicity: per,
                        comment,
                    });
                }
                Section::Nonbond | Section::None => { /* skip or extend */ }
            }
        }
        Ok(result)
    }

    /// Write to file
    pub fn save_frcmod(&self, path: &Path) -> io::Result<()> {
        let mut f = File::create(path)?;

        for r in &self.remarks {
            writeln!(f, "{}", r)?;
        }

        writeln!(f)?;

        writeln!(f, "MASS")?;
        for m in &self.mass {
            if let Some(c) = &m.comment {
                writeln!(f, "{} {:>10.4} ! {}", m.ff_type, m.mass, c)?;
            } else {
                writeln!(f, "{} {:>10.4}", m.ff_type, m.mass)?;
            }
        }
        writeln!(f)?;

        writeln!(f, "BOND")?;
        for b in &self.bond {
            if let Some(c) = &b.comment {
                writeln!(
                    f,
                    "{}-{} {:>8.3} {:>8.3} {}",
                    b.ff_types.0, b.ff_types.1, b.k, b.r_0, c
                )?;
            } else {
                writeln!(
                    f,
                    "{}-{} {:>8.3} {:>8.3}",
                    b.ff_types.0, b.ff_types.1, b.k, b.r_0
                )?;
            }
        }
        writeln!(f)?;

        writeln!(f, "ANGLE")?;
        for a in &self.angle {
            if let Some(c) = &a.comment {
                writeln!(
                    f,
                    "{}-{}-{} {:>8.3} {:>8.3} {}",
                    a.ff_types.0, a.ff_types.1, a.ff_types.2, a.k, a.angle.to_degrees(), c
                )?;
            } else {
                writeln!(
                    f,
                    "{}-{}-{} {:>8.3} {:>8.3}",
                    a.ff_types.0, a.ff_types.1, a.ff_types.2, a.k, a.angle.to_degrees()
                )?;
            }
        }
        writeln!(f)?;

        writeln!(f, "DIHE")?;
        for d in &self.dihedral {
            let names = format!(
                "{}-{}-{}-{}",
                d.ff_types.0, d.ff_types.1, d.ff_types.2, d.ff_types.3
            );
            let mut line = format!(
                "{} {:>3} {:>8.3} {:>8.3} {:>8.3}",
                names, d.scaling_factor, d.barrier_height_vn, d.gamma.to_degrees(), d.periodicity
            );
            if let Some(n) = &d.notes {
                line.push_str(&format!("  {}", n));
            }
            if d.penalty_score != 0.0 {
                line.push_str(&format!(", penalty score={:>5.2}", d.penalty_score));
            }
            writeln!(f, "{}", line)?;
        }
        writeln!(f)?;

        writeln!(f, "IMPROPER")?;
        for imp in &self.improper {
            let names = format!(
                "{}-{}-{}-{}",
                imp.ff_types.0, imp.ff_types.1, imp.ff_types.2, imp.ff_types.3
            );
            if let Some(c) = &imp.comment {
                writeln!(
                    f,
                    "{} {:>8.3} {:>8.3} {:>8.3} {}",
                    names, imp.barrier_height_vn, imp.gamma.to_degrees(), imp.periodicity, c
                )?;
            } else {
                writeln!(
                    f,
                    "{} {:>8.3} {:>8.3} {:>8.3}",
                    names, imp.barrier_height_vn, imp.gamma.to_degrees(), imp.periodicity
                )?;
            }
        }
        writeln!(f)?;

        // todo: Placeholder. A/R.
        writeln!(f, "NONBON")?;

        Ok(())
    }

    /// todo: Sort out the syntax for loading from different sources.
    pub fn load_frcmod(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        let data_str: String = String::from_utf8(buffer)
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid UTF8"))?;

        Self::from_frcmod(&data_str)
    }
}
