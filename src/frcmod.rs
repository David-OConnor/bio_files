use std::{
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
    ptr::write,
};

use crate::Mol2;

///! For operating on frcmod files, which describe Amber force fields for small molecules.

/// Data for a MASS entry: e.g. "CT 12.01100" with optional comment
#[derive(Debug, Clone)]
pub struct MassData {
    pub atom_type: String,
    pub mass: f32,
    pub comment: Option<String>,
}

/// Data for a BOND entry: e.g. "CT-CT  310.0    1.526" with optional comment
#[derive(Debug, Clone)]
pub struct BondData {
    pub pair: (String, String),
    pub k: f32,
    pub length: f32,
    pub comment: Option<String>,
}

/// Data for an ANGLE entry: e.g. "CT-CT-CT  63.0    109.5" with optional comment
#[derive(Debug, Clone)]
pub struct AngleData {
    pub triple: (String, String, String),
    pub k: f32,
    pub angle: f32,
    pub comment: Option<String>,
}

/// Data for a DIHEDRAL (proper) entry
#[derive(Debug, Clone)]
pub struct DihedralData {
    pub atom_names: (String, String, String, String),
    pub periodicity: u8,
    pub phase: f32,
    pub k: f32,
    pub _unused: f32,
    pub notes: Option<String>,
    pub penalty_score: f32,
}

/// Data for an IMPROPER entry
#[derive(Debug, Clone)]
pub struct ImproperData {
    pub atom_names: (String, String, String, String),
    pub k: f32,
    pub phase: f32,
    pub periodicity: f32,
    pub comment: Option<String>,
}

/// Top-level frcmod data
#[derive(Debug, Default)]
pub struct FrcmodData {
    pub remarks: Vec<String>,
    pub mass: Vec<MassData>,
    pub bond: Vec<BondData>,
    pub angle: Vec<AngleData>,
    pub dihedral: Vec<DihedralData>,
    pub improper: Vec<ImproperData>,
    // pub nonbond: Vec<NonbondData>, // implement as needed
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

impl FrcmodData {
    /// From a string of a CIF or PDB text file.
    pub fn new(text: &str) -> io::Result<Self> {
        let mut out = FrcmodData::default();
        let mut section = Section::Remark;

        for raw in text.lines() {
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
                    out.remarks.push(line.to_owned());
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
                    out.mass.push(MassData {
                        atom_type: atom.to_string(),
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
                    out.bond.push(BondData {
                        pair: (a1, a2),
                        k,
                        length,
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
                        })?;
                    let comment = parts.next().map(|s| s.to_string());
                    out.angle.push(AngleData {
                        triple: (a1, a2, a3),
                        k,
                        angle,
                        comment,
                    });
                }
                Section::Dihedral => {
                    let tokens: Vec<&str> = raw
                        .splitn(6, char::is_whitespace)
                        .filter(|s| !s.is_empty())
                        .collect();
                    if tokens.len() < 6 {
                        continue;
                    }
                    let names: Vec<&str> = tokens[0].split('-').collect();
                    let periodicity = tokens[1].parse::<u8>().unwrap_or(1);
                    let phase = tokens[2].parse::<f32>().unwrap_or(0.0);
                    let k = tokens[3].parse::<f32>().unwrap_or(0.0);
                    let unused = tokens[4].parse::<f32>().unwrap_or(0.0);
                    let mut notes_raw = tokens[5].trim().to_string();
                    let mut penalty = 0.0;
                    if let Some(idx) = notes_raw.to_lowercase().find("penalty score=") {
                        let note_part = notes_raw[..idx].trim().to_string();
                        let score_part = notes_raw[idx..]
                            .trim_start_matches(|c: char| !c.is_digit(10) && c != '-' && c != '.')
                            .trim();
                        penalty = score_part.parse::<f32>().unwrap_or(0.0);
                        notes_raw = note_part;
                    }
                    let notes = if notes_raw.is_empty() {
                        None
                    } else {
                        Some(notes_raw)
                    };
                    out.dihedral.push(DihedralData {
                        atom_names: (
                            names.get(0).unwrap_or(&"").to_string(),
                            names.get(1).unwrap_or(&"").to_string(),
                            names.get(2).unwrap_or(&"").to_string(),
                            names.get(3).unwrap_or(&"").to_string(),
                        ),
                        periodicity,
                        phase,
                        k,
                        _unused: unused,
                        notes,
                        penalty_score: penalty,
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
                    let names: Vec<&str> = tokens[0].split('-').collect();
                    let k = tokens[1].parse::<f32>().unwrap_or(0.0);
                    let phase = tokens[2].parse::<f32>().unwrap_or(0.0);
                    let per = tokens[3].parse::<f32>().unwrap_or(0.0);
                    let comment = tokens.get(4).map(|s| s.to_string());
                    out.improper.push(ImproperData {
                        atom_names: (
                            names.get(0).unwrap_or(&"").to_string(),
                            names.get(1).unwrap_or(&"").to_string(),
                            names.get(2).unwrap_or(&"").to_string(),
                            names.get(3).unwrap_or(&"").to_string(),
                        ),
                        k,
                        phase,
                        periodicity: per,
                        comment,
                    });
                }
                Section::Nonbond | Section::None => { /* skip or extend */ }
            }
        }
        Ok(out)
    }

    /// Write back to file
    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut f = File::create(path)?;

        // remarks
        for r in &self.remarks {
            writeln!(f, "{}", r)?;
        }
        writeln!(f)?;
        // MASS
        writeln!(f, "MASS")?;
        for m in &self.mass {
            if let Some(c) = &m.comment {
                writeln!(f, "{} {:>10.4} ! {}", m.atom_type, m.mass, c)?;
            } else {
                writeln!(f, "{} {:>10.4}", m.atom_type, m.mass)?;
            }
        }
        writeln!(f)?;
        // BOND
        writeln!(f, "BOND")?;
        for b in &self.bond {
            if let Some(c) = &b.comment {
                writeln!(
                    f,
                    "{}-{} {:>8.3} {:>8.3} {}",
                    b.pair.0, b.pair.1, b.k, b.length, c
                )?;
            } else {
                writeln!(
                    f,
                    "{}-{} {:>8.3} {:>8.3}",
                    b.pair.0, b.pair.1, b.k, b.length
                )?;
            }
        }
        writeln!(f)?;
        // ANGLE
        writeln!(f, "ANGLE")?;
        for a in &self.angle {
            if let Some(c) = &a.comment {
                writeln!(
                    f,
                    "{}-{}-{} {:>8.3} {:>8.3} {}",
                    a.triple.0, a.triple.1, a.triple.2, a.k, a.angle, c
                )?;
            } else {
                writeln!(
                    f,
                    "{}-{}-{} {:>8.3} {:>8.3}",
                    a.triple.0, a.triple.1, a.triple.2, a.k, a.angle
                )?;
            }
        }
        writeln!(f)?;
        // DIHE
        writeln!(f, "DIHE")?;
        for d in &self.dihedral {
            let names = format!(
                "{}-{}-{}-{}",
                d.atom_names.0, d.atom_names.1, d.atom_names.2, d.atom_names.3
            );
            let mut line = format!(
                "{} {:>3} {:>8.3} {:>8.3} {:>8.3}",
                names, d.periodicity, d.phase, d.k, d._unused
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
        // IMPROPER
        writeln!(f, "IMPROPER")?;
        for imp in &self.improper {
            let names = format!(
                "{}-{}-{}-{}",
                imp.atom_names.0, imp.atom_names.1, imp.atom_names.2, imp.atom_names.3
            );
            if let Some(c) = &imp.comment {
                writeln!(
                    f,
                    "{} {:>8.3} {:>8.3} {:>8.3} {}",
                    names, imp.k, imp.phase, imp.periodicity, c
                )?;
            } else {
                writeln!(
                    f,
                    "{} {:>8.3} {:>8.3} {:>8.3}",
                    names, imp.k, imp.phase, imp.periodicity
                )?;
            }
        }
        writeln!(f)?;
        // NONBON placeholder
        writeln!(f, "NONBON")?;
        // extend when nonbond data present
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
