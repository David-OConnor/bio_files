//! For operating on frcmod files, which describe Amber force fields for small molecules.
use std::{
    fs::File,
    io::{self, ErrorKind, Read, Write},
    path::Path,
};

use crate::md_params::{
    AngleBendingParams, BondStretchingParams, DihedralParams, ForceFieldParams,
    ForceFieldParamsVec, LjParams, MassParams,
};

#[derive(Debug, PartialEq)]
enum Section {
    Remark,
    Mass,
    Bond,
    Angle,
    Dihedral,
    Improper,
    Nonbond,
    /// Sections we don't read, e.g. `CMAP` (as in ff19SB), `HBON`, `IPOL`, and `LJEDIT`.
    Other,
}

impl ForceFieldParamsVec {
    /// From a string of a FRCMOD text file.
    pub fn from_frcmod(text: &str) -> io::Result<Self> {
        let mut result = Self::default();

        let lines: Vec<&str> = text.lines().collect();

        let mut section = Section::Remark;

        for line in lines {
            let line = line.trim();
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
                "ANGLE" | "ANGL" => {
                    section = Section::Angle;
                    continue;
                }
                "DIHE" | "DIHEDRAL" => {
                    section = Section::Dihedral;
                    continue;
                }
                "IMPROPER" | "IMPR" => {
                    section = Section::Improper;
                    continue;
                }
                "NONBON" | "NONB" | "NONBOND" => {
                    section = Section::Nonbond;
                    continue;
                }
                "CMAP" | "HBON" | "HBOND" | "IPOL" | "LJEDIT" => {
                    section = Section::Other;
                    continue;
                }
                _ => {}
            }

            match section {
                Section::Remark => {
                    result.remarks.push(line.to_owned());
                }
                Section::Mass => {
                    result.mass.push(MassParams::from_line(line)?);
                }
                Section::Bond => {
                    result.bond.push(BondStretchingParams::from_line(line)?);
                }
                Section::Angle => {
                    result.angle.push(AngleBendingParams::from_line(line)?);
                }
                Section::Dihedral => {
                    result.dihedral.push(DihedralParams::from_line(line)?.0);
                }
                Section::Improper => {
                    result.improper.push(DihedralParams::from_line(line)?.0);
                }
                // Lines in the same format as a `.dat`'s MOD4 RE section: Type, R_min / 2, ε.
                Section::Nonbond => {
                    if !line.starts_with("MOD4") {
                        result.lennard_jones.push(LjParams::from_line(line)?);
                    }
                }
                Section::Other => {}
            }
        }

        // for r in &result.dihedral {
        //     println!("Dihe: {:?}", r);
        // }
        // println!("\n\n\n");
        // for r in &result.improper {
        //     println!("Imp: {:?}", r);
        // }

        Ok(result)
    }

    /// Write to file
    pub fn save_frcmod(&self, path: &Path) -> io::Result<()> {
        let mut f = File::create(path)?;

        for r in &self.remarks {
            writeln!(f, "{r}")?;
        }

        writeln!(f)?;

        writeln!(f, "MASS")?;
        for m in &self.mass {
            if let Some(c) = &m.comment {
                writeln!(f, "{} {:>10.4} ! {}", m.atom_type, m.mass, c)?;
            } else {
                writeln!(f, "{} {:>10.4}", m.atom_type, m.mass)?;
            }
        }
        writeln!(f)?;

        writeln!(f, "BOND")?;
        for b in &self.bond {
            if let Some(c) = &b.comment {
                writeln!(
                    f,
                    "{}-{} {:>8.3} {:>8.3} {}",
                    b.atom_types.0, b.atom_types.1, b.k_b, b.r_0, c
                )?;
            } else {
                writeln!(
                    f,
                    "{}-{} {:>8.3} {:>8.3}",
                    b.atom_types.0, b.atom_types.1, b.k_b, b.r_0
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
                    a.atom_types.0,
                    a.atom_types.1,
                    a.atom_types.2,
                    a.k,
                    a.theta_0.to_degrees(),
                    c
                )?;
            } else {
                writeln!(
                    f,
                    "{}-{}-{} {:>8.3} {:>8.3}",
                    a.atom_types.0,
                    a.atom_types.1,
                    a.atom_types.2,
                    a.k,
                    a.theta_0.to_degrees()
                )?;
            }
        }
        writeln!(f)?;

        writeln!(f, "DIHE")?;
        for d in &self.dihedral {
            let names = format!(
                "{}-{}-{}-{}",
                d.atom_types.0, d.atom_types.1, d.atom_types.2, d.atom_types.3
            );
            let mut line = format!(
                "{} {:>3} {:>8.3} {:>8.3} {:>8.3}",
                names,
                d.divider,
                d.barrier_height,
                d.phase.to_degrees(),
                d.periodicity
            );
            if let Some(n) = &d.comment {
                line.push_str(&format!("  {n}"));
            }
            writeln!(f, "{line}")?;
        }
        writeln!(f)?;

        writeln!(f, "IMPROPER")?;
        for imp in &self.improper {
            let names = format!(
                "{}-{}-{}-{}",
                imp.atom_types.0, imp.atom_types.1, imp.atom_types.2, imp.atom_types.3
            );
            if let Some(c) = &imp.comment {
                writeln!(
                    f,
                    "{} {:>8.3} {:>8.3} {:>8.3} {}",
                    names,
                    imp.barrier_height,
                    imp.phase.to_degrees(),
                    imp.periodicity,
                    c
                )?;
            } else {
                writeln!(
                    f,
                    "{} {:>8.3} {:>8.3} {:>8.3}",
                    names,
                    imp.barrier_height,
                    imp.phase.to_degrees(),
                    imp.periodicity
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

impl ForceFieldParams {
    /// Write as an Amber `.frcmod` file, including Lennard-Jones parameters in its NONBON section.
    pub fn to_frcmod(&self, remark: &str) -> String {
        let lines = self.amber_lines();

        let mut s = String::new();
        for r in remark.lines() {
            s.push_str(r);
            s.push('\n');
        }

        for (name, section) in [
            ("MASS", &lines.mass),
            ("BOND", &lines.bond),
            ("ANGLE", &lines.angle),
            ("DIHE", &lines.dihedral),
            ("IMPROPER", &lines.improper),
            ("NONBON", &lines.lj),
        ] {
            s.push_str(name);
            s.push('\n');

            for line in section {
                s.push_str(line);
                s.push('\n');
            }
            s.push('\n');
        }

        s
    }

    /// See `to_frcmod`.
    pub fn save_frcmod(&self, path: &Path, remark: &str) -> io::Result<()> {
        std::fs::write(path, self.to_frcmod(remark))
    }
}
