//! Residue topology databases (`.rtp`), as used by `pdb2gmx`. E.g. `aminoacids.rtp` in GROMACS
//! force field directories. Each residue lists its atoms, with their types and charges, and its
//! bonds, dihedrals, and impropers. Atom names prefixed with '-' or '+' refer to the previous or
//! next residue.
//!
//! We don't read `[ exclusions ]`, `[ angles ]` or `[ cmap ]`; residues written omit them.

use std::{collections::HashMap, fmt::Write, fs, io, path::Path};

/// Subsections of a residue. Any other section name starts a new residue.
const RES_SECTIONS: [&str; 7] = [
    "atoms",
    "bonds",
    "angles",
    "dihedrals",
    "impropers",
    "exclusions",
    "cmap",
];

#[derive(Clone, Debug, PartialEq)]
pub struct RtpAtom {
    /// Atom name in the residue, e.g. "CA".
    pub name: String,
    pub ff_type: String,
    /// Elementary charge.
    pub charge: f32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct RtpResidue {
    pub atoms: Vec<RtpAtom>,
    pub bonds: Vec<(String, String)>,
    /// 4 atom names, then optional parameters, or a macro name. E.g. OPLS' `dih_LEU_chi1_N_C_C_C`.
    pub dihedrals: Vec<Vec<String>>,
    /// As `dihedrals`.
    pub impropers: Vec<Vec<String>>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Rtp {
    /// The `[ bondedtypes ]` line: default functional forms, and generation settings.
    pub bonded_types: Option<String>,
    pub residues: HashMap<String, RtpResidue>,
}

fn err(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("RTP: {msg}"))
}

impl Rtp {
    pub fn new(text: &str) -> io::Result<Self> {
        let mut result = Self::default();
        let mut residue = String::new();
        let mut section = String::new();

        for line in text.lines() {
            let line = line.split(';').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }

            if line.starts_with('[') {
                let tag = line.trim_matches(['[', ']']).trim();

                if RES_SECTIONS.contains(&tag) {
                    section = tag.to_owned();
                } else {
                    residue = tag.to_owned();
                    section.clear();
                }
                continue;
            }

            let cols: Vec<_> = line.split_whitespace().collect();

            if residue == "bondedtypes" {
                result.bonded_types = Some(cols.join(" "));
                continue;
            }

            if residue.is_empty() || section.is_empty() {
                continue;
            }

            let res = result.residues.entry(residue.clone()).or_default();

            match section.as_str() {
                "atoms" => {
                    if cols.len() < 3 {
                        return Err(err(&format!("Too few columns in atom: {line}")));
                    }

                    res.atoms.push(RtpAtom {
                        name: cols[0].to_owned(),
                        ff_type: cols[1].to_owned(),
                        charge: cols[2]
                            .parse()
                            .map_err(|_| err(&format!("Invalid charge in {residue}: {line}")))?,
                    });
                }
                "bonds" => {
                    if cols.len() < 2 {
                        return Err(err(&format!("Too few columns in bond: {line}")));
                    }
                    res.bonds.push((cols[0].to_owned(), cols[1].to_owned()));
                }
                "dihedrals" => res
                    .dihedrals
                    .push(cols.iter().map(|c| c.to_string()).collect()),
                "impropers" => res
                    .impropers
                    .push(cols.iter().map(|c| c.to_string()).collect()),
                _ => (),
            }
        }

        Ok(result)
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        Self::new(&fs::read_to_string(path)?)
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        fs::write(path, self.to_rtp())
    }

    /// Residues are written in alphabetical order. Each atom is its own charge group.
    pub fn to_rtp(&self) -> String {
        let mut s = String::new();

        if let Some(b) = &self.bonded_types {
            s.push_str("[ bondedtypes ]\n");
            s.push_str(
                "; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih\n",
            );
            writeln!(s, " {b}\n").unwrap();
        }

        let mut names: Vec<_> = self.residues.keys().collect();
        names.sort_unstable();

        for name in names {
            let res = &self.residues[name];

            writeln!(s, "[ {name} ]").unwrap();

            s.push_str(" [ atoms ]\n");
            for (i, a) in res.atoms.iter().enumerate() {
                writeln!(
                    s,
                    "  {:>6} {:<10} {:>9.5} {:>4}",
                    a.name,
                    a.ff_type,
                    a.charge,
                    i + 1
                )
                .unwrap();
            }

            if !res.bonds.is_empty() {
                s.push_str(" [ bonds ]\n");
                for (a, b) in &res.bonds {
                    writeln!(s, "  {a:>6} {b:>6}").unwrap();
                }
            }

            for (section, rows) in [("dihedrals", &res.dihedrals), ("impropers", &res.impropers)] {
                if rows.is_empty() {
                    continue;
                }

                writeln!(s, " [ {section} ]").unwrap();
                for row in rows {
                    writeln!(s, "  {}", row.join("  ")).unwrap();
                }
            }

            s.push('\n');
        }

        s
    }
}
