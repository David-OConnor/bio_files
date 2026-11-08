//! Force field parameter tables, as in GROMACS force field directories: `[ atomtypes ]` (e.g.
//! `ffnonbonded.itp`), and `[ bondtypes ]`, `[ constrainttypes ]`, `[ angletypes ]`,
//! `[ dihedraltypes ]`, and `[ pairtypes ]` (e.g. `ffbonded.itp`), along with `#define`d parameter
//! macros.
//!
//! Unlike `top_parse`, which resolves each term of a full topology, this keeps the tables as
//! written, in file order, so they can be viewed, edited, and written back. (GROMACS' dihedral
//! matching depends on the order.) Units are GROMACS': nm, kJ/mol, degrees, amu, and elementary
//! charge.
//!
//! We evaluate conditional blocks (`#ifdef` etc.) with only the macros defined in the files
//! themselves, and write the result. E.g. OPLS' `HEAVY_H` variants are dropped.

use std::{collections::HashSet, fmt::Write, fs, io, path::Path};

use super::{
    top_parse::parse_atom_type,
    topology::{ANG_TO_NM, BOND_K_FACTOR, KCAL_TO_KJ, atomic_number_from_mass},
};
use crate::md_params::ForceFieldParams;

/// We don't follow `#include`s deeper than this. Guards against cycles.
const MAX_INCLUDE_DEPTH: usize = 8;

/// An `[ atomtypes ]` entry.
#[derive(Clone, Debug, PartialEq)]
pub struct ItpAtomType {
    pub name: String,
    /// The type used to look up bonded parameters. Usually the same as `name`, but not in OPLS.
    pub bond_type: String,
    pub atomic_number: Option<u32>,
    /// amu
    pub mass: f64,
    /// Elementary charge
    pub charge: f64,
    /// A: atom. D, S, or V: virtual site/dummy.
    pub ptype: String,
    /// σ (nm) for combination rules 2 and 3; C6 for rule 1.
    pub v: f64,
    /// ε (kJ/mol) for combination rules 2 and 3; C12 for rule 1.
    pub w: f64,
    pub comment: Option<String>,
}

/// A line of a `[ *types ]` table, other than `[ atomtypes ]`. Parameter meanings depend on the
/// section and `funct`, as in the GROMACS reference manual. E.g.:
/// - Harmonic bonds (funct 1): b0 (nm), kb (kJ/mol/nm²)
/// - Harmonic angles (funct 1): θ0 (degrees), kθ (kJ/mol/rad²)
/// - Periodic dihedrals (funct 1, 4, 9): φs (degrees), k (kJ/mol), multiplicity
/// - Ryckaert-Bellemans dihedrals (funct 3): C0 - C5 (kJ/mol)
#[derive(Clone, Debug, PartialEq)]
pub struct ItpTypeEntry {
    /// Bond types; atom types for `[ pairtypes ]`. "X" is a wildcard. Dihedral types may list 2:
    /// the central atoms of a proper, or the outer atoms of an improper.
    pub types: Vec<String>,
    pub funct: u8,
    pub params: Vec<f64>,
    pub comment: Option<String>,
}

/// A `#define`d macro, which GROMACS substitutes token-wise. E.g. OPLS' residue-specific
/// dihedrals, which its `.rtp` residue templates reference by name.
#[derive(Clone, Debug, PartialEq)]
pub struct ItpDefine {
    pub name: String,
    pub value: String,
    pub comment: Option<String>,
}

/// The parameter tables of a GROMACS force field. See the module docs.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GromacsFfParams {
    pub atom_types: Vec<ItpAtomType>,
    pub bond_types: Vec<ItpTypeEntry>,
    pub constraint_types: Vec<ItpTypeEntry>,
    pub angle_types: Vec<ItpTypeEntry>,
    /// Propers and impropers, in file order.
    pub dihedral_types: Vec<ItpTypeEntry>,
    pub pair_types: Vec<ItpTypeEntry>,
    pub defines: Vec<ItpDefine>,
}

fn err(msg: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("GROMACS parameters: {msg}"),
    )
}

/// A `[ *types ]` line: `n_types` types, funct, then parameters.
fn type_entry(cols: &[&str], n_types: usize, comment: Option<String>) -> io::Result<ItpTypeEntry> {
    if cols.len() < n_types + 1 {
        return Err(err(&format!("Too few columns: {}", cols.join(" "))));
    }

    let funct = cols[n_types]
        .parse()
        .map_err(|_| err(&format!("Invalid funct: {}", cols.join(" "))))?;

    let mut params = Vec::with_capacity(cols.len() - n_types - 1);
    for v in &cols[n_types + 1..] {
        params.push(
            v.parse()
                .map_err(|_| err(&format!("Invalid value {v:?} in: {}", cols.join(" "))))?,
        );
    }

    Ok(ItpTypeEntry {
        types: cols[..n_types].iter().map(|s| (*s).to_owned()).collect(),
        funct,
        params,
        comment,
    })
}

/// Dihedral types list either 4 types, or 2. This matches `top_parse`'s rule.
fn dihedral_n_types(cols: &[&str]) -> usize {
    let is_funct = |c: Option<&&str>| c.is_some_and(|c| c.parse::<u8>().is_ok());

    if cols.len() >= 5 && is_funct(cols.get(4)) && !is_funct(cols.get(2)) {
        4
    } else {
        2
    }
}

/// Writes `line`, with its comment, if present.
fn push_line(s: &mut String, line: &str, comment: &Option<String>) {
    s.push_str(line.trim_end());

    if let Some(c) = comment
        && !c.trim().is_empty()
    {
        s.push_str(" ; ");
        s.push_str(c.trim());
    }
    s.push('\n');
}

/// Amber parameters are f32. We convert using their shortest representation, so these don't
/// pick up spurious digits.
fn widen(v: f32) -> f64 {
    v.to_string().parse().unwrap_or(v as f64)
}

/// Degrees, converted from radians, to 4 decimal places.
fn widen_deg(rad: f32) -> f64 {
    widen((rad.to_degrees() * 1e4).round() / 1e4)
}

impl GromacsFfParams {
    /// Parse from text. `#include`s are skipped; use `load` to follow them.
    pub fn new(text: &str) -> io::Result<Self> {
        let mut result = Self::default();
        result.parse(text, None, 0)?;

        Ok(result)
    }

    /// Load from a file, following `#include`s relative to it. E.g. a force field's
    /// `forcefield.itp`, or `ffbonded.itp`.
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut result = Self::default();
        result.parse_file(path, 0)?;

        Ok(result)
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        fs::write(path, self.to_itp())
    }

    pub fn define(&self, name: &str) -> Option<&ItpDefine> {
        self.defines.iter().find(|d| d.name == name)
    }

    fn parse_file(&mut self, path: &Path, depth: usize) -> io::Result<()> {
        if depth > MAX_INCLUDE_DEPTH {
            return Err(err(&format!(
                "Includes nested too deeply at {}",
                path.display()
            )));
        }

        let text = fs::read_to_string(path)
            .map_err(|e| err(&format!("Unable to open {}: {e}", path.display())))?;

        self.parse(&text, path.parent(), depth)
    }

    fn parse(&mut self, text: &str, dir: Option<&Path>, depth: usize) -> io::Result<()> {
        // Each level: (this branch is active, some branch of this #if has been taken).
        let mut conditions: Vec<(bool, bool)> = Vec::new();
        let mut section = String::new();
        // A comment on its own line, directly above. We attach these to `#define`s, which
        // are commented this way, e.g. in OPLS.
        let mut comment_above = None;

        let mut pending = String::new();
        for raw in text.lines() {
            // Join continuation lines.
            if let Some(stripped) = raw.strip_suffix('\\') {
                pending.push_str(stripped);
                pending.push(' ');
                continue;
            }
            pending.push_str(raw);
            let line = std::mem::take(&mut pending);

            let (content, comment) = match line.split_once(';') {
                Some((content, comment)) => {
                    let comment = comment.trim();
                    (
                        content.trim(),
                        (!comment.is_empty()).then(|| comment.to_owned()),
                    )
                }
                None => (line.trim(), None),
            };

            let active = conditions.iter().all(|(a, _)| *a);

            if content.is_empty() {
                comment_above = if active { comment } else { None };
                continue;
            }

            if let Some(directive) = content.strip_prefix('#') {
                let directive = directive.trim();
                let (name, arg) = match directive.split_once(char::is_whitespace) {
                    Some((name, arg)) => (name, arg.trim()),
                    None => (directive, ""),
                };

                match name {
                    "ifdef" | "ifndef" => {
                        let defined = self.define(arg).is_some();
                        let take = if name == "ifdef" { defined } else { !defined };
                        conditions.push((take, take));
                    }
                    "else" => {
                        let (_, taken) = conditions
                            .pop()
                            .ok_or_else(|| err("#else without #ifdef"))?;
                        conditions.push((!taken, true));
                    }
                    "endif" => {
                        conditions
                            .pop()
                            .ok_or_else(|| err("#endif without #ifdef"))?;
                    }
                    _ if !active => {}
                    "define" => {
                        let (key, value) = match arg.split_once(char::is_whitespace) {
                            Some((key, value)) => (key, value.trim()),
                            None => (arg, ""),
                        };

                        let define = ItpDefine {
                            name: key.to_owned(),
                            value: value.to_owned(),
                            comment: comment.or(comment_above.take()),
                        };

                        match self.defines.iter_mut().find(|d| d.name == key) {
                            Some(existing) => *existing = define,
                            None => self.defines.push(define),
                        }
                    }
                    "undef" => self.defines.retain(|d| d.name != arg),
                    "include" => {
                        // Without a directory, e.g. when parsing text, we skip these.
                        if let Some(dir) = dir {
                            let file = arg.trim_matches(|c| c == '"' || c == '<' || c == '>');
                            self.parse_file(&dir.join(file), depth + 1)?;
                        }
                    }
                    _ => {}
                }

                comment_above = None;
                continue;
            }

            comment_above = None;
            if !active {
                continue;
            }

            if content.starts_with('[') {
                section = content
                    .trim_start_matches('[')
                    .trim_end_matches(']')
                    .trim()
                    .to_lowercase();
                continue;
            }

            // Substitute macros, as GROMACS does; e.g. parameters given by name.
            let expanded: Vec<String> = content
                .split_whitespace()
                .map(|tok| match self.define(tok) {
                    Some(d) => d.value.clone(),
                    None => tok.to_owned(),
                })
                .collect();
            let expanded = expanded.join(" ");
            let cols: Vec<&str> = expanded.split_whitespace().collect();

            match section.as_str() {
                "atomtypes" => {
                    let t = parse_atom_type(&cols)?;
                    let entry = ItpAtomType {
                        name: t.name,
                        bond_type: t.bond_type,
                        atomic_number: t.atomic_number,
                        mass: t.mass,
                        charge: t.charge,
                        ptype: t.ptype,
                        v: t.v,
                        w: t.w,
                        comment,
                    };

                    // As in GROMACS, a repeated type overrides the earlier one.
                    match self.atom_types.iter_mut().find(|a| a.name == entry.name) {
                        Some(existing) => *existing = entry,
                        None => self.atom_types.push(entry),
                    }
                }
                "bondtypes" => self.bond_types.push(type_entry(&cols, 2, comment)?),
                "constrainttypes" => self.constraint_types.push(type_entry(&cols, 2, comment)?),
                "angletypes" => self.angle_types.push(type_entry(&cols, 3, comment)?),
                "pairtypes" => self.pair_types.push(type_entry(&cols, 2, comment)?),
                "dihedraltypes" => {
                    let n = dihedral_n_types(&cols);
                    self.dihedral_types.push(type_entry(&cols, n, comment)?);
                }
                // E.g. `[ defaults ]`, or the molecule sections of a full topology.
                _ => {}
            }
        }

        if !conditions.is_empty() {
            return Err(err("Unterminated #ifdef"));
        }

        Ok(())
    }

    /// Write as GROMACS `.itp` text: macros first, then each non-empty table.
    pub fn to_itp(&self) -> String {
        let mut s = String::new();

        if !self.defines.is_empty() {
            for d in &self.defines {
                if let Some(c) = &d.comment {
                    writeln!(s, "; {c}").unwrap();
                }
                writeln!(s, "#define {:<32} {}", d.name, d.value).unwrap();
            }
            s.push('\n');
        }

        if !self.atom_types.is_empty() {
            s.push_str("[ atomtypes ]\n");
            s.push_str(
                "; name  bond_type  at.num  mass  charge  ptype  sigma (nm)  epsilon (kJ/mol)\n",
            );

            for t in &self.atom_types {
                let at_num = t.atomic_number.map(|n| n.to_string()).unwrap_or_default();

                let line = format!(
                    " {:<10} {:<6} {:>3} {:>10} {:>8} {:>2} {:>12} {:>12}",
                    t.name, t.bond_type, at_num, t.mass, t.charge, t.ptype, t.v, t.w,
                );
                push_line(&mut s, &line, &t.comment);
            }
            s.push('\n');
        }

        for (name, entries) in [
            ("bondtypes", &self.bond_types),
            ("constrainttypes", &self.constraint_types),
            ("angletypes", &self.angle_types),
            ("dihedraltypes", &self.dihedral_types),
            ("pairtypes", &self.pair_types),
        ] {
            if entries.is_empty() {
                continue;
            }

            writeln!(s, "[ {name} ]").unwrap();

            for e in entries {
                let mut line = String::from(" ");
                for t in &e.types {
                    write!(line, "{t:<6} ").unwrap();
                }

                write!(line, "{:>2}", e.funct).unwrap();
                for p in &e.params {
                    write!(line, " {p:>12}").unwrap();
                }

                push_line(&mut s, &line, &e.comment);
            }
            s.push('\n');
        }

        s
    }

    /// Convert Amber-format parameters, e.g. from a `.dat` or `.frcmod` file. Each Amber type
    /// becomes both an atom type, and a bond type. Charges are 0; Amber assigns them per atom.
    /// Use these with `[ defaults ]` of `1 2 yes 0.5 0.8333`.
    pub fn from_amber(params: &ForceFieldParams) -> Self {
        let mut result = Self::default();

        let names: HashSet<&String> = params
            .mass
            .keys()
            .chain(params.lennard_jones.keys())
            .collect();
        let mut names: Vec<_> = names.into_iter().collect();
        names.sort_unstable();

        for name in names {
            let mass = params.mass.get(name);
            let lj = params.lennard_jones.get(name);
            let m = mass.map_or(0., |m| m.mass);

            // Virtual sites, e.g. water's extra points, have no mass.
            let atomic_number = if m < 0.5 {
                0
            } else {
                atomic_number_from_mass(m) as u32
            };

            result.atom_types.push(ItpAtomType {
                name: name.clone(),
                bond_type: name.clone(),
                atomic_number: Some(atomic_number),
                mass: widen(m),
                charge: 0.,
                ptype: if m < 0.5 { "D" } else { "A" }.to_owned(),
                v: lj.map_or(0., |lj| widen(lj.sigma * ANG_TO_NM)),
                w: lj.map_or(0., |lj| widen(lj.eps * KCAL_TO_KJ)),
                comment: mass.and_then(|m| m.comment.clone()),
            });
        }

        // Note: Amber's harmonic force constants omit the ½ that GROMACS' include.
        let mut bonds: Vec<_> = params.bond.iter().collect();
        bonds.sort_unstable_by(|a, b| a.0.cmp(b.0));
        for ((t0, t1), p) in bonds {
            result.bond_types.push(ItpTypeEntry {
                types: vec![t0.clone(), t1.clone()],
                funct: 1,
                params: vec![widen(p.r_0 * ANG_TO_NM), widen(2. * p.k_b * BOND_K_FACTOR)],
                comment: p.comment.clone(),
            });
        }

        let mut angles: Vec<_> = params.angle.iter().collect();
        angles.sort_unstable_by(|a, b| a.0.cmp(b.0));
        for ((t0, t1, t2), p) in angles {
            result.angle_types.push(ItpTypeEntry {
                types: vec![t0.clone(), t1.clone(), t2.clone()],
                funct: 1,
                params: vec![widen_deg(p.theta_0), widen(2. * p.k * KCAL_TO_KJ)],
                comment: p.comment.clone(),
            });
        }

        // Propers are funct 9, so several terms for the same types (consecutive lines) all
        // apply. Impropers are funct 4.
        for (map, funct) in [(&params.dihedral, 9), (&params.improper, 4)] {
            let mut dihedrals: Vec<_> = map.iter().collect();
            dihedrals.sort_unstable_by(|a, b| a.0.cmp(b.0));

            for ((t0, t1, t2, t3), terms) in dihedrals {
                for p in terms {
                    result.dihedral_types.push(ItpTypeEntry {
                        types: vec![t0.clone(), t1.clone(), t2.clone(), t3.clone()],
                        funct,
                        params: vec![
                            widen_deg(p.phase),
                            widen(p.barrier_height / p.divider.max(1) as f32 * KCAL_TO_KJ),
                            p.periodicity as f64,
                        ],
                        comment: p.comment.clone(),
                    });
                }
            }
        }

        result
    }
}
