//! Reads GROMACS topologies (`.top`, `.itp`), as written by `pdb2gmx`, ParmEd, ACPYPE,
//! CHARMM-GUI, and OpenFF Interchange.
//!
//! We run the preprocessor (`#include`, `#define`, `#ifdef` etc.), then resolve every bonded
//! term's parameters, either from the term itself, or from the `[ *types ]` sections, so that
//! consumers see explicit parameters per term. We read every functional form, and leave it to
//! consumers to decide which they support.
//!
//! Units are GROMACS': nm, kJ/mol, degrees, amu, and elementary charge.
//! [File format reference](https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html)

use std::{
    collections::{HashMap, HashSet},
    fs, io,
    path::{Path, PathBuf},
};

/// `[ defaults ]`
#[derive(Clone, Debug, PartialEq)]
pub struct TopDefaults {
    /// 1: Lennard-Jones. 2: Buckingham.
    pub nbfunc: u8,
    /// 1: V and W are C6 and C12, combined geometrically. 2: σ and ε; Lorentz-Berthelot.
    /// 3: σ and ε; geometric.
    pub comb_rule: u8,
    /// If true, 1-4 pairs without `[ pairtypes ]` entries are generated from the atom types, with
    /// LJ scaled by `fudge_lj`.
    pub gen_pairs: bool,
    pub fudge_lj: f64,
    pub fudge_qq: f64,
}

impl Default for TopDefaults {
    fn default() -> Self {
        Self {
            nbfunc: 1,
            comb_rule: 1,
            gen_pairs: false,
            fudge_lj: 1.,
            fudge_qq: 1.,
        }
    }
}

/// `[ atomtypes ]`
#[derive(Clone, Debug, PartialEq)]
pub struct TopAtomType {
    pub name: String,
    /// The type used to look up bonded parameters. Usually the same as `name`.
    pub bond_type: String,
    pub atomic_number: Option<u32>,
    pub mass: f64,
    pub charge: f64,
    /// A: atom. D, S, or V: virtual site/dummy.
    pub ptype: String,
    /// σ (nm) for combination rules 2 and 3; C6 for rule 1.
    pub v: f64,
    /// ε (kJ/mol) for combination rules 2 and 3; C12 for rule 1.
    pub w: f64,
}

/// An atom in a `[ moleculetype ]`.
#[derive(Clone, Debug, PartialEq)]
pub struct TopAtom {
    pub atom_type: String,
    pub resnr: i32,
    pub residue: String,
    pub name: String,
    /// Elementary charge.
    pub charge: f64,
    /// amu
    pub mass: f64,
}

/// A bonded interaction (bond, pair, angle, dihedral, or constraint), with its A-state parameters
/// resolved. Atom indices are 0-based, within the molecule. Parameter meanings depend on the
/// section and `funct`, as in the GROMACS reference manual. E.g. for harmonic bonds (funct 1):
/// b0 (nm), kb (kJ/mol/nm²), with V = ½ kb (b − b0)². For periodic dihedrals (funct 1, 4, 9):
/// φs (degrees), k (kJ/mol), multiplicity; V = k (1 + cos(n φ − φs)).
#[derive(Clone, Debug, PartialEq)]
pub struct TopInteraction {
    pub atoms: Vec<usize>,
    pub funct: u8,
    pub params: Vec<f64>,
}

/// A rigid water molecule, constrained with SETTLE. `o` is the oxygen's index; the hydrogens
/// follow it.
#[derive(Clone, Debug, PartialEq)]
pub struct TopSettle {
    pub o: usize,
    /// nm
    pub d_oh: f64,
    /// nm
    pub d_hh: f64,
}

/// `[ virtual_sites3 ]`: A massless site constructed from three atoms. For funct 1,
/// r_site = r_i + a (r_j − r_i) + b (r_k − r_i).
#[derive(Clone, Debug, PartialEq)]
pub struct TopVirtualSite3 {
    pub site: usize,
    pub from: [usize; 3],
    pub funct: u8,
    pub params: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TopMoleculeType {
    pub name: String,
    /// Exclude non-bonded interactions between atoms within this many bonds.
    pub nrexcl: usize,
    pub atoms: Vec<TopAtom>,
    pub bonds: Vec<TopInteraction>,
    pub pairs: Vec<TopInteraction>,
    pub angles: Vec<TopInteraction>,
    /// Propers and impropers. Periodic terms from `[ dihedraltypes ]` with several lines for the
    /// same atoms (funct 9) appear as separate entries.
    pub dihedrals: Vec<TopInteraction>,
    pub constraints: Vec<TopInteraction>,
    /// Explicit `[ exclusions ]`, in addition to those implied by `nrexcl`. (i < j)
    pub exclusions: Vec<(usize, usize)>,
    pub settles: Vec<TopSettle>,
    pub virtual_sites: Vec<TopVirtualSite3>,
    /// Number of `[ cmap ]` terms.
    pub n_cmap: usize,
    /// Other sections present in this molecule type, e.g. "virtual_sites2". Position restraints
    /// aren't included.
    pub other_sections: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GromacsTopology {
    pub defaults: TopDefaults,
    pub atom_types: HashMap<String, TopAtomType>,
    /// Pair-specific non-bonded parameters, which override the combining rule. (`[ nonbond_params ]`)
    pub nonbond_params: Vec<TopInteraction>,
    /// Pair-specific non-bonded parameters, by atom type name. Atom indices aren't meaningful.
    pub nonbond_param_types: Vec<(String, String)>,
    pub molecule_types: Vec<TopMoleculeType>,
    /// (Molecule type name, count), in the order the atoms appear in coordinate files.
    pub molecules: Vec<(String, usize)>,
    pub system_name: String,
    /// Sections outside molecule types that this reader doesn't interpret, e.g.
    /// "intermolecular_interactions", or "cmaptypes".
    pub other_sections: Vec<String>,
}

impl GromacsTopology {
    /// Load a topology, resolving `#include`s relative to the including file.
    pub fn load(path: &Path) -> io::Result<Self> {
        Self::load_with(path, &[], &[])
    }

    /// `include_dirs` are searched after the including file's directory, like GROMACS' `GMXLIB`
    /// and `-I`, e.g. for force field directories. `defines` are set before reading, like
    /// `define = -DPOSRES` in an MDP file.
    pub fn load_with(
        path: &Path,
        include_dirs: &[PathBuf],
        defines: &[(&str, &str)],
    ) -> io::Result<Self> {
        let mut pre = Preprocessor::new(include_dirs, defines);
        pre.process_file(path)?;
        Self::from_lines(&pre.lines)
    }

    /// Parse a topology from text. `#include` directives resolve relative to the working
    /// directory.
    pub fn new(text: &str) -> io::Result<Self> {
        let mut pre = Preprocessor::new(&[], &[]);
        pre.process_text(text, Path::new("."))?;
        Self::from_lines(&pre.lines)
    }

    pub fn molecule_type(&self, name: &str) -> Option<&TopMoleculeType> {
        self.molecule_types.iter().find(|m| m.name == name)
    }

    fn from_lines(lines: &[String]) -> io::Result<Self> {
        let mut result = Self::default();
        let mut tables = TypeTables::default();
        let mut section = String::new();
        let mut defaults_seen = false;

        for line in lines {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            if line.starts_with('[') {
                section = line
                    .trim_start_matches('[')
                    .trim_end_matches(']')
                    .trim()
                    .to_lowercase();

                if section == "moleculetype" {
                    result.molecule_types.push(TopMoleculeType::default());
                }
                continue;
            }

            let cols: Vec<&str> = line.split_whitespace().collect();
            let mol = result.molecule_types.last_mut();

            match section.as_str() {
                "defaults" => {
                    if defaults_seen {
                        return Err(err("Multiple [ defaults ] sections"));
                    }
                    defaults_seen = true;
                    result.defaults = TopDefaults {
                        nbfunc: parse(cols[0])?,
                        comb_rule: parse(cols.get(1).copied().unwrap_or("1"))?,
                        gen_pairs: cols
                            .get(2)
                            .is_some_and(|v| v.to_lowercase().starts_with('y')),
                        fudge_lj: parse(cols.get(3).copied().unwrap_or("1"))?,
                        fudge_qq: parse(cols.get(4).copied().unwrap_or("1"))?,
                    };
                }
                "atomtypes" => {
                    let t = parse_atom_type(&cols)?;
                    result.atom_types.insert(t.name.clone(), t);
                }
                "bondtypes" | "constrainttypes" => {
                    let entry = parse_type_entry(&cols, 2)?;
                    if section == "bondtypes" {
                        tables.bonds.push(entry);
                    }
                }
                "angletypes" => tables.angles.push(parse_type_entry(&cols, 3)?),
                "dihedraltypes" => tables.dihedrals.push(parse_dihedral_type(&cols)?),
                "pairtypes" => tables.pairs.push(parse_type_entry(&cols, 2)?),
                "nonbond_params" => {
                    let entry = parse_type_entry(&cols, 2)?;
                    result
                        .nonbond_param_types
                        .push((entry.types[0].clone(), entry.types[1].clone()));
                    result.nonbond_params.push(TopInteraction {
                        atoms: Vec::new(),
                        funct: entry.funct,
                        params: entry.params,
                    });
                }
                "moleculetype" => {
                    let mol = mol.unwrap();
                    mol.name = cols[0].to_owned();
                    mol.nrexcl = parse(cols.get(1).copied().unwrap_or("3"))?;
                }
                "atoms" => {
                    let mol = in_molecule(mol, &section)?;
                    let atom_type = cols[1].to_owned();
                    let at = result.atom_types.get(&atom_type);
                    mol.atoms.push(TopAtom {
                        resnr: parse(cols[2])?,
                        residue: cols[3].to_owned(),
                        name: cols[4].to_owned(),
                        charge: match cols.get(6) {
                            Some(q) => parse(q)?,
                            None => at.map(|t| t.charge).unwrap_or(0.),
                        },
                        mass: match cols.get(7) {
                            Some(m) => parse(m)?,
                            None => at
                                .map(|t| t.mass)
                                .ok_or_else(|| err(&format!("Unknown atom type {atom_type}")))?,
                        },
                        atom_type,
                    });
                }
                "bonds" => in_molecule(mol, &section)?
                    .bonds
                    .push(parse_interaction(&cols, 2)?),
                "pairs" => in_molecule(mol, &section)?
                    .pairs
                    .push(parse_interaction(&cols, 2)?),
                "angles" => in_molecule(mol, &section)?
                    .angles
                    .push(parse_interaction(&cols, 3)?),
                "dihedrals" => in_molecule(mol, &section)?
                    .dihedrals
                    .push(parse_interaction(&cols, 4)?),
                "constraints" => in_molecule(mol, &section)?
                    .constraints
                    .push(parse_interaction(&cols, 2)?),
                "exclusions" => {
                    let mol = in_molecule(mol, &section)?;
                    let i = index(cols[0])?;
                    for c in &cols[1..] {
                        let j = index(c)?;
                        mol.exclusions.push((i.min(j), i.max(j)));
                    }
                }
                "settles" => in_molecule(mol, &section)?.settles.push(TopSettle {
                    o: index(cols[0])?,
                    d_oh: parse(cols[2])?,
                    d_hh: parse(cols[3])?,
                }),
                "virtual_sites3" => {
                    let mol = in_molecule(mol, &section)?;
                    mol.virtual_sites.push(TopVirtualSite3 {
                        site: index(cols[0])?,
                        from: [index(cols[1])?, index(cols[2])?, index(cols[3])?],
                        funct: parse(cols[4])?,
                        params: cols[5..]
                            .iter()
                            .map(|v| parse(v))
                            .collect::<io::Result<_>>()?,
                    });
                }
                "cmap" => in_molecule(mol, &section)?.n_cmap += 1,
                "position_restraints"
                | "distance_restraints"
                | "dihedral_restraints"
                | "orientation_restraints"
                | "angle_restraints"
                | "angle_restraints_z" => {}
                "system" => {
                    if !result.system_name.is_empty() {
                        result.system_name.push(' ');
                    }
                    result.system_name.push_str(line);
                }
                "molecules" => {
                    result.molecules.push((
                        cols[0].to_owned(),
                        parse(cols.get(1).copied().unwrap_or("1"))?,
                    ));
                }
                "implicit_genborn_params" => {}
                other => match mol {
                    Some(mol) if !is_global_section(other) => {
                        if !mol.other_sections.iter().any(|s| s == other) {
                            mol.other_sections.push(other.to_owned());
                        }
                    }
                    _ => {
                        if !result.other_sections.iter().any(|s| s == other) {
                            result.other_sections.push(other.to_owned());
                        }
                    }
                },
            }
        }

        for mol in &mut result.molecule_types {
            resolve_params(mol, &result.atom_types, &result.defaults, &tables)?;
        }

        Ok(result)
    }
}

fn is_global_section(name: &str) -> bool {
    matches!(
        name,
        "cmaptypes" | "intermolecular_interactions" | "nonbond_params" | "implicit_genborn_params"
    )
}

fn err(msg: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("GROMACS topology: {msg}"),
    )
}

fn parse<T: std::str::FromStr>(v: &str) -> io::Result<T> {
    v.parse().map_err(|_| err(&format!("Invalid value {v:?}")))
}

/// 1-based in the file; 0-based here.
fn index(v: &str) -> io::Result<usize> {
    let i: usize = parse(v)?;
    if i == 0 {
        return Err(err("Atom indices start at 1"));
    }
    Ok(i - 1)
}

fn in_molecule<'a>(
    mol: Option<&'a mut TopMoleculeType>,
    section: &str,
) -> io::Result<&'a mut TopMoleculeType> {
    mol.ok_or_else(|| err(&format!("[ {section} ] outside a [ moleculetype ]")))
}

/// `name [bond_type] [at.num] mass charge ptype V W`. `ptype` is always third from last.
fn parse_atom_type(cols: &[&str]) -> io::Result<TopAtomType> {
    let n = cols.len();
    if n < 6 {
        return Err(err(&format!(
            "Too few columns in atom type: {}",
            cols.join(" ")
        )));
    }

    let name = cols[0].to_owned();
    let mut bond_type = name.clone();
    let mut atomic_number = None;
    for extra in &cols[1..n - 5] {
        match extra.parse::<u32>() {
            Ok(v) => atomic_number = Some(v),
            Err(_) => bond_type = (*extra).to_owned(),
        }
    }

    Ok(TopAtomType {
        name,
        bond_type,
        atomic_number,
        mass: parse(cols[n - 5])?,
        charge: parse(cols[n - 4])?,
        ptype: cols[n - 3].to_owned(),
        v: parse(cols[n - 2])?,
        w: parse(cols[n - 1])?,
    })
}

/// An interaction line in a molecule: `n_atoms` indices, then funct, then optional parameters.
fn parse_interaction(cols: &[&str], n_atoms: usize) -> io::Result<TopInteraction> {
    if cols.len() < n_atoms {
        return Err(err(&format!("Too few columns: {}", cols.join(" "))));
    }
    let atoms = cols[..n_atoms]
        .iter()
        .map(|c| index(c))
        .collect::<io::Result<_>>()?;
    // funct defaults to 1 when absent. (Rare, but allowed for some sections.)
    let funct = match cols.get(n_atoms) {
        Some(f) => parse(f)?,
        None => 1,
    };
    let params = cols
        .iter()
        .skip(n_atoms + 1)
        .map(|v| parse(v))
        .collect::<io::Result<_>>()?;

    Ok(TopInteraction {
        atoms,
        funct,
        params,
    })
}

/// A `[ *types ]` line: atom types, funct, parameters.
#[derive(Clone, Debug)]
struct TypeEntry {
    types: Vec<String>,
    funct: u8,
    params: Vec<f64>,
}

fn parse_type_entry(cols: &[&str], n_types: usize) -> io::Result<TypeEntry> {
    if cols.len() < n_types + 1 {
        return Err(err(&format!("Too few columns: {}", cols.join(" "))));
    }
    Ok(TypeEntry {
        types: cols[..n_types].iter().map(|s| (*s).to_owned()).collect(),
        funct: parse(cols[n_types])?,
        params: cols[n_types + 1..]
            .iter()
            .map(|v| parse(v))
            .collect::<io::Result<_>>()?,
    })
}

/// Dihedral types list either 4 atom types, or 2. With 2, they're the central atoms of a proper
/// dihedral, or the outer atoms of an improper (funct 2 and 4); we expand these with wildcards.
fn parse_dihedral_type(cols: &[&str]) -> io::Result<TypeEntry> {
    let is_funct = |c: Option<&&str>| c.is_some_and(|c| c.parse::<u8>().is_ok());
    if cols.len() >= 5 && is_funct(cols.get(4)) && !is_funct(cols.get(2)) {
        return parse_type_entry(cols, 4);
    }

    let mut entry = parse_type_entry(cols, 2)?;
    let (a, b) = (entry.types[0].clone(), entry.types[1].clone());
    let x = "X".to_owned();
    entry.types = if matches!(entry.funct, 2 | 4) {
        vec![a, x.clone(), x, b]
    } else {
        vec![x.clone(), a, b, x]
    };
    Ok(entry)
}

#[derive(Default)]
struct TypeTables {
    bonds: Vec<TypeEntry>,
    angles: Vec<TypeEntry>,
    dihedrals: Vec<TypeEntry>,
    pairs: Vec<TypeEntry>,
}

/// Number of A-state parameters for each section and funct. Lines may also have B-state
/// (perturbed) parameters; we drop those. None if unknown; we keep all parameters then.
fn n_params(section: &str, funct: u8) -> Option<usize> {
    Some(match (section, funct) {
        ("bonds", 1 | 2 | 6 | 7 | 8 | 9) => 2,
        ("bonds", 3) => 3,
        ("bonds", 4) => 3,
        ("bonds", 5) => 0,
        ("pairs", 1) => 2,
        ("pairs", 2) => 5,
        ("angles", 1 | 2 | 8) => 2,
        ("angles", 5) => 4,
        ("dihedrals", 1 | 4 | 9) => 3,
        ("dihedrals", 2) => 2,
        ("dihedrals", 3 | 5) => 6,
        ("constraints", 1 | 2) => 1,
        _ => return None,
    })
}

/// Fill in each interaction's parameters, from the type tables if not given inline.
fn resolve_params(
    mol: &mut TopMoleculeType,
    atom_types: &HashMap<String, TopAtomType>,
    defaults: &TopDefaults,
    tables: &TypeTables,
) -> io::Result<()> {
    let bond_types: Vec<String> = mol
        .atoms
        .iter()
        .map(|a| {
            atom_types
                .get(&a.atom_type)
                .map(|t| t.bond_type.clone())
                .unwrap_or_else(|| a.atom_type.clone())
        })
        .collect();
    let types_of =
        |atoms: &[usize]| -> Vec<&str> { atoms.iter().map(|&i| bond_types[i].as_str()).collect() };
    let mol_name = mol.name.clone();
    let label = |section: &str, atoms: &[usize]| -> String {
        let names: Vec<String> = atoms.iter().map(|i| (i + 1).to_string()).collect();
        format!(
            "Missing parameters for {section} {} ({}) in molecule {mol_name}",
            names.join("-"),
            types_of(atoms).join("-"),
        )
    };

    let truncate = |section: &str, term: &mut TopInteraction| {
        if let Some(n) = n_params(section, term.funct)
            && term.params.len() > n
        {
            term.params.truncate(n);
        }
    };

    for (section, terms, table) in [
        ("bonds", &mut mol.bonds, &tables.bonds),
        ("angles", &mut mol.angles, &tables.angles),
    ] {
        for term in terms.iter_mut() {
            truncate(section, term);
            let needed = n_params(section, term.funct).unwrap_or(0);
            if term.params.len() >= needed {
                continue;
            }
            let types = types_of(&term.atoms);
            let hit = table.iter().find(|e| {
                e.funct == term.funct
                    && (types_match(&e.types, &types, false) || {
                        let rev: Vec<&str> = types.iter().rev().copied().collect();
                        types_match(&e.types, &rev, false)
                    })
            });
            match hit {
                Some(e) => {
                    term.params = e.params.clone();
                    truncate(section, term);
                }
                None => return Err(err(&label(section, &term.atoms))),
            }
        }
    }

    // Dihedrals: Terms without parameters may expand to several, from multiple matching lines.
    let mut dihedrals = Vec::with_capacity(mol.dihedrals.len());
    for mut term in std::mem::take(&mut mol.dihedrals) {
        truncate("dihedrals", &mut term);
        let needed = n_params("dihedrals", term.funct).unwrap_or(0);
        if term.params.len() >= needed {
            dihedrals.push(term);
            continue;
        }

        let types = types_of(&term.atoms);
        let entries = match_dihedral_types(&tables.dihedrals, &types, term.funct);
        if entries.is_empty() {
            return Err(err(&label("dihedral", &term.atoms)));
        }
        for e in entries {
            let mut t = TopInteraction {
                atoms: term.atoms.clone(),
                funct: term.funct,
                params: e.params.clone(),
            };
            truncate("dihedrals", &mut t);
            dihedrals.push(t);
        }
    }
    mol.dihedrals = dihedrals;

    // Pairs: From `[ pairtypes ]`, or generated from atom types.
    for term in &mut mol.pairs {
        truncate("pairs", term);
        if term.funct != 1 || term.params.len() >= 2 {
            continue;
        }
        // Pair types use atom type names, not bond types.
        let types: Vec<&str> = term
            .atoms
            .iter()
            .map(|&i| mol.atoms[i].atom_type.as_str())
            .collect();
        let hit = tables.pairs.iter().find(|e| {
            e.funct == 1
                && ((e.types[0] == types[0] && e.types[1] == types[1])
                    || (e.types[0] == types[1] && e.types[1] == types[0]))
        });
        if let Some(e) = hit {
            term.params = e.params.clone();
            truncate("pairs", term);
        } else if defaults.gen_pairs {
            let (a, b) = match (atom_types.get(types[0]), atom_types.get(types[1])) {
                (Some(a), Some(b)) => (a, b),
                _ => return Err(err(&label("pair", &term.atoms))),
            };
            let (v, w) = combine(defaults.comb_rule, a.v, a.w, b.v, b.w);
            // For combination rule 1, both C6 and C12 scale. For 2 and 3, only ε does.
            term.params = if defaults.comb_rule == 1 {
                vec![v * defaults.fudge_lj, w * defaults.fudge_lj]
            } else {
                vec![v, w * defaults.fudge_lj]
            };
        } else {
            return Err(err(&label("pair", &term.atoms)));
        }
    }

    for term in &mut mol.constraints {
        truncate("constraints", term);
    }

    Ok(())
}

/// Combine per-type non-bonded parameters. See `TopDefaults::comb_rule`.
pub fn combine(comb_rule: u8, v0: f64, w0: f64, v1: f64, w1: f64) -> (f64, f64) {
    match comb_rule {
        2 => (0.5 * (v0 + v1), (w0 * w1).sqrt()),
        _ => ((v0 * v1).sqrt(), (w0 * w1).sqrt()),
    }
}

/// Does a type-table entry match these atom types? "X" is a wildcard. With `exact_only`, we
/// don't allow wildcards.
fn types_match(entry: &[String], types: &[&str], exact_only: bool) -> bool {
    entry.len() == types.len()
        && entry
            .iter()
            .zip(types)
            .all(|(e, t)| e == t || (!exact_only && e == "X"))
}

/// GROMACS' rule: The matching lines with the fewest wildcards win; of those, the first in the
/// file, and any lines directly following it with the same types (multiple terms, funct 9).
fn match_dihedral_types<'a>(
    table: &'a [TypeEntry],
    types: &[&str],
    funct: u8,
) -> Vec<&'a TypeEntry> {
    let compatible = |f: u8| f == funct || (matches!(funct, 1 | 9) && matches!(f, 1 | 9));
    let rev: Vec<&str> = types.iter().rev().copied().collect();

    let mut best: Option<(usize, usize)> = None; // (wildcards, index)
    for (i, e) in table.iter().enumerate() {
        if !compatible(e.funct) {
            continue;
        }
        if types_match(&e.types, types, false) || types_match(&e.types, &rev, false) {
            let wildcards = e.types.iter().filter(|t| *t == "X").count();
            if best.is_none_or(|(w, _)| wildcards < w) {
                best = Some((wildcards, i));
            }
        }
    }

    let Some((_, first)) = best else {
        return Vec::new();
    };
    let key = &table[first].types;
    table[first..]
        .iter()
        .take_while(|e| &e.types == key && compatible(e.funct))
        .collect()
}

/// Runs the C-style preprocessor GROMACS uses, producing logical lines with comments removed.
struct Preprocessor {
    include_dirs: Vec<PathBuf>,
    defines: HashMap<String, String>,
    lines: Vec<String>,
    /// Guards against include cycles.
    active_files: HashSet<PathBuf>,
}

impl Preprocessor {
    fn new(include_dirs: &[PathBuf], defines: &[(&str, &str)]) -> Self {
        Self {
            include_dirs: include_dirs.to_vec(),
            defines: defines
                .iter()
                .map(|(k, v)| ((*k).to_owned(), (*v).to_owned()))
                .collect(),
            lines: Vec::new(),
            active_files: HashSet::new(),
        }
    }

    fn process_file(&mut self, path: &Path) -> io::Result<()> {
        let canonical = fs::canonicalize(path)
            .map_err(|e| err(&format!("Unable to open {}: {e}", path.display())))?;
        if !self.active_files.insert(canonical.clone()) {
            return Err(err(&format!("Include cycle at {}", path.display())));
        }

        let text = fs::read_to_string(path)?;
        let dir = path.parent().unwrap_or(Path::new(".")).to_path_buf();
        self.process_text(&text, &dir)?;

        self.active_files.remove(&canonical);
        Ok(())
    }

    fn process_text(&mut self, text: &str, dir: &Path) -> io::Result<()> {
        // Each level: (this branch is active, some branch of this #if has been taken).
        let mut conditions: Vec<(bool, bool)> = Vec::new();
        let active = |c: &[(bool, bool)]| c.iter().all(|(a, _)| *a);

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

            let content = match line.find(';') {
                Some(i) => &line[..i],
                None => &line,
            };
            let trimmed = content.trim();

            if let Some(directive) = trimmed.strip_prefix('#') {
                let mut parts = directive.trim().splitn(2, char::is_whitespace);
                let name = parts.next().unwrap_or_default();
                let arg = parts.next().unwrap_or_default().trim();

                match name {
                    "ifdef" | "ifndef" => {
                        let defined = self.defines.contains_key(arg);
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
                    _ if !active(&conditions) => {}
                    "define" => {
                        let mut kv = arg.splitn(2, char::is_whitespace);
                        let key = kv.next().unwrap_or_default().to_owned();
                        let value = kv.next().unwrap_or_default().trim().to_owned();
                        self.defines.insert(key, value);
                    }
                    "undef" => {
                        self.defines.remove(arg);
                    }
                    "include" => {
                        let file = arg.trim_matches(|c| c == '"' || c == '<' || c == '>');
                        let path = std::iter::once(dir.to_path_buf())
                            .chain(self.include_dirs.iter().cloned())
                            .map(|d| d.join(file))
                            .find(|p| p.is_file())
                            .ok_or_else(|| {
                                err(&format!(
                                    "Unable to find included file {file}. Pass its directory \
                                     (e.g. the force field's) as an include directory."
                                ))
                            })?;
                        self.process_file(&path)?;
                    }
                    "error" => return Err(err(&format!("#error {arg}"))),
                    _ => {}
                }
                continue;
            }

            if !active(&conditions) || trimmed.is_empty() {
                continue;
            }

            // Substitute macros, token by token. (E.g. GROMOS and CHARMM parameter macros.)
            let expanded: Vec<&str> = trimmed
                .split_whitespace()
                .map(|tok| self.defines.get(tok).map(String::as_str).unwrap_or(tok))
                .collect();
            self.lines.push(expanded.join(" "));
        }

        if !conditions.is_empty() {
            return Err(err("Unterminated #ifdef"));
        }
        Ok(())
    }
}
