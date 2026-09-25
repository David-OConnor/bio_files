//! CHARMM force field files: Topologies (`.rtf`), parameters (`.prm`), and stream files (`.str`),
//! which may contain both. [Format reference](https://academiccharmm.org/documentation/latest/rtop)
//!
//! Units are CHARMM's: Å, kcal/mol, degrees, amu, and elementary charge. Energy conventions:
//! Bonds: Kb (b − b0)². Angles: Kθ (θ − θ0)². Urey-Bradley: Kub (S − S0)². Dihedrals:
//! Kχ (1 + cos(n χ − δ)). Impropers: Kψ (ψ − ψ0)². LJ: ε [(Rmin/r)¹² − 2 (Rmin/r)⁶], with
//! ε_ij = √(ε_i ε_j) and Rmin_ij = Rmin/2_i + Rmin/2_j.

use std::{collections::HashMap, io};

fn err(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("CHARMM: {msg}"))
}

fn parse<T: std::str::FromStr>(v: &str) -> io::Result<T> {
    v.parse().map_err(|_| err(&format!("Invalid value {v:?}")))
}

/// The part of a line before any `!` comment, split into tokens.
fn tokens(line: &str) -> Vec<&str> {
    let content = match line.find('!') {
        Some(i) => &line[..i],
        None => line,
    };
    content.split_whitespace().collect()
}

/// CHARMM keywords may be abbreviated to 4 characters, and are case-insensitive.
fn keyword(tok: &str) -> String {
    tok.chars().take(4).collect::<String>().to_uppercase()
}

/// An atom in a residue or patch.
#[derive(Clone, Debug, PartialEq)]
pub struct RtfAtom {
    pub name: String,
    pub atom_type: String,
    /// Elementary charge.
    pub charge: f32,
}

/// A residue (`RESI`) or patch (`PRES`). Atom names in bonded terms may have prefixes: `-` and `+`
/// for the previous and next residue, or in patches that span residues, `1`, `2` etc.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RtfResidue {
    pub name: String,
    pub charge: f32,
    pub is_patch: bool,
    pub atoms: Vec<RtfAtom>,
    /// Includes `DOUBLE` and `TRIPLE` bonds.
    pub bonds: Vec<[String; 2]>,
    /// Explicit angles and dihedrals. Usually these are generated from bonds instead.
    pub angles: Vec<[String; 3]>,
    pub dihedrals: Vec<[String; 4]>,
    /// The first atom is the central one.
    pub impropers: Vec<[String; 4]>,
    /// Two dihedrals: φ then ψ.
    pub cmaps: Vec<[String; 8]>,
    /// Atoms a patch removes.
    pub delete_atoms: Vec<String>,
    /// Patches to apply when this residue is first or last in a chain, overriding the defaults.
    /// "NONE" for no patch.
    pub first_patch: Option<String>,
    pub last_patch: Option<String>,
    /// `NOANG` / `NODIH`: Don't generate angles or dihedrals from bonds. (E.g. rigid water)
    pub no_auto_angles: bool,
    pub no_auto_dihedrals: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CharmmTopology {
    /// Atom type → (mass in amu, element symbol)
    pub masses: HashMap<String, (f32, String)>,
    /// Residues and patches.
    pub residues: Vec<RtfResidue>,
    /// Default patches for the first and last residues of chains.
    pub default_first: Option<String>,
    pub default_last: Option<String>,
}

impl CharmmTopology {
    pub fn residue(&self, name: &str) -> Option<&RtfResidue> {
        self.residues
            .iter()
            .find(|r| r.name.eq_ignore_ascii_case(name))
    }

    /// Parse an RTF file's content. Later residues with the same name replace earlier ones.
    pub fn new(text: &str) -> io::Result<Self> {
        let mut result = Self::default();
        result.add_rtf(text)?;
        Ok(result)
    }

    /// Add the content of an RTF file, or an RTF block from a stream file.
    pub fn add_rtf(&mut self, text: &str) -> io::Result<()> {
        let mut current: Option<RtfResidue> = None;

        for line in text.lines() {
            if line.trim_start().starts_with('*') {
                continue; // Title
            }
            let toks = tokens(line);
            if toks.is_empty() {
                continue;
            }
            let kw = keyword(toks[0]);
            let args = &toks[1..];

            // Pairs, triples etc. of atom names following a keyword.
            let groups = |n: usize| -> Vec<Vec<String>> {
                args.chunks_exact(n)
                    .map(|c| c.iter().map(|s| s.to_uppercase()).collect())
                    .collect()
            };

            match kw.as_str() {
                "MASS" => {
                    // MASS <number> <type> <mass> [element]
                    if args.len() >= 3 {
                        let element = args.get(3).map(|e| e.to_string()).unwrap_or_default();
                        self.masses
                            .insert(args[1].to_uppercase(), (parse(args[2])?, element));
                    }
                }
                "DEFA" => {
                    // DEFA FIRS <patch> LAST <patch>
                    for pair in args.chunks_exact(2) {
                        match keyword(pair[0]).as_str() {
                            "FIRS" => self.default_first = Some(pair[1].to_uppercase()),
                            "LAST" => self.default_last = Some(pair[1].to_uppercase()),
                            _ => {}
                        }
                    }
                }
                "RESI" | "PRES" => {
                    if let Some(r) = current.take() {
                        self.push(r);
                    }
                    let upper: Vec<String> = args.iter().map(|a| a.to_uppercase()).collect();
                    current = Some(RtfResidue {
                        name: upper.first().cloned().unwrap_or_default(),
                        charge: args.get(1).and_then(|c| c.parse().ok()).unwrap_or(0.),
                        is_patch: kw == "PRES",
                        no_auto_angles: upper.iter().any(|a| a.starts_with("NOAN")),
                        no_auto_dihedrals: upper.iter().any(|a| a.starts_with("NODI")),
                        ..Default::default()
                    });
                }
                "END" => {
                    if let Some(r) = current.take() {
                        self.push(r);
                    }
                }
                _ => {
                    let Some(r) = current.as_mut() else {
                        continue;
                    };
                    match kw.as_str() {
                        "ATOM" => {
                            if args.len() >= 3 {
                                let atom = RtfAtom {
                                    name: args[0].to_uppercase(),
                                    atom_type: args[1].to_uppercase(),
                                    charge: parse(args[2])?,
                                };
                                // A patch may redefine an atom.
                                match r.atoms.iter_mut().find(|a| a.name == atom.name) {
                                    Some(a) => *a = atom,
                                    None => r.atoms.push(atom),
                                }
                            }
                        }
                        "BOND" | "DOUB" | "TRIP" | "AROM" => {
                            for g in groups(2) {
                                r.bonds.push([g[0].clone(), g[1].clone()]);
                            }
                        }
                        "ANGL" | "THET" => {
                            for g in groups(3) {
                                r.angles.push([g[0].clone(), g[1].clone(), g[2].clone()]);
                            }
                        }
                        "DIHE" | "PHI" => {
                            for g in groups(4) {
                                r.dihedrals.push(std::array::from_fn(|i| g[i].clone()));
                            }
                        }
                        "IMPR" | "IMPH" => {
                            for g in groups(4) {
                                r.impropers.push(std::array::from_fn(|i| g[i].clone()));
                            }
                        }
                        "CMAP" => {
                            for g in groups(8) {
                                r.cmaps.push(std::array::from_fn(|i| g[i].clone()));
                            }
                        }
                        "DELE" => {
                            // DELETE ATOM <name>. (Other deletions, e.g. of acceptors, don't
                            // affect the energy.)
                            if args.len() >= 2 && keyword(args[0]) == "ATOM" {
                                r.delete_atoms.push(args[1].to_uppercase());
                            }
                        }
                        "PATC" => {
                            for pair in args.chunks_exact(2) {
                                match keyword(pair[0]).as_str() {
                                    "FIRS" => r.first_patch = Some(pair[1].to_uppercase()),
                                    "LAST" => r.last_patch = Some(pair[1].to_uppercase()),
                                    _ => {}
                                }
                            }
                        }
                        // GROUP, DONOR, ACCEPTOR, IC, DECL, AUTO etc.: Not relevant to energies.
                        _ => {}
                    }
                }
            }
        }

        if let Some(r) = current.take() {
            self.push(r);
        }
        Ok(())
    }

    fn push(&mut self, res: RtfResidue) {
        self.residues.retain(|r| r.name != res.name);
        self.residues.push(res);
    }
}

/// Angle parameters, with an optional Urey-Bradley term between the outer atoms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CharmmAngle {
    /// kcal/mol/rad²
    pub k: f32,
    /// Degrees
    pub theta0: f32,
    /// (Kub kcal/mol/Å², S0 Å)
    pub urey_bradley: Option<(f32, f32)>,
}

/// One term of a dihedral. Several can share the same types.
#[derive(Clone, Debug, PartialEq)]
pub struct CharmmDihedral {
    /// May include "X" wildcards.
    pub types: [String; 4],
    /// kcal/mol
    pub k: f32,
    pub n: u8,
    /// Degrees
    pub delta: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CharmmImproper {
    /// May include "X" wildcards. The first atom is the central one.
    pub types: [String; 4],
    /// kcal/mol/rad²
    pub k: f32,
    /// Degrees
    pub psi0: f32,
}

/// A CMAP correction map for a pair of dihedrals, φ and ψ.
#[derive(Clone, Debug, PartialEq)]
pub struct CharmmCmap {
    /// Types of φ's four atoms, then ψ's.
    pub types: [String; 8],
    /// Grid points per dimension. Spacing is 360° / size, starting at −180°.
    pub size: usize,
    /// kcal/mol. Row-major: `energies[i * size + j]` is at φ = −180° + i Δ, ψ = −180° + j Δ.
    pub energies: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CharmmNonbonded {
    /// kcal/mol. Positive; the file stores it negated.
    pub eps: f32,
    /// Å
    pub rmin_half: f32,
    /// Special parameters for 1-4 pairs, if different.
    pub eps_14: Option<f32>,
    pub rmin_half_14: Option<f32>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CharmmParams {
    /// From `MASS` lines in the `ATOMS` section.
    pub masses: HashMap<String, f32>,
    /// (Kb kcal/mol/Å², b0 Å), keyed by type pair in file order.
    pub bonds: HashMap<(String, String), (f32, f32)>,
    pub angles: HashMap<(String, String, String), CharmmAngle>,
    pub dihedrals: Vec<CharmmDihedral>,
    pub impropers: Vec<CharmmImproper>,
    pub cmaps: Vec<CharmmCmap>,
    pub nonbonded: HashMap<String, CharmmNonbonded>,
    /// Pair-specific LJ: (ε kcal/mol, positive; Rmin Å), keyed by type pair in file order.
    pub nbfix: HashMap<(String, String), (f32, f32)>,
}

impl CharmmParams {
    pub fn new(text: &str) -> io::Result<Self> {
        let mut result = Self::default();
        result.add_prm(text)?;
        Ok(result)
    }

    /// Add the content of a parameter file, or a parameter block from a stream file. Entries
    /// replace earlier ones with the same types; for dihedrals, all terms with those types.
    pub fn add_prm(&mut self, text: &str) -> io::Result<()> {
        let mut section = String::new();
        let mut skip_continuation = false;
        // Dihedral types whose terms this block has started to replace.
        let mut replaced_dihedrals: Vec<[String; 4]> = Vec::new();
        // An in-progress CMAP: (types, size, energies)
        let mut cmap: Option<CharmmCmap> = None;

        for line in text.lines() {
            if line.trim_start().starts_with('*') {
                continue;
            }
            let toks = tokens(line);
            if toks.is_empty() {
                continue;
            }
            if skip_continuation {
                skip_continuation = false;
                continue;
            }

            let kw = keyword(toks[0]);
            let is_header = matches!(
                kw.as_str(),
                "ATOM"
                    | "BOND"
                    | "ANGL"
                    | "THET"
                    | "DIHE"
                    | "PHI"
                    | "IMPR"
                    | "IMPH"
                    | "CMAP"
                    | "NONB"
                    | "NBON"
                    | "NBFI"
                    | "HBON"
                    | "END"
            ) && (toks.len() == 1 || kw == "NONB" || kw == "NBON" || kw == "HBON");

            if is_header {
                if let Some(c) = cmap.take() {
                    self.finish_cmap(c)?;
                }
                section = kw;
                // E.g. "NONBONDED nbxmod 5 ... -", continued on the next line.
                if toks.last() == Some(&"-") {
                    skip_continuation = true;
                }
                continue;
            }

            let upper = |i: usize| toks[i].to_uppercase();

            match section.as_str() {
                "ATOM" => {
                    if keyword(toks[0]) == "MASS" && toks.len() >= 4 {
                        self.masses.insert(upper(2), parse(toks[3])?);
                    }
                }
                "BOND" if toks.len() >= 4 => {
                    self.bonds
                        .insert((upper(0), upper(1)), (parse(toks[2])?, parse(toks[3])?));
                }
                "ANGL" | "THET" if toks.len() >= 5 => {
                    let urey_bradley = if toks.len() >= 7 {
                        Some((parse(toks[5])?, parse(toks[6])?))
                    } else {
                        None
                    };
                    self.angles.insert(
                        (upper(0), upper(1), upper(2)),
                        CharmmAngle {
                            k: parse(toks[3])?,
                            theta0: parse(toks[4])?,
                            urey_bradley,
                        },
                    );
                }
                "DIHE" | "PHI" if toks.len() >= 7 => {
                    let types: [String; 4] = std::array::from_fn(upper);
                    if !replaced_dihedrals.contains(&types) {
                        self.dihedrals.retain(|d| d.types != types);
                        replaced_dihedrals.push(types.clone());
                    }
                    self.dihedrals.push(CharmmDihedral {
                        types,
                        k: parse(toks[4])?,
                        n: parse(toks[5])?,
                        delta: parse(toks[6])?,
                    });
                }
                "IMPR" | "IMPH" if toks.len() >= 7 => {
                    let types: [String; 4] = std::array::from_fn(upper);
                    self.impropers.retain(|d| d.types != types);
                    self.impropers.push(CharmmImproper {
                        types,
                        k: parse(toks[4])?,
                        psi0: parse(toks[6])?,
                    });
                }
                "CMAP" => {
                    // A header: 8 types and the grid size. Otherwise, grid values.
                    if toks.len() == 9
                        && toks[8].parse::<usize>().is_ok()
                        && toks[0].parse::<f64>().is_err()
                    {
                        if let Some(c) = cmap.take() {
                            self.finish_cmap(c)?;
                        }
                        cmap = Some(CharmmCmap {
                            types: std::array::from_fn(upper),
                            size: parse(toks[8])?,
                            energies: Vec::new(),
                        });
                    } else if let Some(c) = cmap.as_mut() {
                        for t in &toks {
                            c.energies.push(parse(t)?);
                        }
                    }
                }
                "NONB" | "NBON" if toks.len() >= 4 => {
                    let (eps_14, rmin_half_14) = if toks.len() >= 7 {
                        (Some(parse::<f32>(toks[5])?.abs()), Some(parse(toks[6])?))
                    } else {
                        (None, None)
                    };
                    self.nonbonded.insert(
                        upper(0),
                        CharmmNonbonded {
                            eps: parse::<f32>(toks[2])?.abs(),
                            rmin_half: parse(toks[3])?,
                            eps_14,
                            rmin_half_14,
                        },
                    );
                }
                "NBFI" if toks.len() >= 4 => {
                    self.nbfix.insert(
                        (upper(0), upper(1)),
                        (parse::<f32>(toks[2])?.abs(), parse(toks[3])?),
                    );
                }
                _ => {}
            }
        }

        if let Some(c) = cmap.take() {
            self.finish_cmap(c)?;
        }
        Ok(())
    }

    fn finish_cmap(&mut self, c: CharmmCmap) -> io::Result<()> {
        if c.energies.len() != c.size * c.size {
            return Err(err(&format!(
                "CMAP {} has {} values; expected {}",
                c.types.join(" "),
                c.energies.len(),
                c.size * c.size
            )));
        }
        self.cmaps.retain(|m| m.types != c.types);
        self.cmaps.push(c);
        Ok(())
    }

    /// Bond parameters for a type pair, in either order.
    pub fn bond(&self, t0: &str, t1: &str) -> Option<(f32, f32)> {
        self.bonds
            .get(&(t0.to_owned(), t1.to_owned()))
            .or_else(|| self.bonds.get(&(t1.to_owned(), t0.to_owned())))
            .copied()
    }

    /// Angle parameters for types, in either direction.
    pub fn angle(&self, t0: &str, t1: &str, t2: &str) -> Option<CharmmAngle> {
        self.angles
            .get(&(t0.to_owned(), t1.to_owned(), t2.to_owned()))
            .or_else(|| {
                self.angles
                    .get(&(t2.to_owned(), t1.to_owned(), t0.to_owned()))
            })
            .copied()
    }

    /// Dihedral terms for types, as CHARMM matches them: All terms of an exact match, in either
    /// direction; otherwise those of a match with wildcard outer atoms (X-B-C-X).
    pub fn dihedral(&self, t: [&str; 4]) -> Vec<&CharmmDihedral> {
        let matches = |d: &CharmmDihedral, pattern: [&str; 4]| {
            let rev = [pattern[3], pattern[2], pattern[1], pattern[0]];
            [pattern, rev]
                .iter()
                .any(|p| d.types.iter().zip(p).all(|(dt, pt)| dt == pt))
        };

        let exact: Vec<_> = self.dihedrals.iter().filter(|d| matches(d, t)).collect();
        if !exact.is_empty() {
            return exact;
        }
        self.dihedrals
            .iter()
            .filter(|d| matches(d, ["X", t[1], t[2], "X"]))
            .collect()
    }

    /// Improper parameters for types (central atom first), as CHARMM matches them: Exact, then
    /// A-X-X-D, X-B-C-D, X-B-C-X, and X-X-C-D; each in either direction.
    pub fn improper(&self, t: [&str; 4]) -> Option<&CharmmImproper> {
        let patterns = [
            t,
            [t[0], "X", "X", t[3]],
            ["X", t[1], t[2], t[3]],
            ["X", t[1], t[2], "X"],
            ["X", "X", t[2], t[3]],
        ];
        for p in patterns {
            let rev = [p[3], p[2], p[1], p[0]];
            if let Some(hit) = self.impropers.iter().find(|imp| {
                [p, rev]
                    .iter()
                    .any(|q| imp.types.iter().zip(q).all(|(it, qt)| it == qt))
            }) {
                return Some(hit);
            }
        }
        None
    }

    /// The CMAP for the types of φ's and ψ's atoms.
    pub fn cmap(&self, t: [&str; 8]) -> Option<&CharmmCmap> {
        self.cmaps
            .iter()
            .find(|m| m.types.iter().zip(&t).all(|(mt, tt)| mt == tt))
    }

    /// Pair-specific LJ (ε, Rmin), in either order.
    pub fn nbfix(&self, t0: &str, t1: &str) -> Option<(f32, f32)> {
        self.nbfix
            .get(&(t0.to_owned(), t1.to_owned()))
            .or_else(|| self.nbfix.get(&(t1.to_owned(), t0.to_owned())))
            .copied()
    }
}

/// Split a stream file (`.str`) into its topology and parameter blocks: The lines between
/// `read rtf card` / `read para card` and the following `END`. Other content is CHARMM script.
pub fn split_stream(text: &str) -> (Vec<String>, Vec<String>) {
    let mut rtf = Vec::new();
    let mut prm = Vec::new();
    let mut current: Option<(bool, String)> = None; // (is_rtf, content)

    for line in text.lines() {
        let toks = tokens(line);
        let lower: Vec<String> = toks.iter().map(|t| t.to_lowercase()).collect();

        if current.is_none() {
            if lower.len() >= 2 && lower[0] == "read" {
                if lower[1].starts_with("rtf") {
                    current = Some((true, String::new()));
                } else if lower[1].starts_with("para") {
                    current = Some((false, String::new()));
                }
            }
            continue;
        }

        let (is_rtf, content) = current.as_mut().unwrap();
        content.push_str(line);
        content.push('\n');
        if lower.first().is_some_and(|t| t == "end") {
            let block = std::mem::take(content);
            if *is_rtf {
                rtf.push(block);
            } else {
                prm.push(block);
            }
            current = None;
        }
    }

    (rtf, prm)
}
