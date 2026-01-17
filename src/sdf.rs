//! For opening Structure Data Format (SDF) files. These are common molecular descriptions for ligands. It's a simpler format
//! than PDB.

use std::{
    collections::HashMap,
    fs,
    fs::File,
    io,
    io::{ErrorKind, Write},
    path::Path,
    str::FromStr,
};

use bio_apis::{drugbank, pdbe, pubchem, pubchem::StructureSearchNamespace};
use lin_alg::f64::Vec3;
use na_seq::Element;

use crate::{
    AtomGeneric, BondGeneric, BondType, ChainGeneric, Mol2, ResidueEnd, ResidueGeneric, ResidueType,
};

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum PharmacophoreType {
    Acceptor,
    Donor,
    Cation,
    Rings,
    Hydrophobe,
    Anion,
}

impl PharmacophoreType {
    fn from_pubchem_str(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "acceptor" => Some(Self::Acceptor),
            "donor" => Some(Self::Donor),
            "cation" => Some(Self::Cation),
            "rings" | "ring" => Some(Self::Rings),
            "hydrophobe" => Some(Self::Hydrophobe),
            "anion" => Some(Self::Anion),
            _ => None,
        }
    }

    fn to_pubchem_str(self) -> &'static str {
        match self {
            Self::Acceptor => "acceptor",
            Self::Donor => "donor",
            Self::Cation => "cation",
            Self::Rings => "rings",
            Self::Hydrophobe => "hydrophobe",
            Self::Anion => "anion",
        }
    }
}

#[derive(Clone, Debug)]
pub struct PharmacaphoreFeatures {
    pub atom_sns: Vec<u32>,       // 1-based atom indices (SDF serial numbers)
    pub type_: PharmacophoreType, // e.g. "acceptor", "cation", "rings"
}

/// It's a format used for small organic molecules, and is a common format on online databases
/// like PubChem and Drugbank. This struct will likely
/// be used as an intermediate format, and converted to something application-specific.
#[derive(Clone, Debug)]
pub struct Sdf {
    pub ident: String,
    pub metadata: HashMap<String, String>,
    pub atoms: Vec<AtomGeneric>,
    pub bonds: Vec<BondGeneric>,
    pub chains: Vec<ChainGeneric>,
    pub residues: Vec<ResidueGeneric>,
    pub pharmacophore_features: Vec<PharmacaphoreFeatures>,
}

/// Workaround for some SDFs we've seen on PubChem. Applies to both counts, and bond
/// rows. This variant is safer for the atom count where atom and bond counts are likely
/// to have the same number of digits; things are trickier with bonds.
///
/// This still causes trouble when, for example, atom count is 3-digit and bond count is 2-digit,
/// for example.
fn split_atom_bond_count_col(col: &str, split: usize) -> io::Result<(usize, usize)> {
    if col.len() < split {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            format!("Atom/bond count column too short to split: {col}"),
        ));
    }

    // // todo temp
    // println!(
    //     "Split: {split}, part 0: {} part 1: {}",
    //     &col[..split],
    //     &col[split..]
    // );

    let n_atoms = col[..split]
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse number of atoms"))?;

    let n_bonds = col[split..]
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse number of bonds"))?;

    Ok((n_atoms, n_bonds))
}

/// Similar to our atom_bond_count split function, but takes into account the nearby
/// lines to make an educated guess; otherwise
///
/// Both of these split functions are from observed problems on PubChem where spaces are omitted
/// in bond and atom/bond count lines.
fn split_bond_indices_line(
    col: &str,
    prev: Option<&str>,
    _next: Option<&str>, // We don't strictly need 'next' if we trust 'prev' + validity
    max_atom_index: usize,
) -> io::Result<(usize, usize)> {
    // 1. Identify all Chemically Valid Splits
    // A split is valid ONLY if both resulting numbers are <= max_atom_index.
    let mut valid_splits = Vec::new();

    for i in 1..col.len() {
        let s1 = &col[..i];
        let s2 = &col[i..];

        if let (Ok(n1), Ok(n2)) = (s1.parse::<usize>(), s2.parse::<usize>()) {
            // Hard Constraint: Atoms must exist
            if n1 <= max_atom_index && n2 <= max_atom_index && n1 > 0 && n2 > 0 {
                valid_splits.push((n1, n2));
            }
        }
    }

    // 2. Decision Logic
    match valid_splits.len() {
        // Case A: No valid splits found. The data is corrupt or the column is garbage.
        0 => Err(io::Error::new(
            ErrorKind::InvalidData,
            format!("No valid bond split found for '{col}' with max atom {max_atom_index}"),
        )),

        // Case B: Exactly one valid split.
        // We don't care about sorting/context here. This is the only chemical possibility.
        // This fixes your '55114' case: (55, 114) is valid, (551, 14) is not.
        1 => Ok(valid_splits[0]),

        // Case C: Ambiguity! Multiple splits are chemically valid.
        // Example: "1234" with 200 atoms could be (1, 234-No), (12, 34-Yes), (123, 4-Yes).
        // NOW we use the context (prev) to guess which one is correct.
        _ => {
            let p_val = prev.and_then(|s| s.parse::<usize>().ok()).unwrap_or(0);

            // Find the split where the first atom (n1) is closest to the previous line's atom.
            // In a sorted list, n1 should be >= prev and close to it.
            valid_splits.sort_by_key(|(n1, _)| {
                let diff = if *n1 >= p_val {
                    n1 - p_val
                } else {
                    // If n1 < prev, it's out of order. Penalize it heavily so we prefer sorted options.
                    // But don't make it impossible, just unlikely.
                    (p_val - n1) + 10000
                };
                diff
            });

            // Return the "closest" valid split
            Ok(valid_splits[0])
        }
    }
}

impl Sdf {
    /// From a string of an SDF text file.
    pub fn new(text: &str) -> io::Result<Self> {
        let lines: Vec<&str> = text.lines().collect();

        // SDF files typically have at least 4 lines before the atom block:
        //   1) A title or identifier
        //   2) Usually blank or comments
        //   3) Often blank or comments
        //   4) "counts" line: e.g. " 50  50  0  ..." for V2000
        if lines.len() < 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines to parse an SDF header",
            ));
        }

        // todo: Incorporate more cols A/R.
        // After element:
        // Mass difference (0, unless an isotope)
        // Charge (+1 for cation etc)
        // Stereo, valence, other flags

        // todo: Do bonds too
        // first atom index
        // second atom index
        // 1 for single, 2 for double etc
        // 0 for no stereochemistry, 1=up, 6=down etc
        // Other properties: Bond topology, reaction center flags etc. Usually 0

        // This is the "counts" line, e.g. " 50 50  0  0  0  0  0  0  0999 V2000"
        let counts_line = lines[3];
        let counts_cols: Vec<&str> = counts_line.split_whitespace().collect();

        if counts_cols.len() < 2 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Counts line doesn't have enough fields",
            ));
        }

        // Typically, the first number is the number of atoms (natoms)
        // and the second number is the number of bonds (nbonds).
        let mut n_atoms = counts_cols[0].parse::<usize>().map_err(|_| {
            io::Error::new(ErrorKind::InvalidData, "Could not parse number of atoms")
        })?;

        let mut n_bonds = counts_cols[1].parse::<usize>().map_err(|_| {
            io::Error::new(ErrorKind::InvalidData, "Could not parse number of bonds")
        })?;

        // Now read the next 'natoms' lines as the atom block.
        // Each line usually looks like:
        //   X Y Z Element ??? ??? ...
        //   e.g. "    1.4386   -0.8054   -0.4963 O   0  0  0  0  0  0  0  0  0  0  0  0"
        //

        let first_atom_line = 4;
        let mut last_atom_line = first_atom_line + n_atoms;
        let mut first_bond_line = last_atom_line;
        let mut last_bond_line = first_bond_line + n_bonds;

        if lines.len() < last_atom_line {
            // We have observed that sometimes the space between atom and bond counts is omitted,
            // from PubChem molecules that have triple-digit counts, so we will assume that here.
            let atom_bond_count = counts_cols[0];

            // If atom 0 is 3 digits, there is no space in these cases. So could be 4-6 total len.
            if atom_bond_count.len() >= 4 {
                (n_atoms, n_bonds) = split_atom_bond_count_col(atom_bond_count, 3)?;
                // Note: We have observed cases where we need to split at atom 2, i.e. if
                // the count is 5 digits long. So it will still fail if the space is ommitted.

                last_atom_line = first_atom_line + n_atoms;
                first_bond_line = last_atom_line;
                last_bond_line = first_bond_line + n_bonds;

                if lines.len() < last_atom_line {
                    return Err(io::Error::new(
                        ErrorKind::InvalidData,
                        format!(
                            "Not enough lines for the declared atom block (after 3/3 split.) \
                        Lines: {}, expected: {} atoms",
                            lines.len(),
                            last_atom_line - 4
                        ),
                    ));
                }
            } else {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Not enough lines for the declared atom block",
                ));
            }
        }

        let mut atoms = Vec::with_capacity(n_atoms);

        let mut pharmacophore_features: Vec<PharmacaphoreFeatures> = Vec::new();

        for i in first_atom_line..last_atom_line {
            let line = lines[i];
            let cols: Vec<&str> = line.split_whitespace().collect();

            if cols.len() < 4 {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Atom line {i} does not have enough columns"),
                ));
            }

            let x = cols[0].parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse X coordinate")
            })?;
            let y = cols[1].parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse Y coordinate")
            })?;
            let z = cols[2].parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse Z coordinate")
            })?;
            let element = cols[3];

            atoms.push(AtomGeneric {
                // SDF doesn't explicitly include incices.
                serial_number: (i - first_atom_line) as u32 + 1,
                posit: Vec3 { x, y, z }, // or however you store coordinates
                element: Element::from_letter(element)?,
                hetero: true,
                ..Default::default()
            });
        }

        let mut bonds = Vec::with_capacity(n_bonds);
        for i in first_bond_line..last_bond_line {
            let line = lines[i];
            let cols: Vec<&str> = line.split_whitespace().collect();

            if cols.len() < 3 {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Bond line {i} does not have enough columns"),
                ));
            }

            let mut atom_0_sn = cols[0].parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse bond atom 0")
            })?;
            let mut atom_1_sn = cols[1].parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse bond atom 1")
            })?;

            // See note above on atom/bond counts. We observe the same thing from PubChem in some
            // cases for the atom SNs listed in bond lines; space omitted if len of atom0 is 3.
            let bond_type = if cols[2] == "0" {
                let prev = if i != first_bond_line {
                    let line_p = lines[i - 1];
                    let cols_p: Vec<&str> = line_p.split_whitespace().collect();
                    Some(cols_p[0])
                } else {
                    None
                };
                let next = if i != last_bond_line - 1 {
                    let line_p = lines[i + 1];
                    let cols_p: Vec<&str> = line_p.split_whitespace().collect();
                    Some(cols_p[0])
                } else {
                    None
                };

                let (atom_0, atom_1) = split_bond_indices_line(cols[0], prev, next, atoms.len())?;
                atom_0_sn = atom_0 as u32;
                atom_1_sn = atom_1 as u32;

                if atom_0_sn as usize > atoms.len() || atom_1_sn as usize > atoms.len() {
                    return Err(io::Error::new(
                        ErrorKind::InvalidData,
                        format!(
                            "Bond indices out of bounds, during split: {}. Atom len: {}. A0: {atom_0_sn} A1: {atom_1_sn}",
                            cols[0],
                            atoms.len(),
                        ),
                    ));
                }
                // }

                BondType::from_str(cols[1])?
            } else {
                BondType::from_str(cols[2])?
            };

            bonds.push(BondGeneric {
                atom_0_sn,
                atom_1_sn,
                bond_type,
            })
        }

        // Look for a molecule identifier in the file. Check for either
        // "> <PUBCHEM_COMPOUND_CID>" or "> <DATABASE_ID>" and take the next nonempty line.
        let mut _pubchem_cid = None;
        let mut _drugbank_id = None;

        for (i, line) in lines.iter().enumerate() {
            if line.contains("> <PUBCHEM_COMPOUND_CID>")
                && let Some(value_line) = lines.get(i + 1)
            {
                let value = value_line.trim();
                if let Ok(v) = value.parse::<u32>() {
                    _pubchem_cid = Some(v);
                }
            }
            if line.contains("> <DATABASE_ID>")
                && let Some(value_line) = lines.get(i + 1)
            {
                let value = value_line.trim();
                if !value.is_empty() {
                    _drugbank_id = Some(value.to_string());
                }
            }
        }

        let ident = lines[0].trim().to_string();
        // We observe that on at least some DrugBank files, this line
        // is the PubChem ID, even if the PUBCHEM_COMPOUND_CID line is omitted.
        if let Ok(v) = lines[0].parse::<u32>() {
            _pubchem_cid = Some(v);
        }

        // We could now skip over the bond lines if we want:
        //   let first_bond_line = last_atom_ line;
        //   let last_bond_line = first_bond_line + nbonds;
        // etc.
        // Then we look for "M  END" or the data fields, etc.

        // For now, just return the Sdf with the atoms we parsed:

        let mut chains = Vec::new();
        let mut residues = Vec::new();

        // let atom_indices: Vec<usize> = (0..atoms.len()).collect();
        let atom_sns: Vec<_> = atoms.iter().map(|a| a.serial_number).collect();

        residues.push(ResidueGeneric {
            serial_number: 0,
            res_type: ResidueType::Other("Unknown".to_string()),
            atom_sns: atom_sns.clone(),
            end: ResidueEnd::Hetero,
        });

        chains.push(ChainGeneric {
            id: "A".to_string(),
            residue_sns: vec![0],
            atom_sns,
        });

        // Load metadata. We use a separate pass for simplicity, although this is a bit slower.
        let metadata = {
            let mut md: HashMap<String, String> = HashMap::new();

            let mut idx = if let Some(m_end) = lines.iter().position(|l| l.trim() == "M  END") {
                m_end + 1
            } else {
                last_bond_line
            };

            while idx < lines.len() {
                let line = lines[idx].trim();
                if line == "$$$$" {
                    break;
                }
                if line.starts_with('>')
                    && let (Some(l), Some(r)) = (line.find('<'), line.rfind('>'))
                    && r > l + 1
                {
                    let key = &line[l + 1..r];
                    idx += 1;

                    let mut vals: Vec<&str> = Vec::new();
                    while idx < lines.len() {
                        let v = lines[idx];
                        let v_trim = v.trim_end();
                        if v_trim.is_empty() || v_trim == "$$$$" || v_trim.starts_with("> <") {
                            break;
                        }
                        vals.push(v_trim);
                        idx += 1;
                    }
                    if key == "PUBCHEM_PHARMACOPHORE_FEATURES" {
                        match parse_pubchem_pharmacophore_features(&vals) {
                            Ok(v) => pharmacophore_features = v,
                            Err(e) => {
                                md.insert(key.to_string(), vals.join("\n"));
                                eprintln!("Failed to parse PUBCHEM_PHARMACOPHORE_FEATURES: {e}");
                            }
                        }
                    } else {
                        md.insert(key.to_string(), vals.join("\n"));
                    }

                    // OpenFF format.
                    if key == "atom.dprop.PartialCharge" {
                        let joined = vals.join(" ");
                        let charges: Vec<&str> = joined.split_whitespace().collect();

                        for (i, q) in charges.into_iter().enumerate() {
                            if i < atoms.len() {
                                atoms[i].partial_charge = Some(q.parse().unwrap_or(0.));
                            }
                        }
                    }

                    // Pubchem format.
                    if key == "PUBCHEM_MMFF94_PARTIAL_CHARGES" {
                        if vals.is_empty() {
                            eprintln!("No values for PUBCHEM_MMFF94_PARTIAL_CHARGES");
                        } else {
                            let n = vals[0].trim().parse::<usize>().unwrap_or(0);
                            if vals.len().saturating_sub(1) != n {
                                eprintln!(
                                    "Charge count mismatch: expected {}, got {}",
                                    n,
                                    vals.len().saturating_sub(1)
                                );
                            }
                            for line in vals.iter().skip(1).take(n) {
                                let mut it = line.split_whitespace();
                                let i1 = it.next().and_then(|s| s.parse::<usize>().ok());
                                let q = it.next().and_then(|s| s.parse::<f32>().ok());

                                if let (Some(i1), Some(q)) = (i1, q) {
                                    if (1..=atoms.len()).contains(&i1) {
                                        atoms[i1 - 1].partial_charge = Some(q); // 1-based -> 0-based
                                    } else {
                                        eprintln!(
                                            "Atom index {} out of range (n_atoms={})",
                                            i1,
                                            atoms.len()
                                        );
                                    }
                                }
                            }
                        }
                    }

                    continue;
                }
                idx += 1;
            }
            md
        };

        Ok(Self {
            ident,
            metadata,
            atoms,
            chains,
            residues,
            bonds,
            pharmacophore_features,
        })
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut file = File::create(path)?;

        // 1) Title line (often the first line in SDF).
        //    We use the molecule's name/identifier here:
        // todo: There is a subtlety here. Add that to your parser as well. There are two values
        // todo in the files we have; this top ident is not the DB id.
        writeln!(file, "{}", self.ident)?;

        // 2) Write two blank lines:
        writeln!(file)?;
        writeln!(file)?;

        let natoms = self.atoms.len();
        let nbonds = self.bonds.len();

        // Format the counts line. We loosely mimic typical spacing,
        // though it's not strictly required to line up exactly.
        writeln!(
            file,
            "{:>3}{:>3}  0  0  0  0           0999 V2000",
            natoms, nbonds
        )?;

        for atom in &self.atoms {
            let x = atom.posit.x;
            let y = atom.posit.y;
            let z = atom.posit.z;
            let symbol = atom.element.to_letter();

            // MDL v2000 format often uses fixed-width fields,
            // but for simplicity we use whitespace separation:
            writeln!(
                file,
                "{:>10.4}{:>10.4}{:>10.4} {:<2}  0  0  0  0  0  0  0  0  0  0",
                x, y, z, symbol
            )?;
        }

        for bond in &self.bonds {
            writeln!(
                file,
                "{:>3}{:>3}{:>3}  0  0  0  0",
                bond.atom_0_sn,
                bond.atom_1_sn,
                bond.bond_type.to_str_sdf()
            )?;
        }

        writeln!(file, "M  END")?;

        for m in &self.metadata {
            write_metadata(m.0, m.1, &mut file)?;
        }

        // If partial charges are available, write them to metadata. This is an OpenFF convention.
        // todo: Should we use the Pubhcem format instead? e.g.
        // > <PUBCHEM_MMFF94_PARTIAL_CHARGES>
        // 16
        // 1 -0.53
        // 10 0.63
        // 11 0.15
        // 12 0.15
        // 13 0.15
        // 14 0.15
        // 15 0.45
        // 16 0.5
        // 2 -0.65
        // 3 -0.57
        // 4 0.09
        // 5 -0.15
        // 6 -0.15
        // 7 0.08
        // 8 -0.15
        // 9 -0.15

        let mut partial_charges = Vec::new();
        let mut all_partial_charges_present = true;
        for atom in &self.atoms {
            match atom.partial_charge {
                Some(q) => partial_charges.push(q),
                None => {
                    all_partial_charges_present = false;
                    break;
                }
            }
        }

        // Pubchem format.
        if all_partial_charges_present {
            let charges_formated: Vec<_> =
                partial_charges.iter().map(|q| format!("{q:.8}")).collect();
            let charge_str = charges_formated.join(" ");
            write_metadata("atom.dprop.PartialCharge", &charge_str, &mut file)?;
        }

        if !self.pharmacophore_features.is_empty() {
            let v = format_pubchem_pharmacophore_features(&self.pharmacophore_features);
            write_metadata("PUBCHEM_PHARMACOPHORE_FEATURES", &v, &mut file)?;
        }

        // End of this molecule record in SDF
        writeln!(file, "$$$$")?;

        Ok(())
    }

    // todo: Generic fn for this and save, among all text-based types.
    pub fn load(path: &Path) -> io::Result<Self> {
        let data_str = fs::read_to_string(path)?;
        Self::new(&data_str)
    }

    /// Download from DrugBank from a Drugbank ID.
    pub fn load_drugbank(ident: &str) -> io::Result<Self> {
        let data_str = drugbank::load_sdf(ident)
            .map_err(|e| io::Error::other(format!("Error loading: {e:?}")))?;
        Self::new(&data_str)
    }

    /// Download from PubChem from a CID.
    pub fn load_pubchem(id_type: StructureSearchNamespace, id: &str) -> io::Result<Self> {
        let data_str = pubchem::load_sdf(id_type, id)
            .map_err(|e| io::Error::other(format!("Error loading: {e:?}")))?;
        Self::new(&data_str)
    }

    /// Download from PDBe from a PDBe ID.
    pub fn load_pdbe(ident: &str) -> io::Result<Self> {
        let data_str =
            pdbe::load_sdf(ident).map_err(|e| io::Error::other(format!("Error loading: {e:?}")))?;
        Self::new(&data_str)
    }
}

impl From<Mol2> for Sdf {
    fn from(m: Mol2) -> Self {
        Self {
            ident: m.ident.clone(),
            metadata: m.metadata.clone(),
            atoms: m.atoms.clone(),
            bonds: m.bonds.clone(),
            chains: Vec::new(),
            residues: Vec::new(),
            pharmacophore_features: Vec::new(),
        }
    }
}

fn write_metadata(key: &str, val: &str, file: &mut File) -> io::Result<()> {
    writeln!(file, "> <{key}>")?;
    writeln!(file, "{val}")?;
    writeln!(file)?; // blank line

    Ok(())
}

fn parse_pubchem_pharmacophore_features(vals: &[&str]) -> io::Result<Vec<PharmacaphoreFeatures>> {
    if vals.is_empty() {
        return Ok(Vec::new());
    }

    let n = vals[0].trim().parse::<usize>().map_err(|_| {
        io::Error::new(
            ErrorKind::InvalidData,
            "Bad PUBCHEM_PHARMACOPHORE_FEATURES count",
        )
    })?;

    let mut out = Vec::with_capacity(n);

    for line in vals.iter().skip(1).take(n) {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 2 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Bad PUBCHEM_PHARMACOPHORE_FEATURES line",
            ));
        }

        let k = toks[0].parse::<usize>().map_err(|_| {
            io::Error::new(
                ErrorKind::InvalidData,
                "Bad PUBCHEM_PHARMACOPHORE_FEATURES atom-count",
            )
        })?;

        if toks.len() < 1 + k + 1 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "PUBCHEM_PHARMACOPHORE_FEATURES line missing atoms/type",
            ));
        }

        let mut atom_sns = Vec::with_capacity(k);
        for t in toks.iter().skip(1).take(k) {
            let sn = t.parse::<u32>().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    "Bad PUBCHEM_PHARMACOPHORE_FEATURES atom index",
                )
            })?;
            atom_sns.push(sn);
        }

        let type_str = toks[1 + k..].join(" ");
        let type_ = PharmacophoreType::from_pubchem_str(&type_str).ok_or_else(|| {
            io::Error::new(
                ErrorKind::InvalidData,
                format!("Unknown pharmacophore type: {type_str}"),
            )
        })?;

        out.push(PharmacaphoreFeatures { atom_sns, type_ });
    }

    Ok(out)
}

fn format_pubchem_pharmacophore_features(features: &[PharmacaphoreFeatures]) -> String {
    let mut s = String::new();
    s.push_str(&format!("{}\n", features.len()));

    for f in features {
        s.push_str(&format!("{}", f.atom_sns.len()));
        for sn in &f.atom_sns {
            s.push_str(&format!(" {sn}"));
        }
        s.push(' ');
        s.push_str(f.type_.to_pubchem_str());
        s.push('\n');
    }

    s
}
