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
    AtomGeneric, BondGeneric, BondType, ChainGeneric, Mol2, PharmacophoreFeatureGeneric,
    PharmacophoreTypeGeneric, ResidueEnd, ResidueGeneric, ResidueType,
};

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
    pub pharmacophore_features: Vec<PharmacophoreFeatureGeneric>,
}

/// https://en.wikipedia.org/wiki/Chemical_table_file
#[derive(Clone, Copy, Debug, Default)]
pub enum SdfFormat {
    #[default]
    /// Example: DrugBank, PubChem
    V2000,
    /// Example: Chemdraw
    V3000,
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
                if *n1 >= p_val {
                    n1 - p_val
                } else {
                    // If n1 < prev, it's out of order. Penalize it heavily so we prefer sorted options.
                    // But don't make it impossible, just unlikely.
                    (p_val - n1) + 10000
                }
            });

            // Return the "closest" valid split
            Ok(valid_splits[0])
        }
    }
}

/// Parse V2000 atom and bond blocks. Returns `(atoms, bonds, last_bond_line)`.
/// `last_bond_line` is used as a fallback start index for metadata if `M  END` is absent.
fn parse_v2000_ctab(lines: &[&str]) -> io::Result<(Vec<AtomGeneric>, Vec<BondGeneric>, usize)> {
    let counts_line = lines[3];
    let counts_cols: Vec<&str> = counts_line.split_whitespace().collect();

    if counts_cols.len() < 2 {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "Counts line doesn't have enough fields",
        ));
    }

    let mut n_atoms = counts_cols[0]
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse number of atoms"))?;

    let mut n_bonds = counts_cols[1]
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse number of bonds"))?;

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

    #[allow(clippy::needless_range_loop)]
    for i in first_atom_line..last_atom_line {
        let line = lines[i];
        let cols: Vec<&str> = line.split_whitespace().collect();

        if cols.len() < 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Atom line {i} does not have enough columns"),
            ));
        }

        let x = cols[0]
            .parse::<f64>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse X coordinate"))?;
        let y = cols[1]
            .parse::<f64>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse Y coordinate"))?;
        let z = cols[2]
            .parse::<f64>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse Z coordinate"))?;
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

        let mut atom_0_sn = cols[0]
            .parse::<u32>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse bond atom 0"))?;
        let mut atom_1_sn = cols[1]
            .parse::<u32>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse bond atom 1"))?;

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

    Ok((atoms, bonds, last_bond_line))
}

/// Parse V3000 atom and bond blocks from the CTAB section.
fn parse_v3000_ctab(lines: &[&str]) -> io::Result<(Vec<AtomGeneric>, Vec<BondGeneric>)> {
    // Find BEGIN CTAB
    let ctab_begin = lines
        .iter()
        .position(|l| l.trim() == "M  V30 BEGIN CTAB")
        .ok_or_else(|| {
            io::Error::new(ErrorKind::InvalidData, "V3000: no M  V30 BEGIN CTAB found")
        })?;

    // Find COUNTS line: "M  V30 COUNTS na nb nsg n3d chiral"
    let counts_idx = lines[ctab_begin..]
        .iter()
        .position(|l| {
            let t = l.trim();
            t.starts_with("M  V30 COUNTS") || t.starts_with("M V30 COUNTS")
        })
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "V3000: no COUNTS line"))?
        + ctab_begin;

    let counts_cols: Vec<&str> = lines[counts_idx].split_whitespace().collect();
    // ["M", "V30", "COUNTS", na, nb, nsg, n3d, chiral, ...]
    if counts_cols.len() < 5 {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "V3000 COUNTS line too short",
        ));
    }
    let n_atoms = counts_cols[3]
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: could not parse atom count"))?;
    let n_bonds = counts_cols[4]
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: could not parse bond count"))?;

    // Find ATOM block
    let atom_begin = lines[ctab_begin..]
        .iter()
        .position(|l| l.trim() == "M  V30 BEGIN ATOM")
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "V3000: no BEGIN ATOM"))?
        + ctab_begin
        + 1; // +1: skip the BEGIN ATOM line itself

    let atom_end = lines[atom_begin..]
        .iter()
        .position(|l| l.trim() == "M  V30 END ATOM")
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "V3000: no END ATOM"))?
        + atom_begin;

    let mut atoms = Vec::with_capacity(n_atoms);
    for i in atom_begin..atom_end {
        let line = lines[i];
        // Format: "M  V30 idx elem x y z map_no [KEY=VAL ...]"
        // After split_whitespace: ["M", "V30", idx, elem, x, y, z, map_no, ...]
        //                           0    1      2    3     4  5  6    7
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 7 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("V3000 atom line {i} too short (got {} cols)", cols.len()),
            ));
        }

        let serial_number = cols[2]
            .parse::<u32>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: bad atom serial number"))?;
        let element_str = cols[3];
        let x = cols[4]
            .parse::<f64>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: bad X coordinate"))?;
        let y = cols[5]
            .parse::<f64>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: bad Y coordinate"))?;
        let z = cols[6]
            .parse::<f64>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: bad Z coordinate"))?;

        atoms.push(AtomGeneric {
            serial_number,
            posit: Vec3 { x, y, z },
            element: Element::from_letter(element_str)?,
            hetero: true,
            ..Default::default()
        });
    }

    // Find BOND block
    let bond_begin = lines[ctab_begin..]
        .iter()
        .position(|l| l.trim() == "M  V30 BEGIN BOND")
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "V3000: no BEGIN BOND"))?
        + ctab_begin
        + 1;

    let bond_end = lines[bond_begin..]
        .iter()
        .position(|l| l.trim() == "M  V30 END BOND")
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "V3000: no END BOND"))?
        + bond_begin;

    let mut bonds = Vec::with_capacity(n_bonds);
    for i in bond_begin..bond_end {
        let line = lines[i];
        // Format: "M  V30 idx bond_type atom1 atom2 [KEY=VAL ...]"
        // After split_whitespace: ["M", "V30", idx, bond_type, atom1, atom2, ...]
        //                           0    1      2       3         4      5
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 6 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("V3000 bond line {i} too short (got {} cols)", cols.len()),
            ));
        }

        let bond_type = BondType::from_str(cols[3])?;
        let atom_0_sn = cols[4]
            .parse::<u32>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: bad bond atom 0"))?;
        let atom_1_sn = cols[5]
            .parse::<u32>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "V3000: bad bond atom 1"))?;

        bonds.push(BondGeneric {
            atom_0_sn,
            atom_1_sn,
            bond_type,
        });
    }

    Ok((atoms, bonds))
}

/// Parse metadata data fields (the `> <KEY>` sections after `M  END`) and apply any
/// partial-charge fields to `atoms` in place.
///
/// `fallback_start` is the line index to begin scanning from if `M  END` is not found
/// (used as a V2000 compatibility measure; pass `lines.len()` for V3000).
fn parse_metadata_section(
    lines: &[&str],
    atoms: &mut Vec<AtomGeneric>,
    fallback_start: usize,
) -> io::Result<(HashMap<String, String>, Vec<PharmacophoreFeatureGeneric>)> {
    let mut metadata: HashMap<String, String> = HashMap::new();
    let mut pharmacophore_features: Vec<PharmacophoreFeatureGeneric> = Vec::new();

    // Look for molecule identifiers in the data fields.
    for (i, line) in lines.iter().enumerate() {
        if line.contains("> <PUBCHEM_COMPOUND_CID>")
            && let Some(value_line) = lines.get(i + 1)
        {
            let value = value_line.trim();
            let _ = value.parse::<u32>(); // retain for future use
        }
        if line.contains("> <DATABASE_ID>")
            && let Some(value_line) = lines.get(i + 1)
        {
            let value = value_line.trim();
            let _ = (!value.is_empty()).then(|| value.to_string()); // retain for future use
        }
    }

    let mut idx = if let Some(m_end) = lines.iter().position(|l| l.trim() == "M  END") {
        m_end + 1
    } else {
        fallback_start
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

            let mut rows_pharm: Vec<&str> = Vec::new();
            while idx < lines.len() {
                let v = lines[idx];
                let v_trim = v.trim_end();
                if v_trim.is_empty() || v_trim == "$$$$" || v_trim.starts_with("> <") {
                    break;
                }
                rows_pharm.push(v_trim);
                idx += 1;
            }
            if key == "PUBCHEM_PHARMACOPHORE_FEATURES" {
                match parse_pharmacophore_features(&rows_pharm) {
                    Ok(v) => pharmacophore_features = v,
                    Err(e) => {
                        metadata.insert(key.to_string(), rows_pharm.join("\n"));
                        eprintln!("Failed to parse PUBCHEM_PHARMACOPHORE_FEATURES: {e}");
                    }
                }
            } else {
                metadata.insert(key.to_string(), rows_pharm.join("\n"));
            }

            // OpenFF format.
            if key == "atom.dprop.PartialCharge" {
                let joined = rows_pharm.join(" ");
                let charges: Vec<&str> = joined.split_whitespace().collect();

                for (i, q) in charges.into_iter().enumerate() {
                    if i < atoms.len() {
                        atoms[i].partial_charge = Some(q.parse().unwrap_or(0.));
                    }
                }
            }

            // Pubchem format.
            if key == "PUBCHEM_MMFF94_PARTIAL_CHARGES" {
                if rows_pharm.is_empty() {
                    eprintln!("No values for PUBCHEM_MMFF94_PARTIAL_CHARGES");
                } else {
                    let n = rows_pharm[0].trim().parse::<usize>().unwrap_or(0);
                    if rows_pharm.len().saturating_sub(1) != n {
                        eprintln!(
                            "Charge count mismatch: expected {}, got {}",
                            n,
                            rows_pharm.len().saturating_sub(1)
                        );
                    }
                    for line in rows_pharm.iter().skip(1).take(n) {
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

    Ok((metadata, pharmacophore_features))
}

impl Sdf {
    /// From a string of an SDF text file. Automatically detects V2000 or V3000 format.
    pub fn new(text: &str) -> io::Result<Self> {
        let lines: Vec<&str> = text.lines().collect();

        // SDF files typically have at least 4 lines before the atom block:
        //   1) A title or identifier
        //   2) Usually blank or comments
        //   3) Often blank or comments
        //   4) "counts" line: e.g. " 50  50  0  ..." for V2000,
        //      or "  0  0  0     0  0              0 V3000" for V3000
        if lines.len() < 4 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines to parse an SDF header",
            ));
        }

        let ident = lines[0].trim().to_string();

        // Detect format from the counts line (line index 3).
        let format = if lines[3].contains("V3000") {
            SdfFormat::V3000
        } else {
            SdfFormat::V2000
        };

        let (mut atoms, bonds, fallback_start) = match format {
            SdfFormat::V2000 => parse_v2000_ctab(&lines)?,
            SdfFormat::V3000 => {
                let (a, b) = parse_v3000_ctab(&lines)?;
                // For V3000, M  END is always present; pass lines.len() as an unreachable fallback.
                (a, b, lines.len())
            }
        };

        let (metadata, pharmacophore_features) =
            parse_metadata_section(&lines, &mut atoms, fallback_start)?;

        let atom_sns: Vec<_> = atoms.iter().map(|a| a.serial_number).collect();

        let mut chains = Vec::new();
        let mut residues = Vec::new();

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

    pub fn save(&self, path: &Path, format: SdfFormat) -> io::Result<()> {
        let mut file = File::create(path)?;

        // 1) Title line (often the first line in SDF).
        writeln!(file, "{}", self.ident)?;

        // 2) Two blank lines (program/timestamp and comment lines).
        writeln!(file)?;
        writeln!(file)?;

        let natoms = self.atoms.len();
        let nbonds = self.bonds.len();

        match format {
            SdfFormat::V2000 => {
                // Counts line with V2000 stamp.
                // Format the counts line. We loosely mimic typical spacing,
                // though it's not strictly required to line up exactly.
                let w = natoms.max(nbonds).to_string().len().max(3) + 1; // +1 keeps a guaranteed gap between the two fields

                writeln!(
                    file,
                    "{:>w$}{:>w$}  0  0  0  0           0999 V2000",
                    natoms,
                    nbonds,
                    w = w,
                )?;

                for atom in &self.atoms {
                    let x = atom.posit.x;
                    let y = atom.posit.y;
                    let z = atom.posit.z;
                    let symbol = atom.element.to_letter();

                    writeln!(
                        file,
                        "{:>10.4}{:>10.4}{:>10.4} {:<2}  0  0  0  0  0  0  0  0  0  0",
                        x, y, z, symbol
                    )?;
                }

                // For formatting the bond atom SNs correctly
                let max_sn = self
                    .bonds
                    .iter()
                    .flat_map(|b| [b.atom_0_sn, b.atom_1_sn])
                    .max()
                    .unwrap_or(0);

                let w = max_sn.to_string().len().max(3) + 1; // +1 ensures a gap even at max width

                for bond in &self.bonds {
                    writeln!(
                        file,
                        "{:>w$}{:>w$}{:>3}  0  0  0  0",
                        bond.atom_0_sn,
                        bond.atom_1_sn,
                        bond.bond_type.to_str_sdf(),
                        w = w,
                    )?;
                }

                writeln!(file, "M  END")?;
            }

            SdfFormat::V3000 => {
                // V2000-compatibility counts line with V3000 stamp.
                // Counts are always 0 0 0 ... for V3000; atoms/bonds are in the CTAB block.
                writeln!(file, "  0  0  0     0  0              0 V3000")?;

                writeln!(file, "M  V30 BEGIN CTAB")?;
                writeln!(file, "M  V30 COUNTS {} {} 0 0 0", natoms, nbonds)?;

                writeln!(file, "M  V30 BEGIN ATOM")?;
                for (i, atom) in self.atoms.iter().enumerate() {
                    // Format: "M  V30 idx elem x y z map_no"
                    writeln!(
                        file,
                        "M  V30 {} {} {:.6} {:.6} {:.6} 0",
                        i + 1,
                        atom.element.to_letter(),
                        atom.posit.x,
                        atom.posit.y,
                        atom.posit.z,
                    )?;
                }
                writeln!(file, "M  V30 END ATOM")?;

                writeln!(file, "M  V30 BEGIN BOND")?;
                for (i, bond) in self.bonds.iter().enumerate() {
                    // Format: "M  V30 idx bond_type atom1 atom2"
                    writeln!(
                        file,
                        "M  V30 {} {} {} {}",
                        i + 1,
                        bond.bond_type.to_str_sdf(),
                        bond.atom_0_sn,
                        bond.atom_1_sn,
                    )?;
                }
                writeln!(file, "M  V30 END BOND")?;

                writeln!(file, "M  V30 END CTAB")?;
                writeln!(file, "M  END")?;
            }
        }

        // Metadata data fields are format-agnostic — they follow M  END in both V2000 and V3000.
        for m in &self.metadata {
            write_metadata(m.0, m.1, &mut file)?;
        }

        // If partial charges are available, write them to metadata. This is an OpenFF convention.
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

        if all_partial_charges_present {
            let charges_formated: Vec<_> =
                partial_charges.iter().map(|q| format!("{q:.8}")).collect();
            let charge_str = charges_formated.join(" ");
            write_metadata("atom.dprop.PartialCharge", &charge_str, &mut file)?;
        }

        if !self.pharmacophore_features.is_empty() {
            let v = format_pharmacophore_features(&self.pharmacophore_features);
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

pub(crate) fn parse_pharmacophore_features(
    rows: &[&str],
) -> io::Result<Vec<PharmacophoreFeatureGeneric>> {
    if rows.is_empty() {
        return Ok(Vec::new());
    }

    let n = rows[0].trim().parse::<usize>().map_err(|_| {
        io::Error::new(
            ErrorKind::InvalidData,
            "Bad PUBCHEM_PHARMACOPHORE_FEATURES count",
        )
    })?;

    let mut out = Vec::with_capacity(n);

    for line in rows.iter().skip(1).take(n) {
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 2 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Bad PUBCHEM_PHARMACOPHORE_FEATURES line",
            ));
        }

        let k = cols[0].parse::<usize>().map_err(|_| {
            io::Error::new(
                ErrorKind::InvalidData,
                "Bad PUBCHEM_PHARMACOPHORE_FEATURES atom-count",
            )
        })?;

        if cols.len() < 1 + k + 1 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "PUBCHEM_PHARMACOPHORE_FEATURES line missing atoms/type",
            ));
        }

        let mut atom_sns = Vec::with_capacity(k);
        for t in cols.iter().skip(1).take(k) {
            let sn = t.parse::<u32>().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    "Bad PUBCHEM_PHARMACOPHORE_FEATURES atom index",
                )
            })?;
            atom_sns.push(sn);
        }

        let type_str = cols[1 + k..].join(" ");
        let type_ = PharmacophoreTypeGeneric::from_pubchem_str(&type_str).ok_or_else(|| {
            io::Error::new(
                ErrorKind::InvalidData,
                format!("Unknown pharmacophore type: {type_str}"),
            )
        })?;

        out.push(PharmacophoreFeatureGeneric { atom_sns, type_ });
    }

    Ok(out)
}

pub(crate) fn format_pharmacophore_features(features: &[PharmacophoreFeatureGeneric]) -> String {
    let mut s = String::new();
    s.push_str(&format!("{}\n", features.len()));

    for f in features {
        s.push_str(&format!("{}", f.atom_sns.len()));
        for sn in &f.atom_sns {
            s.push_str(&format!(" {sn}"));
        }
        s.push(' ');
        s.push_str(&f.type_.clone().to_pubchem_str());
        s.push('\n');
    }

    s
}
