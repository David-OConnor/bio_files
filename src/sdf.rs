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

const PERIODIC_TABLE_SYMBOLS: &[&str] = &[
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

/// Preserve records containing a real element that `na_seq::Element` does not model yet. The
/// generic atom representation cannot retain the original symbol in that case, so it uses
/// `Element::Other`; malformed and SDF query-atom symbols still produce an error.
fn parse_sdf_element(symbol: &str) -> io::Result<Element> {
    match Element::from_letter(symbol) {
        Ok(element) => Ok(element),
        Err(_)
            if PERIODIC_TABLE_SYMBOLS
                .iter()
                .any(|candidate| candidate.eq_ignore_ascii_case(symbol)) =>
        {
            Ok(Element::Other)
        }
        Err(error) => Err(error),
    }
}

/// MDL V2000 also defines query bond types 5-8. `BondType` cannot retain those distinctions, so
/// keep the molecule and represent the ambiguous bond as unknown instead of rejecting the record.
fn parse_sdf_bond_type(value: &str) -> io::Result<BondType> {
    match value.trim() {
        "5" | "6" | "7" | "8" => Ok(BondType::Unknown),
        _ => BondType::from_str(value),
    }
}

fn parse_v2000_counts(line: &str) -> io::Result<(usize, usize)> {
    let cols: Vec<&str> = line.split_whitespace().collect();

    // Accommodate relaxed writers (including our historical writer) that put explicit whitespace
    // between the count fields. A joined token larger than V2000's three-column maximum indicates
    // that the two standard fixed-width fields touch, so parse by columns below.
    if cols.len() >= 2
        && let (Ok(n_atoms), Ok(n_bonds)) = (cols[0].parse::<usize>(), cols[1].parse::<usize>())
        && n_atoms <= 999
        && n_bonds <= 999
    {
        return Ok((n_atoms, n_bonds));
    }

    let n_atoms = line
        .get(0..3)
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Counts line is too short"))?
        .trim()
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse number of atoms"))?;
    let n_bonds = line
        .get(3..6)
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Counts line is too short"))?
        .trim()
        .parse::<usize>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse number of bonds"))?;

    Ok((n_atoms, n_bonds))
}

fn parse_v2000_atom_fields(line: &str) -> io::Result<(f64, f64, f64, &str)> {
    let cols: Vec<&str> = line.split_whitespace().collect();
    if cols.len() >= 4
        && let (Ok(x), Ok(y), Ok(z)) = (
            cols[0].parse::<f64>(),
            cols[1].parse::<f64>(),
            cols[2].parse::<f64>(),
        )
    {
        return Ok((x, y, z, cols[3]));
    }

    // Coordinates occupy three adjacent ten-character fields in V2000. Large negative values can
    // fill a field completely, leaving no whitespace for `split_whitespace` to find.
    let x = line
        .get(0..10)
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Atom line is too short"))?
        .trim()
        .parse::<f64>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse X coordinate"))?;
    let y = line
        .get(10..20)
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Atom line is too short"))?
        .trim()
        .parse::<f64>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse Y coordinate"))?;
    let z = line
        .get(20..30)
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Atom line is too short"))?
        .trim()
        .parse::<f64>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse Z coordinate"))?;
    let element = line
        .get(31..34)
        .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Atom line is too short"))?
        .trim();
    if element.is_empty() {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "Atom line has no element",
        ));
    }

    Ok((x, y, z, element))
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

fn parse_v2000_bond_fields(
    line: &str,
    prev: Option<&str>,
    next: Option<&str>,
    max_atom_index: usize,
) -> io::Result<(u32, u32, BondType)> {
    let cols: Vec<&str> = line.split_whitespace().collect();

    // Prefer ordinary whitespace-delimited rows so relaxed writers remain compatible.
    if cols.len() >= 3
        && let (Ok(atom_0_sn), Ok(atom_1_sn), Ok(bond_type)) = (
            cols[0].parse::<u32>(),
            cols[1].parse::<u32>(),
            parse_sdf_bond_type(cols[2]),
        )
        && atom_0_sn > 0
        && atom_1_sn > 0
        && atom_0_sn as usize <= max_atom_index
        && atom_1_sn as usize <= max_atom_index
    {
        return Ok((atom_0_sn, atom_1_sn, bond_type));
    }

    // The two atom indices and bond type are adjacent three-character fields in V2000. When an
    // index reaches three digits there may be no separating whitespace.
    if let (Some(atom_0), Some(atom_1), Some(bond_type)) =
        (line.get(0..3), line.get(3..6), line.get(6..9))
        && let (Ok(atom_0_sn), Ok(atom_1_sn), Ok(bond_type)) = (
            atom_0.trim().parse::<u32>(),
            atom_1.trim().parse::<u32>(),
            parse_sdf_bond_type(bond_type),
        )
        && atom_0_sn > 0
        && atom_1_sn > 0
        && atom_0_sn as usize <= max_atom_index
        && atom_1_sn as usize <= max_atom_index
    {
        return Ok((atom_0_sn, atom_1_sn, bond_type));
    }

    // Retain compatibility with observed non-fixed-width PubChem rows whose two indices are joined.
    if cols.len() >= 2 {
        let prev_col = prev.and_then(|row| row.split_whitespace().next());
        let next_col = next.and_then(|row| row.split_whitespace().next());
        let (atom_0_sn, atom_1_sn) =
            split_bond_indices_line(cols[0], prev_col, next_col, max_atom_index)?;
        let bond_type = parse_sdf_bond_type(cols[1])?;
        return Ok((atom_0_sn as u32, atom_1_sn as u32, bond_type));
    }

    Err(io::Error::new(
        ErrorKind::InvalidData,
        "Bond line does not have enough fields",
    ))
}

/// Parse V2000 atom and bond blocks. Returns `(atoms, bonds, last_bond_line)`.
/// `last_bond_line` is used as a fallback start index for metadata if `M  END` is absent.
fn parse_v2000_ctab(
    lines: &[&str],
    counts_line_index: usize,
) -> io::Result<(Vec<AtomGeneric>, Vec<BondGeneric>, usize)> {
    let counts_line = lines[counts_line_index];
    let (n_atoms, n_bonds) = parse_v2000_counts(counts_line)?;

    let first_atom_line = counts_line_index + 1;
    let last_atom_line = first_atom_line + n_atoms;
    let first_bond_line = last_atom_line;
    let last_bond_line = first_bond_line + n_bonds;

    if lines.len() < last_atom_line {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            format!(
                "Not enough lines for the declared atom block: found {}, expected {n_atoms} atoms",
                lines.len().saturating_sub(first_atom_line),
            ),
        ));
    }
    if lines.len() < last_bond_line {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            format!(
                "Not enough lines for the declared bond block: found {}, expected {n_bonds} bonds",
                lines.len().saturating_sub(first_bond_line),
            ),
        ));
    }

    let mut atoms = Vec::with_capacity(n_atoms);
    for (atom_index, line) in lines[first_atom_line..last_atom_line].iter().enumerate() {
        let (x, y, z, element) = parse_v2000_atom_fields(line)?;
        atoms.push(AtomGeneric {
            serial_number: atom_index as u32 + 1,
            posit: Vec3 { x, y, z },
            element: parse_sdf_element(element)?,
            hetero: true,
            ..Default::default()
        });
    }

    let mut bonds = Vec::with_capacity(n_bonds);
    for i in first_bond_line..last_bond_line {
        let prev = (i != first_bond_line).then(|| lines[i - 1]);
        let next = (i + 1 != last_bond_line).then(|| lines[i + 1]);
        let (atom_0_sn, atom_1_sn, bond_type) =
            parse_v2000_bond_fields(lines[i], prev, next, atoms.len())?;

        bonds.push(BondGeneric {
            atom_0_sn,
            atom_1_sn,
            bond_type,
        });
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
            element: parse_sdf_element(element_str)?,
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

        let bond_type = parse_sdf_bond_type(cols[3])?;
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

        if lines.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines to parse an SDF header",
            ));
        }

        let ident = lines[0].trim().to_string();

        // A conforming molfile uses three header lines, but some catalog producers omit the title
        // or comment line. Find the stamped counts line in the small header region while retaining
        // the historical line-four fallback for unstamped V2000 input.
        let (counts_line_index, format) = lines
            .iter()
            .take(8)
            .enumerate()
            .find_map(|(index, line)| {
                if line.contains("V3000") {
                    Some((index, SdfFormat::V3000))
                } else if line.contains("V2000") {
                    Some((index, SdfFormat::V2000))
                } else {
                    None
                }
            })
            .or_else(|| (lines.len() > 3).then_some((3, SdfFormat::V2000)))
            .ok_or_else(|| {
                io::Error::new(ErrorKind::InvalidData, "Could not find an SDF counts line")
            })?;

        let (mut atoms, bonds, fallback_start) = match format {
            SdfFormat::V2000 => parse_v2000_ctab(&lines, counts_line_index)?,
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
        self.write(&mut file, format)
    }

    /// Write one molecule record, including its trailing `$$$$`. Factored out of [`Sdf::save`]
    /// so that multiple molecules can share a file.
    fn write(&self, file: &mut File, format: SdfFormat) -> io::Result<()> {
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
            write_metadata(m.0, m.1, file)?;
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
            write_metadata("atom.dprop.PartialCharge", &charge_str, file)?;
        }

        if !self.pharmacophore_features.is_empty() {
            let v = format_pharmacophore_features(&self.pharmacophore_features);
            write_metadata("PUBCHEM_PHARMACOPHORE_FEATURES", &v, file)?;
        }

        // End of this molecule record in SDF
        writeln!(file, "$$$$")?;

        Ok(())
    }

    /// From a string that may hold any number of molecules, each terminated by a `$$$$` line.
    /// A single-molecule SDF parses as a one-element `Vec`, with or without the terminator.
    ///
    /// Records that fail to parse are skipped with a warning, so that one malformed entry in a
    /// large catalog file doesn't discard the rest.
    pub fn new_multi(text: &str) -> io::Result<Vec<Self>> {
        let mut result = Vec::new();
        let mut record: Vec<&str> = Vec::new();
        let mut record_number = 1;

        for line in text.lines() {
            if line.trim_end() == "$$$$" {
                push_record(&mut result, &record, record_number);
                record.clear();
                record_number += 1;
            } else {
                record.push(line);
            }
        }

        // A trailing record without the `$$$$` terminator.
        push_record(&mut result, &record, record_number);

        if result.is_empty() {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "No molecules parsed from the SDF",
            ));
        }

        Ok(result)
    }

    // todo: Generic fn for this and save, among all text-based types.
    pub fn load(path: &Path) -> io::Result<Self> {
        let data_str = fs::read_to_string(path)?;
        Self::new(&data_str)
    }

    /// Load an SDF file that may hold more than one molecule. See [`Sdf::new_multi`].
    pub fn load_multi(path: &Path) -> io::Result<Vec<Self>> {
        let data_str = fs::read_to_string(path)?;
        Self::new_multi(&data_str)
    }

    /// Save any number of molecules to a single file, each terminated by `$$$$`.
    pub fn save_multi(mols: &[Self], path: &Path, format: SdfFormat) -> io::Result<()> {
        let mut file = File::create(path)?;
        for mol in mols {
            mol.write(&mut file, format)?;
        }

        Ok(())
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

/// Parse one `$$$$`-delimited record and append it. Blank trailing sections (the text after the
/// final `$$$$`) are ignored; a record that fails to parse is reported and skipped.
fn push_record(out: &mut Vec<Sdf>, record: &[&str], record_number: usize) {
    if record.iter().all(|l| l.trim().is_empty()) {
        return;
    }

    match Sdf::new(&record.join("\n")) {
        Ok(mol) => out.push(mol),
        Err(e) => eprintln!(
            "Warning: Skipping SDF record {} ({:?}): {e}",
            record_number,
            record.first().unwrap_or(&"").trim(),
        ),
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
