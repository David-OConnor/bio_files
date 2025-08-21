//! For opening Structure Data Format (SDF) files. These are common molecular descriptions for ligands. It's a simpler format
//! than PDB.

use std::{
    collections::HashMap,
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
    str::FromStr,
};

use lin_alg::f64::Vec3;
use na_seq::Element;

use crate::{AtomGeneric, BondGeneric, BondType, ChainGeneric, ResidueGeneric, ResidueType};

#[derive(Debug)]
pub struct Sdf {
    /// These fields aren't universal to the format.
    pub ident: String,
    pub metadata: HashMap<String, String>,
    pub atoms: Vec<AtomGeneric>,
    pub bonds: Vec<BondGeneric>,
    pub chains: Vec<ChainGeneric>,
    pub residues: Vec<ResidueGeneric>,
    pub pubchem_cid: Option<u32>,
    pub drugbank_id: Option<String>,
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
        let n_atoms = counts_cols[0].parse::<usize>().map_err(|_| {
            io::Error::new(ErrorKind::InvalidData, "Could not parse number of atoms")
        })?;
        let n_bonds = counts_cols[1].parse::<usize>().map_err(|_| {
            io::Error::new(ErrorKind::InvalidData, "Could not parse number of bonds")
        })?;

        // Now read the next 'natoms' lines as the atom block.
        // Each line usually looks like:
        //   X Y Z Element ??? ??? ...
        //   e.g. "    1.4386   -0.8054   -0.4963 O   0  0  0  0  0  0  0  0  0  0  0  0"
        //

        let first_atom_line = 4;
        let last_atom_line = first_atom_line + n_atoms;
        let first_bond_line = last_atom_line;
        let last_bond_line = first_bond_line + n_bonds;

        if lines.len() < last_atom_line {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines for the declared atom block",
            ));
        }

        let mut atoms = Vec::with_capacity(n_atoms);

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
                type_in_res: None,
                posit: Vec3 { x, y, z }, // or however you store coordinates
                element: Element::from_letter(element)?,
                occupancy: None,
                partial_charge: None,
                force_field_type: None,
                hetero: true,
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

            let atom_0_sn = cols[0].parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse bond atom 0")
            })?;
            let atom_1_sn = cols[1].parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse bond atom 1")
            })?;

            let bond_type = BondType::from_str(cols[2])?;

            bonds.push(BondGeneric {
                atom_0_sn,
                atom_1_sn,
                bond_type,
            })
        }

        // Look for a molecule identifier in the file. Check for either
        // "> <PUBCHEM_COMPOUND_CID>" or "> <DATABASE_ID>" and take the next nonempty line.
        let mut pubchem_cid = None;
        let mut drugbank_id = None;

        // todo: Handle more metadata?

        for (i, line) in lines.iter().enumerate() {
            if line.contains("> <PUBCHEM_COMPOUND_CID>")
                && let Some(value_line) = lines.get(i + 1)
            {
                let value = value_line.trim();
                if let Ok(v) = value.parse::<u32>() {
                    pubchem_cid = Some(v);
                }
            }
            if line.contains("> <DATABASE_ID>") {
                if let Some(value_line) = lines.get(i + 1) {
                    let value = value_line.trim();
                    if !value.is_empty() {
                        drugbank_id = Some(value.to_string());
                    }
                }
            }
        }

        let ident = lines[0].trim().to_string();
        // We observe that on at least some DrugBank files, this line
        // is the PubChem ID, even if the PUBCHEM_COMPOUND_CID line is omitted.
        if let Ok(v) = lines[0].parse::<u32>() {
            pubchem_cid = Some(v);
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
        });

        chains.push(ChainGeneric {
            id: "A".to_string(),
            residue_sns: vec![0],
            atom_sns,
        });

        Ok(Self {
            ident,
            atoms,
            chains,
            residues,
            pubchem_cid,
            drugbank_id,
            metadata: HashMap::new(), // todo: A/R
            bonds,
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
                bond.atom_0_sn, bond.atom_1_sn, &bond.bond_type
            )?;
        }

        writeln!(file, "M  END")?;

        // Metadata
        if let Some(cid) = self.pubchem_cid {
            writeln!(file, "> <PUBCHEM_COMPOUND_CID>")?;
            writeln!(file, "{cid}")?;
            writeln!(file)?; // blank line
        }
        if let Some(ref dbid) = self.drugbank_id {
            writeln!(file, "> <DATABASE_ID>")?;
            writeln!(file, "{dbid}")?;
            writeln!(file)?; // blank line
            writeln!(file, "> <DATABASE_NAME>")?;
            writeln!(file, "drugbank")?;
            writeln!(file)?; // blank line
        }

        // If you have a general metadata HashMap, you could do:
        // for (key, value) in &self.metadata {
        //     writeln!(file, "> <{}>", key)?;
        //     writeln!(file, "{}", value)?;
        //     writeln!(file)?;
        // }

        // 8) End of this molecule record in SDF
        writeln!(file, "$$$$")?;

        Ok(())
    }

    // todo: Generic fn for this and save, among all text-based types.
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        let data_str: String = String::from_utf8(buffer)
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid UTF8"))?;

        Self::new(&data_str)
    }
}
