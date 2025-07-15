//! Load atoms and related data (e.g. secondary structure) from mmCIF files.
//! These are the coordinate files that come from the RCSB PDB.
//!
//! Note: This does not parse bond data. This is usually not available. Bonds can
//! be inferred from computations. (Send in a Github issue or PR if you want bond support, and
//! include an example mmCIF that has them.).

use std::{
    collections::HashMap,
    fs::File,
    io,
    io::{Cursor, ErrorKind, Read},
    path::Path,
    str::FromStr,
};

use lin_alg::f64::Vec3;
use na_seq::{AtomTypeInRes, Element};

use crate::{
    AtomGeneric, BackboneSS, ChainGeneric, ExperimentalMethod, ResidueGeneric, ResidueType,
    mmcif_aux::load_ss_method,
};

pub struct MmCif {
    pub ident: String,
    pub metadata: HashMap<String, String>,
    pub atoms: Vec<AtomGeneric>,
    // pub bonds: Vec<BondGeneric>,
    pub chains: Vec<ChainGeneric>,
    pub residues: Vec<ResidueGeneric>,
    pub secondary_structure: Vec<BackboneSS>,
    pub experimental_method: Option<ExperimentalMethod>,
}

#[derive(Clone, Copy, PartialEq)]
enum Section {
    None,
    Atoms,
    Helix, // todo etc
}

impl MmCif {
    pub fn new(text: &str) -> io::Result<Self> {
        println!("Loading mmCIF data...");
        // todo: For these `new` methods in general that take a &str param: Should we use
        // todo R: Reed + Seek instead, and pass a Cursor or File object? Probably doesn't matter.
        // todo Either way, we should keep it consistent between the files.

        // todo: This is far too slow.

        let mut metadata = HashMap::<String, String>::new();
        let mut atoms = Vec::<AtomGeneric>::new();
        let mut residues = Vec::<ResidueGeneric>::new();
        let mut chains = Vec::<ChainGeneric>::new();
        let mut res_idx = HashMap::<(String, u32), usize>::new();
        let mut chain_idx = HashMap::<String, usize>::new();

        let lines: Vec<&str> = text.lines().collect();
        let mut i = 0;
        let n = lines.len();

        while i < n {
            let mut line = lines[i].trim();
            if line.is_empty() {
                i += 1;
                continue;
            }

            if line == "loop_" {
                i += 1;
                let mut headers = Vec::<&str>::new();
                while i < n {
                    line = lines[i].trim();
                    if line.starts_with('_') {
                        headers.push(line);
                        i += 1;
                    } else {
                        break;
                    }
                }

                // Not an atom loop? → skip its rows.
                if !headers
                    .first()
                    .map_or(false, |h| h.starts_with("_atom_site."))
                {
                    while i < n {
                        line = lines[i].trim();
                        if line == "#" || line == "loop_" || line.starts_with('_') {
                            break;
                        }
                        i += 1;
                    }
                    continue;
                }

                let col = |tag: &str| -> io::Result<usize> {
                    headers.iter().position(|h| *h == tag).ok_or_else(|| {
                        io::Error::new(ErrorKind::InvalidData, format!("mmCIF missing {}", tag))
                    })
                };
                let c_id = col("_atom_site.id")?;
                let c_x = col("_atom_site.Cartn_x")?;
                let c_y = col("_atom_site.Cartn_y")?;
                let c_z = col("_atom_site.Cartn_z")?;
                let c_el = col("_atom_site.type_symbol")?;
                let c_name = col("_atom_site.label_atom_id")?;
                let c_res = col("_atom_site.label_comp_id")?;
                let c_chain = col("_atom_site.label_asym_id")?;
                let c_res_sn = col("_atom_site.label_seq_id")?;
                let c_occ = col("_atom_site.occupancy")?;

                while i < n {
                    line = lines[i].trim();
                    if line.is_empty() || line == "#" || line == "loop_" || line.starts_with('_') {
                        break;
                    }
                    let fields: Vec<&str> = line.split_whitespace().collect();
                    if fields.len() < headers.len() {
                        i += 1;
                        continue;
                    }

                    // Atom lines.
                    let serial_number = fields[c_id].parse::<u32>().unwrap_or(0);
                    let x = fields[c_x].parse::<f64>().unwrap_or(0.0);
                    let y = fields[c_y].parse::<f64>().unwrap_or(0.0);
                    let z = fields[c_z].parse::<f64>().unwrap_or(0.0);

                    let element = Element::from_letter(fields[c_el])?;
                    let atom_name = fields[c_name];
                    let type_in_res = AtomTypeInRes::from_str(atom_name).ok();
                    let occ = match fields[c_occ] {
                        "?" | "." => None,
                        v => v.parse().ok(),
                    };

                    atoms.push(AtomGeneric {
                        serial_number,
                        posit: Vec3::new(x, y, z),
                        element,
                        type_in_res,
                        force_field_type: None,
                        occupancy: occ,
                        partial_charge: None,
                    });

                    // --------- Residue / Chain bookkeeping -----------
                    let res_sn = fields[c_res_sn].parse::<u32>().unwrap_or(0);
                    let chain_id = fields[c_chain];
                    let res_key = (chain_id.to_string(), res_sn);

                    // Residues
                    let r_i = *res_idx.entry(res_key.clone()).or_insert_with(|| {
                        let idx = residues.len();
                        residues.push(ResidueGeneric {
                            serial_number: res_sn,
                            res_type: ResidueType::from_str(fields[c_res]),
                            atom_sns: Vec::new(),
                        });
                        idx
                    });
                    residues[r_i].atom_sns.push(serial_number);

                    // Chains
                    let c_i = *chain_idx.entry(chain_id.to_string()).or_insert_with(|| {
                        let idx = chains.len();
                        chains.push(ChainGeneric {
                            id: chain_id.to_string(),
                            residue_sns: Vec::new(),
                            atom_sns: Vec::new(),
                        });
                        idx
                    });
                    chains[c_i].atom_sns.push(serial_number);
                    if !chains[c_i].residue_sns.contains(&res_sn) {
                        chains[c_i].residue_sns.push(res_sn);
                    }

                    i += 1;
                }
                continue; // outer while will handle terminator line
            }

            if line.starts_with('_') {
                if let Some((tag, val)) = line.split_once(char::is_whitespace) {
                    metadata.insert(tag.to_string(), val.trim_matches('\'').to_string());
                } else {
                    metadata.insert(line.to_string(), String::new());
                }
            }

            i += 1; // advance to next top-level line
        }

        let ident = metadata
            .get("_struct.entry_id")
            .or_else(|| metadata.get("_entry.id"))
            .cloned()
            .unwrap_or_else(|| "UNKNOWN".to_string());

        let mut cursor = Cursor::new(text);

        // todo: Integraet this so it's not taking a second line loop through the whole file.
        // todo: It'll be faster this way.
        let (secondary_structure, experimental_method) = load_ss_method(&mut cursor)?;

        println!("CIF load complete");
        Ok(Self {
            ident,
            metadata,
            atoms,
            chains,
            residues,
            secondary_structure,
            experimental_method,
        })
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        //todo: Fix this so it outputs mol2 instead of sdf.
        let mut file = File::create(path)?;

        // todo: Implement this once loading works.
        //
        // // There is a subtlety here. Add that to your parser as well. There are two values
        // // todo in the files we have; this top ident is not the DB id.
        // writeln!(file, "@<TRIPOS>MOLECULE")?;
        // writeln!(file, "{}", self.ident)?;
        // writeln!(file, "{} {}", self.atoms.len(), self.bonds.len())?;
        // writeln!(file, "{}", self.mol_type.to_str())?;
        // writeln!(file, "{}", self.charge_type)?;

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
