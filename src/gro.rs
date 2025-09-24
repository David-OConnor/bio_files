//! For opening GROMACS .gro files (GRO). These are common molecular descriptions for small molecules and proteins.
//! Positions in GRO are in nm; we convert to Å to match the rest of the codebase.

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

use crate::{
    AtomGeneric, BondGeneric, BondType, ChainGeneric, Mol2, ResidueEnd, ResidueGeneric, ResidueType,
};

#[derive(Clone, Debug)]
pub struct Gro {
    pub ident: String,
    pub metadata: HashMap<String, String>,
    pub atoms: Vec<AtomGeneric>,
    pub bonds: Vec<BondGeneric>,
    pub chains: Vec<ChainGeneric>,
    pub residues: Vec<ResidueGeneric>,
    pub box_nm: Option<(f64, f64, f64)>,
}

impl Gro {
    pub fn new(text: &str) -> io::Result<Self> {
        let lines: Vec<&str> = text.lines().collect();

        if lines.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines to parse a GRO header",
            ));
        }

        let ident = lines[0].trim().to_string();

        // Line 2: number of atoms (optionally with time; we only read the leading integer)
        let natoms = lines[1]
            .split_whitespace()
            .next()
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Missing atom count"))?
            .parse::<usize>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse atom count"))?;

        let first_atom_line = 2;
        let last_atom_line = first_atom_line + natoms;

        if lines.len() < last_atom_line + 1 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines for declared atom block",
            ));
        }

        let mut atoms = Vec::with_capacity(natoms);
        let mut res_map: HashMap<i32, (usize, String)> = HashMap::new(); // resid -> (residue index, resname)
        let mut residues: Vec<ResidueGeneric> = Vec::new();

        fn slice(s: &str, a: usize, b: usize) -> &str {
            let len = s.len();
            let a = a.min(len);
            let b = b.min(len);
            &s[a..b]
        }

        fn infer_element(name: &str) -> io::Result<Element> {
            // GRO atom name often starts with element symbol; e.g. " C ", "CA", "CL", "HG11"
            let trimmed = name.trim();
            if trimmed.is_empty() {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Empty atom name; cannot infer element",
                ));
            }
            let chars: Vec<char> = trimmed.chars().collect();
            // Try two-letter then one-letter
            if chars.len() >= 2 {
                let two = format!(
                    "{}{}",
                    chars[0].to_ascii_uppercase(),
                    chars[1].to_ascii_lowercase()
                );
                if let Ok(el) = Element::from_letter(&two) {
                    return Ok(el);
                }
            }
            let one = chars[0].to_ascii_uppercase().to_string();
            Element::from_letter(&one)
        }

        for i in first_atom_line..last_atom_line {
            let line = lines[i];

            // Fixed-width GRO fields (0-based, end-exclusive):
            //  0-5   resid
            //  5-10  resname
            // 10-15  atom name
            // 15-20  atom serial
            // 20-28  x (nm)
            // 28-36  y (nm)
            // 36-44  z (nm)
            let resid_str = slice(line, 0, 5).trim();
            let resname = slice(line, 5, 10).trim().to_string();
            let atom_name = slice(line, 10, 15).trim().to_string();
            let serial_str = slice(line, 15, 20).trim();
            let x_nm_str = slice(line, 20, 28).trim();
            let y_nm_str = slice(line, 28, 36).trim();
            let z_nm_str = slice(line, 36, 44).trim();

            if resid_str.is_empty()
                || serial_str.is_empty()
                || x_nm_str.is_empty()
                || y_nm_str.is_empty()
                || z_nm_str.is_empty()
            {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Atom line {i} does not have enough columns"),
                ));
            }

            let resid = resid_str.parse::<i32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse residue id")
            })?;
            let serial_number = serial_str.parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse atom serial")
            })?;

            let x_nm = x_nm_str
                .parse::<f64>()
                .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse X (nm)"))?;
            let y_nm = y_nm_str
                .parse::<f64>()
                .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse Y (nm)"))?;
            let z_nm = z_nm_str
                .parse::<f64>()
                .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse Z (nm)"))?;

            let element = infer_element(&atom_name)?;

            // Convert nm -> Å
            let posit = Vec3 {
                x: x_nm * 10.0,
                y: y_nm * 10.0,
                z: z_nm * 10.0,
            };

            // Track residues; group by resid
            let (res_idx, _name) = if let Some((idx, _)) = res_map.get(&resid).cloned() {
                (idx, resname.clone())
            } else {
                let idx = residues.len();
                res_map.insert(resid, (idx, resname.clone()));

                residues.push(ResidueGeneric {
                    serial_number: resid as u32,
                    res_type: ResidueType::Other(resname.clone()),
                    atom_sns: Vec::new(),
                    end: ResidueEnd::Hetero,
                });
                (idx, resname.clone())
            };

            atoms.push(AtomGeneric {
                serial_number,
                type_in_res: None, // todo: A/R.
                posit,
                element,
                occupancy: None,
                partial_charge: None,
                force_field_type: None,
                hetero: true,
                alt_conformation_id: None,
            });

            residues[res_idx].atom_sns.push(serial_number);
        }

        // Box line (nm). If absent, set None.
        let mut box_nm = None;
        if let Some(box_line) = lines.get(last_atom_line) {
            let parts: Vec<&str> = box_line.split_whitespace().collect();
            if !parts.is_empty() {
                // Typical GRO has 3 floats (triclinic can have 9, we only take 3)
                let mut nums = parts
                    .iter()
                    .filter_map(|s| s.parse::<f64>().ok())
                    .collect::<Vec<_>>();
                if nums.len() >= 3 {
                    box_nm = Some((nums[0], nums[1], nums[2]));
                }
            }
        }

        // One chain that includes all residues
        let mut residue_sns: Vec<u32> = residues.iter().map(|r| r.serial_number).collect();
        if residue_sns.is_empty() {
            // If there were no residue lines (shouldn't happen), make a single residue with all atoms
            let atom_sns: Vec<_> = atoms.iter().map(|a| a.serial_number).collect();
            residues.push(ResidueGeneric {
                serial_number: 1,
                res_type: ResidueType::Other("RES".to_string()),
                atom_sns: atom_sns.clone(),
                end: ResidueEnd::Hetero,
            });
            residue_sns = vec![1];
        }

        let atom_sns: Vec<_> = atoms.iter().map(|a| a.serial_number).collect();

        let chains = vec![ChainGeneric {
            id: "A".to_string(),
            residue_sns,
            atom_sns,
        }];

        Ok(Self {
            ident,
            metadata: HashMap::new(),
            atoms,
            bonds: Vec::new(),
            chains,
            residues,
            box_nm,
        })
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut file = File::create(path)?;

        writeln!(file, "{}", self.ident)?;
        writeln!(file, "{:>5}", self.atoms.len())?;

        // Map atom -> (resid, resname, atom_name)
        // If residues exist, use them; otherwise default to resid 1 / RES.
        let mut atom_to_resid: HashMap<u32, (i32, String)> = HashMap::new();
        for r in &self.residues {
            let resid = r.serial_number as i32;
            let resname = match &r.res_type {
                ResidueType::Other(s) => s.clone(),
                _ => "RES".to_string(),
            };
            for &sn in &r.atom_sns {
                atom_to_resid.insert(sn, (resid, resname.clone()));
            }
        }

        for atom in &self.atoms {
            let (resid, resname) = atom_to_resid
                .get(&atom.serial_number)
                .cloned()
                .unwrap_or((1, "RES".to_string()));

            let atom_name = match &atom.type_in_res {
                Some(tir) => tir.to_string(),
                None => String::new(),
            };

            // Convert Å -> nm for GRO
            let x_nm = atom.posit.x / 10.0;
            let y_nm = atom.posit.y / 10.0;
            let z_nm = atom.posit.z / 10.0;

            // GRO fixed-width formatting
            // resid(5) resname(5) atom(5) atom_id(5) x(8.3) y(8.3) z(8.3)
            writeln!(
                file,
                "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                resid,
                truncate_to(&resname, 5),
                right_truncate_to(&atom_name, 5),
                atom.serial_number,
                x_nm,
                y_nm,
                z_nm
            )?;
        }

        if let Some((bx, by, bz)) = self.box_nm {
            writeln!(file, "{:>10.5}{:>10.5}{:>10.5}", bx, by, bz)?;
        } else {
            writeln!(file, "{:>10.5}{:>10.5}{:>10.5}", 0.0, 0.0, 0.0)?;
        }

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

impl From<Mol2> for Gro {
    fn from(m: Mol2) -> Self {
        Self {
            ident: m.ident.clone(),
            metadata: m.metadata.clone(),
            atoms: m.atoms.clone(),
            bonds: m.bonds.clone(),
            chains: Vec::new(),
            residues: Vec::new(),
            box_nm: None,
        }
    }
}

fn truncate_to(s: &str, n: usize) -> String {
    let mut out = s.chars().take(n).collect::<String>();
    if out.len() < n {
        out.push_str(&" ".repeat(n - out.len()));
    }
    out
}

fn right_truncate_to(s: &str, n: usize) -> String {
    // Right-aligned within width n; truncate if longer
    let trimmed = s
        .chars()
        .rev()
        .take(n)
        .collect::<String>()
        .chars()
        .rev()
        .collect::<String>();
    format!("{trimmed:>width$}", width = n)
}
