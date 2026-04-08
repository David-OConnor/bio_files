//! For reading and writing .GRO files

use std::{
    fmt::Write as _,
    fs,
    fs::File,
    io,
    io::{BufWriter, ErrorKind, Write},
    path::Path,
};

use lin_alg::f64::Vec3;
use na_seq::Element;

use crate::{el_from_atom_name, gromacs::MoleculeInput};

#[derive(Clone, Debug, PartialEq)]
pub struct AtomGro {
    /// Residue sequence number (resid), columns 0–4. Uniquely identifies each residue/molecule instance.
    pub mol_id: u32,
    /// Residue name (resname), columns 5–9.
    pub mol_name: String,
    pub element: Element,
    pub atom_type: String,
    pub serial_number: u32,
    pub posit: Vec3,
    pub velocity: Option<Vec3>,
}

#[derive(Clone, Debug, PartialEq)]
/// Contains data on atoms, e.g. as the initial state of a MD sim. Similar to Mol2 and SDF in some ways.
/// Can also be used as MD snapshots.
pub struct Gro {
    pub atoms: Vec<AtomGro>,
    pub head_text: String,
    /// The simulation box vectors (diagonal components, in nm).
    pub box_vec: Vec3,
}

impl Gro {
    /// Parse from a GRO-format string.
    pub fn new(text: &str) -> io::Result<Self> {
        let lines: Vec<&str> = text.lines().collect();

        if lines.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines to parse a GRO file",
            ));
        }

        let head_text = lines[0].to_owned();

        let num_atoms = lines[1]
            .trim()
            .parse::<usize>()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse atom count"))?;

        if lines.len() < 2 + num_atoms + 1 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "File has fewer lines than declared atom count",
            ));
        }

        let mut atoms = Vec::with_capacity(num_atoms);

        for line in &lines[2..2 + num_atoms] {
            // GRO uses fixed-width columns. When numbers exceed their field width the
            // spaces between fields disappear, so split_whitespace() fails. Slice directly:
            //   resid    [0 .. 5]   (5 chars, right-justified integer)
            //   resname  [5 ..10]   (5 chars, left-justified string)  → mol_name
            //   atomname [10..15]   (5 chars)                         → atom_type
            //   atomser  [15..20]   (5 chars, right-justified integer) → serial_number
            //   x        [20..28]   (8.3f nm)
            //   y        [28..36]
            //   z        [36..44]
            let gro_slice = |a: usize, b: usize| -> &str {
                let len = line.len();
                line[a.min(len)..b.min(len)].trim()
            };

            if line.len() < 44 {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Atom line is shorter than 44 chars: {line}"),
                ));
            }

            let mol_id = gro_slice(0, 5).parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse residue ID")
            })?;
            let mol_name = gro_slice(5, 10).to_owned(); // resname
            let atom_type = gro_slice(10, 15).to_owned(); // atomname

            let serial_number = gro_slice(15, 20).parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse atom serial number")
            })?;
            let x = gro_slice(20, 28).parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse X coordinate")
            })?;
            let y = gro_slice(28, 36).parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse Y coordinate")
            })?;
            let z = gro_slice(36, 44).parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse Z coordinate")
            })?;

            // Optional velocity columns: vx[44..52] vy[52..60] vz[60..68] (8.4f, nm/ps)
            let velocity = if line.len() >= 68 {
                let vx = gro_slice(44, 52).parse::<f64>().ok();
                let vy = gro_slice(52, 60).parse::<f64>().ok();
                let vz = gro_slice(60, 68).parse::<f64>().ok();
                match (vx, vy, vz) {
                    (Some(vx), Some(vy), Some(vz)) => Some(Vec3 {
                        x: vx,
                        y: vy,
                        z: vz,
                    }),
                    _ => None,
                }
            } else {
                None
            };

            atoms.push(AtomGro {
                mol_id,
                mol_name,
                element: el_from_atom_name(&atom_type),
                atom_type,
                serial_number,
                posit: Vec3 { x, y, z },
                velocity,
            });
        }

        // Last line: box vectors (nm). Take the first three floats.
        let box_line = lines[2 + num_atoms];
        let box_cols: Vec<f64> = box_line
            .split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();

        let box_vec = if box_cols.len() >= 3 {
            Vec3 {
                x: box_cols[0],
                y: box_cols[1],
                z: box_cols[2],
            }
        } else {
            Vec3 {
                x: 0.,
                y: 0.,
                z: 0.,
            }
        };

        Ok(Self {
            atoms,
            head_text,
            box_vec,
        })
    }

    pub fn write_to(&self, w: &mut impl Write) -> io::Result<()> {
        writeln!(w, "{}", self.head_text)?;
        writeln!(w, "{:>5}", self.atoms.len())?;

        for atom in &self.atoms {
            // Standard GRO fixed-width: resid(5) resname(5) atomname(5) atomserial(5) x(8.3) y(8.3) z(8.3)
            // Optional velocities: vx(8.4) vy(8.4) vz(8.4)
            match atom.velocity {
                Some(v) => writeln!(
                    w,
                    "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}",
                    atom.mol_id % 100_000,
                    &atom.mol_name[..atom.mol_name.len().min(5)],
                    &atom.atom_type[..atom.atom_type.len().min(5)],
                    atom.serial_number % 100_000,
                    atom.posit.x,
                    atom.posit.y,
                    atom.posit.z,
                    v.x,
                    v.y,
                    v.z,
                )?,
                None => writeln!(
                    w,
                    "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                    atom.mol_id % 100_000,
                    &atom.mol_name[..atom.mol_name.len().min(5)],
                    &atom.atom_type[..atom.atom_type.len().min(5)],
                    atom.serial_number % 100_000,
                    atom.posit.x,
                    atom.posit.y,
                    atom.posit.z,
                )?,
            }
        }

        writeln!(
            w,
            "{:>10.5}{:>10.5}{:>10.5}",
            self.box_vec.x, self.box_vec.y, self.box_vec.z,
        )?;

        Ok(())
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);
        self.write_to(&mut w)
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let data_str = fs::read_to_string(path)?;
        Self::new(&data_str)
    }
}

/// Generate a GROMACS `.gro` structure file from all molecule atoms.
/// [Example format](https://manual.gromacs.org/2026.1/reference-manual/file-formats.html#gro)
///
/// Positions are converted from Å (internal) to nm (GROMACS). This has similar data to that stored
/// in various molecule formats (SDF, Mol2, mmCIF etc), but also includes velocities,
pub fn make_gro(mols: &[MoleculeInput], box_nm: &Option<(f64, f64, f64)>) -> String {
    let total_atoms: usize = mols.iter().map(|m| m.atoms.len() * m.count).sum();

    // Shift all atoms so their centroid lands at the box center.
    // This ensures `gmx solvate` places water around the solute rather than in a corner.
    let (shift_x, shift_y, shift_z) = if let Some((bx, by, bz)) = box_nm {
        let (mut sx, mut sy, mut sz) = (0.0, 0.0, 0.0);
        let mut n = 0usize;
        for mol in mols {
            for _ in 0..mol.count {
                for atom in &mol.atoms {
                    sx += atom.posit.x;
                    sy += atom.posit.y;
                    sz += atom.posit.z;
                    n += 1;
                }
            }
        }
        if n > 0 {
            // centroid in nm; box_nm already in nm
            let cx = sx / n as f64 / 10.0;
            let cy = sy / n as f64 / 10.0;
            let cz = sz / n as f64 / 10.0;
            (bx / 2.0 - cx, by / 2.0 - cy, bz / 2.0 - cz)
        } else {
            (0.0, 0.0, 0.0)
        }
    } else {
        (0.0, 0.0, 0.0)
    };

    let mut s = String::from("Generated by Bio Files\n");
    let _ = writeln!(s, "{}", total_atoms);

    let mut atom_serial = 1;
    let mut res_serial = 1;

    for mol in mols {
        for _copy in 0..mol.count {
            for atom in &mol.atoms {
                // atom name priority:
                // - hetero / small molecule: FF type (GAFF, e.g. "oh", "c3") first,
                //   then mol2 atom name, then element+index
                // - protein / standard residue: residue atom name ("CA", "CB") first,
                //   then FF type, then element+index
                let atom_name = if atom.hetero {
                    atom.force_field_type
                        .clone()
                        .or_else(|| atom.type_in_res.as_ref().map(|t| t.to_string()))
                        .or_else(|| atom.type_in_res_general.clone())
                        .unwrap_or_else(|| format!("{}{}", atom.element.to_letter(), atom_serial))
                } else {
                    atom.type_in_res
                        .as_ref()
                        .map(|t| t.to_string())
                        .or_else(|| atom.force_field_type.clone())
                        .or_else(|| atom.type_in_res_general.clone())
                        .unwrap_or_else(|| format!("{}{}", atom.element.to_letter(), atom_serial))
                };

                let x_nm = atom.posit.x / 10.0 + shift_x;
                let y_nm = atom.posit.y / 10.0 + shift_y;
                let z_nm = atom.posit.z / 10.0 + shift_z;

                // GRO fixed-width format:
                // resid(5) resname(5) atom(5) serial(5) x(8.3) y(8.3) z(8.3)
                let _ = writeln!(
                    s,
                    "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                    res_serial,
                    &mol.name[..mol.name.len().min(5)],
                    &atom_name[..atom_name.len().min(5)],
                    atom_serial % 100_000,
                    x_nm,
                    y_nm,
                    z_nm,
                );
                atom_serial += 1;
            }
            res_serial += 1;
        }
    }

    // Box vector line (nm)
    let (bx, by, bz) = box_nm.unwrap_or((0.0, 0.0, 0.0));
    let _ = writeln!(s, "{:>10.5}{:>10.5}{:>10.5}", bx, by, bz);

    s
}
