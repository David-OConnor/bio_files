//! For reading and writing .GRO files

use std::{
    fs,
    fs::File,
    io,
    io::{BufWriter, ErrorKind, Write},
    path::Path,
};

use lin_alg::f64::Vec3;
use na_seq::Element;

use crate::el_from_atom_name;

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

            atoms.push(AtomGro {
                mol_id,
                mol_name,
                element: el_from_atom_name(&atom_type),
                atom_type,
                serial_number,
                posit: Vec3 { x, y, z },
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
            writeln!(
                w,
                "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                atom.mol_id % 100_000,
                &atom.mol_name[..atom.mol_name.len().min(5)],
                &atom.atom_type[..atom.atom_type.len().min(5)],
                atom.serial_number % 100_000,
                atom.posit.x,
                atom.posit.y,
                atom.posit.z,
            )?;
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
