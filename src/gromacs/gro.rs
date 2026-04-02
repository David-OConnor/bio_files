//! For reading and writing .GRO files

use std::{
    fs,
    fs::File,
    io,
    io::{BufWriter, ErrorKind, Write},
    path::Path,
};

use lin_alg::f64::Vec3;

#[derive(Clone, Debug, PartialEq)]
pub struct AtomGro {
    pub mol_name: String,
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
            let cols: Vec<&str> = line.split_whitespace().collect();

            if cols.len() < 6 {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Atom line has fewer than 6 columns: {line}"),
                ));
            }

            let mol_name = cols[0].to_owned();
            let atom_type = cols[1].to_owned();

            let serial_number = cols[2].parse::<u32>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse atom serial number")
            })?;
            let x = cols[3].parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse X coordinate")
            })?;
            let y = cols[4].parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse Y coordinate")
            })?;
            let z = cols[5].parse::<f64>().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Could not parse Z coordinate")
            })?;

            atoms.push(AtomGro {
                mol_name,
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
            // mol_name(9) atom_type(6) serial(5) x(8.3) y(8.3) z(8.3)
            writeln!(
                w,
                "{:>9}{:>6}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                atom.mol_name,
                atom.atom_type,
                atom.serial_number,
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
