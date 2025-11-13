use std::{
    fs,
    fs::File,
    io::{self, ErrorKind, Write},
    path::Path,
};

use lin_alg::f64::Vec3;

use crate::AtomGeneric;

// todo: Support trajectory files, e.g. ones with multiple segments of  XYZ.

#[derive(Clone, Debug)]
pub struct Xyz {
    pub atoms: Vec<AtomGeneric>,
    pub comment: String,
}

impl Xyz {
    pub fn new(text: &str) -> io::Result<Self> {
        let lines: Vec<&str> = text.lines().collect();

        if lines.len() < 3 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Not enough lines to parse an XYZ file. Must be at least 3.",
            ));
        }

        let comment = lines[1].to_string();

        let mut atoms = Vec::new();
        for (i, line) in lines.iter().enumerate() {
            if i < 2 || line.trim().is_empty() {
                continue;
            }

            let mut parts = line.split_whitespace();

            let el_str = parts.next().ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Missing element symbol on atom line {}", i),
                )
            })?;

            let x_str = parts.next().ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Missing x coordinate on atom line {}", i),
                )
            })?;
            let y_str = parts.next().ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Missing y coordinate on atom line {}", i),
                )
            })?;
            let z_str = parts.next().ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Missing z coordinate on atom line {}", i),
                )
            })?;

            let x: f64 = x_str.parse().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Invalid x coordinate on atom line {}", i),
                )
            })?;
            let y: f64 = y_str.parse().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Invalid y coordinate on atom line {}", i),
                )
            })?;
            let z: f64 = z_str.parse().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Invalid z coordinate on atom line {}", i),
                )
            })?;

            // Adjust these two lines to match your actual types if needed:
            let element = crate::Element::from_letter(el_str)?;
            let posit = Vec3::new(x, y, z);

            atoms.push(AtomGeneric {
                element,
                posit,
                ..Default::default()
            });
        }

        Ok(Self { atoms, comment })
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let data_str = fs::read_to_string(path)?;
        Self::new(&data_str)
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut file = File::create(path)?;

        writeln!(file, "{}", self.atoms.len())?;
        writeln!(file, "{}", self.comment)?;

        // Note: I'm not sure if there are standards regarding coordinate precision,
        // or indentation. For example, have seen variants with a 2-space indent, and ones with none.
        // I believe 6 spaces between digits not including - is the move though.
        for atom in &self.atoms {
            writeln!(
                file,
                "{:<2} {:>17.10} {:>17.10} {:>17.10}",
                atom.element.to_letter(),
                atom.posit.x,
                atom.posit.y,
                atom.posit.z
            )?;
        }

        Ok(())
    }
}
