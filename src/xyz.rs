use std::{
    fs,
    fs::File,
    io::{self, ErrorKind, Write},
    path::Path,
};

use lin_alg::f64::Vec3;

use crate::AtomGeneric;

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

/// xyz files can contain multiple sets, e.g. in a molecular dynamics
/// trajectory.
pub fn new_xyz_trajectory(text: &str) -> io::Result<Vec<Xyz>> {
    let lines: Vec<&str> = text.lines().collect();

    let mut items: Vec<Xyz> = Vec::new();
    let mut i: usize = 0;

    while i < lines.len() {
        while i < lines.len() && lines[i].trim().is_empty() {
            i += 1;
        }
        if i >= lines.len() {
            break;
        }

        let n_atoms: usize = lines[i].trim().parse().map_err(|_| {
            io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid atom count on line {}", i + 1),
            )
        })?;
        if n_atoms == 0 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid atom count (0) on line {}", i + 1),
            ));
        }

        let comment_line = lines.get(i + 1).ok_or_else(|| {
            io::Error::new(
                ErrorKind::InvalidData,
                format!("Missing comment line after atom count on line {}", i + 1),
            )
        })?;

        let mut atom_lines: Vec<&str> = Vec::with_capacity(n_atoms);
        let mut j = i + 2;
        while atom_lines.len() < n_atoms && j < lines.len() {
            let l = lines[j];
            if !l.trim().is_empty() {
                atom_lines.push(l);
            }
            j += 1;
        }

        if atom_lines.len() != n_atoms {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!(
                    "Unexpected EOF while reading atoms for frame starting at line {} (expected {}, got {})",
                    i + 1,
                    n_atoms,
                    atom_lines.len()
                ),
            ));
        }

        let frame_text = format!(
            "{n_atoms}\n{comment}\n{atoms}\n",
            comment = comment_line,
            atoms = atom_lines.join("\n")
        );

        items.push(Xyz::new(&frame_text)?);
        i = j;
    }

    Ok(items)
}

pub fn load_xyz_trajectory(path: &Path) -> io::Result<Vec<Xyz>> {
    let data_str = fs::read_to_string(path)?;
    new_xyz_trajectory(&data_str)
}

pub fn save_xyz_trajectory(items: &[Xyz], path: &Path) -> io::Result<()> {
    let mut file = File::create(path)?;

    for (idx, item) in items.iter().enumerate() {
        writeln!(file, "{}", item.atoms.len())?;
        writeln!(file, "{}", item.comment)?;

        for atom in &item.atoms {
            writeln!(
                file,
                "{:<2} {:>17.10} {:>17.10} {:>17.10}",
                atom.element.to_letter(),
                atom.posit.x,
                atom.posit.y,
                atom.posit.z
            )?;
        }

        if idx + 1 != items.len() {
            writeln!(file)?;
        }
    }

    Ok(())
}
