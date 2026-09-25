//! Amber ASCII coordinate files (`.inpcrd`, `.rst7`, `.crd`), as written by tleap, ParmEd, and
//! `sander`/`pmemd` restarts. NetCDF restarts aren't supported.
//!
//! Layout: A title line; the atom count (and optionally the time); coordinates in `6F12.7` (two
//! atoms per line); optionally velocities, in the same layout; and optionally a box line with
//! lengths and angles.

use std::{fs, io, path::Path};

use lin_alg::f64::Vec3;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AmberCoords {
    pub title: String,
    /// Å
    pub posits: Vec<Vec3>,
    /// Å per 1/20.455 ps (Amber's internal time unit), as stored. Present in restart files.
    pub velocities: Option<Vec<Vec3>>,
    /// Box lengths (Å) and angles (degrees), if periodic.
    pub box_dims: Option<[f64; 6]>,
}

impl AmberCoords {
    pub fn load(path: &Path) -> io::Result<Self> {
        Self::new(&fs::read_to_string(path)?)
    }

    pub fn new(text: &str) -> io::Result<Self> {
        let invalid =
            |msg: &str| io::Error::new(io::ErrorKind::InvalidData, format!("inpcrd: {msg}"));

        let mut lines = text.lines();
        let title = lines.next().unwrap_or_default().trim().to_owned();
        let n_atoms: usize = lines
            .next()
            .and_then(|l| l.split_whitespace().next())
            .and_then(|v| v.parse().ok())
            .ok_or_else(|| invalid("Missing atom count"))?;

        // Fixed-width 12-character fields. (Values can run together, so we don't split on
        // whitespace.)
        let mut values = Vec::new();
        for line in lines {
            let chars: Vec<char> = line.chars().collect();
            for chunk in chars.chunks(12) {
                let field: String = chunk.iter().collect();
                let field = field.trim();
                if field.is_empty() {
                    continue;
                }
                values.push(
                    field
                        .parse::<f64>()
                        .map_err(|_| invalid(&format!("Invalid number {field:?}")))?,
                );
            }
        }

        let n = 3 * n_atoms;
        if values.len() < n {
            return Err(invalid(&format!(
                "{} values for {n_atoms} atoms",
                values.len()
            )));
        }

        let to_vecs = |v: &[f64]| -> Vec<Vec3> {
            v.chunks_exact(3)
                .map(|c| Vec3::new(c[0], c[1], c[2]))
                .collect()
        };

        let posits = to_vecs(&values[..n]);
        let rest = &values[n..];

        // What follows is velocities (3N values), a box (6 values), or both.
        let (velocities, box_values) = if rest.len() >= n + 6 || (rest.len() == n && n != 6) {
            (Some(to_vecs(&rest[..n])), &rest[n..])
        } else {
            (None, rest)
        };
        let box_dims = match box_values.len() {
            0 => None,
            6 => Some([
                box_values[0],
                box_values[1],
                box_values[2],
                box_values[3],
                box_values[4],
                box_values[5],
            ]),
            len => return Err(invalid(&format!("{len} unexpected trailing values"))),
        };

        Ok(Self {
            title,
            posits,
            velocities,
            box_dims,
        })
    }
}
