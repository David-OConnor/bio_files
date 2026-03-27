//! Parsing GROMACS simulation output.
//!
//! After an `mdrun`, we use `gmx trjconv` to export the trajectory as a
//! multi-model `.gro` file, then parse each frame here into [`GromacsFrame`].
//! Energy and log data are extracted from the `.log` file when present.

use std::io::{self, ErrorKind};

use lin_alg::f64::Vec3;

/// A single frame of a GROMACS trajectory.
#[derive(Clone, Debug)]
pub struct GromacsFrame {
    /// Simulation time in **ps**.
    pub time: f64,
    /// Atom positions in **Å** (converted from nm on parse).
    pub atom_posits: Vec<Vec3>,
}

/// The collected output of a `gmx mdrun` run.
#[derive(Clone, Debug)]
pub struct GromacsOutput {
    /// Full text of the GROMACS `.log` file.
    pub log_text: String,
    /// Parsed trajectory frames.
    pub trajectory: Vec<GromacsFrame>,
    /// Number of solute (non-water) atoms per frame. Atoms beyond this index
    /// in each frame's `atom_posits` are water molecules.
    pub solute_atom_count: usize,
}

impl GromacsOutput {
    /// Parse a GROMACS run's output given the raw log text, a multi-frame
    /// `.gro` string produced by `gmx trjconv`, and the number of solute atoms.
    ///
    /// todo: Would it be faster to, instead, parse an output binary (TRR? XTC?)
    pub fn new(log_text: String, traj_gro: &str, solute_atom_count: usize) -> io::Result<Self> {
        let trajectory = parse_multi_gro(traj_gro)?;
        Ok(Self {
            log_text,
            trajectory,
            solute_atom_count,
        })
    }
}

/// Parse a multi-model `.gro` file (as produced by `gmx trjconv`) into a
/// sequence of [`GromacsFrame`] values.
///
/// Each model block has the layout:
/// ```text
/// <title — may contain "t= <time_ps>">
/// <natoms>
/// <atom lines…>
/// <box vector line>
/// ```
pub fn parse_multi_gro(text: &str) -> io::Result<Vec<GromacsFrame>> {
    let mut frames = Vec::new();
    let mut lines = text.lines().peekable();

    while lines.peek().is_some() {
        // Title line
        let title = match lines.next() {
            Some(l) => l,
            None => break,
        };

        // Extract time from title if present, e.g. "MD of system, t= 0.00000"
        let time_ps = extract_time(title);

        // Atom count
        let natoms_str = lines.next().ok_or_else(|| {
            io::Error::new(ErrorKind::UnexpectedEof, "Expected atom count after title")
        })?;
        let natoms: usize = natoms_str
            .trim()
            .split_whitespace()
            .next()
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Empty atom count line"))?
            .parse()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse atom count"))?;

        let mut posits = Vec::with_capacity(natoms);

        for _ in 0..natoms {
            let line = lines.next().ok_or_else(|| {
                io::Error::new(ErrorKind::UnexpectedEof, "Unexpected end of atom block")
            })?;

            // GRO fixed-width: x at cols 20-28, y at 28-36, z at 36-44 (nm)
            let x_nm = parse_col(line, 20, 28)?;
            let y_nm = parse_col(line, 28, 36)?;
            let z_nm = parse_col(line, 36, 44)?;

            posits.push(Vec3 {
                x: x_nm * 10.0, // nm → Å
                y: y_nm * 10.0,
                z: z_nm * 10.0,
            });
        }

        // Box vector line — consume it even if we don't use the value.
        lines.next();

        frames.push(GromacsFrame {
            time: time_ps,
            atom_posits: posits,
        });
    }

    Ok(frames)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Extract simulation time (ps) from a GRO title line such as
/// `"Protein in water, t= 12.00000"`. Returns `0.0` if not found.
fn extract_time(title: &str) -> f64 {
    // Look for "t=" or "t =" anywhere in the title.
    for sep in ["t=", "t ="] {
        if let Some(pos) = title.find(sep) {
            let rest = title[pos + sep.len()..].trim();
            // Take the first whitespace-delimited token.
            if let Some(tok) = rest.split_whitespace().next() {
                if let Ok(v) = tok.parse::<f64>() {
                    return v;
                }
            }
        }
    }
    0.0
}

/// Parse a fixed-width column slice as f64. Returns 0.0 on an empty slice
/// rather than an error, to tolerate short lines at end of file.
fn parse_col(line: &str, start: usize, end: usize) -> io::Result<f64> {
    let len = line.len();
    if start >= len {
        return Ok(0.0);
    }
    let end = end.min(len);
    let s = line[start..end].trim();
    if s.is_empty() {
        return Ok(0.0);
    }
    s.parse::<f64>()
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, format!("Cannot parse '{s}' as f64")))
}
