//! For parsing [Amber force field](https://ambermd.org/AmberModels.php) dat files, like
//! `gaff2.dat` (small molecules/ligands). These include parameters used in molecular dynamics.
//! See also: `frcmod`, which patches these parameters for specific molecules.

// todo: For now, this has much overlap with loading frcmod data.
// todo: Reoncosider the API, how this and the frcmod modules are related, and which
// todo: to export the main FF struct from

use std::{
    fs::File,
    io,
    io::{ErrorKind, Read},
    path::Path,
};

use crate::frcmod::{
    AngleData, BondData, DihedralData, ForceFieldParams, ImproperDihedralData, MassData, VdwData,
};

impl ForceFieldParams {
    /// From a string of a dat text file, from Amber.
    pub fn from_dat(text: &str) -> io::Result<Self> {
        let mut result = Self::default();

        let mut in_mod4 = false;

        // These dat text-based files are tabular data, and don't have clear delineations bewteen sections.
        // we parse each line based on its content. Notably, the first column alone is a good indicator
        // based on the number of dashes in it. Three atom names separated by dashes, e.g. `pf-p2-s` is angle data.
        // Two, e.g. `ca-s6` is linear bond data. (e.g. springiness of the covalent bond). FOur names indicates
        // a dihedral (aka torsion) angle.
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() {
                // blank line ends MOD4 block
                in_mod4 = false;
                continue;
            }
            // header that *starts* the block
            if line.starts_with("MOD4") {
                in_mod4 = true;
                continue; // nothing else to parse on this header line
            }

            // 1. Break the line into whitespace‐split tokens.
            let tokens: Vec<&str> = line.split_whitespace().collect();

            // 2. Find the first token that _is_ a number (we assume f32 here).
            let num_idx = tokens
                .iter()
                .position(|t| t.parse::<f32>().is_ok())
                .unwrap_or(tokens.len());

            // 3. Everything before that is our "atom field"
            let atom_field = tokens[..num_idx].join(" ");

            // 4. Everything from num_idx onward are the numeric + comment tokens
            let data_tokens = &tokens[num_idx..];
            let mut cols = data_tokens.iter();

            // 5. Helper to grab a trailing comment after consuming N numeric tokens
            let remainder_as_comment = |consumed: usize| -> Option<String> {
                if data_tokens.len() > consumed {
                    Some(data_tokens[consumed..].join(" "))
                } else {
                    None
                }
            };

            // 6. Now split atom_field on hyphens or whitespace to get the atom names
            let atoms: Vec<&str> = atom_field
                .split(|c: char| c == '-' || c.is_whitespace())
                .filter(|s| !s.is_empty())
                .collect();

            match atoms.len() {
                1 => {
                    if in_mod4 {
                        // ---------- VdW line (R*, ε) ---------------------------
                        let (Ok(r_star), Ok(eps)) = (
                            cols.next().unwrap_or(&"0").parse(),
                            cols.next().unwrap_or(&"0").parse(),
                        ) else {
                            result.remarks.push(line.to_string());
                            continue;
                        };
                        let comment = remainder_as_comment(2);
                        result.van_der_waals.push(VdwData {
                            ff_type: atoms[0].to_string(),
                            r_star,
                            eps,
                        });
                    } else {
                        // ---------- normal MASS line ---------------------------
                        let Ok(mass_val) = cols.next().unwrap_or(&"0").parse() else {
                            result.remarks.push(line.to_string());
                            continue;
                        };
                        let _ = cols.next(); // skip polarizability
                        let comment = remainder_as_comment(2);
                        result.mass.push(MassData {
                            ff_type: atoms[0].to_string(),
                            mass: mass_val,
                            comment,
                        });
                    }
                }

                2 => {
                    // (linear) Bond data.
                    let (Ok(k), Ok(r0)) = (
                        cols.next().unwrap_or(&"0").parse(),
                        cols.next().unwrap_or(&"0").parse(),
                    ) else {
                        result.remarks.push(line.to_string());
                        continue;
                    };
                    let comment = remainder_as_comment(2);
                    result.bond.push(BondData {
                        ff_types: (atoms[0].to_string(), atoms[1].to_string()),
                        k,
                        r_0: r0,
                        comment,
                    });
                }

                3 => {
                    // Valence angle data between 3 atoms.
                    let (Ok(k), Ok(angle)) = (
                        cols.next().unwrap_or(&"0").parse(),
                        cols.next().unwrap_or(&"0").parse::<f32>(),
                    ) else {
                        result.remarks.push(line.to_string());
                        continue;
                    };

                    let comment = remainder_as_comment(2);

                    result.angle.push(AngleData {
                        ff_types: (
                            atoms[0].to_string(),
                            atoms[1].to_string(),
                            atoms[2].to_string(),
                        ),
                        k,
                        angle: angle.to_radians(),
                        comment,
                    });
                }

                4 => {
                    // Either DIHEDRAL or IMPROPER.  We decide by counting how
                    // many numeric tokens follow.
                    let numeric: Vec<f32> = cols
                        .clone() // don’t consume yet
                        .filter_map(|t| t.parse().ok())
                        .collect();

                    match numeric.len() {
                        3 => {
                            // IMPROPER: k  phase  periodicity
                            let barrier_height_vn = numeric[0];
                            let gamma = numeric[1].to_radians();
                            let periodicity = numeric[2] as i8;

                            let comment = remainder_as_comment(3);
                            result.improper.push(ImproperDihedralData {
                                ff_types: (
                                    atoms[0].to_string(),
                                    atoms[1].to_string(),
                                    atoms[2].to_string(),
                                    atoms[3].to_string(),
                                ),
                                barrier_height_vn,
                                gamma,
                                periodicity,
                                comment,
                            });
                        }
                        4.. => {
                            // DIHEDRAL: scaling  Vn  γ  n  [notes…]
                            let scaling_factor = cols.next().unwrap().parse::<u8>().unwrap_or(1);
                            let barrier_height_vn = cols.next().unwrap().parse().unwrap_or(0.0);
                            let gamma = cols
                                .next()
                                .unwrap()
                                .parse::<f32>()
                                .unwrap_or(0.0)
                                .to_radians();
                            let periodicity = cols.next().unwrap().parse().unwrap_or(1.0) as i8;

                            let rest_of_line: String = cols
                                .copied() // turn &&str → &str
                                .collect::<Vec<&str>>()
                                .join(" ");

                            let (notes, penalty_score) = if rest_of_line.is_empty() {
                                (None, 0.0)
                            } else {
                                let lower = rest_of_line.to_lowercase();
                                if let Some(idx) = lower.find("penalty score=") {
                                    let after = &rest_of_line[(idx + 14)..];
                                    let penalty_val = after
                                        .split_whitespace()
                                        .next()
                                        .and_then(|s| s.parse().ok())
                                        .unwrap_or(0.0);
                                    (Some(rest_of_line.clone()), penalty_val)
                                } else {
                                    (Some(rest_of_line.clone()), 0.0)
                                }
                            };

                            result.dihedral.push(DihedralData {
                                ff_types: (
                                    atoms[0].to_string(),
                                    atoms[1].to_string(),
                                    atoms[2].to_string(),
                                    atoms[3].to_string(),
                                ),
                                scaling_factor,
                                barrier_height_vn,
                                gamma,
                                periodicity,
                                notes,
                                penalty_score,
                            });
                        }
                        _ => result.remarks.push(line.to_string()),
                    }
                }

                _ => {
                    // anything else
                    result.remarks.push(line.to_string());
                }
            }
        }

        println!("Loaded mass data: {:?}", result.mass);
        println!("Loaded LJ data: {:?}", result.van_der_waals);

        Ok(result)
    }

    /// Write to file
    pub fn save_dat(&self, path: &Path) -> io::Result<()> {
        let mut f = File::create(path)?;

        Ok(())
    }

    /// todo: Sort out the syntax for loading from different sources.
    pub fn load_dat(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        let data_str: String = String::from_utf8(buffer)
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid UTF8"))?;

        Self::from_dat(&data_str)
    }
}
