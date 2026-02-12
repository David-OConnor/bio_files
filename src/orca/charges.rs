//! For computing atom-centered charges, e.g. MBIS, CHELPG, and RESP.
//! [Docs](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html?q=mbis&n=0#mbis-charges)

// todo: Support CHELPG and RESP.

use std::io;

use lin_alg::f64::Vec3;

use crate::orca::make_inp_block;

/// Note: these fields currently do not do anything; all are returned.
#[derive(Debug, Clone, Default)]
pub struct MbisChargesCfg {
    pub dipole: bool,
    pub quadrupole: bool,
    pub octopole: bool,
}

impl MbisChargesCfg {
    pub fn make_inp(&self) -> String {
        let contents = vec![("MBIS_LARGEPRINT", "true".to_string())];
        make_inp_block("method", &contents, &[])
    }
}

#[derive(Debug, Clone)]
pub struct AtomChargeData {
    pub charge: f64,
    pub population: f64,
    pub spin: f64,
}

#[derive(Debug, Clone)]
pub struct Quadrupole {
    pub xx: f64,
    pub yy: f64,
    pub zz: f64,
    pub xy: f64,
    pub xz: f64,
    pub yz: f64,
}

#[derive(Debug, Clone)]
pub struct Octopole {
    pub xxx: f64,
    pub yyy: f64,
    pub zzz: f64,
    pub xxy: f64,
    pub xxz: f64,
    pub xyy: f64,
    pub xyz: f64,
    pub xzz: f64,
    pub yyz: f64,
    pub yzz: f64,
}

/// [Charge tutorial](https://www.faccts.de/docs/orca/5.0/tutorials/prop/charges.html)
/// [MBIS Charges](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html?q=mbis&n=0#mbis-charges)
#[derive(Debug, Clone)]
pub struct ChargesOutput {
    pub text: String,
    pub convergence_thresh: f64,
    pub num_iters: u32,
    pub total_integrated_alpha_density: f64,
    pub total_integrated_beta_density: f64,
    pub charges: Vec<AtomChargeData>,
    pub dipole: Vec<Vec3>,
    pub quadrupole: Vec<Quadrupole>,
    pub octopole: Vec<Octopole>,
}

impl ChargesOutput {
    pub fn new(text: String) -> io::Result<Self> {
        let start = text.rfind("MBIS ANALYSIS").ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "MBIS ANALYSIS section not found",
            )
        })?;

        let mut convergence_thresh: Option<f64> = None;
        let mut num_iters: Option<u32> = None;
        let mut total_alpha: Option<f64> = None;
        let mut total_beta: Option<f64> = None;

        let mut charges: Vec<AtomChargeData> = Vec::new();
        let mut dipole: Vec<Vec3> = Vec::new();
        let mut quadrupole: Vec<Quadrupole> = Vec::new();
        let mut octopole: Vec<Octopole> = Vec::new();

        let mut lines = text[start..].lines();

        let parse_f64 = |s: &str| -> io::Result<f64> {
            s.parse::<f64>()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        };

        let parse_u32 = |s: &str| -> io::Result<u32> {
            s.parse::<u32>()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        };

        // Scan until we hit the header line for the atomic charges table.
        for line in lines.by_ref() {
            let t = line.trim();

            if t.starts_with("Convergence threshold (charges)") {
                if let Some(last) = t.split_whitespace().last() {
                    convergence_thresh = Some(parse_f64(last)?);
                }
            } else if t.starts_with("Number of iterations") {
                if let Some(last) = t.split_whitespace().last() {
                    num_iters = Some(parse_u32(last)?);
                }
            } else if t.starts_with("Total integrated alpha density") {
                if let Some(last) = t.split_whitespace().last() {
                    total_alpha = Some(parse_f64(last)?);
                }
            } else if t.starts_with("Total integrated beta density") {
                if let Some(last) = t.split_whitespace().last() {
                    total_beta = Some(parse_f64(last)?);
                }
            } else if t.starts_with("ATOM") && t.contains("CHARGE") && t.contains("POPULATION") {
                break;
            }
        }

        // Now parse the atomic charge rows until we hit TOTAL or the valence-shell section.
        for line in lines {
            let t = line.trim();
            if t.is_empty() {
                continue;
            }

            let t_no_leading = t.trim_start();

            if t_no_leading.starts_with("TOTAL")
                || t_no_leading.starts_with("MBIS VALENCE-SHELL DATA")
            {
                break;
            }

            // Expect lines like: "0 C    0.208633    5.791367    0.000000"
            let parts: Vec<_> = t.split_whitespace().collect();
            if parts.len() < 5 {
                continue;
            }

            let charge = parse_f64(parts[2])?;
            let population = parse_f64(parts[3])?;
            let spin = parse_f64(parts[4])?;

            charges.push(AtomChargeData {
                charge,
                population,
                spin,
            });
        }

        let convergence_thresh = convergence_thresh.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Convergence threshold not found",
            )
        })?;

        let num_iters = num_iters.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Number of iterations not found")
        })?;

        let total_integrated_alpha_density = total_alpha.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Total integrated alpha density not found",
            )
        })?;

        let total_integrated_beta_density = total_beta.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Total integrated beta density not found",
            )
        })?;

        if charges.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "No atomic MBIS charges found",
            ));
        }

        let mbis_text = &text[start..];

        // Dipoles
        if let Some(after) = section_after(mbis_text, "MBIS ATOMIC DIPOLE MOMENT (A.U.):")
            && let Some(mut rows) = parse_table_rows(after, |t| {
                t.starts_with("ATOM") && t.contains("X") && t.contains("Y") && t.contains("Z")
            })
        {
            for line in rows.by_ref() {
                let t = line.trim();
                if t.is_empty() {
                    break;
                }
                // "0 O 0.000000 -0.000000 -0.165797"
                let parts: Vec<_> = t.split_whitespace().collect();
                if parts.len() < 5 {
                    continue;
                }
                let x = parse_f64(parts[2])?;
                let y = parse_f64(parts[3])?;
                let z = parse_f64(parts[4])?;

                // If Vec3::new doesn't exist in your lin_alg, replace with your constructor.
                dipole.push(Vec3::new(x, y, z));
            }
        }

        // Quadrupoles
        if let Some(after) = section_after(mbis_text, "MBIS ATOMIC QUADRUPOLE MOMENT (A.U.):")
            && let Some(mut rows) = parse_table_rows(after, |t| {
                t.starts_with("ATOM")
                    && t.contains("XX")
                    && t.contains("YY")
                    && t.contains("ZZ")
                    && t.contains("XY")
                    && t.contains("XZ")
                    && t.contains("YZ")
            })
        {
            for line in rows.by_ref() {
                let t = line.trim();
                if t.is_empty() {
                    break;
                }
                // "0 O -4.750163 -5.099543 -4.972574 -0.000000 -0.000000 -0.000000"
                let parts: Vec<_> = t.split_whitespace().collect();
                if parts.len() < 8 {
                    continue;
                }
                let xx = parse_f64(parts[2])?;
                let yy = parse_f64(parts[3])?;
                let zz = parse_f64(parts[4])?;
                let xy = parse_f64(parts[5])?;
                let xz = parse_f64(parts[6])?;
                let yz = parse_f64(parts[7])?;
                quadrupole.push(Quadrupole {
                    xx,
                    yy,
                    zz,
                    xy,
                    xz,
                    yz,
                });
            }
        }

        // Octopoles
        if let Some(after) = section_after(mbis_text, "MBIS ATOMIC OCTUPOLE MOMENT (A.U.):")
            && let Some(mut rows) = parse_table_rows(after, |t| {
                t.starts_with("ATOM")
                    && t.contains("XXX")
                    && t.contains("YYY")
                    && t.contains("ZZZ")
                    && t.contains("XXY")
                    && t.contains("XXZ")
                    && t.contains("XYY")
                    && t.contains("XYZ")
                    && t.contains("XZZ")
                    && t.contains("YYZ")
                    && t.contains("YZZ")
            })
        {
            for line in rows.by_ref() {
                let t = line.trim();
                if t.is_empty() {
                    break;
                }
                // "0 O xxx yyy zzz xxy xxz xyy xyz xzz yyz yzz"
                let parts: Vec<_> = t.split_whitespace().collect();
                if parts.len() < 12 {
                    continue;
                }
                let xxx = parse_f64(parts[2])?;
                let yyy = parse_f64(parts[3])?;
                let zzz = parse_f64(parts[4])?;
                let xxy = parse_f64(parts[5])?;
                let xxz = parse_f64(parts[6])?;
                let xyy = parse_f64(parts[7])?;
                let xyz = parse_f64(parts[8])?;
                let xzz = parse_f64(parts[9])?;
                let yyz = parse_f64(parts[10])?;
                let yzz = parse_f64(parts[11])?;
                octopole.push(Octopole {
                    xxx,
                    yyy,
                    zzz,
                    xxy,
                    xxz,
                    xyy,
                    xyz,
                    xzz,
                    yyz,
                    yzz,
                });
            }
        }

        // Optional: if present, enforce counts match atoms (comment out if you prefer “best-effort”)
        let n = charges.len();
        if !dipole.is_empty() && dipole.len() != n {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "MBIS dipole count does not match atom count",
            ));
        }
        if !quadrupole.is_empty() && quadrupole.len() != n {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "MBIS quadrupole count does not match atom count",
            ));
        }
        if !octopole.is_empty() && octopole.len() != n {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "MBIS octopole count does not match atom count",
            ));
        }

        Ok(Self {
            text,
            convergence_thresh,
            num_iters,
            total_integrated_alpha_density,
            total_integrated_beta_density,
            charges,
            dipole,
            quadrupole,
            octopole,
        })
    }
}

fn section_after<'a>(haystack: &'a str, needle: &str) -> Option<&'a str> {
    let i = haystack.find(needle)?;
    Some(&haystack[i + needle.len()..])
}

fn parse_table_rows(
    s: &str,
    header_predicate: impl Fn(&str) -> bool,
) -> Option<impl Iterator<Item = &str>> {
    let mut it = s.lines();
    while let Some(line) = it.next() {
        let t = line.trim();
        if header_predicate(t) {
            return Some(it);
        }
    }
    None
}
