//! For computing atom-centered charges, e.g. MBIS, CHELPG, and RESP.
//! [Docs](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html?q=mbis&n=0#mbis-charges)

// todo: Support CHELPG and RESP.

use std::io;

#[derive(Debug, Clone)]
pub struct AtomChargeData {
    pub charge: f64,
    pub population: f64,
    pub spin: f64,
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

        Ok(Self {
            text,
            convergence_thresh,
            num_iters,
            total_integrated_alpha_density,
            total_integrated_beta_density,
            charges,
        })
    }
}
