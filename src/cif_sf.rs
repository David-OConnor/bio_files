//! Loads 2fo-fc mmCIF data. This contains electron density. Does not perform a Fourier transform
//! to convert this data into a density grid; returns that file contents in a struct that
//! can then be converted to densities.

use std::{fmt::Display, fs, io, path::Path};

use crate::{DensityHeaderInner, UnitCell};

/// Data on lattice planes. `h`, `k`, and `l` are Miller Indices.
#[derive(Clone, Debug)]
pub struct MillerIndices {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    pub amp: Option<f32>, // amplitude of 2Fo-Fc
    pub phase: Option<f32>,
    // todo: Complex struct for this? Or can they be present one or the other?
    pub re: Option<f32>,
    pub im: Option<f32>,
}

impl Display for MillerIndices {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {}", self.h, self.k, self.l)?;

        if let Some(v) = self.amp {
            write!(f, "  Amp: {v:.2}")?;
        }
        if let Some(v) = self.phase {
            write!(f, "  Phase: {v:.2}")?;
        }
        if let Some(v) = self.re {
            write!(f, "  Re: {v:.2}")?;
        }
        if let Some(v) = self.im {
            write!(f, "  Im: {v:.2}")?;
        }

        Ok(())
    }
}

/// Note: This has much overlap with the `MapHeader`.
#[derive(Clone, Debug)]
pub struct CifStructureFactors {
    pub header: DensityHeaderInner,
    pub miller_indices: Vec<MillerIndices>,
}

impl Display for CifStructureFactors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Header: {:.3?}", self.header)?;
        writeln!(f, "\nMiller indices:")?;

        for index in &self.miller_indices {
            writeln!(f, "{index}")?;
        }

        Ok(())
    }
}

impl CifStructureFactors {
    pub fn new(cif_data: &str) -> io::Result<Self> {
        let mut cell_a = None;
        let mut cell_b = None;
        let mut cell_c = None;
        let mut cell_alpha = None;
        let mut cell_beta = None;
        let mut cell_gamma = None;
        let mut ispg = 0;

        let mut loops = Vec::new();
        parse_loops(cif_data, &mut loops);

        for (tags, _rows) in &loops {
            // cell constants (non-loop case handled below too)
            if tags.iter().any(|t| t.starts_with("_cell.")) {
                // ignore; single-value keys are parsed later
            }
            // reflections
            if tags
                .iter()
                .any(|t| t.eq_ignore_ascii_case("_refln.index_h"))
                || tags.iter().any(|t| t.eq_ignore_ascii_case("_refln.h"))
            {
                // handled below after gathering all loops
            }
        }

        // single-value tags
        let kv = parse_keyvals(cif_data);
        if let Some(v) = kv.get_ci("_cell.length_a") {
            cell_a = some_f(v);
        }
        if let Some(v) = kv.get_ci("_cell.length_b") {
            cell_b = some_f(v);
        }
        if let Some(v) = kv.get_ci("_cell.length_c") {
            cell_c = some_f(v);
        }
        if let Some(v) = kv.get_ci("_cell.angle_alpha") {
            cell_alpha = some_f(v);
        }
        if let Some(v) = kv.get_ci("_cell.angle_beta") {
            cell_beta = some_f(v);
        }
        if let Some(v) = kv.get_ci("_cell.angle_gamma") {
            cell_gamma = some_f(v);
        }

        if ispg == 0 {
            if let Some(v) = kv.get_ci("_space_group.it_number") {
                ispg = some_i(v).unwrap_or(0);
            } else if let Some(v) = kv.get_ci("_symmetry.Int_Tables_number") {
                ispg = some_i(v).unwrap_or(0);
            }
        }

        let (a, b, c) = (
            cell_a.ok_or_else(err("missing _cell.length_a"))?,
            cell_b.ok_or_else(err("missing _cell.length_b"))?,
            cell_c.ok_or_else(err("missing _cell.length_c"))?,
        );
        let (alpha, beta, gamma) = (
            cell_alpha.ok_or_else(err("missing _cell.angle_alpha"))?,
            cell_beta.ok_or_else(err("missing _cell.angle_beta"))?,
            cell_gamma.ok_or_else(err("missing _cell.angle_gamma"))?,
        );
        let cell = UnitCell::new(a, b, c, alpha, beta, gamma);

        // pick the reflections loop
        let mut milller = Vec::<MillerIndices>::new();

        for (tags, rows) in &loops {
            let idx_h = find_tag(tags, &["_refln.index_h", "_refln.h"]);
            let idx_k = find_tag(tags, &["_refln.index_k", "_refln.k"]);
            let idx_l = find_tag(tags, &["_refln.index_l", "_refln.l"]);
            if idx_h.is_none() || idx_k.is_none() || idx_l.is_none() {
                continue;
            }
            let ih = idx_h.unwrap();
            let ik = idx_k.unwrap();
            let il = idx_l.unwrap();

            // 2Fo-Fc amplitude/phase (FWT/PHWT or synonyms)
            let ia = find_tag_relaxed(
                tags,
                &["FWT", "F_2FOFC", "F_2FOFCWT", "F_2FO_FC", "FWT_2FOFC"],
            );
            let ip = find_tag_relaxed(
                tags,
                &["PHWT", "PH_2FOFC", "PH_2FOFCWT", "PH_2FO_FC", "PHWT_2FOFC"],
            );

            // or complex re/im
            let ire = find_tag_relaxed(tags, &["C_2FOFC_RE", "FWT_RE", "MAPC_RE", "RE_2FOFC"]);
            let iim = find_tag_relaxed(tags, &["C_2FOFC_IM", "FWT_IM", "MAPC_IM", "IM_2FOFC"]);

            for row in rows {
                let h = some_i(&row[ih]).unwrap_or(0);
                let k = some_i(&row[ik]).unwrap_or(0);
                let l = some_i(&row[il]).unwrap_or(0);

                if let (Some(ia), Some(ip)) = (ia, ip) {
                    let amp = some_f(&row[ia]);
                    let ph = some_f(&row[ip]);

                    if let (Some(amp), Some(ph)) = (amp, ph) {
                        milller.push(MillerIndices {
                            h,
                            k,
                            l,
                            amp: Some(amp as f32),
                            phase: Some((ph as f32).to_radians()),
                            re: None,
                            im: None,
                        });
                        continue;
                    }
                }

                if let (Some(ire), Some(iim)) = (ire, iim) {
                    let re = some_f(&row[ire]);
                    let im = some_f(&row[iim]);
                    if let (Some(re), Some(im)) = (re, im) {
                        milller.push(MillerIndices {
                            h,
                            k,
                            l,
                            amp: None,
                            phase: None,
                            re: Some(re as f32),
                            im: Some(im as f32),
                        });
                        continue;
                    }
                }
            }
            if !milller.is_empty() {
                break;
            }
        }

        if milller.is_empty() {
            return Err(io_err("no reflections with 2Fo-Fc coefficients found"));
        }

        // grid from Miller span
        let (mut max_h, mut max_k, mut max_l) = (0i32, 0i32, 0i32);
        for r in &milller {
            max_h = max_h.max(r.h.abs());
            max_k = max_k.max(r.k.abs());
            max_l = max_l.max(r.l.abs());
        }
        let mx = next_good_fft_len((2 * max_h + 1).unsigned_abs() as usize).max(64);
        let my = next_good_fft_len((2 * max_k + 1).unsigned_abs() as usize).max(64);
        let mz = next_good_fft_len((2 * max_l + 1).unsigned_abs() as usize).max(64);

        let header = DensityHeaderInner {
            cell,
            // Is this always true? Lots of hard-coded values here. I don't
            // see this values from observing these files, so it's probably OK.
            mapc: 1,
            mapr: 2,
            maps: 3,
            nxstart: 0,
            nystart: 0,
            nzstart: 0,
            mx: mx as i32,
            my: my as i32,
            mz: mz as i32,
            ispg,
            nsymbt: 0,
            version: 20140,
            xorigin: None,
            yorigin: None,
            zorigin: None,
        };

        Ok(Self {
            header,
            miller_indices: milller,
        })
    }

    pub fn new_from_path(path: &Path) -> io::Result<Self> {
        let cif_str = fs::read_to_string(path)?;
        Self::new(&cif_str)
    }
}

// --- helpers ---

fn err(msg: &'static str) -> impl FnOnce() -> io::Error {
    move || io_err(msg)
}
fn io_err(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

trait GetCI {
    fn get_ci(&self, k: &str) -> Option<&String>;
}
impl GetCI for std::collections::HashMap<String, String> {
    fn get_ci(&self, k: &str) -> Option<&String> {
        self.get(&k.to_ascii_lowercase())
    }
}

fn some_f(s: &str) -> Option<f64> {
    let t = s.trim_matches(|c| c == '\'' || c == '"');
    if t.eq("?") || t.eq(".") {
        return None;
    }
    t.parse::<f64>().ok()
}
fn some_i(s: &str) -> Option<i32> {
    let t = s.trim_matches(|c| c == '\'' || c == '"');
    if t.eq("?") || t.eq(".") {
        return None;
    }
    t.parse::<i32>().ok()
}

fn next_good_fft_len(n: usize) -> usize {
    // power of two for simplicity
    n.next_power_of_two()
}

fn normalize_tag(t: &str) -> String {
    let mut s = String::with_capacity(t.len());
    for ch in t.chars() {
        if ch.is_alphanumeric() {
            s.push(ch.to_ascii_uppercase());
        } else if ch == '_' { /* drop */
        }
    }
    s
}

fn find_tag(tags: &[String], keys: &[&str]) -> Option<usize> {
    for (i, t) in tags.iter().enumerate() {
        for k in keys {
            if t.eq_ignore_ascii_case(k) {
                return Some(i);
            }
        }
    }
    None
}

fn find_tag_relaxed(tags: &[String], aliases: &[&str]) -> Option<usize> {
    let norm_tags: Vec<String> = tags.iter().map(|t| normalize_tag(t)).collect();
    for (i, nt) in norm_tags.iter().enumerate() {
        for a in aliases {
            let na = normalize_tag(a);
            if nt.contains(&na) || nt.ends_with(&na) {
                return Some(i);
            }
        }
    }
    None
}

fn parse_keyvals(input: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for line in input.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("loop_") {
            continue;
        }
        if line.starts_with('_') {
            let mut it = line.splitn(2, char::is_whitespace);
            let key = it.next().unwrap().to_ascii_lowercase();
            let val = it.next().unwrap_or("").trim().to_string();
            if !val.is_empty() {
                map.insert(key, val);
            }
        }
    }
    map
}

fn parse_loops(input: &str, out: &mut Vec<(Vec<String>, Vec<Vec<String>>)>) {
    let mut it = input.lines().peekable();
    while let Some(line) = it.next() {
        let l = line.trim();

        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        if l == "loop_" {
            let mut tags = Vec::new();
            while let Some(nl) = it.peek() {
                let s = nl.trim();
                if s.starts_with('_') {
                    tags.push(s.to_string());
                    it.next();
                } else {
                    break;
                }
            }
            let mut rows = Vec::new();
            let mut row = Vec::new();
            while let Some(nl) = it.peek() {
                let s = nl.trim();
                if s.is_empty() || s.starts_with('#') {
                    it.next();
                    continue;
                }
                if s == "loop_" || s.starts_with('_') {
                    break;
                }
                let toks = cif_tokenize(s);
                if toks.is_empty() {
                    it.next();
                    continue;
                }
                row.extend(toks);
                if row.len() >= tags.len() {
                    if row.len() == tags.len() {
                        rows.push(row.clone());
                        row.clear();
                    } else {
                        // rare multiline; flush by chunks of tags.len()
                        while row.len() >= tags.len() {
                            rows.push(row.drain(..tags.len()).collect());
                        }
                    }
                }
                it.next();
            }
            out.push((tags, rows));
        }
    }
}

fn cif_tokenize(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0;
    let b = line.as_bytes();
    while i < b.len() {
        while i < b.len() && b[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= b.len() {
            break;
        }
        let c = b[i] as char;
        if c == '\'' || c == '"' {
            let q = c;
            i += 1;
            let start = i;
            while i < b.len() && (b[i] as char) != q {
                i += 1;
            }
            out.push(line[start..i].to_string());
            i += 1;
        } else {
            let start = i;
            while i < b.len() && !b[i].is_ascii_whitespace() {
                i += 1;
            }
            out.push(line[start..i].to_string());
        }
    }
    out
}
