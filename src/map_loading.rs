use std::{collections::HashMap, f32::consts::PI, io};

use num_complex::Complex32;
use rustfft::{FftPlanner, num_traits::Zero};

use crate::{DensityMap, MapHeader, UnitCell, map::get_origin_frac};
// ----------------- mmCIF parsing helpers -----------------

#[derive(Clone, Debug)]
struct SymOp {
    // rotation acts on fractional x as x' = R x + t (R has integer entries)
    r: [[i32; 3]; 3],
    t: [f64; 3], // fractional
}

fn parse_frac_token(tok: &str) -> f64 {
    // handles "x", "-y", "z+1/3", "1/2", "0", etc.
    // returns coefficient for x/y/z plus constant shift; we parse per axis below.
    // This function is used by parse_symop() helper.
    match tok.trim() {
        "" => 0.0,
        _ => tok.trim().parse::<f64>().unwrap_or(0.0),
    }
}

fn parse_symop_xyz(expr: &str) -> SymOp {
    // expr like "x,y,z" or "-y,x-y,z+1/3"
    let mut r = [[0i32; 3]; 3];
    let mut t = [0f64; 3];

    for (axis_i, term) in expr.trim().split(',').enumerate() {
        let s = term.trim().replace(" ", "");
        // parse linear part: coefficients of x,y,z are -1,0,+1
        let mut coeff = [0i32; 3];
        let mut shift = 0f64;

        // tokenize like "-x+y+1/2"
        let mut i = 0usize;
        let bytes = s.as_bytes();
        let mut cur = String::new();
        let mut parts: Vec<String> = Vec::new();
        while i < s.len() {
            let c = bytes[i] as char;
            if c == '+' || c == '-' {
                if !cur.is_empty() {
                    parts.push(cur.clone());
                    cur.clear();
                }
                cur.push(c);
            } else {
                cur.push(c);
            }
            i += 1;
        }
        if !cur.is_empty() {
            parts.push(cur);
        }

        for p in parts {
            let p = p.trim();
            if p == "x" || p == "+x" {
                coeff[0] += 1;
            } else if p == "-x" {
                coeff[0] -= 1;
            } else if p == "y" || p == "+y" {
                coeff[1] += 1;
            } else if p == "-y" {
                coeff[1] -= 1;
            } else if p == "z" || p == "+z" {
                coeff[2] += 1;
            } else if p == "-z" {
                coeff[2] -= 1;
            } else {
                // fraction like +1/2 or -2/3 or integer
                let val = if let Some((n, d)) = p.split_once('/') {
                    n.parse::<f64>().unwrap_or(0.0) / d.parse::<f64>().unwrap_or(1.0)
                } else {
                    p.parse::<f64>().unwrap_or(0.0)
                };
                shift += val;
            }
        }

        r[axis_i] = coeff;
        t[axis_i] = shift;
    }

    SymOp { r, t }
}

fn parse_mmcif_symops_and_cell(cif_text: &str) -> (Vec<SymOp>, [f32; 6], i32) {
    // Returns (symops, cell[a,b,c,alpha,beta,gamma], ispg_guess)
    // Symops from either _space_group_symop.operation_xyz or _symmetry_equiv.pos_as_xyz
    let mut ops: Vec<SymOp> = Vec::new();
    let mut in_loop = false;
    let mut in_ops = false;

    // cell
    let mut cell = [0f32; 6];
    let mut ispg = 1;

    // quick-and-dirty tag scrapes:
    for line in cif_text.lines() {
        let l = line.trim();

        if l.starts_with("_cell.length_a") {
            cell[0] = l.split_whitespace().last().unwrap().parse().unwrap_or(0.0);
        }
        if l.starts_with("_cell.length_b") {
            cell[1] = l.split_whitespace().last().unwrap().parse().unwrap_or(0.0);
        }
        if l.starts_with("_cell.length_c") {
            cell[2] = l.split_whitespace().last().unwrap().parse().unwrap_or(0.0);
        }
        if l.starts_with("_cell.angle_alpha") {
            cell[3] = l.split_whitespace().last().unwrap().parse().unwrap_or(0.0);
        }
        if l.starts_with("_cell.angle_beta") {
            cell[4] = l.split_whitespace().last().unwrap().parse().unwrap_or(0.0);
        }
        if l.starts_with("_cell.angle_gamma") {
            cell[5] = l.split_whitespace().last().unwrap().parse().unwrap_or(0.0);
        }

        if l.starts_with("_space_group.IT_number") || l.starts_with("_symmetry.Int_Tables_number") {
            if let Some(tok) = l.split_whitespace().last() {
                ispg = tok.parse::<i32>().unwrap_or(1);
            }
        }

        if l.starts_with("loop_") {
            in_loop = true;
            in_ops = false;
            continue;
        }
        if in_loop
            && (l.starts_with("_space_group_symop.operation_xyz")
                || l.starts_with("_symmetry_equiv.pos_as_xyz"))
        {
            in_ops = true;
            continue;
        }
        if in_loop
            && l.starts_with('_')
            && !l.contains("operation_xyz")
            && !l.contains("pos_as_xyz")
        {
            // another category starts
            in_ops = false;
        }

        if in_loop && in_ops {
            if l.is_empty() {
                continue;
            }
            // usually quoted: 'x,y,z'
            let trimmed = l.trim_matches('\'').trim_matches('"');
            if trimmed.contains('x') || trimmed.contains('y') || trimmed.contains('z') {
                ops.push(parse_symop_xyz(trimmed));
            }
        }
    }

    if ops.is_empty() {
        // fallback identity if CIF did not include ops
        ops.push(SymOp {
            r: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            t: [0.0, 0.0, 0.0],
        });
    }

    (ops, cell, ispg)
}

// LCM of denominators needed by symmetry translations so Nx is compatible
fn lcm(a: usize, b: usize) -> usize {
    fn gcd(mut x: usize, mut y: usize) -> usize {
        while y != 0 {
            let t = y;
            y = x % y;
            x = t;
        }
        x
    }
    a / gcd(a, b) * b
}
fn denom_of_fraction(x: f64) -> usize {
    // crude: detect 1/2, 1/3, 1/4, 1/6, etc. up to 12
    let mut best = 1usize;
    let mut err_best = f64::MAX;
    for d in 1..=12 {
        let n = (x * d as f64).round();
        let err = ((x * d as f64) - n).abs();
        if err < err_best {
            err_best = err;
            best = d;
        }
    }
    best
}

// FFT-friendly next size with factors 2,3,5
fn next_good_fft(n: usize) -> usize {
    fn is_good(mut x: usize) -> bool {
        for p in [2, 3, 5] {
            while x % p == 0 {
                x /= p;
            }
        }
        x == 1
    }
    let mut k = n.max(1);
    while !is_good(k) {
        k += 1;
    }
    k
}

// --------------- Main conversion ----------------

/// Converts the CIF 2fo-fc data to a map data. Similar to Gemmi's `sf2map` functionality.
pub fn sf_cif_to_map(txt: &str) -> io::Result<DensityMap> {
    // 2) parse symmetry ops, unit cell, ispg
    let (symops, cell, ispg) = parse_mmcif_symops_and_cell(&txt);

    // 3) parse reflections (h,k,l,FWT,PHWT) from the refln loop
    // Accept both modern mmCIF names (pdbx_FWT/pdbx_PHWT) and CCP4 aliases (FWT/PHWT)
    let mut col_idx: HashMap<String, usize> = HashMap::new();
    let mut in_loop = false;
    let mut in_refln = false;
    let mut rows: Vec<Vec<String>> = Vec::new();

    for line in txt.lines() {
        let l = line.trim();
        if l.starts_with("loop_") {
            in_loop = true;
            in_refln = false;
            rows.clear();
            continue;
        }
        if in_loop && l.starts_with("_refln.") {
            in_refln = true;
            let name = l.split_whitespace().next().unwrap().to_string();
            col_idx.insert(name, col_idx.len());
            continue;
        }
        if in_loop && in_refln {
            if l.starts_with('_') {
                // next loop starts
                in_loop = false;
                in_refln = false;
                continue;
            }
            if l.is_empty() {
                continue;
            }
            // accumulate row
            // split on whitespace, but handle quoted items
            let mut row: Vec<String> = Vec::new();
            let mut buf = String::new();
            let mut in_quote = false;
            for tok in l.split_whitespace() {
                let starts_q = tok.starts_with('\'') || tok.starts_with('"');
                let ends_q = tok.ends_with('\'') || tok.ends_with('"');
                if starts_q && !ends_q {
                    in_quote = true;
                    buf.clear();
                    buf.push_str(tok);
                } else if in_quote {
                    buf.push(' ');
                    buf.push_str(tok);
                    if ends_q {
                        in_quote = false;
                        row.push(buf.trim_matches('\'').trim_matches('"').to_string());
                        buf.clear();
                    }
                } else {
                    row.push(tok.trim_matches('\'').trim_matches('"').to_string());
                }
            }
            if !row.is_empty() {
                rows.push(row);
            }
        }
    }

    // required columns
    let hcol = col_idx
        .get("_refln.index_h")
        .copied()
        .or_else(|| col_idx.get("_refln.h").copied())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "CIF missing _refln.index_h"))?;
    let kcol = col_idx
        .get("_refln.index_k")
        .copied()
        .or_else(|| col_idx.get("_refln.k").copied())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "CIF missing _refln.index_k"))?;
    let lcol = col_idx
        .get("_refln.index_l")
        .copied()
        .or_else(|| col_idx.get("_refln.l").copied())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "CIF missing _refln.index_l"))?;

    let fcol = col_idx
        .get("_refln.pdbx_FWT")
        .copied()
        .or_else(|| col_idx.get("_refln.FWT").copied())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "CIF missing 2Fo-Fc amplitudes (pdbx_FWT/FWT)",
            )
        })?;

    let phcol = col_idx
        .get("_refln.pdbx_PHWT")
        .copied()
        .or_else(|| col_idx.get("_refln.PHWT").copied())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "CIF missing 2Fo-Fc phases (pdbx_PHWT/PHWT)",
            )
        })?;

    // 4) collect ASU reflections
    #[derive(Clone)]
    struct Ref {
        h: i32,
        k: i32,
        l: i32,
        amp: f32,
        phase_deg: f32,
    }
    let mut refls: Vec<Ref> = Vec::new();
    let mut hmax = 0i32;
    let mut kmax = 0i32;
    let mut lmax = 0i32;

    for row in &rows {
        let h = row
            .get(hcol)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0);
        let k = row
            .get(kcol)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0);
        let l = row
            .get(lcol)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0);

        let amp = row
            .get(fcol)
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(0.0);
        let phi = row
            .get(phcol)
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(0.0);

        if amp.is_finite() && phi.is_finite() {
            refls.push(Ref {
                h,
                k,
                l,
                amp,
                phase_deg: phi,
            });
            hmax = hmax.max(h.abs());
            kmax = kmax.max(k.abs());
            lmax = lmax.max(l.abs());
        }
    }

    // 5) choose reciprocal grid size (at least to hold all h,k,l; make FFT-friendly; SG-compatible)
    let mut nx = next_good_fft((hmax as usize) * 2 + 1);
    let mut ny = next_good_fft((kmax as usize) * 2 + 1);
    let mut nz = next_good_fft((lmax as usize) * 2 + 1);

    // SG translation denominators
    let mut dx = 1usize;
    let mut dy = 1usize;
    let mut dz = 1usize;
    for op in &symops {
        dx = lcm(dx, denom_of_fraction(op.t[0]));
        dy = lcm(dy, denom_of_fraction(op.t[1]));
        dz = lcm(dz, denom_of_fraction(op.t[2]));
    }
    if nx % dx != 0 {
        nx = next_good_fft(((nx + dx - 1) / dx) * dx);
    }
    if ny % dy != 0 {
        ny = next_good_fft(((ny + dy - 1) / dy) * dy);
    }
    if nz % dz != 0 {
        nz = next_good_fft(((nz + dz - 1) / dz) * dz);
    }

    // 6) put F(hkl) onto full reciprocal grid with symmetry expansion and proper phase
    let mut grid = vec![Complex32::new(0.0, 0.0); nx * ny * nz];
    let mut count = vec![0u32; nx * ny * nz];

    let idx = |h: i32, k: i32, l: i32| -> usize {
        let wrap = |q: i32, n: usize| -> usize {
            let m = n as i32;
            ((q % m) + m) as usize % n
        };
        let x = wrap(h, nx);
        let y = wrap(k, ny);
        let z = wrap(l, nz);
        (z * ny + y) * nx + x // X-fastest
    };

    let two_pi = 2.0f32 * std::f32::consts::PI;

    for r in &refls {
        let fh = Complex32::from_polar(r.amp, r.phase_deg.to_radians());

        for op in &symops {
            // ---- reciprocal mapping: h' = R^T h  (NOT R h) ----
            // R is in x' = R x + t, with rows giving x',y',z' in terms of x,y,z.
            // So R^T has columns from those rows.
            let hp = op.r[0][0] * r.h + op.r[1][0] * r.k + op.r[2][0] * r.l;
            let kp = op.r[0][1] * r.h + op.r[1][1] * r.k + op.r[2][1] * r.l;
            let lp = op.r[0][2] * r.h + op.r[1][2] * r.k + op.r[2][2] * r.l;

            // ---- phase factor: exp(-2π i h·t)  (minus sign!) ----
            let phase_shift = -two_pi
                * (r.h as f32 * op.t[0] as f32
                    + r.k as f32 * op.t[1] as f32
                    + r.l as f32 * op.t[2] as f32);
            let fhp = fh * Complex32::from_polar(1.0, phase_shift);

            // write F(h') and its Friedel mate to enforce real density
            let p = idx(hp, kp, lp);
            let p_neg = idx(-hp, -kp, -lp);

            grid[p] += fhp;
            grid[p_neg] += fhp.conj();

            count[p] += 1;
            count[p_neg] += 1;
        }
    }

    // average duplicates from different symops
    for i in 0..grid.len() {
        if count[i] > 0 {
            grid[i] /= count[i] as f32;
        }
    }

    // 7) inverse 3D FFT of the **conjugated** array (Gemmi uses ifftn(full.conj()))
    //    We'll do separable 1D FFTs: x -> y -> z
    //    rustfft does not normalize; we divide by N afterwards.
    {
        let mut planner = FftPlanner::new();
        // conj in-place
        for v in &mut grid {
            *v = v.conj();
        }

        // X dimension batched FFTs for each (y,z)
        let fft_x = planner.plan_fft_inverse(nx);
        for z in 0..nz {
            for y in 0..ny {
                let start = (z * ny + y) * nx;
                let slice = &mut grid[start..start + nx];
                fft_x.process(slice);
            }
        }
        // Y dimension: we need strided transforms
        let fft_y = planner.plan_fft_inverse(ny);
        let mut work_y = vec![Complex32::ZERO; ny];
        for z in 0..nz {
            for x in 0..nx {
                // gather
                for (iy, w) in work_y.iter_mut().enumerate() {
                    *w = grid[(z * ny + iy) * nx + x];
                }
                // fft
                fft_y.process(&mut work_y);
                // scatter back
                for (iy, w) in work_y.iter().enumerate() {
                    grid[(z * ny + iy) * nx + x] = *w;
                }
            }
        }
        // Z dimension
        let fft_z = planner.plan_fft_inverse(nz);
        let mut work_z = vec![Complex32::ZERO; nz];
        for y in 0..ny {
            for x in 0..nx {
                for (iz, w) in work_z.iter_mut().enumerate() {
                    *w = grid[(iz * ny + y) * nx + x];
                }
                fft_z.process(&mut work_z);
                for (iz, w) in work_z.iter().enumerate() {
                    grid[(iz * ny + y) * nx + x] = *w;
                }
            }
        }
    }

    // 8) take real part and scale
    let nxyz = (nx * ny * nz) as f32;
    let vol = {
        let (a, b, c, al, be, ga) = (
            cell[0] as f64,
            cell[1] as f64,
            cell[2] as f64,
            (cell[3] as f64).to_radians(),
            (cell[4] as f64).to_radians(),
            (cell[5] as f64).to_radians(),
        );
        let v = a
            * b
            * c
            * (1.0 + 2.0 * (al.cos() * be.cos() * ga.cos())
                - al.cos().powi(2)
                - be.cos().powi(2)
                - ga.cos().powi(2))
            .sqrt();
        v as f32
    };

    let mut data = Vec::with_capacity(nx * ny * nz);
    for v in &grid {
        // rustfft inverse has no 1/N; divide here. If you want exact e·Å⁻³, also divide by volume:
        data.push(v.re / nxyz /* / vol */);
    }

    // 9) build header + DensityMap
    let dmin = data.iter().cloned().fold(f32::INFINITY, f32::min);
    let dmax = data.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let dmean = data.iter().sum::<f32>() / data.len() as f32;

    let hdr = MapHeader {
        nx: nx as i32,
        ny: ny as i32,
        nz: nz as i32,
        mode: 2,
        nxstart: 0,
        nystart: 0,
        nzstart: 0,
        mx: nx as i32,
        my: ny as i32,
        mz: nz as i32,
        cell,
        mapc: 1,
        mapr: 2,
        maps: 3,
        dmin,
        dmax,
        dmean,
        ispg,
        nsymbt: 0,
        version: 20140,
        xorigin: Some(0.0),
        yorigin: Some(0.0),
        zorigin: Some(0.0),
    };

    let cell_uc = UnitCell::new(
        hdr.cell[0] as f64,
        hdr.cell[1] as f64,
        hdr.cell[2] as f64,
        hdr.cell[3] as f64,
        hdr.cell[4] as f64,
        hdr.cell[5] as f64,
    );

    let perm_f2c = [
        hdr.mapc as usize - 1,
        hdr.mapr as usize - 1,
        hdr.maps as usize - 1,
    ];
    let mut perm_c2f = [0; 3];
    for (f, c) in perm_f2c.iter().enumerate() {
        perm_c2f[*c] = f;
    }

    let origin_frac = get_origin_frac(&hdr, &cell_uc);

    // compute mean/sigma for display normalization
    let n = data.len() as f32;
    let mean = dmean;
    let variance: f32 = data.iter().map(|v| (*v - mean).powi(2)).sum::<f32>() / n;
    let sigma = variance.sqrt().max(1e-6);
    let inv_sigma = 1.0 / sigma;

    Ok(DensityMap {
        hdr,
        cell: cell_uc,
        origin_frac,
        perm_f2c,
        perm_c2f,
        data,
        mean,
        inv_sigma,
    })
}
