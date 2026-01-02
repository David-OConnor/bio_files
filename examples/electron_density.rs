//! This example demonstrates how to open and save electron density data. For example, map, MTZ or
//! 2fo-fc mmCIF

use std::{io, path::Path, time::Instant};

use bio_files::{DensityMap, MapHeader, cif_sf::CifStructureFactors};
use ewald::fft3d_c2r;
use rustfft::{FftPlanner, num_complex::Complex};

fn wrap_idx(i: i32, n: usize) -> usize {
    let n_i32 = n as i32;
    let m = i % n_i32;
    if m < 0 {
        (m + n_i32) as usize
    } else {
        m as usize
    }
}

pub fn density_map_from_sf(
    sf: &CifStructureFactors,
    planner: &mut FftPlanner<f32>,
) -> io::Result<DensityMap> {
    println!("Computing electron density from mmCIF 2fo-fc data...");
    let start = Instant::now();

    let (nx, ny, nz) = (
        sf.header.mx as usize,
        sf.header.my as usize,
        sf.header.mz as usize,
    );

    let perm_f2c = [
        (sf.header.mapc - 1) as usize,
        (sf.header.mapr - 1) as usize,
        (sf.header.maps - 1) as usize,
    ];

    // Reciprocal grid in CRYSTAL order with X-fast layout ---
    let mut data_k = vec![Complex::<f32>::new(0.0, 0.0); nx * ny * nz];

    for r in &sf.miller_indices {
        let c = if let (Some(re), Some(im)) = (r.re, r.im) {
            Complex::new(re, im)
        } else {
            let Some(amp) = r.amp else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Miller index out of bounds",
                ));
            };
            let Some(phase) = r.phase else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Miller index out of bounds",
                ));
            };
            Complex::from_polar(amp, phase)
        };

        // Crystal grid indices
        let u = wrap_idx(r.h, nx);
        let v = wrap_idx(r.k, ny);
        let w = wrap_idx(r.l, nz);

        // place at (u,v,w) and its conjugate at (-u,-v,-w) in X-fast layout
        let i0 = idx_xfast(u, v, w, nx, ny, nz);
        data_k[i0] = c;

        let u2 = wrap_idx(-r.h, nx);
        let v2 = wrap_idx(-r.k, ny);
        let w2 = wrap_idx(-r.l, nz);

        let i1 = idx_xfast(u2, v2, w2, nx, ny, nz);

        if i1 != i0 && data_k[i1] == Complex::new(0.0, 0.0) {
            data_k[i1] = c.conj();
        }
    }

    // --- inverse FFT on X-fast layout; dims are crystal (mx,my,mz) ---
    let mut rho_crystal = fft3d_c2r(&mut data_k, (nx, ny, nz), planner);

    // If your FFT is unscaled, divide by N
    let nvox = (nx * ny * nz) as f32;
    for v in &mut rho_crystal {
        *v /= nvox;
    }

    // --- remap to FILE order buffer if you want DensityMap.data in file order ---
    // let mut rho_file = vec![0f32; mx * my * mz];
    let mut density_data = vec![0f32; nz * ny * nx];

    // file indices → crystal indices → pull from rho_crystal (Z-fast) → store in file buffer
    // todo: Experiment with order here too

    for i_f in 0..nx {
        for j_f in 0..ny {
            for k_f in 0..nz {
                let mut ic = [0usize; 3];
                ic[perm_f2c[0]] = i_f; // crystal X index
                ic[perm_f2c[1]] = j_f; // crystal Y index
                ic[perm_f2c[2]] = k_f; // crystal Z index

                let src_idx = idx_xfast(ic[0], ic[1], ic[2], nx, ny, nz);
                let dst_idx = idx_file(i_f, j_f, k_f, nx, ny);

                density_data[dst_idx] = rho_crystal[src_idx];
            }
        }
    }

    // stats on rho_file
    let mut sum = 0.0f64;
    let mut sum2 = 0.0f64;
    let mut min_v = f32::INFINITY;
    let mut max_v = f32::NEG_INFINITY;

    for &x in &density_data {
        let xd = x as f64;
        sum += xd;
        sum2 += xd * xd;
        if x < min_v {
            min_v = x;
        }
        if x > max_v {
            max_v = x;
        }
    }
    let mean = (sum / (nvox as f64)) as f32;

    let hdr = MapHeader {
        inner: sf.header.clone(),
        nx: nx as i32,
        ny: ny as i32,
        nz: nz as i32,
        mode: 2, // f32

        dmin: min_v,
        dmax: max_v,
        dmean: mean,
    };

    let elapsed = start.elapsed().as_millis();
    println!("Complete in {elapsed} ms");

    // let density_data = xfast_to_zfast(&density_data, nx, ny, nz);

    DensityMap::new(hdr, density_data)
}

fn main() {
    // Load electron density structure factors data, to be processed with a FFT:
    let path = Path::new("8s6p_validation_2fo-fc_map_coef.cif");
    let data = CifStructureFactors::new_from_path(path).unwrap();

    // These functions aren't included; an example of turning loaded structure factor data
    // into a density map.
    let mut fft_planner = FftPlanner::new();
    let dm = density_map_from_sf(&data, &mut fft_planner).unwrap();

    // For MTZ files, or 2fo-fc:
    let dm = DensityMap::load_sf_or_mtz(path, None).unwrap();

    // Or if you have a Map file:
    let path = Path::new("8s6p.map");
    let dm = DensityMap::load(path).unwrap();
}
