//! This example demonstrates how to open and save electron density data. For example, map, MTZ or
//! 2fo-fc mmCIF

use std::{io, path::Path, time::Instant};

use bio_files::{DensityMap, MapHeader, cif_sf::CifStructureFactors};
use rustfft::{FftPlanner, num_complex::Complex};

fn wrap_idx(i: i32, n: usize) -> usize {
    i.rem_euclid(n as i32) as usize
}

/// Crystal order (X, Y, Z), X contiguous.
fn idx_xfast(ix: usize, iy: usize, iz: usize, nx: usize, ny: usize) -> usize {
    ix + nx * (iy + ny * iz)
}

/// In-place inverse 3D FFT on an X-fast grid. Unscaled.
fn ifft3d(
    data: &mut [Complex<f32>],
    (nx, ny, nz): (usize, usize, usize),
    planner: &mut FftPlanner<f32>,
) {
    let ifft_x = planner.plan_fft_inverse(nx);
    for row in data.chunks_exact_mut(nx) {
        ifft_x.process(row);
    }

    let ifft_y = planner.plan_fft_inverse(ny);
    let mut buf = vec![Complex::default(); ny];
    for iz in 0..nz {
        for ix in 0..nx {
            for iy in 0..ny {
                buf[iy] = data[idx_xfast(ix, iy, iz, nx, ny)];
            }
            ifft_y.process(&mut buf);
            for iy in 0..ny {
                data[idx_xfast(ix, iy, iz, nx, ny)] = buf[iy];
            }
        }
    }

    let ifft_z = planner.plan_fft_inverse(nz);
    let mut buf = vec![Complex::default(); nz];
    for iy in 0..ny {
        for ix in 0..nx {
            for iz in 0..nz {
                buf[iz] = data[idx_xfast(ix, iy, iz, nx, ny)];
            }
            ifft_z.process(&mut buf);
            for iz in 0..nz {
                data[idx_xfast(ix, iy, iz, nx, ny)] = buf[iz];
            }
        }
    }
}

/// An example of turning structure factor data into a density map with an FFT.
/// `DensityMap::load_sf_or_mtz` does this using Gemmi instead.
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

    // Reciprocal grid in crystal order, X-fast.
    let mut data = vec![Complex::<f32>::new(0.0, 0.0); nx * ny * nz];

    for r in &sf.miller_indices {
        let c = match (r.re, r.im, r.amp, r.phase) {
            (Some(re), Some(im), _, _) => Complex::new(re, im),
            (_, _, Some(amp), Some(phase)) => Complex::from_polar(amp, phase),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Miller index missing amplitude or phase",
                ));
            }
        };

        // Place at (h, k, l), and its conjugate at (-h, -k, -l), so the result is real.
        let i0 = idx_xfast(
            wrap_idx(r.h, nx),
            wrap_idx(r.k, ny),
            wrap_idx(r.l, nz),
            nx,
            ny,
        );
        data[i0] = c;

        let i1 = idx_xfast(
            wrap_idx(-r.h, nx),
            wrap_idx(-r.k, ny),
            wrap_idx(-r.l, nz),
            nx,
            ny,
        );
        if i1 != i0 && data[i1] == Complex::new(0.0, 0.0) {
            data[i1] = c.conj();
        }
    }

    ifft3d(&mut data, (nx, ny, nz), planner);

    // The FFT is unscaled; divide by N.
    let n_vox = (nx * ny * nz) as f32;

    // Remap from crystal order to file order. File order is also fast-axis contiguous.
    let mut density = vec![0f32; nx * ny * nz];
    for k_f in 0..nz {
        for j_f in 0..ny {
            for i_f in 0..nx {
                let mut ic = [0usize; 3];
                ic[perm_f2c[0]] = i_f;
                ic[perm_f2c[1]] = j_f;
                ic[perm_f2c[2]] = k_f;

                density[idx_xfast(i_f, j_f, k_f, nx, ny)] =
                    data[idx_xfast(ic[0], ic[1], ic[2], nx, ny)].re / n_vox;
            }
        }
    }

    let dmin = density.iter().copied().fold(f32::INFINITY, f32::min);
    let dmax = density.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let dmean = (density.iter().map(|&v| v as f64).sum::<f64>() / n_vox as f64) as f32;

    let hdr = MapHeader {
        inner: sf.header.clone(),
        nx: nx as i32,
        ny: ny as i32,
        nz: nz as i32,
        mode: 2, // f32
        dmin,
        dmax,
        dmean,
    };

    println!("Complete in {} ms", start.elapsed().as_millis());

    DensityMap::new(hdr, density)
}

fn main() {
    // Load electron density structure factors data, to be processed with a FFT:
    let path = Path::new("8s6p_validation_2fo-fc_map_coef.cif");
    let data = CifStructureFactors::new_from_path(path).unwrap();

    let mut fft_planner = FftPlanner::new();
    let _dm = density_map_from_sf(&data, &mut fft_planner).unwrap();

    // For MTZ files, or 2fo-fc. Uses Gemmi; pass its path if it's not on PATH.
    let _dm = DensityMap::load_sf_or_mtz(path, None).unwrap();

    // Or if you have a Map file:
    let path = Path::new("8s6p.map");
    let dm = DensityMap::load(path).unwrap();

    // Saving a density map:
    dm.save(Path::new("8s6p.map")).unwrap();
    // Uses the file extension `.mmcif` or `.mtz` to determine which format to save as.
    dm.save_sf_or_mtz(Path::new("8s6p.mtz"), None).unwrap();
}
