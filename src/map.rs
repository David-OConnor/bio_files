//! For parsing CCP4/MRC .map files. These contain electron density; they are computed
//! from reflection data, using a fourier transform.
//!
//! This [source file from Gemmi](https://github.com/project-gemmi/gemmi/blob/master/include/gemmi/ccp4.hpp)
//! is a decent ref of the spec. See for example, the `prepare_ccp4_header_except_mode_and_stats` function
//! for header fields.

use std::{
    fs,
    fs::File,
    io::{self, ErrorKind, Read, Seek, SeekFrom, Write},
    path::Path,
    process::Command,
};

use bio_apis::rcsb;
use byteorder::{LittleEndian, ReadBytesExt};
use lin_alg::f64::{Mat3, Vec3};

const HEADER_SIZE: u64 = 1_024;

#[derive(Clone, Debug)]
pub struct Density {
    /// In Å
    pub coords: Vec3,
    /// Normalized, using the unit cell volume, as reported in the reflection data.
    pub density: f64,
}

/// Minimal subset of the 1024-byte CCP4/MRC header
#[allow(unused)]
#[derive(Clone, Debug)]
pub struct MapHeader {
    /// Numbemr of grid points along each axis. nx × ny × nz is the number of voxels.
    pub nx: i32,
    pub ny: i32,
    pub nz: i32,
    pub mode: i32, // data type (0=int8, 1=int16, 2=float32, …)
    /// the `n[xyz]start` values specify the grid starting point (offset) in each dimension.
    /// They are 0 in most cryo-EM maps.
    pub nxstart: i32,
    pub nystart: i32,
    pub nzstart: i32,
    /// the m values give number of grid-sampling intervals along each axis. This is usually equal
    /// to the n values. In this case, each voxel corresponds to one spacing unit. They will differ
    /// in the case of over and under sampling.
    pub mx: i32,
    pub my: i32,
    pub mz: i32,
    /// Unit cell dimensions. [XYZ] length. Then α: Angle between Y and Z, β: Angle
    /// between X and Z, and γ: ANgle between X and Y. Distances are in Å, and angles are in degrees.
    pub cell: [f32; 6], // cell dimensions: a, b, c, alpha, beta, gamma
    pub mapc: i32,  // which axis is fast (1=X, 2=Y, 3=Z)
    pub mapr: i32,  // which axis is medium
    pub maps: i32,  // which axis is slow
    pub dmin: f32,  // minimum density value
    pub dmax: f32,  // maximum density value
    pub dmean: f32, // mean density
    pub ispg: i32,  // space group number
    /// Number of bytes used by the symmetry block. (Usually 0 for cry-EM, and 80xn for crystallography.)
    pub nsymbt: i32,
    pub version: i32,
    /// origin in Å (MRC-2014) **or** derived from *\*START if absent**
    pub xorigin: Option<f32>,
    pub yorigin: Option<f32>,
    pub zorigin: Option<f32>, // todo: More header items A/R.
}

fn read_map_header<R: Read + Seek>(mut r: R) -> io::Result<MapHeader> {
    // ensure we are at the start
    r.seek(SeekFrom::Start(0))?;

    // Classic part of the header.
    let nx = r.read_i32::<LittleEndian>()?;
    let ny = r.read_i32::<LittleEndian>()?;
    let nz = r.read_i32::<LittleEndian>()?;

    let mode = r.read_i32::<LittleEndian>()?;

    // Word 5
    let nxstart = r.read_i32::<LittleEndian>()?;
    let nystart = r.read_i32::<LittleEndian>()?;
    let nzstart = r.read_i32::<LittleEndian>()?;

    let mx = r.read_i32::<LittleEndian>()?;
    let my = r.read_i32::<LittleEndian>()?;
    let mz = r.read_i32::<LittleEndian>()?;

    // Word 11
    let mut cell = [0f32; 6];
    for c in &mut cell {
        *c = r.read_f32::<LittleEndian>()?;
    }

    // Word 17
    let mapc = r.read_i32::<LittleEndian>()?;
    let mapr = r.read_i32::<LittleEndian>()?;
    let maps = r.read_i32::<LittleEndian>()?;

    let dmin = r.read_f32::<LittleEndian>()?;
    let dmax = r.read_f32::<LittleEndian>()?;
    let dmean = r.read_f32::<LittleEndian>()?;

    // Word 23
    let ispg = r.read_i32::<LittleEndian>()?;
    let nsymbt = r.read_i32::<LittleEndian>()?;

    r.seek(SeekFrom::Start(27 * 4))?;
    let version = r.read_i32::<LittleEndian>()?; // e.g. 20140

    // offset notes:
    // nystart / my * cell_b. e.g. 4/60 * 85.142 = 5.6

    // let mut test_str = String::with_capacity(4);
    // r.read_to_string(&mut test_str)?;
    // println!("TEST STR: {:?}", test_str);

    // Words 25-49 are “extra”; skip straight to word 50 (49 × 4 bytes)
    r.seek(SeekFrom::Start(49 * 4))?;

    // words 50-52 = XORIGIN, YORIGIN, ZORIGIN   (MRC-2014)
    let mut xorigin_ = r.read_f32::<LittleEndian>()?;
    let mut yorigin_ = r.read_f32::<LittleEndian>()?;
    let mut zorigin_ = r.read_f32::<LittleEndian>()?;

    let mut xorigin = None;
    let mut yorigin = None;
    let mut zorigin = None;

    let mut tag = [0u8; 4];
    r.seek(SeekFrom::Start(52 * 4))?;
    r.read_exact(&mut tag)?;

    if &tag != b"MAP " {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "Invalid MAP tag in header.",
        ));
    }

    const EPS: f32 = 0.0001;

    if xorigin_.abs() > EPS {
        xorigin = Some(xorigin_);
    }

    if yorigin_.abs() > EPS {
        yorigin = Some(yorigin_);
    }

    if zorigin_.abs() > EPS {
        zorigin = Some(zorigin_);
    }

    // skip ahead to end of 1024-byte header
    r.seek(SeekFrom::Start(HEADER_SIZE))?;

    Ok(MapHeader {
        nx,
        ny,
        nz,
        mode,
        nxstart,
        nystart,
        nzstart,
        mx,
        my,
        mz,
        cell,
        mapc,
        mapr,
        maps,
        dmin,
        dmax,
        dmean,
        ispg,
        nsymbt,
        version,
        xorigin,
        yorigin,
        zorigin,
    })
}

// todo: Good candidate to make generic? Used in both reflections/electron map, and docking.
#[derive(Clone, Debug)]
pub struct UnitCell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    /// For converting fractional to cartesian coordinates
    pub ortho: Mat3,
    /// For converting cartesian to fractional coordinates.
    pub ortho_inv: Mat3,
}

impl UnitCell {
    pub fn new(a: f64, b: f64, c: f64, alpha_deg: f64, beta_deg: f64, gamma_deg: f64) -> Self {
        let (α, β, γ) = (
            alpha_deg.to_radians(),
            beta_deg.to_radians(),
            gamma_deg.to_radians(),
        );

        // components of the three cell vectors in Cartesian space
        let v_a = Vec3::new(a, 0.0, 0.0);
        let v_b = Vec3::new(b * γ.cos(), b * γ.sin(), 0.0);

        let cx = c * β.cos();
        let cy = c * (α.cos() - β.cos() * γ.cos()) / γ.sin();
        let cz = c * (1.0 - β.cos().powi(2) - cy.powi(2) / c.powi(2)).sqrt();
        let v_c = Vec3::new(cx, cy, cz);

        // 3×3 matrix whose columns are a,b,c
        let ortho = Mat3::from_cols(v_a, v_b, v_c);
        let ortho_inv = ortho.inverse().expect("unit-cell matrix is singular");

        Self {
            a,
            b,
            c,
            alpha: α,
            beta: β,
            gamma: γ,
            ortho,
            ortho_inv,
        }
    }

    #[inline]
    pub fn fractional_to_cartesian(&self, f: Vec3) -> Vec3 {
        // todo: Don't clone!
        self.ortho.clone() * f
    }

    #[inline]
    pub fn cartesian_to_fractional(&self, c: Vec3) -> Vec3 {
        // todo: Don't clone!
        self.ortho_inv.clone() * c
    }
}

/// Load the header, and density from Map data.
fn read_header_dens<R: Read + Seek>(data: &mut R) -> io::Result<(MapHeader, Vec<f32>)> {
    let hdr = read_map_header(&mut *data)?;

    if hdr.mode != 2 {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            format!("Unsupported mode: {}", hdr.mode),
        ));
    }

    data.seek(SeekFrom::Start(HEADER_SIZE + hdr.nsymbt as u64))?;

    let n = (hdr.nx * hdr.ny * hdr.nz) as usize;
    let mut dens = Vec::with_capacity(n);

    for _ in 0..n {
        dens.push(data.read_f32::<LittleEndian>()?);
    }

    Ok((hdr, dens))
}

fn get_origin_frac(hdr: &MapHeader, cell: &UnitCell) -> Vec3 {
    if let (Some(ox), Some(oy), Some(oz)) = (hdr.xorigin, hdr.yorigin, hdr.zorigin) {
        cell.cartesian_to_fractional(Vec3::new(ox as f64, oy as f64, oz as f64))
    } else {
        Vec3::new(
            hdr.nxstart as f64 / hdr.mx as f64,
            hdr.nystart as f64 / hdr.my as f64,
            hdr.nzstart as f64 / hdr.mz as f64,
        )
    }
}

#[derive(Clone, Debug)]
/// Represents electron density data. Set up in a flexible way that can be overlayed over
/// atom coordinates, with symmetry considerations taken into account.
pub struct DensityMap {
    pub hdr: MapHeader,
    pub cell: UnitCell,
    /// header origin, already converted to fractional
    pub origin_frac: Vec3,
    /// file-axis → cryst-axis (from MAPC/MAPR/MAPS)
    pub perm_f2c: [usize; 3],
    /// cryst-axis → file-axis (inverse permutation)
    pub perm_c2f: [usize; 3],
    pub data: Vec<f32>,
}

impl DensityMap {
    /// Create a new density map, e.g. from a File or byte array.
    pub fn new<R: Read + Seek>(data: &mut R) -> io::Result<Self> {
        let (hdr, data) = read_header_dens(data)?;

        let cell = UnitCell::new(
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

        let origin_frac = get_origin_frac(&hdr, &cell);

        Ok(Self {
            hdr,
            cell,
            origin_frac,
            perm_f2c,
            perm_c2f,
            data,
        })
    }

    /// Uses nearest-neighbour lookup to calculate density at a point.
    pub fn density_at_point(&self, cart: Vec3) -> f32 {
        // Cartesian to fractional (wrap to [0,1) )
        let mut frac = self.cell.cartesian_to_fractional(cart);
        frac.x -= frac.x.floor();
        frac.y -= frac.y.floor();
        frac.z -= frac.z.floor();

        // Fractional to crystallographic voxel index.
        let ic = [
            (frac.x * self.hdr.mx as f64 - 0.5).round() as isize,
            (frac.y * self.hdr.my as f64 - 0.5).round() as isize,
            (frac.z * self.hdr.mz as f64 - 0.5).round() as isize,
        ];

        // Crystallographic ➜ file order
        let ifile = [
            pmod(ic[self.perm_c2f[0]], self.hdr.nx as usize),
            pmod(ic[self.perm_c2f[1]], self.hdr.ny as usize),
            pmod(ic[self.perm_c2f[2]], self.hdr.nz as usize),
        ];

        // Linear offset in file order (x is fastest dimesion, from the experimental data)
        let offset = (ifile[2] * self.hdr.ny as usize + ifile[1]) * self.hdr.nx as usize + ifile[0];

        self.data[offset]
    }

    /// Load a map from file.
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;

        Self::new(&mut file)
    }

    // todo: Implement
    // pub fn save(&self, path: &Path) -> io::Result<()> {
    //
    // }
}

/// Assumes `gemmi` is available on the path.
pub fn gemmi_cif_to_map(cif_path: &str) -> io::Result<DensityMap> {
    let _status = Command::new("gemmi")
        .args(["sf2map", cif_path, "temp_map.map"])
        .status()?;

    let map = DensityMap::load(Path::new("temp_map.map"))?;

    fs::remove_file(Path::new(cif_path))?;
    fs::remove_file(Path::new("temp_map.map"))?;

    Ok(map)
}

/// Assumes `gemmi` is available on the path.
pub fn density_from_2fo_fc_rcsb_gemmi(ident: &str) -> io::Result<DensityMap> {
    println!("Downloading Map data for {ident}...");

    let map_2fo_fc = rcsb::load_validation_2fo_fc_cif(ident)
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Problem loading 2fo-fc from RCSB"))?;

    fs::write("temp_map.cif", map_2fo_fc)?;

    let _status = Command::new("gemmi")
        .args(["sf2map", "temp_map.cif", "temp_map.map"])
        .status()?;

    let map = DensityMap::load(Path::new("temp_map.map"))?;

    fs::remove_file(Path::new("temp_map.cif"))?;
    fs::remove_file(Path::new("temp_map.map"))?;

    Ok(map)
}

/// Positive modulus that always lands in 0..n-1
fn pmod(i: isize, n: usize) -> usize {
    ((i % n as isize) + n as isize) as usize % n
}
