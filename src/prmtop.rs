use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{self, Read, Write},
    path::Path,
};

use crate::{
    AtomGeneric,
    md_params::{ForceFieldParams, LjParams, MassParams},
};

const AMBER_CHARGE_SCALE: f32 = 18.2223; // prmtop stores q * 18.2223 (real q = stored/18.2223)
const INT_WIDTH: usize = 8;
const STR4_PER_LINE: usize = 20;
const INT_PER_LINE: usize = 10;
const FLO_PER_LINE: usize = 5;

fn wline(file: &mut File, s: &str) -> io::Result<()> {
    file.write_all(s.as_bytes())?;
    file.write_all(b"\n")
}

fn fmt_i(v: i32) -> String {
    format!("{v:>width$}", width = INT_WIDTH)
}

fn fmt_e(v: f32) -> String {
    // Amber uses E16.8; this matches width/precision.
    format!("{:>16.8E}", v as f64)
}

fn fmt_a4(s: &str) -> String {
    let mut t = s.chars().take(4).collect::<String>();
    while t.len() < 4 {
        t.push(' ');
    }
    t
}

fn tri_index(i: usize, j: usize) -> usize {
    // 0-based upper-tri (i<=j): idx = i*nt - i*(i-1)/2 + (j-i), but easier with formula below
    if j < i {
        return tri_index(j, i);
    }
    // number of elements in rows < i + offset in row i
    i * (i + 1) / 2 + (j - i)
}

fn write_flag<T: Fn(usize) -> String>(
    file: &mut File,
    name: &str,
    fmt: &str,
    n_per_line: usize,
    n_items: usize,
    item: T,
) -> io::Result<()> {
    wline(file, &format!("%FLAG {name}"))?;
    wline(file, &format!("%FORMAT({fmt})"))?;
    if n_items == 0 {
        return Ok(());
    }
    let mut line = String::new();
    for i in 0..n_items {
        if i > 0 && i % n_per_line == 0 {
            wline(file, &line)?;
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(&item(i));
    }
    if !line.is_empty() {
        wline(file, &line)?;
    }
    Ok(())
}

pub fn save_prmtop(
    atoms: &[AtomGeneric],
    params: &ForceFieldParams,
    path: &Path,
) -> io::Result<()> {
    // Collect LJ types actually used by atoms (stable order).
    let mut used_types: BTreeMap<String, ()> = BTreeMap::new();
    for a in atoms {
        let t = a.force_field_type.as_ref().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "Atom missing force_field_type")
        })?;
        used_types.insert(t.clone(), ());
    }
    let type_names: Vec<String> = used_types.keys().cloned().collect();
    let nt = type_names.len();
    if nt == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No atom types found",
        ));
    }

    // Per-type LJ (require present).
    let mut sigma: Vec<f32> = Vec::with_capacity(nt);
    let mut eps: Vec<f32> = Vec::with_capacity(nt);
    for t in &type_names {
        let lj = params.lennard_jones.get(t).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Missing LJ params for type {t}"),
            )
        })?;
        sigma.push(lj.sigma);
        eps.push(lj.eps);
    }

    // Map type -> 1-based LJ index
    let mut type_to_idx: HashMap<&str, usize> = HashMap::new();
    for (i, t) in type_names.iter().enumerate() {
        type_to_idx.insert(t.as_str(), i + 1);
    }

    // ATOM_TYPE_INDEX (1-based)
    let mut atom_type_index: Vec<i32> = Vec::with_capacity(atoms.len());
    for a in atoms {
        let t = a.force_field_type.as_ref().unwrap();
        let idx = *type_to_idx.get(t.as_str()).unwrap();
        atom_type_index.push(idx as i32);
    }

    // NONBONDED_PARM_INDEX (nt * nt): 1-based pointer into triangular A/B
    let mut nb_index: Vec<i32> = Vec::with_capacity(nt * nt);
    for i in 0..nt {
        for j in 0..nt {
            let k = tri_index(i, j) + 1;
            nb_index.push(k as i32);
        }
    }

    // Triangular LENNARD_JONES_{A,B} from Lorentz–Berthelot on (sigma, eps)
    let ntri = nt * (nt + 1) / 2;
    let mut acoef = vec![0f32; ntri];
    let mut bcoef = vec![0f32; ntri];
    for i in 0..nt {
        for j in i..nt {
            let sij = 0.5 * (sigma[i] + sigma[j]);
            let eij = (eps[i] * eps[j]).sqrt();
            let s2 = sij * sij;
            let s3 = s2 * sij;
            let s6 = s3 * s3;
            let s12 = s6 * s6;
            let a = 4.0 * eij * s12;
            let b = 4.0 * eij * s6;
            let k = tri_index(i, j);
            acoef[k] = a;
            bcoef[k] = b;
        }
    }

    // Per-atom arrays
    let mut amber_types_4: Vec<String> = Vec::with_capacity(atoms.len());
    let mut atom_names_4: Vec<String> = Vec::with_capacity(atoms.len());
    let mut charges_stored: Vec<f32> = Vec::with_capacity(atoms.len());
    let mut masses: Vec<f32> = Vec::with_capacity(atoms.len());

    for a in atoms {
        let t = a.force_field_type.as_ref().unwrap();
        amber_types_4.push(fmt_a4(t));

        // Use type_in_res if present for ATOM_NAME, otherwise reuse ff type.
        let nm = a
            .type_in_res
            .as_ref()
            .map(|x| format!("{x:?}")) // fallback; adjust if your AtomTypeInRes has Display
            .unwrap_or_else(|| t.clone());
        atom_names_4.push(fmt_a4(&nm));

        let q = a.partial_charge.unwrap_or(0.0) * AMBER_CHARGE_SCALE;
        charges_stored.push(q);

        let m = params.mass.get(t).map(|x| x.mass).unwrap_or(0.0);
        masses.push(m);
    }

    // Minimal RESIDUE_* (single residue spanning all atoms)
    let nres = 1_i32;
    let residue_labels = [fmt_a4("SYS")];
    let residue_ptr = [1_i32]; // 1-based start of residue

    // POINTERS (31 ints; NCOPY present and set to 1)
    // Order per Amber spec.
    let mut pointers: [i32; 31] = [0; 31];
    pointers[0] = atoms.len() as i32; // NATOM
    pointers[1] = nt as i32; // NTYPES
    pointers[11] = nres; // NRES
    pointers[18] = 0; // NPHB
    pointers[20] = 0; // NBPER
    pointers[21] = 0; // NGPER
    pointers[22] = 0; // NDPER
    pointers[23] = 0; // MBPER
    pointers[24] = 0; // MGPER
    pointers[25] = 0; // MDPER
    pointers[26] = 0; // IFBOX
    pointers[27] = 0; // NMXRS
    pointers[28] = 0; // IFCAP
    pointers[29] = 0; // NUMEXTRA
    pointers[30] = 1; // NCOPY

    let mut f = File::create(path)?;

    // VERSION and TITLE (minimal)
    wline(
        &mut f,
        "%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/01  00:00:00",
    )?;
    write_flag(&mut f, "TITLE", "20a4", STR4_PER_LINE, 1, |_i| {
        fmt_a4("GENERATED")
    })?;

    // POINTERS
    write_flag(
        &mut f,
        "POINTERS",
        "10I8",
        INT_PER_LINE,
        pointers.len(),
        |i| fmt_i(pointers[i]),
    )?;

    // AMBER_ATOM_TYPE / ATOM_NAME
    write_flag(
        &mut f,
        "AMBER_ATOM_TYPE",
        "20a4",
        STR4_PER_LINE,
        amber_types_4.len(),
        |i| amber_types_4[i].clone(),
    )?;
    write_flag(
        &mut f,
        "ATOM_NAME",
        "20a4",
        STR4_PER_LINE,
        atom_names_4.len(),
        |i| atom_names_4[i].clone(),
    )?;

    // CHARGE, MASS
    write_flag(
        &mut f,
        "CHARGE",
        "5E16.8",
        FLO_PER_LINE,
        charges_stored.len(),
        |i| fmt_e(charges_stored[i]),
    )?;
    write_flag(&mut f, "MASS", "5E16.8", FLO_PER_LINE, masses.len(), |i| {
        fmt_e(masses[i])
    })?;

    // Type indices and NB tables
    write_flag(
        &mut f,
        "ATOM_TYPE_INDEX",
        "10I8",
        INT_PER_LINE,
        atom_type_index.len(),
        |i| fmt_i(atom_type_index[i]),
    )?;
    write_flag(
        &mut f,
        "NONBONDED_PARM_INDEX",
        "10I8",
        INT_PER_LINE,
        nb_index.len(),
        |i| fmt_i(nb_index[i]),
    )?;

    // LJ coefficients (triangular)
    write_flag(
        &mut f,
        "LENNARD_JONES_ACOEF",
        "5E16.8",
        FLO_PER_LINE,
        acoef.len(),
        |i| fmt_e(acoef[i]),
    )?;
    write_flag(
        &mut f,
        "LENNARD_JONES_BCOEF",
        "5E16.8",
        FLO_PER_LINE,
        bcoef.len(),
        |i| fmt_e(bcoef[i]),
    )?;

    // Minimal residues
    write_flag(
        &mut f,
        "RESIDUE_LABEL",
        "20a4",
        STR4_PER_LINE,
        residue_labels.len(),
        |i| residue_labels[i].clone(),
    )?;
    write_flag(
        &mut f,
        "RESIDUE_POINTER",
        "10I8",
        INT_PER_LINE,
        residue_ptr.len(),
        |i| fmt_i(residue_ptr[i]),
    )?;

    Ok(())
}

pub fn load_prmtop(path: &Path) -> io::Result<(Vec<AtomGeneric>, ForceFieldParams)> {
    let mut file = File::open(path)?;
    let mut buf = String::new();
    file.read_to_string(&mut buf)?;

    #[derive(Default)]
    struct Block {
        fmt: String,
        data: Vec<String>,
    }
    let mut blocks: HashMap<String, Block> = HashMap::new();

    let mut cur: Option<String> = None;
    for line in buf.lines() {
        let line = line.trim_end();
        if line.starts_with("%FLAG") {
            let name = line.split_whitespace().nth(1).unwrap().to_string();
            blocks.entry(name.clone()).or_default();
            cur = Some(name);
        } else if line.starts_with("%FORMAT") {
            if let Some(k) = &cur {
                blocks.get_mut(k).unwrap().fmt = line["%FORMAT(".len()..line.len() - 1].to_string();
            }
        } else if let Some(k) = &cur {
            let b = blocks.get_mut(k).unwrap();
            if !line.is_empty() {
                // Split into tokens, preserving 4-char fields for 20a4 by whitespace split (OK).
                b.data
                    .extend(line.split_whitespace().map(|s| s.to_string()));
            }
        }
    }

    let get_i = |s: &str| -> io::Result<i32> {
        s.parse::<i32>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    };
    let get_f = |s: &str| -> io::Result<f32> {
        s.parse::<f64>()
            .map(|x| x as f32)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    };

    // POINTERS (order per spec)
    let p = blocks
        .get("POINTERS")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing POINTERS"))?;
    if p.data.len() < 31 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "POINTERS too short",
        ));
    }
    let mut pointers = [0i32; 31];
    for (i, pointer) in pointers.iter_mut().enumerate() {
        *pointer = get_i(&p.data[i])?;
    }

    let natom = pointers[0] as usize;
    let ntypes = pointers[1] as usize;
    let _nres = pointers[11] as usize;

    // Arrays we use
    let atm_types = blocks
        .get("AMBER_ATOM_TYPE")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing AMBER_ATOM_TYPE"))?;
    if atm_types.data.len() < natom {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "AMBER_ATOM_TYPE too short",
        ));
    }
    let mut type_names: Vec<String> = Vec::with_capacity(natom);
    for i in 0..natom {
        type_names.push(atm_types.data[i].trim().to_string());
    }

    let charge_b = blocks
        .get("CHARGE")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing CHARGE"))?;
    if charge_b.data.len() < natom {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "CHARGE too short",
        ));
    }
    let mut charges: Vec<f32> = Vec::with_capacity(natom);
    for i in 0..natom {
        charges.push(get_f(&charge_b.data[i])? / AMBER_CHARGE_SCALE);
    }

    let mass_b = blocks
        .get("MASS")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing MASS"))?;
    if mass_b.data.len() < natom {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "MASS too short"));
    }
    let mut masses: Vec<f32> = Vec::with_capacity(natom);
    for i in 0..natom {
        masses.push(get_f(&mass_b.data[i])?);
    }

    let ati_b = blocks
        .get("ATOM_TYPE_INDEX")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing ATOM_TYPE_INDEX"))?;
    if ati_b.data.len() < natom {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "ATOM_TYPE_INDEX too short",
        ));
    }
    let mut atom_type_index: Vec<usize> = Vec::with_capacity(natom);
    for i in 0..natom {
        atom_type_index.push((get_i(&ati_b.data[i])? as usize).max(1) - 1);
    }

    let nb_b = blocks.get("NONBONDED_PARM_INDEX").ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "Missing NONBONDED_PARM_INDEX")
    })?;
    if nb_b.data.len() < ntypes * ntypes {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "NONBONDED_PARM_INDEX too short",
        ));
    }
    let mut nb_index: Vec<usize> = Vec::with_capacity(ntypes * ntypes);
    for i in 0..ntypes * ntypes {
        nb_index.push((get_i(&nb_b.data[i])? as usize).max(1) - 1);
    }

    let a_b = blocks
        .get("LENNARD_JONES_ACOEF")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing LENNARD_JONES_ACOEF"))?;
    let b_b = blocks
        .get("LENNARD_JONES_BCOEF")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing LENNARD_JONES_BCOEF"))?;
    let ntri = ntypes * (ntypes + 1) / 2;
    if a_b.data.len() < ntri || b_b.data.len() < ntri {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "LJ coef arrays too short",
        ));
    }
    let mut acoef: Vec<f32> = Vec::with_capacity(ntri);
    let mut bcoef: Vec<f32> = Vec::with_capacity(ntri);
    for i in 0..ntri {
        acoef.push(get_f(&a_b.data[i])?);
        bcoef.push(get_f(&b_b.data[i])?);
    }

    // Recover per-type sigma, eps from diagonal A/B:
    // For pair V(r)=A/r^12 - B/r^6, diagonal (i,i): Rmin_ii = (2A/B)^(1/6); eps_ii = B^2/(4A).
    // For atom-type parameters compatible with LB on (sigma, eps):
    // sigma_i = (A/B)^(1/6); eps_i = B^2/(4A) on the diagonal.
    let mut lj_by_type: Vec<(f32, f32)> = vec![(0.0, 0.0); ntypes];
    for i in 0..ntypes {
        let k = nb_index[i * ntypes + i]; // pointer into triangular
        let a = acoef[k];
        let b = bcoef[k];
        if a <= 0.0 || b <= 0.0 {
            lj_by_type[i] = (0.0, 0.0);
        } else {
            let sigma_i = (a / b).powf(1.0 / 6.0);
            let eps_i = (b * b) / (4.0 * a);
            lj_by_type[i] = (sigma_i, eps_i);
        }
    }

    // Build outputs
    let mut atoms_out = Vec::with_capacity(natom);

    for i in 0..natom {
        let a = AtomGeneric {
            serial_number: (i + 1) as u32,
            force_field_type: Some(type_names[i].clone()),
            partial_charge: Some(charges[i]),
            ..Default::default()
        };

        // element/posit/occupancy left as defaults
        atoms_out.push(a);
    }

    let mut ff = ForceFieldParams::default();
    // mass per type: first occurrence wins
    let mut seen_mass: HashMap<String, f32> = HashMap::new();
    for i in 0..natom {
        let tname = &type_names[i];
        let ti = atom_type_index[i];
        let (sigma_i, eps_i) = lj_by_type[ti];
        ff.lennard_jones.entry(tname.clone()).or_insert(LjParams {
            atom_type: tname.clone(),
            sigma: sigma_i,
            eps: eps_i,
        });
        if !seen_mass.contains_key(tname) {
            seen_mass.insert(tname.clone(), masses[i]);
            ff.mass.insert(
                tname.clone(),
                MassParams {
                    atom_type: tname.clone(),
                    mass: masses[i],
                    comment: None,
                },
            );
        }
    }

    Ok((atoms_out, ff))
}

// /// Create an Amber PRMTOP file from atom and forcefield data.
// pub fn save_prmtop(
//     atoms: &[AtomGeneric],
//     params: &ForceFieldParams,
//     path: &Path,
// ) -> io::Result<()> {
//     Ok(())
// }
//
// /// Load atom and forcefield data from an AMBER PRMTOP file.
// pub fn load_prmtop(path: &Path) -> io::Result<(Vec<AtomGeneric>, ForceFieldParams)> {
//     let mut file = File::open(path)?;
//     let mut buffer = Vec::new();
//     file.read_to_end(&mut buffer)?;
// }

// ---------------------------------------------------------------------------------------------
// Full topology reader. `load_prmtop` above extracts types, charges, masses, and LJ parameters
// only; `AmberPrmtop` reads every term, per atom and per interaction.
// ---------------------------------------------------------------------------------------------

/// Charges in prmtop files are stored as q * 18.2223 (the square root of Amber's electrostatic
/// constant, 332.0522 kcal·Å/(mol·e²)); we store them in elementary charge.
/// [Format reference](https://ambermd.org/prmtop.pdf)
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AmberPrmtop {
    pub title: String,
    pub atom_names: Vec<String>,
    /// E.g. "CT", "HC". Called `AMBER_ATOM_TYPE` in the file.
    pub atom_types: Vec<String>,
    /// Elementary charge.
    pub charges: Vec<f32>,
    /// amu
    pub masses: Vec<f32>,
    /// Absent in older files.
    pub atomic_numbers: Option<Vec<i32>>,
    /// Per atom: A 0-based index into the LJ type tables.
    pub lj_type_index: Vec<usize>,
    pub n_lj_types: usize,
    /// Per LJ type pair (row-major, `n_lj_types`²): A 1-based index into `lj_acoef` and
    /// `lj_bcoef`. Negative values index 10-12 H-bond terms instead.
    pub nonbonded_parm_index: Vec<i32>,
    /// E = A/r¹² − B/r⁶, in kcal/mol and Å.
    pub lj_acoef: Vec<f64>,
    pub lj_bcoef: Vec<f64>,
    pub residue_labels: Vec<String>,
    /// The 0-based index of each residue's first atom.
    pub residue_starts: Vec<usize>,
    pub bonds: Vec<PrmtopBond>,
    pub angles: Vec<PrmtopAngle>,
    /// Includes impropers; see `PrmtopDihedral::improper`.
    pub dihedrals: Vec<PrmtopDihedral>,
    /// Pairs excluded from normal non-bonded interactions: 1-2, 1-3, and 1-4 pairs. 1-4 pairs
    /// interact through the scaled terms of their dihedrals instead. (i < j)
    pub excluded_pairs: Vec<(usize, usize)>,
    /// Box lengths (Å) and angle β (degrees), if periodic.
    pub box_dims: Option<PrmtopBox>,
    /// Number of extra points (e.g. TIP4P's massless charge sites).
    pub num_extra_points: usize,
    /// Terms present in the file whose energy contributions aren't represented by the fields
    /// above. E.g. CHARMM Urey-Bradley and improper terms, CMAP, 1-4 LJ tables, and polarizability.
    pub unsupported_terms: Vec<String>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrmtopBox {
    /// Degrees
    pub beta: f32,
    /// Å
    pub lengths: [f32; 3],
}

/// E = k (r − r0)²
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrmtopBond {
    pub atoms: (usize, usize),
    /// kcal/mol/Å²
    pub k: f32,
    /// Å
    pub r0: f32,
}

/// E = k (θ − θ0)²
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrmtopAngle {
    pub atoms: (usize, usize, usize),
    /// kcal/mol/rad²
    pub k: f32,
    /// Radians
    pub theta0: f32,
}

/// E = k (1 + cos(n φ − phase)). Several terms can share the same atoms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrmtopDihedral {
    pub atoms: [usize; 4],
    /// kcal/mol
    pub k: f32,
    pub periodicity: f32,
    /// Radians
    pub phase: f32,
    /// The 1-4 electrostatic interaction is divided by this. Defaults to 1.2 in files that
    /// predate the per-dihedral scale factors.
    pub scee: f32,
    /// The 1-4 LJ interaction is divided by this. Defaults to 2.0.
    pub scnb: f32,
    /// Impropers conventionally list the central atom third.
    pub improper: bool,
    /// If true, this term doesn't contribute a 1-4 interaction for its end atoms, e.g. because
    /// another term with the same end atoms does, or because they're also 1-2 or 1-3 in a ring.
    pub skip_14: bool,
}

impl AmberPrmtop {
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut text = String::new();
        File::open(path)?.read_to_string(&mut text)?;
        Self::new(&text)
    }

    pub fn new(text: &str) -> io::Result<Self> {
        let sections = PrmtopSections::parse(text)?;

        let pointers = sections.ints("POINTERS")?;
        if pointers.len() < 30 {
            return Err(invalid("POINTERS section is too short"));
        }
        let ptr = |i: usize| pointers[i].max(0) as usize;
        let n_atoms = ptr(0);
        let n_lj_types = ptr(1);
        let (n_bonds_h, n_bonds) = (ptr(2), ptr(12));
        let (n_angles_h, n_angles) = (ptr(4), ptr(13));
        let (n_dihedrals_h, n_dihedrals) = (ptr(6), ptr(14));
        let n_residues = ptr(11);
        let has_box = ptr(27) > 0;
        let num_extra_points = pointers.get(30).map(|v| (*v).max(0) as usize).unwrap_or(0);

        let atom_names = sections.strings_n("ATOM_NAME", n_atoms)?;
        let atom_types = sections.strings_n("AMBER_ATOM_TYPE", n_atoms)?;
        let charges = sections
            .floats_n("CHARGE", n_atoms)?
            .into_iter()
            .map(|q| (q / AMBER_CHARGE_SCALE as f64) as f32)
            .collect();
        let masses = sections
            .floats_n("MASS", n_atoms)?
            .into_iter()
            .map(|m| m as f32)
            .collect();
        let atomic_numbers = if sections.has("ATOMIC_NUMBER") {
            Some(sections.ints_n("ATOMIC_NUMBER", n_atoms)?)
        } else {
            None
        };

        let lj_type_index = sections
            .ints_n("ATOM_TYPE_INDEX", n_atoms)?
            .into_iter()
            .map(|i| (i.max(1) - 1) as usize)
            .collect();
        let nonbonded_parm_index =
            sections.ints_n("NONBONDED_PARM_INDEX", n_lj_types * n_lj_types)?;
        let n_lj_pairs = n_lj_types * (n_lj_types + 1) / 2;
        let lj_acoef = sections.floats_n("LENNARD_JONES_ACOEF", n_lj_pairs)?;
        let lj_bcoef = sections.floats_n("LENNARD_JONES_BCOEF", n_lj_pairs)?;

        let residue_labels = sections.strings_n("RESIDUE_LABEL", n_residues)?;
        let residue_starts = sections
            .ints_n("RESIDUE_POINTER", n_residues)?
            .into_iter()
            .map(|i| (i.max(1) - 1) as usize)
            .collect();

        // Bonded parameter tables, indexed by the terms below (1-based).
        let bond_k = sections.floats("BOND_FORCE_CONSTANT")?;
        let bond_r0 = sections.floats("BOND_EQUIL_VALUE")?;
        let angle_k = sections.floats("ANGLE_FORCE_CONSTANT")?;
        let angle_theta0 = sections.floats("ANGLE_EQUIL_VALUE")?;
        let dihe_k = sections.floats("DIHEDRAL_FORCE_CONSTANT")?;
        let dihe_n = sections.floats("DIHEDRAL_PERIODICITY")?;
        let dihe_phase = sections.floats("DIHEDRAL_PHASE")?;
        let scee = if sections.has("SCEE_SCALE_FACTOR") {
            sections.floats("SCEE_SCALE_FACTOR")?
        } else {
            vec![1.2; dihe_k.len()]
        };
        let scnb = if sections.has("SCNB_SCALE_FACTOR") {
            sections.floats("SCNB_SCALE_FACTOR")?
        } else {
            vec![2.0; dihe_k.len()]
        };

        // Atom indices in term lists are stored as 3 × (0-based index), i.e. coordinate-array
        // offsets. Signs carry flags for dihedrals.
        let atom_i = |v: i32| (v.unsigned_abs() / 3) as usize;
        let param = |table: &[f64], i: i32, label: &str| -> io::Result<f32> {
            table
                .get((i.max(1) - 1) as usize)
                .map(|v| *v as f32)
                .ok_or_else(|| invalid(&format!("{label} parameter index {i} is out of range")))
        };

        let mut bonds = Vec::with_capacity(n_bonds_h + n_bonds);
        for (flag, n) in [
            ("BONDS_INC_HYDROGEN", n_bonds_h),
            ("BONDS_WITHOUT_HYDROGEN", n_bonds),
        ] {
            for t in sections.ints_n(flag, 3 * n)?.chunks_exact(3) {
                bonds.push(PrmtopBond {
                    atoms: (atom_i(t[0]), atom_i(t[1])),
                    k: param(&bond_k, t[2], "Bond")?,
                    r0: param(&bond_r0, t[2], "Bond")?,
                });
            }
        }

        let mut angles = Vec::with_capacity(n_angles_h + n_angles);
        for (flag, n) in [
            ("ANGLES_INC_HYDROGEN", n_angles_h),
            ("ANGLES_WITHOUT_HYDROGEN", n_angles),
        ] {
            for t in sections.ints_n(flag, 4 * n)?.chunks_exact(4) {
                angles.push(PrmtopAngle {
                    atoms: (atom_i(t[0]), atom_i(t[1]), atom_i(t[2])),
                    k: param(&angle_k, t[3], "Angle")?,
                    theta0: param(&angle_theta0, t[3], "Angle")?,
                });
            }
        }

        let mut dihedrals = Vec::with_capacity(n_dihedrals_h + n_dihedrals);
        for (flag, n) in [
            ("DIHEDRALS_INC_HYDROGEN", n_dihedrals_h),
            ("DIHEDRALS_WITHOUT_HYDROGEN", n_dihedrals),
        ] {
            for t in sections.ints_n(flag, 5 * n)?.chunks_exact(5) {
                dihedrals.push(PrmtopDihedral {
                    atoms: [atom_i(t[0]), atom_i(t[1]), atom_i(t[2]), atom_i(t[3])],
                    k: param(&dihe_k, t[4], "Dihedral")?,
                    periodicity: param(&dihe_n, t[4], "Dihedral")?,
                    phase: param(&dihe_phase, t[4], "Dihedral")?,
                    scee: param(&scee, t[4], "SCEE")?,
                    scnb: param(&scnb, t[4], "SCNB")?,
                    improper: t[3] < 0,
                    skip_14: t[2] < 0,
                });
            }
        }

        // Per-atom exclusion counts, then a flat list of 1-based partner indices. Atoms without
        // exclusions have a single placeholder 0 entry.
        let n_excluded = sections.ints_n("NUMBER_EXCLUDED_ATOMS", n_atoms)?;
        let excluded_list = sections.ints("EXCLUDED_ATOMS_LIST")?;
        let mut excluded_pairs = Vec::new();
        let mut offset = 0;
        for (i, count) in n_excluded.iter().enumerate() {
            let count = (*count).max(0) as usize;
            let partners = excluded_list
                .get(offset..offset + count)
                .ok_or_else(|| invalid("EXCLUDED_ATOMS_LIST is too short"))?;
            for &j in partners {
                if j > 0 {
                    let j = (j - 1) as usize;
                    excluded_pairs.push((i.min(j), i.max(j)));
                }
            }
            offset += count;
        }

        let box_dims = if has_box && sections.has("BOX_DIMENSIONS") {
            let b = sections.floats_n("BOX_DIMENSIONS", 4)?;
            Some(PrmtopBox {
                beta: b[0] as f32,
                lengths: [b[1] as f32, b[2] as f32, b[3] as f32],
            })
        } else {
            None
        };

        // Terms that change the energy function beyond what the fields above describe.
        let mut unsupported_terms = Vec::new();
        for (flag, label) in [
            ("CHARMM_UREY_BRADLEY_COUNT", "CHARMM Urey-Bradley terms"),
            ("CHARMM_NUM_IMPROPERS", "CHARMM harmonic impropers"),
            ("CMAP_COUNT", "CMAP backbone corrections"),
            ("CHARMM_CMAP_COUNT", "CMAP backbone corrections"),
        ] {
            if sections.has(flag) && sections.ints(flag)?.first().is_some_and(|n| *n > 0) {
                unsupported_terms.push(label.to_owned());
            }
        }
        if sections.has("LENNARD_JONES_14_ACOEF") {
            unsupported_terms.push("Separate 1-4 LJ parameters (CHARMM)".to_owned());
        }
        if sections.has("IPOL") && sections.ints("IPOL")?.first().is_some_and(|n| *n > 0) {
            unsupported_terms.push("Polarizability".to_owned());
        }
        let hbond_pairs = nonbonded_parm_index.iter().any(|i| *i < 0);
        if hbond_pairs {
            let a = sections.floats("HBOND_ACOEF").unwrap_or_default();
            let b = sections.floats("HBOND_BCOEF").unwrap_or_default();
            if a.iter().chain(&b).any(|v| *v != 0.) {
                unsupported_terms.push("10-12 H-bond terms".to_owned());
            }
        }
        unsupported_terms.dedup();

        Ok(Self {
            title: sections.strings("TITLE").unwrap_or_default().join(""),
            atom_names,
            atom_types,
            charges,
            masses,
            atomic_numbers,
            lj_type_index,
            n_lj_types,
            nonbonded_parm_index,
            lj_acoef,
            lj_bcoef,
            residue_labels,
            residue_starts,
            bonds,
            angles,
            dihedrals,
            excluded_pairs,
            box_dims,
            num_extra_points,
            unsupported_terms,
        })
    }

    pub fn n_atoms(&self) -> usize {
        self.atom_names.len()
    }

    /// The LJ A and B coefficients between two atoms. (A/r¹² − B/r⁶). None for 10-12 H-bond
    /// pairs.
    pub fn lj_ab(&self, atom_0: usize, atom_1: usize) -> Option<(f64, f64)> {
        let t0 = self.lj_type_index[atom_0];
        let t1 = self.lj_type_index[atom_1];
        let i = self.nonbonded_parm_index[t0 * self.n_lj_types + t1];
        if i <= 0 {
            return None;
        }
        let i = (i - 1) as usize;
        Some((self.lj_acoef[i], self.lj_bcoef[i]))
    }

    /// The residue index of each atom.
    pub fn atom_residues(&self) -> Vec<usize> {
        let mut result = vec![0; self.n_atoms()];
        for (res_i, start) in self.residue_starts.iter().enumerate() {
            let end = self
                .residue_starts
                .get(res_i + 1)
                .copied()
                .unwrap_or(self.n_atoms());
            for r in result.iter_mut().take(end).skip(*start) {
                *r = res_i;
            }
        }
        result
    }
}

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("prmtop: {msg}"))
}

/// Raw `%FLAG` sections, split into fields using each section's Fortran `%FORMAT`.
struct PrmtopSections {
    sections: HashMap<String, Vec<String>>,
}

impl PrmtopSections {
    fn parse(text: &str) -> io::Result<Self> {
        let mut sections: HashMap<String, Vec<String>> = HashMap::new();
        let mut current: Option<(String, usize, usize)> = None; // (flag, per line, width)

        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("%FLAG") {
                let flag = rest.trim().to_owned();
                sections.entry(flag.clone()).or_default();
                current = Some((flag, 0, 0));
            } else if let Some(rest) = line.strip_prefix("%FORMAT") {
                let Some((flag, _, _)) = current.take() else {
                    continue;
                };
                let (per_line, width) = parse_fortran_format(rest)?;
                current = Some((flag, per_line, width));
            } else if line.starts_with('%') {
                // E.g. %VERSION, or %COMMENT
                continue;
            } else if let Some((flag, per_line, width)) = &current {
                if *width == 0 {
                    return Err(invalid(&format!("Section {flag} has no %FORMAT")));
                }
                let fields = sections.get_mut(flag).unwrap();
                let chars: Vec<char> = line.chars().collect();
                for chunk in chars.chunks(*width).take(*per_line) {
                    let field: String = chunk.iter().collect();
                    // Fixed-width numeric fields may be blank-padded; skip all-blank ones.
                    if !field.trim().is_empty() || !flag_is_numeric(flag) {
                        fields.push(field);
                    }
                }
            }
        }

        Ok(Self { sections })
    }

    fn has(&self, flag: &str) -> bool {
        self.sections.contains_key(flag)
    }

    fn get(&self, flag: &str) -> io::Result<&Vec<String>> {
        self.sections
            .get(flag)
            .ok_or_else(|| invalid(&format!("Missing section {flag}")))
    }

    fn strings(&self, flag: &str) -> io::Result<Vec<String>> {
        Ok(self
            .get(flag)?
            .iter()
            .map(|s| s.trim().to_owned())
            .collect())
    }

    fn strings_n(&self, flag: &str, n: usize) -> io::Result<Vec<String>> {
        let v = self.strings(flag)?;
        check_len(flag, v, n)
    }

    fn ints(&self, flag: &str) -> io::Result<Vec<i32>> {
        self.get(flag)?
            .iter()
            .map(|s| {
                s.trim()
                    .parse::<i32>()
                    .map_err(|_| invalid(&format!("Invalid integer {s:?} in {flag}")))
            })
            .collect()
    }

    fn ints_n(&self, flag: &str, n: usize) -> io::Result<Vec<i32>> {
        let v = self.ints(flag)?;
        check_len(flag, v, n)
    }

    fn floats(&self, flag: &str) -> io::Result<Vec<f64>> {
        self.get(flag)?
            .iter()
            .map(|s| {
                s.trim()
                    .replace(['D', 'd'], "E")
                    .parse::<f64>()
                    .map_err(|_| invalid(&format!("Invalid number {s:?} in {flag}")))
            })
            .collect()
    }

    fn floats_n(&self, flag: &str, n: usize) -> io::Result<Vec<f64>> {
        let v = self.floats(flag)?;
        check_len(flag, v, n)
    }
}

/// Numeric sections are parsed as numbers; others (names, labels, the title) are text. We use
/// this to decide whether blank fixed-width fields are meaningful.
fn flag_is_numeric(flag: &str) -> bool {
    !matches!(
        flag,
        "TITLE"
            | "CTITLE"
            | "ATOM_NAME"
            | "AMBER_ATOM_TYPE"
            | "RESIDUE_LABEL"
            | "TREE_CHAIN_CLASSIFICATION"
            | "RADIUS_SET"
            | "FORCE_FIELD_TYPE"
    )
}

fn check_len<T>(flag: &str, mut v: Vec<T>, n: usize) -> io::Result<Vec<T>> {
    if v.len() < n {
        return Err(invalid(&format!(
            "Section {flag} has {} entries; expected {n}",
            v.len()
        )));
    }
    v.truncate(n);
    Ok(v)
}

/// Parses e.g. "(10I8)", "(5E16.8)", or "(20a4)" into (fields per line, field width).
fn parse_fortran_format(spec: &str) -> io::Result<(usize, usize)> {
    let spec = spec.trim().trim_start_matches('(').trim_end_matches(')');
    let letter_i = spec
        .find(|c: char| c.is_ascii_alphabetic())
        .ok_or_else(|| invalid(&format!("Invalid %FORMAT {spec:?}")))?;
    let per_line = if letter_i == 0 {
        1
    } else {
        spec[..letter_i]
            .parse()
            .map_err(|_| invalid(&format!("Invalid %FORMAT {spec:?}")))?
    };
    let width_str: String = spec[letter_i + 1..]
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    let width = width_str
        .parse()
        .map_err(|_| invalid(&format!("Invalid %FORMAT {spec:?}")))?;
    Ok((per_line, width))
}
