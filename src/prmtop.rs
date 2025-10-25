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
    for i in 0..31 {
        pointers[i] = get_i(&p.data[i])?;
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
