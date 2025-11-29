//! We use this to create molecule data (atoms and bonds) from template files for
//! building-block molecules like lipids, nucleic acids, and amino acids. This can be used to make
//! create things like DNA and RNA strands.
//!
//! This assumes Amber format template libraries like `ff-nucleic-OL24.lib`, `RNA.lib`, `amino19.lib`,
//! and `lipid21.lib`.

use std::{collections::HashMap, io};

use lin_alg::f64::Vec3;
use na_seq::Element;

use crate::{AtomGeneric, BondGeneric, BondType};

#[derive(Clone)]
struct AtomRow {
    /// per-residue atom name (eg. "C116")
    name: String,
    /// Amber/GAFF atom type (eg. "cD")
    ff_type: String,
    /// atomic number
    z: u32,
    /// partial charge
    q: f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct UnitConnect {
    pub head: Option<u32>,
    pub tail: Option<u32>,
}

pub type ResidueConnectRow = [u32; 6];

#[derive(Default)]
struct Work {
    atoms: Vec<AtomRow>,
    positions: Vec<Vec3>,
    bonds: Vec<(u32, u32, u32)>, // (a1, a2, flag)

    connect_vals: Vec<u32>,
    connect_present: bool,

    residue_connect: Vec<ResidueConnectRow>,
}

#[derive(Clone, Copy, PartialEq)]
enum Section {
    Atoms,
    Positions,
    Connectivity,
    ConnectArray,
    ResidueConnect,
}

#[derive(Clone, Debug)]
pub struct TemplateData {
    pub atoms: Vec<AtomGeneric>,
    pub bonds: Vec<BondGeneric>,
    pub unit_connect: UnitConnect,
    pub res_connect: Vec<ResidueConnectRow>,
}

impl TemplateData {
    /// Find the indices of a template that join it to previous ones. None for terminal templates.
    pub fn attach_points(&self) -> io::Result<(Option<usize>, Option<usize>)> {
        let head_i = match self.unit_connect.head {
            Some(v) => Some(v as usize - 1),
            None => None,
        };
        let tail_i = match self.unit_connect.tail {
            Some(v) => Some(v as usize - 1),
            None => None,
        };

        Ok((head_i, tail_i))
    }

    pub fn find_atom_i_by_name(&self, name: &str) -> Option<usize> {
        self.atoms
            .iter()
            .enumerate()
            .find(|(_, a)| a.type_in_res_general.as_deref() == Some(name))
            .map(|(i, _)| i)
    }
    pub fn find_atom_by_name(&self, name: &str) -> Option<&AtomGeneric> {
        self.atoms
            .iter()
            .find(|a| a.type_in_res_general.as_deref() == Some(name))
    }
}

/// Creates a set of  atom and bonds for all items in a `.lib` template file, e.g. `lipid21.lib`,
/// or `amino19.lib` from Amber.
/// The hashmap key is the lipid name, e.g. "AR", "CHL" etc.
pub fn load_templates(template_text: &str) -> io::Result<HashMap<String, TemplateData>> {
    let mut out: HashMap<String, TemplateData> = HashMap::new();

    let mut cur_key: Option<String> = None;
    let mut cur_sec: Option<Section> = None;
    let mut cur: Work = Work::default();

    let element_from_z = |z: u32| -> Element {
        match z {
            1 => Element::Hydrogen,
            6 => Element::Carbon,
            7 => Element::Nitrogen,
            8 => Element::Oxygen,
            9 => Element::Fluorine,
            15 => Element::Phosphorus,
            16 => Element::Sulfur,
            17 => Element::Chlorine,
            35 => Element::Bromine,
            53 => Element::Iodine,
            _ => Element::Tellurium,
        }
    };

    let bond_from_flag = |f: u32| -> BondType {
        match f {
            1 => BondType::Single,
            2 => BondType::Double,
            3 => BondType::Triple,
            _ => BondType::Single,
        }
    };

    let mut finalize = |key: Option<String>, work: &mut Work| {
        if let Some(k) = key {
            if !work.atoms.is_empty() {
                let n = work.atoms.len();

                if work.positions.len() < n {
                    work.positions.extend(
                        std::iter::repeat_with(|| Vec3::new(0.0, 0.0, 0.0))
                            .take(n - work.positions.len()),
                    );
                }

                let atoms: Vec<AtomGeneric> = work
                    .atoms
                    .iter()
                    .enumerate()
                    .map(|(i, ar)| AtomGeneric {
                        serial_number: (i as u32) + 1,
                        posit: work.positions[i],
                        element: element_from_z(ar.z),
                        type_in_res: None,
                        type_in_res_general: Some(ar.name.clone()),
                        force_field_type: Some(ar.ff_type.clone()),
                        partial_charge: Some(ar.q as f32),
                        hetero: false,
                        occupancy: None,
                        alt_conformation_id: None,
                    })
                    .collect();

                let bonds: Vec<BondGeneric> = work
                    .bonds
                    .iter()
                    .map(|&(a1, a2, fl)| BondGeneric {
                        bond_type: bond_from_flag(fl),
                        atom_0_sn: a1,
                        atom_1_sn: a2,
                    })
                    .collect();

                let unit_connect = if work.connect_present {
                    let head = work.connect_vals.get(0).copied().unwrap_or(0);
                    let tail = work.connect_vals.get(1).copied().unwrap_or(0);
                    UnitConnect {
                        head: if head == 0 { None } else { Some(head) },
                        tail: if tail == 0 { None } else { Some(tail) },
                    }
                } else {
                    UnitConnect {
                        head: None,
                        tail: None,
                    }
                };

                let res_connect = work.residue_connect.clone();
                out.insert(
                    k,
                    TemplateData {
                        atoms,
                        bonds,
                        unit_connect,
                        res_connect,
                    },
                );
            }
            *work = Work::default();
        }
    };

    for line in template_text.lines() {
        let l = line.trim_start();

        if let Some(rest) = l.strip_prefix("!entry.") {
            if let Some(dot_unit_idx) = rest.find(".unit.") {
                let key = &rest[..dot_unit_idx];
                let after_unit = &rest[dot_unit_idx + ".unit.".len()..];

                let key_changed = match &cur_key {
                    Some(k) => k != key,
                    None => true,
                };
                if key_changed {
                    finalize(cur_key.take(), &mut cur);
                    cur_key = Some(key.to_string());
                }

                if after_unit.starts_with("atoms table") {
                    cur_sec = Some(Section::Atoms);
                    continue;
                } else if after_unit.starts_with("positions table") {
                    cur_sec = Some(Section::Positions);
                    continue;
                } else if after_unit.starts_with("connectivity table") {
                    cur_sec = Some(Section::Connectivity);
                    continue;
                } else if after_unit.starts_with("connect array") {
                    cur_sec = Some(Section::ConnectArray);
                    cur.connect_present = true;
                    cur.connect_vals.clear();
                    continue;
                } else if after_unit.starts_with("residueconnect table") {
                    cur_sec = Some(Section::ResidueConnect);
                    continue;
                } else {
                    cur_sec = None;
                    continue;
                }
            } else {
                cur_sec = None;
                continue;
            }
        }

        match cur_sec {
            Some(Section::Atoms) => {
                let bytes = l.as_bytes();
                let mut qpos = Vec::with_capacity(4);
                for (i, &b) in bytes.iter().enumerate() {
                    if b == b'"' {
                        qpos.push(i);
                    }
                    if qpos.len() == 4 {
                        break;
                    }
                }

                if qpos.len() == 4 {
                    let name = &l[qpos[0] + 1..qpos[1]];
                    let ff_type = &l[qpos[2] + 1..qpos[3]];
                    let tail = &l[qpos[3] + 1..];

                    let nums: Vec<&str> = tail.split_whitespace().collect();
                    if nums.len() >= 6 {
                        let elmnt = nums[4].parse::<u32>().unwrap_or(0);
                        let chg = nums[5].parse::<f64>().unwrap_or(0.0);
                        cur.atoms.push(AtomRow {
                            name: name.to_string(),
                            ff_type: ff_type.to_string(),
                            z: elmnt,
                            q: chg,
                        });
                    }
                }
            }
            Some(Section::Positions) => {
                let mut it = l.split_whitespace();
                if let (Some(x), Some(y), Some(z)) = (it.next(), it.next(), it.next())
                    && let (Ok(xv), Ok(yv), Ok(zv)) =
                        (x.parse::<f64>(), y.parse::<f64>(), z.parse::<f64>())
                {
                    cur.positions.push(Vec3::new(xv, yv, zv));
                }
            }
            Some(Section::Connectivity) => {
                let mut it = l.split_whitespace();
                if let (Some(a1), Some(a2), Some(flg)) = (it.next(), it.next(), it.next())
                    && let (Ok(a1v), Ok(a2v), Ok(fv)) =
                        (a1.parse::<u32>(), a2.parse::<u32>(), flg.parse::<u32>())
                {
                    cur.bonds.push((a1v, a2v, fv));
                }
            }
            Some(Section::ConnectArray) => {
                if let Some(tok) = l.split_whitespace().next()
                    && let Ok(v) = tok.parse::<u32>()
                {
                    cur.connect_vals.push(v);
                }
            }
            Some(Section::ResidueConnect) => {
                let vals: Vec<u32> = l
                    .split_whitespace()
                    .filter_map(|t| t.parse::<u32>().ok())
                    .collect();
                if vals.len() >= 6 {
                    cur.residue_connect
                        .push([vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]]);
                }
            }
            None => {}
        }
    }

    finalize(cur_key.take(), &mut cur);
    Ok(out)
}
