//! Infers protein secondary structure (helices and β-strands) from backbone geometry, using the
//! DSSP algorithm: [Kabsch & Sander, 1983](https://doi.org/10.1002/bip.360221211). We follow the
//! reference implementation (DSSP 2) in the details.
//!
//! This is for structures that don't include secondary structure. For example, mmCIF files from
//! structure prediction tools like Boltz, Chai, and AlphaFold 3.

use std::collections::HashMap;

use lin_alg::f64::Vec3;
use na_seq::{AminoAcid, AtomTypeInRes};

use crate::{
    AtomGeneric, BackboneSS, ChainGeneric, ResidueGeneric, ResidueType, SecondaryStructure,
};

/// kcal/mol. Pairs with an energy below this are H-bonded.
const MAX_HBOND_ENERGY: f64 = -0.5;
/// kcal/mol. H-bond energies are clamped to this.
const MIN_HBOND_ENERGY: f64 = -9.9;
/// kcal·Å/mol. 0.084 * 332: The C, O, N, and H partial charges, and a dimensional factor.
const COUPLING_CONST: f64 = -27.888;
/// Å. Atoms closer than this make the strongest H-bond.
const MIN_DIST: f64 = 0.5;
/// Å. We only look for H-bonds between residues whose Cα atoms are closer than this.
const MAX_CA_DIST: f64 = 9.;
/// Å. Residues whose C and the next one's N are further apart than this aren't bonded.
const MAX_PEPTIDE_BOND_LEN: f64 = 2.5;

#[derive(Clone, Copy, PartialEq, Debug)]
enum Dssp {
    Loop,
    AlphaHelix,
    Helix310,
    PiHelix,
    /// Part of a ladder: Consecutive β-bridges.
    Strand,
    /// An isolated β-bridge.
    Bridge,
}

impl Dssp {
    fn sec_struct(self) -> Option<SecondaryStructure> {
        match self {
            Self::AlphaHelix | Self::Helix310 | Self::PiHelix => Some(SecondaryStructure::Helix),
            Self::Strand => Some(SecondaryStructure::Sheet),
            Self::Loop | Self::Bridge => None,
        }
    }
}

struct BackboneRes {
    n: Vec3,
    ca: Vec3,
    c: Vec3,
    o: Vec3,
    /// The amide hydrogen, placed as DSSP does. None for proline, and at the start of a segment.
    h: Option<Vec3>,
    ca_sn: u32,
    chain: usize,
    /// Residues are bonded to their neighbors in the same segment. A segment ends at the end of a
    /// chain, or at a gap in it.
    segment: usize,
    /// The (up to) two lowest-energy H-bonds from this residue's N-H: (Acceptor index, energy).
    /// An energy of 0 is an empty slot.
    acceptors: [(usize, f64); 2],
}

/// A run of β-bridges between residues `i` and `j`, possibly joined across β-bulges.
struct Ladder {
    parallel: bool,
    i_start: usize,
    i_end: usize,
    j_start: usize,
    j_end: usize,
    bridge_count: usize,
}

/// The electrostatic energy of an H-bond from `acc`'s C=O to `don`'s N-H, in kcal/mol.
fn hbond_energy(don: &BackboneRes, acc: &BackboneRes) -> f64 {
    let Some(h) = don.h else {
        return 0.;
    };

    let d_ho = (h - acc.o).magnitude();
    let d_hc = (h - acc.c).magnitude();
    let d_nc = (don.n - acc.c).magnitude();
    let d_no = (don.n - acc.o).magnitude();

    if d_ho < MIN_DIST || d_hc < MIN_DIST || d_nc < MIN_DIST || d_no < MIN_DIST {
        return MIN_HBOND_ENERGY;
    }

    let e = COUPLING_CONST * (1. / d_ho - 1. / d_hc + 1. / d_nc - 1. / d_no);
    ((e * 1_000.).round() / 1_000.).max(MIN_HBOND_ENERGY)
}

/// Collect the amino acid residues with a complete backbone, in order.
fn backbone_residues(
    atoms: &[AtomGeneric],
    residues: &[ResidueGeneric],
    chains: &[ChainGeneric],
) -> Vec<BackboneRes> {
    let atoms_by_sn: HashMap<u32, &AtomGeneric> =
        atoms.iter().map(|a| (a.serial_number, a)).collect();

    let chain_by_sn: HashMap<u32, usize> = chains
        .iter()
        .enumerate()
        .flat_map(|(i, c)| c.atom_sns.iter().map(move |&sn| (sn, i)))
        .collect();

    let mut result: Vec<BackboneRes> = Vec::new();
    let mut segment = 0;

    for res in residues {
        let ResidueType::AminoAcid(aa) = res.res_type else {
            continue;
        };

        // If there are alternate conformations, this is the first.
        let find = |tir: AtomTypeInRes| {
            res.atom_sns
                .iter()
                .filter_map(|sn| atoms_by_sn.get(sn))
                .find(|a| a.type_in_res.as_ref() == Some(&tir))
        };

        // A residue without a complete backbone can't take part, and breaks its segment.
        let (Some(n), Some(ca), Some(c), Some(o)) = (
            find(AtomTypeInRes::N),
            find(AtomTypeInRes::CA),
            find(AtomTypeInRes::C),
            find(AtomTypeInRes::O),
        ) else {
            continue;
        };

        let chain = chain_by_sn
            .get(&ca.serial_number)
            .copied()
            .unwrap_or_default();

        let prev = result
            .last()
            .filter(|p| p.chain == chain && (p.c - n.posit).magnitude() <= MAX_PEPTIDE_BOND_LEN);

        // DSSP places H 1Å from N, opposite the previous residue's C=O.
        let h = match prev {
            Some(p) if aa != AminoAcid::Pro => Some(n.posit + (p.c - p.o).to_normalized()),
            _ => None,
        };

        if prev.is_none() && !result.is_empty() {
            segment += 1;
        }

        result.push(BackboneRes {
            n: n.posit,
            ca: ca.posit,
            c: c.posit,
            o: o.posit,
            h,
            ca_sn: ca.serial_number,
            chain,
            segment,
            acceptors: [(0, 0.); 2],
        });
    }

    result
}

/// Whether there's an H-bond from `co`'s C=O to `nh`'s N-H.
fn hbond(res: &[BackboneRes], co: usize, nh: usize) -> bool {
    res[nh]
        .acceptors
        .iter()
        .any(|&(i, e)| i == co && e < MAX_HBOND_ENERGY)
}

/// Find each residue's two strongest H-bonds from its N-H.
fn calc_hbonds(res: &mut [BackboneRes]) {
    // Bin Cα atoms into cells, so we only compare nearby residues.
    let cell = |p: Vec3| {
        (
            (p.x / MAX_CA_DIST).floor() as i32,
            (p.y / MAX_CA_DIST).floor() as i32,
            (p.z / MAX_CA_DIST).floor() as i32,
        )
    };

    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    for (i, r) in res.iter().enumerate() {
        grid.entry(cell(r.ca)).or_default().push(i);
    }

    let mut pairs = Vec::new();
    for (i, r) in res.iter().enumerate() {
        let (x, y, z) = cell(r.ca);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(cell) = grid.get(&(x + dx, y + dy, z + dz)) else {
                        continue;
                    };
                    for &j in cell {
                        if j > i && (res[j].ca - r.ca).magnitude() < MAX_CA_DIST {
                            pairs.push((i, j));
                        }
                    }
                }
            }
        }
    }
    pairs.sort_unstable();

    let mut record = |don: usize, acc: usize| {
        let e = hbond_energy(&res[don], &res[acc]);
        let slots = &mut res[don].acceptors;

        if e < slots[0].1 {
            slots[1] = slots[0];
            slots[0] = (acc, e);
        } else if e < slots[1].1 {
            slots[1] = (acc, e);
        }
    };

    for (i, j) in pairs {
        record(i, j);
        if j != i + 1 {
            record(j, i);
        }
    }
}

/// Find β-ladders, including those joined by β-bulges.
fn calc_ladders(res: &[BackboneRes]) -> Vec<Ladder> {
    let len = res.len();

    let hb = |co: usize, nh: usize| hbond(res, co, nh);

    // H-bond partners, in either direction. Used to limit which pairs we test for bridges.
    let mut partners = vec![Vec::new(); len];
    for (nh, r) in res.iter().enumerate() {
        for &(co, e) in &r.acceptors {
            if e < MAX_HBOND_ENERGY {
                partners[nh].push(co);
                partners[co].push(nh);
            }
        }
    }

    let same_segment = |i: usize, j: usize| res[i].segment == res[j].segment;

    // Some(true) for a parallel bridge between `i` and `j`, and Some(false) for antiparallel.
    let bridge = |i: usize, j: usize| {
        if !same_segment(i - 1, i + 1) || !same_segment(j - 1, j + 1) {
            return None;
        }

        if (hb(i - 1, j) && hb(j, i + 1)) || (hb(j - 1, i) && hb(i, j + 1)) {
            Some(true)
        } else if (hb(i, j) && hb(j, i)) || (hb(i - 1, j + 1) && hb(j - 1, i + 1)) {
            Some(false)
        } else {
            None
        }
    };

    let mut ladders: Vec<Ladder> = Vec::new();

    for i in 1..len.saturating_sub(4) {
        // Each bridge pattern includes an H-bond from `i - 1`, `i`, or `i + 1` to `j - 1`, `j`, or
        // `j + 1`.
        let mut candidates: Vec<usize> = (i - 1..=i + 1)
            .flat_map(|k| &partners[k])
            .flat_map(|&p| [p.checked_sub(1), Some(p), Some(p + 1)])
            .flatten()
            .filter(|&j| j >= i + 3 && j + 1 < len)
            .collect();
        candidates.sort_unstable();
        candidates.dedup();

        for j in candidates {
            let Some(parallel) = bridge(i, j) else {
                continue;
            };

            let extends = ladders.iter_mut().find(|l| {
                l.parallel == parallel
                    && l.i_end + 1 == i
                    && if parallel {
                        l.j_end + 1 == j
                    } else {
                        l.j_start == j + 1
                    }
            });

            match extends {
                Some(l) => {
                    l.i_end = i;
                    if parallel {
                        l.j_end = j;
                    } else {
                        l.j_start = j;
                    }
                    l.bridge_count += 1;
                }
                None => ladders.push(Ladder {
                    parallel,
                    i_start: i,
                    i_end: i,
                    j_start: j,
                    j_end: j,
                    bridge_count: 1,
                }),
            }
        }
    }

    // Join ladders separated by a β-bulge: A gap of up to 1 residue on one strand, and 4 on the
    // other.
    ladders.sort_by_key(|l| (res[l.i_start].chain, l.i_start));

    // Whether `to` is after `from`, by less than `max`.
    let within = |from: usize, to: usize, max: usize| to >= from && to - from < max;

    let mut a = 0;
    while a < ladders.len() {
        let mut b = a + 1;

        while b < ladders.len() {
            let (l0, l1) = (&ladders[a], &ladders[b]);

            let joinable = l0.parallel == l1.parallel
                && res[l0.i_start.min(l1.i_start)].chain == res[l0.i_end.max(l1.i_end)].chain
                && res[l0.j_start.min(l1.j_start)].chain == res[l0.j_end.max(l1.j_end)].chain
                && within(l0.i_end, l1.i_start, 6)
                && !(l0.i_end >= l1.i_start && l0.i_start <= l1.i_end);

            let bulge = joinable
                && if l0.parallel {
                    (within(l0.j_end, l1.j_start, 6) && within(l0.i_end, l1.i_start, 3))
                        || within(l0.j_end, l1.j_start, 3)
                } else {
                    (within(l1.j_end, l0.j_start, 6) && within(l0.i_end, l1.i_start, 3))
                        || within(l1.j_end, l0.j_start, 3)
                };

            if bulge {
                let l1 = ladders.remove(b);
                let l0 = &mut ladders[a];

                l0.i_end = l1.i_end;
                if l0.parallel {
                    l0.j_end = l1.j_end;
                } else {
                    l0.j_start = l1.j_start;
                }
                l0.bridge_count += l1.bridge_count;
            } else {
                b += 1;
            }
        }

        a += 1;
    }

    ladders
}

/// Assign each residue's DSSP secondary structure.
fn assign(res: &[BackboneRes]) -> Vec<Dssp> {
    let len = res.len();
    let mut result = vec![Dssp::Loop; len];

    // β-strands and bridges
    for l in calc_ladders(res) {
        let ss = if l.bridge_count > 1 {
            Dssp::Strand
        } else {
            Dssp::Bridge
        };

        for k in (l.i_start..=l.i_end).chain(l.j_start..=l.j_end) {
            if result[k] != Dssp::Strand {
                result[k] = ss;
            }
        }
    }

    // An n-turn at `i`: an H-bond from `i`'s C=O to the N-H `n` residues along.
    let turn = |n: usize, i: usize| {
        i + n < len && res[i].segment == res[i + n].segment && hbond(res, i, i + n)
    };

    // Helices: two consecutive n-turns. α-helices take priority over everything; 3₁₀ and
    // π-helices only fill in loops.
    for (n, ss) in [
        (4, Dssp::AlphaHelix),
        (3, Dssp::Helix310),
        (5, Dssp::PiHelix),
    ] {
        for i in 1..len {
            if !(turn(n, i - 1) && turn(n, i)) {
                continue;
            }

            let range = i..i + n;
            if ss == Dssp::AlphaHelix
                || result[range.clone()]
                    .iter()
                    .all(|&s| s == Dssp::Loop || s == ss)
            {
                result[range].fill(ss);
            }
        }
    }

    result
}

/// Infer secondary structure (helices and β-strands) from the backbone geometry of a protein's
/// amino acid residues, using DSSP. Use this when a structure file doesn't include it.
///
/// 3₁₀ and π-helices are reported as helices. Isolated β-bridges, turns, and bends are not
/// reported.
pub fn infer_secondary_structure(
    atoms: &[AtomGeneric],
    residues: &[ResidueGeneric],
    chains: &[ChainGeneric],
) -> Vec<BackboneSS> {
    let mut res = backbone_residues(atoms, residues, chains);
    calc_hbonds(&mut res);
    let ss = assign(&res);

    // Combine consecutive residues into segments.
    let mut result = Vec::new();
    let mut start = 0;

    for i in 1..=res.len() {
        let same = i < res.len()
            && ss[i].sec_struct() == ss[start].sec_struct()
            && res[i].segment == res[start].segment;
        if same {
            continue;
        }

        if let Some(sec_struct) = ss[start].sec_struct() {
            result.push(BackboneSS {
                start_sn: res[start].ca_sn,
                end_sn: res[i - 1].ca_sn,
                sec_struct,
            });
        }
        start = i;
    }

    result
}
