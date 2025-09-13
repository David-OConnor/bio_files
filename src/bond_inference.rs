//! This module creates bonds between protein components. Most macromolecule PDB/CIF files don't include
//! explicit bond information. Infer bond lengths
//! by comparing each interactomic bond distance, and matching against known amino acid bond lengths.
//!
//! Some info here: https://www.ruppweb.org/Xray/tutorial/protein_structure.htm
//! https://itp.uni-frankfurt.de/~engel/amino.html
//!
//! All lengths are in angstrom (Å)

use std::{collections::HashMap, f64::consts::TAU};

use lin_alg::f64::Vec3;
use na_seq::{
    Element,
    Element::{Carbon, Fluorine, Hydrogen, Nitrogen, Oxygen, Sulfur},
};
use rayon::prelude::*;

use crate::{AtomGeneric, BondGeneric, BondType};
// use crate::{molecule::{Atom, Bond, HydrogenBond}, util::{find_atom, setup_neighbor_pairs}, AtomGeneric, BondGeneric, BondType};

struct BondSpecs {
    len: f64,
    elements: (Element, Element),
    bond_type: BondType,
}

impl BondSpecs {
    pub fn new(len: f64, elements: (Element, Element), bond_type: BondType) -> Self {
        Self {
            len,
            elements,
            bond_type,
        }
    }
}

// If interatomic distance is within this distance of one of our known bond lenghts, consider it to be a bond.
// Relevant to this is both bond variability under various conditions, and measurement precision.
const COV_BOND_LEN_THRESH: f64 = 0.04; // todo: Adjust A/R based on performannce.
const COV_DIST_GRID: f64 = 1.6; // Slightly larger than the largest bond distance + thresh.

#[rustfmt::skip]
fn get_specs() -> Vec<BondSpecs> {
    // Code shorteners
    let single = BondType::Single;
    let amide = BondType::Aromatic;
    let double = BondType::Double;
    let triple = BondType::Triple;
    // todo: Can we identify other bond types, mainly Amide?

    vec![
        // --------------------
        // Carbon–Carbon Bonds
        // --------------------

        // C–C single bond
        // The most frequently encountered bond length for saturated, sp³-hybridized carbons (e.g., in alkanes).
        BondSpecs::new(1.54, (Carbon, Carbon), single),

        // Cα–C′: ~1.50 - 1.52 Å
        BondSpecs::new(1.51, (Carbon, Carbon), single),

        // C–C sp²–sp³ single bond, e.g. connecting Phe's ring to the rest of the atom.
        BondSpecs::new(1.50, (Carbon, Carbon), amide),

        // Workaround for Phe's ring in some cases.
        BondSpecs::new(1.47, (Carbon, Carbon), amide),
        BondSpecs::new(1.44, (Carbon, Carbon), amide),
        BondSpecs::new(1.41, (Carbon, Carbon), amide),

        // C-C phenyl (aromatic) ring bond, or benzene ring.
        // Found in alkynes, where carbons are sp-hybridized (linear). ~1.37-1.40 Å
        BondSpecs::new(1.39, (Carbon, Carbon), amide),

        // C-C Seems to be required for one fo the Trp rings?
        BondSpecs::new(1.36, (Carbon, Carbon), amide),

        // C=C double bond
        // Common in alkenes (sp²-hybridized). Range: ~1.33–1.34 Å
        BondSpecs::new(1.33, (Carbon, Carbon), double),

        // C≡C triple bond
        // Found in alkynes, where carbons are sp-hybridized (linear). ~1.20 Å
        BondSpecs::new(1.20, (Carbon, Carbon), triple),

        // --------------------
        // Carbon–Nitrogen Bonds
        // --------------------

        // C–N single bond
        // Typical for amines or alkyl–amine bonds. ~1.45-1.47 Å
        // Also covers Amide Nitrogen to C-alpha bond in protein backbones.
        BondSpecs::new(1.46, (Carbon, Nitrogen), single),

        // C-N Indole N in 5-member aromatic ring, e.g. Trp. 1.36-1.39
        // BondSpecs::new(1.37, (Carbon, Nitrogen), type_hybrid),
        BondSpecs::new(1.37, (Carbon, Nitrogen), single),

        // todo: Some adjustments here may be required regarding single vs hybrid N-C bonds.

        // C-N (amide). Partial double-bond character due to resonance in the amide.
        // BondSpecs::new(1.33, (Carbon, Nitrogen), type_hybrid),
        BondSpecs::new(1.33, (Carbon, Nitrogen), single),

        // C=N double bond
        // Typical for imines (Schiff bases). ~1.28 Å
        BondSpecs::new(1.28, (Carbon, Nitrogen), double),

        // C≡N triple bond
        // Typical of nitriles (–C≡N). ~1.16 Å
        BondSpecs::new(1.16, (Carbon, Nitrogen), triple),
        // NOTE:
        // In proteins, the amide (peptide) bond between C=O and N has partial double-bond character,
        // and the C–N bond length in an amide is around 1.32–1.33 Å.

        // --------------------
        // Carbon–Oxygen Bonds
        // --------------------

        // C–O single bond
        // Found in alcohols, ethers (sp³–O). ~1.43 Å
        BondSpecs::new(1.43, (Carbon, Oxygen), single),

        // C(phenyl)–O. Phenolic C–O bond often shorter than a typical aliphatic C–O. 1.36-1.38 Å
        BondSpecs::new(1.37, (Carbon, Oxygen), single),

        // C′–O (in –COO⁻). 1.25-1.27 Å
        // BondSpecs::new(1.26, (Carbon, Oxygen), type_singl),
        BondSpecs::new(1.26, (Carbon, Oxygen), double),

        // C=O double bond
        // Typical for carbonyl groups (aldehydes, ketones, carboxylic acids, amides). ~1.21–1.23 Å
        BondSpecs::new(1.22, (Carbon, Oxygen), double),

        // --------------------
        // Carbon–Hydrogen Bonds
        // --------------------

        BondSpecs::new(1.09, (Hydrogen, Carbon), single),

        // 1.01–1.02 Å
        BondSpecs::new(1.01, (Hydrogen, Nitrogen), single),

        // 0.96 – 0.98 Å
        BondSpecs::new(1.01, (Hydrogen, Oxygen), single),
        // BondSpecs::new(1.01, (Hydrogen, Oxygen), single),
        BondSpecs::new(0.95, (Hydrogen, Oxygen), single),


        // Non-protein-backbond bond lengths.

        // 1.34 - 1.35. Example: Cys.
        BondSpecs::new(1.34, (Sulfur, Hydrogen), single),

        // 1.81 - 1.82. Example: Cys.
        BondSpecs::new(1.81, (Sulfur, Carbon), single),
    ]
}

/// Infer bonds from atom distances. Uses spacial partitioning for efficiency.
/// We Check pairs only within nearby bins.
pub fn create_bonds(atoms: &[AtomGeneric]) -> Vec<BondGeneric> {
    let specs = get_specs();

    // We use spacial partitioning, so as not to copmare every pair of atoms.
    let posits: Vec<_> = atoms.iter().map(|a| &a.posit).collect();
    // Indices are all values here.
    let indices: Vec<_> = (0..posits.len()).collect();
    let neighbor_pairs = setup_neighbor_pairs(&posits, &indices, COV_DIST_GRID);

    // todo: Should we create an Vec of neighbors for each atom. (Maybe storeed in a hashmap etc)
    // todo, then iterate over that for neighbors in the j loop? WOuld be more generalizable/extract
    // todo it out from the bus logic.

    neighbor_pairs
        .par_iter()
        .filter_map(|(i, j)| {
            let atom_0 = &atoms[*i];
            let atom_1 = &atoms[*j];
            let dist = (atom_0.posit - atom_1.posit).magnitude();

            specs.iter().find_map(|spec| {
                let matches_elements = (atom_0.element == spec.elements.0
                    && atom_1.element == spec.elements.1)
                    || (atom_0.element == spec.elements.1 && atom_1.element == spec.elements.0);

                // If both the element match and distance-threshold check pass,
                // we create a Bond and stop searching any further specs.
                if matches_elements && (dist - spec.len).abs() < COV_BOND_LEN_THRESH {
                    Some(BondGeneric {
                        bond_type: spec.bond_type,
                        atom_0_sn: atom_0.serial_number,
                        atom_1_sn: atom_1.serial_number,
                        // atom_0: *i,
                        // atom_1: *j,
                        // is_backbone: atom_0.is_backbone() && atom_1.is_backbone(),
                    })
                } else {
                    None
                }
            })
        })
        .collect()
}

/// A helper fn. Maps from a global index, to a local atom from a subset.
fn find_atom<'a>(
    atoms: &'a [AtomGeneric],
    indices: &[usize],
    i_to_find: usize,
) -> Option<&'a AtomGeneric> {
    for (i_set, atom) in atoms.iter().enumerate() {
        if indices[i_set] == i_to_find {
            return Some(atom);
        }
    }

    None
}

/// Creates pairs of all *nearby* positions. Much faster than comparing every combination, if only nearly
/// ones are relevant.
/// The separate `indexes` parameter allows `posits` to be a subset of the array we're indexing into,
/// e.g. a filtered set of atoms.
fn setup_neighbor_pairs(
    posits: &[&Vec3],
    indexes: &[usize],
    grid_size: f64,
) -> Vec<(usize, usize)> {
    // Build a spatial grid for atom indices.
    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();

    for (i, posit) in posits.iter().enumerate() {
        let grid_pos = (
            (posit.x / grid_size).floor() as i32,
            (posit.y / grid_size).floor() as i32,
            (posit.z / grid_size).floor() as i32,
        );

        grid.entry(grid_pos).or_default().push(indexes[i]);
    }

    // Collect candidate atom pairs based on neighboring grid cells.
    let mut result = Vec::new();
    for (&cell, indices) in &grid {
        // Look at this cell and its neighbors.
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor_cell = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                    if let Some(neighbor_indices) = grid.get(&neighbor_cell) {
                        // Attempt to prevent duplicates as we iterate. Note working.
                        for &i in indices {
                            for &j in neighbor_indices {
                                // The ordering prevents duplicates.
                                if i < j {
                                    result.push((i, j));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    result
}
