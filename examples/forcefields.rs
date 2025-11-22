//! This example demonstrates how to open and access Amber-style force field parameters.

use std::path::Path;

use bio_files::md_params::ForceFieldParams;
use bio_files::{MmCif, Mol2};

use dynamics::{FfParamSet, prepare_peptide};

fn load() {
    let param_set = FfParamSet::new_amber().unwrap();

    let mut protein = MmCif::load(Path::new("1c8k.cif")).unwrap();
    let mol = Mol2::load(Path::new("CPB.mol2")).unwrap();
    let mol_specific = ForceFieldParams::load_frcmod(Path::new("CPB.frcmod")).unwrap();

    // Or, instead of loading atoms and mol-specific params separately:
    // let (mol, lig_specific) = load_prmtop("my_mol.prmtop");

    // Or, if you have a small molecule available in Amber Geostd, load it remotely:
    // let data = bio_apis::amber_geostd::load_mol_files("CPB");
    // let mol = Mol2::new(&data.mol2);
    // let mol_specific = ForceFieldParams::from_frcmod(&data.frcmod);

    // Add Hydrogens, force field type, and partial charge to atoms in the protein; these usually aren't
    // included from RSCB PDB. You can also call `populate_hydrogens_dihedrals()`, and
    // `populate_peptide_ff_and_q() separately. Add bonds.
    prepare_peptide(
        &mut protein.atoms,
        &mut protein.bonds,
        &mut protein.residues,
        &mut protein.chains,
        &param_set.peptide_ff_q_map.as_ref().unwrap(),
        7.0,
    )
        .unwrap();
}

fn main() {
    // Loading Force field parameters:
    let p = Path::new("gaff2.dat");
    let params = ForceFieldParams::load_dat(p).unwrap();
}
