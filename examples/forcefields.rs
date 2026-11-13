//! This example demonstrates how to open and access Amber-style force field parameters.

use std::path::Path;

use bio_files::{MmCif, Mol2, md_params::ForceFieldParams};

fn main() {
    // General force field parameters:
    let _params = ForceFieldParams::load_dat(Path::new("gaff2.dat")).unwrap();

    let _protein = MmCif::load(Path::new("1c8k.cif")).unwrap();
    let _mol = Mol2::load(Path::new("CPB.mol2")).unwrap();
    let _mol_specific = ForceFieldParams::load_frcmod(Path::new("CPB.frcmod")).unwrap();

    // Or, instead of loading atoms and mol-specific params separately:
    // let (mol, lig_specific) = load_prmtop("my_mol.prmtop");

    // Or, if you have a small molecule available in Amber Geostd, load it remotely:
    // let data = bio_apis::amber_geostd::load_mol_files("CPB");
    // let mol = Mol2::new(&data.mol2);
    // let mol_specific = ForceFieldParams::from_frcmod(&data.frcmod);

    // RCSB PDB proteins usually lack hydrogens, force field types, partial charges, and bonds.
    // The `dynamics` crate adds them:
    // let param_set = dynamics::FfParamSet::new_amber().unwrap();
    // dynamics::prepare_peptide(&mut atoms, &mut bonds, &mut residues, &mut chains, ...);
}
