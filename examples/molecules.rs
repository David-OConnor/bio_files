//! This example demonstrates how to open and close molecules. For example, proteins and
//! small organic molecules

use std::path::Path;

use bio_files::{MmCif, Mol2, Sdf, Xyz};

fn main() {
    let sdf = Sdf::load(Path::new("./molecules/DB03496.sdf")).unwrap();

    // sdf_data.atoms[0]; // (as above)
    sdf.atoms[0].posit; // (as above, but lin_alg::Vec3))
    sdf.save(Path::new("test.sdf")).unwrap();

    let mut mol2: Mol2 = sdf.into();
    mol2.ident = String::from("ABC");
    mol2.save(Path::new("test.mol2")).unwrap();

    let xyz_data = Xyz::load(Path::new("./atom_posits.xyz")).unwrap();

    // Load molecules from databases using identifiers:
    let mol = Sdf::load_drugbank("DB00198").unwrap();
    let mol = Sdf::load_pubchem(12345).unwrap();
    let mol = Sdf::load_pdbe("CPB").unwrap();
    let mol = Mol2::load_amber_geostd("CPB").unwrap();

    let peptide = MmCif::load_rcsb("8S6P").unwrap();
}
