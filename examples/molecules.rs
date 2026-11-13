//! This example demonstrates how to open and close molecules. For example, proteins and
//! small organic molecules

use std::path::Path;

use bio_apis::pubchem::StructureSearchNamespace;
use bio_files::{MmCif, Mol2, Sdf, Xyz, sdf::SdfFormat};

fn main() {
    let sdf = Sdf::load(Path::new("./molecules/DB03496.sdf")).unwrap();

    let _posit = sdf.atoms[0].posit; // lin_alg::f64::Vec3
    sdf.save(Path::new("test.sdf"), SdfFormat::V2000).unwrap();

    let mut mol2: Mol2 = sdf.into();
    mol2.ident = String::from("ABC");
    mol2.save(Path::new("test.mol2")).unwrap();

    let _xyz_data = Xyz::load(Path::new("./atom_posits.xyz")).unwrap();

    // Load molecules from databases using identifiers:
    let _mol = Sdf::load_drugbank("DB00198").unwrap();
    let _mol = Sdf::load_pubchem(StructureSearchNamespace::Cid, "12345").unwrap();
    let _mol = Sdf::load_pdbe("CPB").unwrap();
    let _mol = Mol2::load_amber_geostd("CPB").unwrap();

    let _peptide = MmCif::load_rcsb("8S6P").unwrap();
}
