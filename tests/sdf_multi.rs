use bio_files::Sdf;

const TWO_MOLS: &str = "\
water


  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0
    0.9500    0.0000    0.0000 H   0  0
   -0.2400    0.9200    0.0000 H   0  0
  1  2  1  0
  1  3  1  0
M  END
> <NAME>
Water

$$$$
dioxygen


  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0
    1.2100    0.0000    0.0000 O   0  0
  1  2  2  0
M  END
> <NAME>
Dioxygen

$$$$
";

#[test]
fn parses_every_record() {
    let mols = Sdf::new_multi(TWO_MOLS).unwrap();
    assert_eq!(mols.len(), 2);

    assert_eq!(mols[0].ident, "water");
    assert_eq!(mols[0].atoms.len(), 3);
    assert_eq!(mols[0].bonds.len(), 2);
    assert_eq!(mols[0].metadata.get("NAME").unwrap(), "Water");

    assert_eq!(mols[1].ident, "dioxygen");
    assert_eq!(mols[1].atoms.len(), 2);
    assert_eq!(mols[1].bonds.len(), 1);
    assert_eq!(mols[1].metadata.get("NAME").unwrap(), "Dioxygen");
}

#[test]
fn final_record_without_terminator() {
    let text = TWO_MOLS.trim_end().trim_end_matches("$$$$");
    let mols = Sdf::new_multi(text).unwrap();
    assert_eq!(mols.len(), 2);
    assert_eq!(mols[1].ident, "dioxygen");
}

#[test]
fn single_mol_file_still_works() {
    let one = TWO_MOLS.split("$$$$").next().unwrap();
    let mols = Sdf::new_multi(one).unwrap();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].ident, "water");
}

#[test]
fn bad_record_is_skipped_not_fatal() {
    let text = format!("{TWO_MOLS}garbage\nnot an sdf\n$$$$\n");
    let mols = Sdf::new_multi(&text).unwrap();
    assert_eq!(mols.len(), 2);
}

#[test]
fn round_trips_through_save_multi() {
    let mols = Sdf::new_multi(TWO_MOLS).unwrap();
    let path = std::env::temp_dir().join("bio_files_sdf_multi_tmp.sdf");
    Sdf::save_multi(&mols, &path, Default::default()).unwrap();

    let reloaded = Sdf::load_multi(&path).unwrap();
    assert_eq!(reloaded.len(), 2);
    assert_eq!(reloaded[0].ident, "water");
    assert_eq!(reloaded[0].atoms.len(), 3);
    assert_eq!(reloaded[1].ident, "dioxygen");
    assert_eq!(reloaded[1].bonds.len(), 1);

    let _ = std::fs::remove_file(&path);
}
