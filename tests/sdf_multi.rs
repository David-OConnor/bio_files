use std::fmt::Write as _;

use bio_files::{BondType, Sdf};

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

#[test]
fn parses_short_header_touching_coordinates_query_bonds_and_valid_elements() {
    let text = "Mrv1652309142106383D

  3  2  0  0  0  0            999 V2000
 1234.5678-1234.5678    0.0000 Gd  0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 Nb  0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 Sm  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  6  0  0  0  0
  2  3  1  0  0  0  0
M  END
";

    let mol = Sdf::new(text).unwrap();
    assert_eq!(mol.atoms.len(), 3);
    assert_eq!(mol.atoms[0].posit.x, 1234.5678);
    assert_eq!(mol.atoms[0].posit.y, -1234.5678);
    assert!(mol.atoms.iter().all(|atom| atom.element.to_letter() == "X"));
    assert_eq!(mol.bonds[0].bond_type, BondType::Unknown);
}

#[test]
fn parses_touching_fixed_width_count_fields() {
    let mut text = String::from("joined counts\n\n\n 96101  0  0  0  0            999 V2000\n");
    for _ in 0..96 {
        writeln!(
            text,
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
        )
        .unwrap();
    }
    for _ in 0..101 {
        writeln!(text, "  1  2  1  0  0  0  0").unwrap();
    }
    writeln!(text, "M  END").unwrap();

    let mol = Sdf::new(&text).unwrap();
    assert_eq!(mol.atoms.len(), 96);
    assert_eq!(mol.bonds.len(), 101);
}

#[test]
fn parses_touching_fixed_width_bond_indices() {
    let mut text =
        String::from("joined bond indices\n\n\n114  1  0  0  0  0            999 V2000\n");
    for _ in 0..114 {
        writeln!(
            text,
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
        )
        .unwrap();
    }
    writeln!(text, " 96114  1  0  0  0  0").unwrap();
    writeln!(text, "M  END").unwrap();

    let mol = Sdf::new(&text).unwrap();
    assert_eq!(mol.bonds[0].atom_0_sn, 96);
    assert_eq!(mol.bonds[0].atom_1_sn, 114);
}
