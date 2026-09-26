use bio_files::{MmCif, SecondaryStructure};
use na_seq::AtomTypeInRes;

/// Backbone atoms of residues 28-80 and 325-351 of a Boltz-2 model of TdT. Like most structure
/// prediction output, it has no secondary structure records. It includes a parallel β-sheet with
/// a β-bulge, an α-helix, an antiparallel β-hairpin, and a gap in the chain.
const EXCERPT: &str = include_str!("data/boltz2_tdt_excerpt.cif");

/// Per-residue secondary structure, from mkdssp 4.4.5. (H, G, and I are shown as H.)
const EXPECTED: &str =
    "----EEEEEEEE------HHHHHHHHHHHHH--EEE---------EEE-------EEEE-HHHH----------EEEE--";

#[test]
fn infers_missing_secondary_structure() {
    let m = MmCif::new(EXCERPT).unwrap();

    let ss: String = m
        .atoms
        .iter()
        .filter(|a| a.type_in_res == Some(AtomTypeInRes::CA))
        .map(|a| {
            let seg = m
                .secondary_structure
                .iter()
                .find(|s| (s.start_sn..=s.end_sn).contains(&a.serial_number));
            match seg.map(|s| s.sec_struct) {
                Some(SecondaryStructure::Helix) => 'H',
                Some(SecondaryStructure::Sheet) => 'E',
                _ => '-',
            }
        })
        .collect();

    assert_eq!(ss, EXPECTED);
}
