//! For parsing [Amber force field](https://ambermd.org/AmberModels.php) dat files, like
//! `gaff2.dat` (small molecules/ligands). These include parameters used in molecular dynamics.
//! See also: `frcmod`, which patches these parameters for specific molecules.

// todo: For now, this has much overlap with loading frcmod data.
// todo: Reoncosider the API, how this and the frcmod modules are related, and which
// todo: to export the main FF struct from

use std::{
    fs::File,
    io,
    io::{ErrorKind, Read},
    path::Path,
};

use crate::amber_params::{
    AngleBendingParams, BondStretchingParams, DihedralParams, ForceFieldParams, MassParams,
    VdwParams, get_atom_types,
};

impl ForceFieldParams {
    /// From a string of a dat text file, from Amber.
    pub fn from_dat(text: &str) -> io::Result<Self> {
        let mut result = Self::default();

        let mut in_mod4 = false;

        // These dat text-based files are tabular data, and don't have clear delineations bewteen sections.
        // we parse each line based on its content. Notably, the first column alone is a good indicator
        // based on the number of dashes in it. Three atom names separated by dashes, e.g. `pf-p2-s` is angle data.
        // Two, e.g. `ca-s6` is linear bond data. (e.g. springiness of the covalent bond). FOur names indicates
        // a dihedral (aka torsion) angle.
        for (i, line) in text.lines().enumerate() {
            if i == 0 {
                continue; // header line.
            }

            let line = line.trim();

            if line.starts_with("hn  ho  hs") || line.starts_with("hw  ow") {
                continue; // Fragile.
            }

            if line.starts_with("END") {
                break;
            }

            if line.is_empty() {
                // blank line ends MOD4 block
                in_mod4 = false;
                continue;
            }
            // header that *starts* the block
            if line.starts_with("MOD4") {
                in_mod4 = true;
                continue; // nothing else to parse on this header line
            }

            // 1. Break the line into whitespace‐split tokens.
            let cols: Vec<_> = line.split_whitespace().collect();

            let (atom_types, _) = get_atom_types(&cols);

            match atom_types.len() {
                1 => {
                    if in_mod4 {
                        result.van_der_waals.push(VdwParams::from_line(line)?);
                    } else {
                        result.mass.push(MassParams::from_line(line)?);
                    }
                }

                2 => {
                    result.bond.push(BondStretchingParams::from_line(line)?);
                }

                3 => {
                    result.angle.push(AngleBendingParams::from_line(line)?);
                }

                4 => {
                    let (dihedral, improper) = DihedralParams::from_line(line)?;

                    if improper {
                        result.improper.push(dihedral);
                    } else {
                        result.dihedral.push(dihedral);
                    }
                }

                _ => {
                    // anything else
                    println!(
                        "Pushing a remark: ff len: {:?}, {}, {:?}",
                        atom_types,
                        atom_types.len(),
                        line
                    );
                    result.remarks.push(line.to_string());
                }
            }
        }

        // for r in &result.van_der_waals {
        //     println!("Vdw: {:?}", r);
        // }
        // for r in &result.bond {
        //     println!("Bond: {:?}", r);
        // }
        //
        //
        // for r in &result.mass {
        //     println!("Mass: {:?}", r);
        // }
        //
        // for r in &result.angle {
        //     println!("Ang: {:?}", r);
        // }
        //
        // for r in &result.dihedral {
        //     println!("Dih: {:?}", r);
        // }
        // for r in &result.improper {
        //     println!("Imp: {:?}", r);
        // }

        Ok(result)
    }

    /// Write to file
    pub fn save_dat(&self, path: &Path) -> io::Result<()> {
        let mut f = File::create(path)?;

        Ok(())
    }

    /// todo: Sort out the syntax for loading from different sources.
    pub fn load_dat(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        let data_str: String = String::from_utf8(buffer)
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid UTF8"))?;

        Self::from_dat(&data_str)
    }
}
