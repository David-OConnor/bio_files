//! For parsing [Amber force field](https://ambermd.org/AmberModels.php) dat files, like
//! `gaff2.dat` (small molecules/ligands). These include parameters used in molecular dynamics.
//! See also: `frcmod`, which patches these parameters for specific molecules.

// todo: For now, this has much overlap with loading frcmod data.
// todo: Reoncosider the API, how this and the frcmod modules are related, and which
// todo: to export the main FF struct from

use std::{
    collections::HashMap,
    fs::File,
    io,
    io::{ErrorKind, Read},
    path::Path,
};

use crate::md_params::{
    AngleBendingParams, BondStretchingParams, DihedralParams, ForceFieldParamsVec, LjParams,
    MassParams, get_atom_types,
};

impl ForceFieldParamsVec {
    /// From a string of a dat text file, from Amber.
    pub fn from_dat(text: &str) -> io::Result<Self> {
        let mut result = Self::default();

        // Handles the lines `N   NA  N2  N*  NC  NB  NT  NY` above vdw data. This, for example,
        // maps NA, N2 etc to the same parameters as N.
        let mut vdw_alias_map: HashMap<String, String> = HashMap::new();

        let mut in_mod4 = false;

        // These dat text-based files are tabular data, and don't have clear delineations bewteen sections.
        // we parse each line based on its content. Notably, the first column alone is a good indicator
        // based on the number of dashes in it. Three atom names separated by dashes, e.g. `pf-p2-s` is angle data.
        // Two, e.g. `ca-s6` is linear bond data. (e.g. springiness of the covalent bond). FOur names indicates
        // a dihedral (aka torsion) angle.
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            // Header or blank – also resets the MOD4 block if we just left it.
            if i == 0 || line.is_empty() {
                in_mod4 = false;
                continue;
            }

            if line.starts_with("hw  ow  0") {
                // gaff2.dat's fast-water line. Will cause a parsing error unless handled.
                continue;
            }

            if line.starts_with("END") {
                break;
            }

            let cols: Vec<&str> = line.split_whitespace().collect();

            // header that *starts* the block
            if line.starts_with("MOD4") {
                in_mod4 = true;
                continue; // nothing else to parse on this header line
            }

            // Handle alias lines
            if !in_mod4
                && cols.len() > 1
                && cols
                    .iter()
                    .all(|t| t.chars().all(|c| c.is_ascii_alphanumeric() || c == '*'))
            {
                let canonical = cols[0].to_string();
                for alias in &cols {
                    vdw_alias_map.insert((*alias).to_string(), canonical.clone());
                }
                continue; // don’t try to parse this line any further
            }

            let (atom_types, _) = get_atom_types(&cols);

            match atom_types.len() {
                1 => {
                    if in_mod4 {
                        let vdw = LjParams::from_line(line)?;

                        // Produce copies for all matching the alias. (The alias line should
                        // be above all individual VDW lines).
                        result.lennard_jones.push(vdw.clone());

                        for (alias, canonical) in vdw_alias_map.iter() {
                            if canonical == &vdw.atom_type && alias != canonical {
                                let mut alias_vdw = vdw.clone();
                                alias_vdw.atom_type = alias.clone();
                                result.lennard_jones.push(alias_vdw);
                            }
                        }
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
