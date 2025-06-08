# Bio Files: Read and write common biology file formats

[![Crate](https://img.shields.io/crates/v/bio_files.svg)](https://crates.io/crates/bio_files)
[![Docs](https://docs.rs/bio_files/badge.svg)](https://docs.rs/bio_files)


This library contains functionality to load and save data in common biology file formats. It operates
on data structures that are specific to each file format; you will need to convert to and from the structures
used by your application. The API docs, and examples below are sufficient to get started.


### Currently supported formats:
- Mol2 (Molecules)
- SDF (Molecules)
- Map (Electron density, e.g. from crystallography, Cryo EM)
- AB1 (Sequence tracing)

### Planned:
- PDBQT
- MTZ
- DNA
- CIF structure formats (2fo-fc etc)


For Genbank, we recommend [GB-IO](https://docs.rs/gb-io/latest/gb_io/). For atom coordinate mmCIF
and PDB, we recommend [PDBTBX](https://docs.rs/pdbtbx/latest/pdbtbx/). We do not plan to support
these in these formats in bio_files due to these high-quality libraries.

The API is straightforward to use. It operates using structs with public fields, which you can explore
using the API docs, or your IDE. These structs generally include these three methods: `new()`, `save()` and `load()`.
`new()` accepts `&str` for text files, and a `R: Read + Seek` for binary. `save()` and `load()` accept `&Path`. 
Example use:

```rust
pub fn open_molecule(&mut self, path: &Path) -> io::Result<()> {
    let binding = path.extension().unwrap_or_default().to_ascii_lowercase();
    let extension = binding;

    let molecule = match extension.to_str().unwrap() {
        "sdf" => Ok(Sdf::load(path)?.into()),
        "mol2" => Ok(Mol2::load(path)?.into()),
        _ => ()
    };
}

pub fn open_map(&mut self, path: &Path) -> io::Result<()> {
    let dm = DensityMap::load(path)?;
    
    // Call dm.density_at_point_trilinear(coord) to get density
    // Run `density_to_sig` to get sigma-normalized density, for uniform display.
    
    self.load_density(dm);

    Ok(())
}

/// A single endpoint to save a number of file types
pub fn save(&mut self, path: &Path) -> io::Result<()> {
    let binding = path.extension().unwrap_or_default().to_ascii_lowercase();
    let extension = binding;

    match extension.to_str().unwrap_or_default() {
        "sdf" => match &self.ligand {
            Some(lig) => {
                lig.molecule.to_sdf().save(path)?;
            }
            None => return Err(io::Error::new(ErrorKind::InvalidData, "No ligand to save")),
        },
        "mol2" => match &self.ligand {
            Some(lig) => {
                lig.molecule.to_mol2().save(path)?;
            }
            None => return Err(io::Error::new(ErrorKind::InvalidData, "No ligand to save")),
        }
    }
}
```

Note that the above examples expect that your application has a struct representing the molecule that has
`From<Mol2>`, and `to_mol2(&self)` (etc) methods. The details of these depend on the application. For example:

```rust
impl From<Sdf> for Molecule {
    fn from(m: Sdf) -> Self {
        // We've implemented `From<AtomGeneric>` and `From<ResidueGeneric>` for our application's `Atom` and
        // `Residue`
        let atoms = m.atoms.iter().map(|a| a.into()).collect();
        let residues = m.residues.iter().map(|r| r.into()).collect();

        Self::new(m.ident, atoms, m.chains.clone(), residues, None, None);
    }
}
```
