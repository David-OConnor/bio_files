# Bio Files: Read and write common biology file formats

[![Crate](https://img.shields.io/crates/v/bio_files.svg)](https://crates.io/crates/bio_files)
[![Docs](https://docs.rs/bio_files/badge.svg)](https://docs.rs/bio_files)


### Currently supported:
- Mol2 (Molecules)
- Map (Electron density, e.g. from crystallography)
- AB1 (Sequence tracing)

### Planned:
- SDF
- PDBQT
- MTZ
- DNA
- CIF auxiliary formats (2fo-fc etc)


For Genbank, we recommend [GB-IO](https://docs.rs/gb-io/latest/gb_io/). For atom coordinate mmCIF
and PDB, we recommend [PDBTBX](https://docs.rs/pdbtbx/latest/pdbtbx/). We do not plan to support
these in these formats in bio_files due to these high-quality libraries.