use std::{collections::HashMap, path::PathBuf};

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

use crate::{AtomGeneric, ChainGeneric, ResidueGeneric};

#[pyclass(module = "bio_files")]
pub struct MmCif {
    inner: bio_files_rs::MmCif,
}

//     pub secondary_structure: Vec<BackboneSS>,
//     pub experimental_method: Option<ExperimentalMethod>,
#[pymethods]
impl MmCif {
    #[getter]
    fn ident(&self) -> &str {
        &self.inner.ident
    }
    #[setter(ident)]
    fn ident_set(&mut self, val: String) {
        self.inner.ident = val;
    }

    #[getter]
    fn metadata(&self) -> &HashMap<String, String> {
        &self.inner.metadata
    }
    #[setter(metadata)]
    fn metadata_set(&mut self, val: HashMap<String, String>) {
        self.inner.metadata = val;
    }

    #[getter]
    fn atoms(&self) -> Vec<AtomGeneric> {
        self.inner
            .atoms
            .iter()
            .map(|a| AtomGeneric { inner: a.clone() })
            .collect()
    }
    #[setter(atoms)]
    fn atoms_set(&mut self, val: Vec<PyRef<'_, AtomGeneric>>) {
        let atoms = val.iter().map(|a| a.inner.clone()).collect();

        self.inner.atoms = atoms;
    }

    // #[getter]
    // fn bonds(&self) -> Vec<BondGeneric> {
    //     self.inner
    //         .bonds
    //         .iter()
    //         .cloned()
    //         .map(|b| BondGeneric { inner: b.clone() })
    //         .collect()
    // }
    // #[setter(bonds)]
    // fn bonds_set(&mut self, val: Vec<PyRef<'_, BondGeneric>>) {
    //     let bonds = val.iter().map(|a| a.inner.clone()).collect();
    //
    //     self.inner.bonds = bonds;
    // }

    #[getter]
    fn chains(&self) -> Vec<ChainGeneric> {
        self.inner
            .chains
            .iter()
            .map(|c| ChainGeneric { inner: c.clone() })
            .collect()
    }

    #[getter]
    fn residues(&self) -> Vec<ResidueGeneric> {
        self.inner
            .residues
            .iter()
            .cloned()
            .map(|r| ResidueGeneric { inner: r.clone() })
            .collect()
    }

    // todo: String for now
    #[getter]
    fn secondary_structure(&self) -> Vec<String> {
        self.inner
            .secondary_structure
            .iter()
            .map(|s| format!("{s:?}"))
            .collect()
    }

    // todo: String for now
    #[getter]
    fn experimental_method(&self) -> Option<String> {
        match self.inner.experimental_method {
            Some(m) => Some(m.to_string()),
            None => None,
        }
    }

    #[new]
    fn new(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::MmCif::new(text)?,
        })
    }

    fn save(&self, path: PathBuf) -> PyResult<()> {
        Ok(self.inner.save(&path)?)
    }

    #[classmethod]
    fn load(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::MmCif::load(&path)?,
        })
    }

    #[classmethod]
    fn load_rcsb(_cls: &Bound<'_, PyType>, ident: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::MmCif::load_rcsb(ident)?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
