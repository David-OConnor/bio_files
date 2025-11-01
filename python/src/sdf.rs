use std::{collections::HashMap, path::PathBuf};

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

use crate::{AtomGeneric, BondGeneric, ChainGeneric, ResidueGeneric, mol2::Mol2};

#[pyclass(module = "bio_files")]
pub struct Sdf {
    pub inner: bio_files_rs::Sdf,
}

#[pymethods]
impl Sdf {
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

    #[getter]
    fn bonds(&self) -> Vec<BondGeneric> {
        self.inner
            .bonds
            .iter()
            .cloned()
            .map(|b| BondGeneric { inner: b.clone() })
            .collect()
    }
    #[setter(bonds)]
    fn bonds_set(&mut self, val: Vec<PyRef<'_, BondGeneric>>) {
        let bonds = val.iter().map(|a| a.inner.clone()).collect();

        self.inner.bonds = bonds;
    }

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

    #[new]
    fn new(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Sdf::new(text)?,
        })
    }

    fn to_mol2(&self) -> Mol2 {
        Mol2 {
            inner: self.inner.clone().into(),
        }
    }

    fn save(&self, path: PathBuf) -> PyResult<()> {
        Ok(self.inner.save(&path)?)
    }

    #[classmethod]
    fn load(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Sdf::load(&path)?,
        })
    }

    #[classmethod]
    fn load_pubchem(_cls: &Bound<'_, PyType>, cid: u32) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Sdf::load_pubchem(cid)?,
        })
    }

    #[classmethod]
    fn load_drugbank(_cls: &Bound<'_, PyType>, ident: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Sdf::load_drugbank(ident)?,
        })
    }

    #[classmethod]
    fn load_pdbe(_cls: &Bound<'_, PyType>, ident: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Sdf::load_pdbe(ident)?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
