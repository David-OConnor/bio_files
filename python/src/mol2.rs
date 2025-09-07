use std::{collections::HashMap, path::PathBuf};

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

use crate::{AtomGeneric, BondGeneric, make_enum, sdf::Sdf};

make_enum!(
    MolType,
    bio_files_rs::mol2::MolType,
    Small,
    Bipolymer,
    Protein,
    NucleicAcid,
    Saccharide,
);

#[pymethods]
impl MolType {
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}

// todo: ChargeType as well.

#[pyclass(module = "bio_files")]
pub struct Mol2 {
    pub inner: bio_files_rs::Mol2,
}

#[pymethods]
impl Mol2 {
    // todo: Blocked by Pyo3 on macros here.
    // field!(ident, String);
    // field!(mol_type, MolType);

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
    fn mol_type(&self) -> MolType {
        MolType::from_native(self.inner.mol_type)
    }
    #[setter(mol_type)]
    fn mol_type_set(&mut self, val: MolType) {
        self.inner.mol_type = val.into();
    }

    // todo: str for now
    #[getter]
    // fn charge_type(&self) -> ChargeType {
    fn charge_type(&self) -> String {
        self.inner.charge_type.to_string()
    }

    #[getter]
    fn comment(&self) -> Option<String> {
        self.inner.comment.clone()
    }
    #[setter(comment)]
    fn comment_set(&mut self, val: String) {
        self.inner.comment = val.into();
    }

    #[new]
    fn new(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Mol2::new(text)?,
        })
    }

    fn to_sdf(&self) -> Sdf {
        Sdf {
            inner: self.inner.clone().into(),
        }
    }

    fn save(&self, path: PathBuf) -> PyResult<()> {
        Ok(self.inner.save(&path)?)
    }

    #[classmethod]
    fn load(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Mol2::load(&path)?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
