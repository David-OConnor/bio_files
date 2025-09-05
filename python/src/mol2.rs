use std::path::PathBuf;

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

use crate::{AtomGeneric, BondGeneric, make_enum};

make_enum!(MolType, bio_files_rs::mol2::MolType,
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
    inner: bio_files_rs::Mol2,
}

#[pymethods]
impl Mol2 {
    #[getter]
    fn ident(&self) -> &str {
        &self.inner.ident
    }

    #[getter]
    fn mol_type(&self) -> MolType {
        MolType::from_native(self.inner.mol_type)
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

    #[getter]
    fn atoms(&self) -> Vec<AtomGeneric> {
        self.inner
            .atoms
            .iter()
            .map(|a| AtomGeneric { inner: a.clone() })
            .collect()
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

    #[new]
    fn new(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Mol2::new(text)?,
        })
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
