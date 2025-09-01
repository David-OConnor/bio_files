use std::path::PathBuf;

use bio_files_rs::Mol2 as Mol2Rs;
use pyo3::{prelude::*, types::PyType};

use crate::map_io;

#[pyclass(module = "bio_files")]
pub struct Mol2 {
    inner: Mol2Rs,
}

#[pymethods]
impl Mol2 {
    #[new]
    fn new(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(Mol2Rs::new(text))?,
        })
    }

    fn save(&self, path: PathBuf) -> PyResult<()> {
        map_io(self.inner.save(&path))
    }

    #[classmethod]
    fn load(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(Mol2Rs::load(&path))?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
