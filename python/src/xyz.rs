use std::{collections::HashMap, path::PathBuf};

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

use crate::AtomGeneric;

#[pyclass(module = "bio_files")]
pub struct Xyz {
    pub inner: bio_files_rs::Xyz,
}

#[pymethods]
impl Xyz {
    #[getter]
    fn atoms(&self) -> Vec<AtomGeneric> {
        self.inner
            .atoms
            .iter()
            .map(|a| AtomGeneric { inner: a.clone() })
            .collect()
    }

    #[getter]
    fn comment(&self) -> String {
        self.inner.comment.clone()
    }
    #[setter(comment)]
    fn comment_set(&mut self, val: String) {
        self.inner.comment = val.into();
    }

    #[new]
    fn new(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Xyz::new(text)?,
        })
    }

    fn save(&self, path: PathBuf) -> PyResult<()> {
        Ok(self.inner.save(&path)?)
    }

    #[classmethod]
    fn load(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::Xyz::load(&path)?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
