use std::{collections::HashMap, path::PathBuf};

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

// todo: Header; shared with map.
// todo: Miller indices?

#[pyclass(module = "bio_files")]
pub struct MillerIndices {
    pub inner: bio_files_rs::cif_sf::MillerIndices,
}

#[pyclass(module = "bio_files")]
pub struct DensityHeaderInner {
    pub inner: bio_files_rs::DensityHeaderInner,
}

#[pyclass(module = "bio_files")]
pub struct CifStructureFactors {
    pub inner: bio_files_rs::cif_sf::CifStructureFactors,
}

#[pymethods]
impl CifStructureFactors {
    #[getter]
    fn header(&self) -> DensityHeaderInner {
        DensityHeaderInner {
            inner: self.inner.header.clone(),
        }
    }

    #[getter]
    fn miller_indices(&self) -> Vec<MillerIndices> {
        self.inner
            .miller_indices
            .clone()
            .into_iter()
            .map(|i| MillerIndices { inner: i })
            .collect()
    }

    #[new]
    fn new(text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::cif_sf::CifStructureFactors::new(text)?,
        })
    }

    #[classmethod]
    fn new_from_path(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::cif_sf::CifStructureFactors::new_from_path(&path)?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
