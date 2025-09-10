use std::path::PathBuf;

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

#[pyclass]
pub struct ForceFieldParams {
    inner: bio_files_rs::md_params::ForceFieldParams,
}

#[pymethods]
impl ForceFieldParams {
    // todo: new if you also impl ForceFieldParams.
    // #[new]
    // fn new(inner: bio_files_rs::ForceFieldParams) -> Self {
    //     Self { inner }
    // }

    #[classmethod]
    fn from_frcmod(_cls: &Bound<'_, PyType>, text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::md_params::ForceFieldParams::from_frcmod(text)?,
        })
    }

    #[classmethod]
    fn from_dat(_cls: &Bound<'_, PyType>, text: &str) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::md_params::ForceFieldParams::from_dat(text)?,
        })
    }

    #[classmethod]
    fn load_frcmod(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::md_params::ForceFieldParams::load_frcmod(&path)?,
        })
    }

    #[classmethod]
    fn load_dat(_cls: &Bound<'_, PyType>, path: PathBuf) -> PyResult<Self> {
        Ok(Self {
            inner: bio_files_rs::md_params::ForceFieldParams::load_dat(&path)?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
