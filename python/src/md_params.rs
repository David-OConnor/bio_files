use std::{
    fs::File,
    io,
    path::{Path, PathBuf},
};

use bio_files_rs;
use pyo3::{prelude::*, types::PyType};

use crate::{AtomGeneric, BondGeneric};

#[derive(Clone)]
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

#[pyfunction]
pub fn save_prmtop(
    atoms: Vec<PyRef<AtomGeneric>>,
    params: PyRef<ForceFieldParams>,
    path: PathBuf,
) -> io::Result<()> {
    // Requires that the inner types implement Clone.
    let atoms: Vec<_> = atoms.into_iter().map(|a| a.inner.clone()).collect();
    bio_files_rs::prmtop::save_prmtop(&atoms, &params.inner, &path)
}

#[pyfunction]
pub fn load_prmtop(path: PathBuf) -> io::Result<(Vec<AtomGeneric>, ForceFieldParams)> {
    let (atoms, params) = bio_files_rs::prmtop::load_prmtop(&path)?;
    Ok((
        atoms.into_iter().map(|a| AtomGeneric { inner: a }).collect(),
        ForceFieldParams { inner: params },
    ))
}