use std::str::FromStr;

use bio_files_rs;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyType};

mod cif_sf;
mod gro;
mod md_params;
mod mmcif;
mod mol2;
mod pdbqt;
mod sdf;

/// Candidate for standalone helper lib.
#[macro_export]
macro_rules! make_enum {
    ($Py:ident, $Native:path, $( $Var:ident ),+ $(,)?) => {
        #[pyclass]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub enum $Py { $( $Var ),+ }

        impl ::core::convert::From<$Py> for $Native {
            fn from(v: $Py) -> Self { match v { $( $Py::$Var => <$Native>::$Var ),+ } }
        }

        impl ::core::convert::From<$Native> for $Py {
            fn from(v: $Native) -> Self { match v { $( <$Native>::$Var => $Py::$Var ),+ } }
        }

        impl $Py {
            pub fn to_native(self) -> $Native {
                self.into()
            }

            pub fn from_native(native: $Native) -> Self {
               native.into()
            }
        }
    };
}

// todo: Blocked due to a restriction in PYO3
/// Candidate for standalone helper lib.
#[macro_export]
macro_rules! field {
    ($name:ident, $ty:ty) => {
        #[getter]
        fn $name(&self) -> $ty {
            self.inner.$name.into()
        }

        // #[setter($name)]
        // // todo: Do we need to use paste! here?
        // fn $name##_set(&mut self, val: $ty) -> $ty {
        //     self.inner.$name = val.into();
        //     val
        // }
    };
}

#[pyclass]
struct AtomGeneric {
    inner: bio_files_rs::AtomGeneric,
}

#[pymethods]
impl AtomGeneric {
    #[getter]
    fn serial_number(&self) -> u32 {
        self.inner.serial_number
    }

    #[getter]
    fn posit(&self) -> [f64; 3] {
        self.inner.posit.to_arr()
    }

    #[getter]
    // todo: String for now
    fn element(&self) -> String {
        self.inner.element.to_string()
    }

    #[getter]
    // todo: String for now
    fn type_in_res(&self) -> Option<String> {
        self.inner.type_in_res.as_ref().map(|v| v.to_string())
    }

    #[getter]
    fn force_field_type(&self) -> Option<String> {
        self.inner.force_field_type.clone()
    }

    #[getter]
    fn partial_charge(&self) -> Option<f32> {
        self.inner.partial_charge
    }

    #[getter]
    fn hetero(&self) -> bool {
        self.inner.hetero
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

make_enum!(
    BondType,
    bio_files_rs::BondType,
    Single,
    Double,
    Triple,
    Aromatic,
    Amide,
    Dummy,
    Unknown,
    NotConnected,
    Quadruple,
    Delocalized,
    PolymericLink
);

#[pymethods]
impl BondType {
    fn to_str_sdf(&self) -> String {
        self.to_native().to_str_sdf()
    }

    #[classmethod]
    fn from_str(_cls: &Bound<'_, PyType>, str: &str) -> PyResult<Self> {
        Ok(bio_files_rs::BondType::from_str(str)?.into())
    }
    fn __str__(&self) -> String {
        self.to_native().to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}

#[pyclass]
struct BondGeneric {
    inner: bio_files_rs::BondGeneric,
}

#[pymethods]
impl BondGeneric {
    #[getter]
    fn bond_type(&self) -> BondType {
        self.bond_type().into()
    }
    #[setter(bond_type)]
    fn bond_type_set(&mut self, val: BondType) {
        self.inner.bond_type = val.into();
    }

    #[getter]
    fn atom_0_sn(&self) -> u32 {
        self.inner.atom_0_sn
    }
    #[setter(atom_0_sn)]
    fn atom_0_sn_set(&mut self, val: u32) {
        self.inner.atom_0_sn = val;
    }

    #[getter]
    fn atom_1_sn(&self) -> u32 {
        self.inner.atom_1_sn
    }
    #[setter(atom_1_sn)]
    fn atom_1_sn_set(&mut self, val: u32) {
        self.inner.atom_1_sn = val;
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass]
struct ResidueType {
    inner: bio_files_rs::ResidueType,
}

#[pymethods]
impl ResidueType {
    #[classmethod]
    fn from_str(_cls: &Bound<'_, PyType>, str: &str) -> Self {
        Self {
            inner: bio_files_rs::ResidueType::from_str(str),
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass]
struct ResidueGeneric {
    inner: bio_files_rs::ResidueGeneric,
}

#[pymethods]
impl ResidueGeneric {
    #[getter]
    fn serial_number(&self) -> u32 {
        self.inner.serial_number
    }
    #[setter(serial_number)]
    fn serial_number_set(&mut self, val: u32) {
        self.inner.serial_number = val;
    }

    #[getter]
    fn res_type<'py>(&self, py: Python<'py>) -> PyResult<Py<ResidueType>> {
        Py::new(
            py,
            ResidueType {
                inner: self.inner.res_type.clone(),
            },
        )
    }

    #[getter]
    fn end(&self) -> ResidueEnd {
        self.inner.end.into()
    }
    #[setter(end)]
    fn end_set(&mut self, val: ResidueEnd) {
        self.inner.end = val.into();
    }

    #[getter]
    fn atom_sns(&self) -> Vec<u32> {
        self.inner.atom_sns.clone()
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

make_enum!(
    ResidueEnd,
    crate::bio_files_rs::ResidueEnd,
    Internal,
    NTerminus,
    CTerminus,
    Hetero
);

#[pymethods]
impl ResidueEnd {
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}

#[pyclass]
struct ChainGeneric {
    inner: bio_files_rs::ChainGeneric,
}

#[pymethods]
impl ChainGeneric {
    #[getter]
    fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[getter]
    fn residue_sns(&self) -> Vec<u32> {
        self.inner.residue_sns.clone()
    }

    #[getter]
    fn atom_sns(&self) -> Vec<u32> {
        self.inner.atom_sns.clone()
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

make_enum!(
    SecondaryStructure,
    bio_files_rs::SecondaryStructure,
    Helix,
    Sheet,
    Coil
);

#[pymethods]
impl SecondaryStructure {
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}

#[pyclass]
struct BackboneSS {
    inner: bio_files_rs::BackboneSS,
}

#[pymethods]
impl BackboneSS {
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

make_enum!(
    ExperimentalMethod,
    bio_files_rs::ExperimentalMethod,
    XRayDiffraction,
    ElectronDiffraction,
    NeutronDiffraction,
    ElectronMicroscopy,
    SolutionNmr
);

#[pymethods]
impl ExperimentalMethod {
    fn to_str_short(&self) -> String {
        self.to_native().to_str_short()
    }

    #[classmethod]
    fn from_str(_cls: &Bound<'_, PyType>, str: &str) -> PyResult<Self> {
        Ok(bio_files_rs::ExperimentalMethod::from_str(str)?.into())
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
    fn __str__(&self) -> String {
        self.to_native().to_string()
    }
}

#[pymodule]
fn biology_files(_py: Python<'_>, m: &Bound<PyModule>) -> PyResult<()> {
    // General
    m.add_class::<AtomGeneric>()?;
    m.add_class::<BondType>()?;
    m.add_class::<BondGeneric>()?;
    m.add_class::<ResidueType>()?;
    m.add_class::<ResidueGeneric>()?;
    m.add_class::<ChainGeneric>()?;
    m.add_class::<SecondaryStructure>()?;
    m.add_class::<BackboneSS>()?;
    m.add_class::<ExperimentalMethod>()?;

    m.add_class::<ResidueEnd>()?;

    // Small molecules
    m.add_class::<mmcif::MmCif>()?;
    m.add_class::<mol2::Mol2>()?;
    m.add_class::<sdf::Sdf>()?;
    m.add_class::<gro::Gro>()?;
    m.add_class::<pdbqt::Pdbqt>()?;

    // Electron density;
    // todo: Map
    m.add_class::<cif_sf::CifStructureFactors>()?;

    m.add_class::<mol2::MolType>()?;

    m.add_class::<md_params::ForceFieldParams>()?;

    m.add_function(wrap_pyfunction!(md_params::load_prmtop, m)?)?;
    m.add_function(wrap_pyfunction!(md_params::save_prmtop, m)?)?;

    Ok(())
}
