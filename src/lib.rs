pub mod map;
pub mod mol2;
mod mtz;
mod sdf;
pub mod ab1;

use lin_alg::f64::Vec3;
use na_seq::Element;
pub use map::*;
pub use mol2::*;
pub use ab1::*;



#[derive(Clone, Debug)]
pub struct AtomGeneric {
    pub serial_number: usize,
    pub posit: Vec3,
    pub element: Element,
    // pub name: String,
    // pub role: Option<AtomRole>,
    // pub residue: Option<usize>,
    pub occupancy: Option<f32>,
    pub partial_charge: Option<f32>,
}

#[derive(Clone, Debug)]
pub struct BondGeneric {
    pub bond_type: String, // todo: Enum
    pub atom_0: usize,
    pub atom_1: usize,
}