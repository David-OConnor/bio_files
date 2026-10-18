//! Embedded PubChem GHS pictograms for dark user interfaces.
//!
//! SVGs use white hazard symbols, red frames and transparent backgrounds. The fish's
//! interior cutout is dark. These display variants require no network retrieval.

pub const GHS01: &[u8] = include_bytes!("pubchem/ghs/GHS01.svg");
pub const GHS02: &[u8] = include_bytes!("pubchem/ghs/GHS02.svg");
pub const GHS03: &[u8] = include_bytes!("pubchem/ghs/GHS03.svg");
pub const GHS04: &[u8] = include_bytes!("pubchem/ghs/GHS04.svg");
pub const GHS05: &[u8] = include_bytes!("pubchem/ghs/GHS05.svg");
pub const GHS06: &[u8] = include_bytes!("pubchem/ghs/GHS06.svg");
pub const GHS07: &[u8] = include_bytes!("pubchem/ghs/GHS07.svg");
pub const GHS08: &[u8] = include_bytes!("pubchem/ghs/GHS08.svg");
pub const GHS09: &[u8] = include_bytes!("pubchem/ghs/GHS09.svg");

/// Look up the bundled SVG by its standard pictogram code (GHS01 through GHS09).
pub fn ghs_pictogram(code: &str) -> Option<&'static [u8]> {
    match code {
        "GHS01" => Some(GHS01),
        "GHS02" => Some(GHS02),
        "GHS03" => Some(GHS03),
        "GHS04" => Some(GHS04),
        "GHS05" => Some(GHS05),
        "GHS06" => Some(GHS06),
        "GHS07" => Some(GHS07),
        "GHS08" => Some(GHS08),
        "GHS09" => Some(GHS09),
        _ => None,
    }
}
