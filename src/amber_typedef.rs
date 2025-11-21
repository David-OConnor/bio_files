//! Parses Amber Antechamber DEF files, which are used to load forcefield names.

use std::{fs, io, path::Path};

use na_seq::Element;

#[derive(Debug)]
pub struct WildAtom {
    pub name: String,          // e.g. "XX", "XA"
    pub elements: Vec<String>, // e.g. ["C", "N", "O", "O2", "P2"]
}

#[derive(Debug)]
pub struct AtomTypeDef {
    /// f2: Atom type. e.g. "c3", "ca"
    pub name: String,
    /// f3: Residue names, e.g. "AA", or None for "*"
    pub residue: Option<String>,
    /// f4: element (Atomic number in files)
    pub element: Option<Element>,
    /// f5: Number of attached atoms
    pub attached_atoms: Option<u8>,
    /// f6: Number of attached hydrogen atoms
    pub attached_h: Option<u8>, // f6
    /// f7: For hydrogen, number of the electron-withdrawal atoms connected
    pub electron_withdrawal_count: Option<u8>,
    /// f8: Atomic property. e.g. "[RG3]" or "[sb,db,AR2]"
    pub atomic_property: Option<String>,
    /// f9: Chemical environment definitions. e.g. "(C3(C3))"
    pub chem_env: Option<String>,
}

fn parse_u8_field(s: &str) -> Option<u8> {
    if s == "*" { None } else { s.parse().ok() }
}

fn parse_str_field(s: Option<&str>) -> Option<String> {
    match s {
        None => None,
        Some("*") => None,
        Some(v) => Some(v.to_string()),
    }
}

#[derive(Debug)]
pub struct AmberDef {
    pub wildatoms: Vec<WildAtom>,
    pub atomtypes: Vec<AtomTypeDef>,
}

impl AmberDef {
    pub fn new(text: &str) -> io::Result<Self> {
        let mut wildatoms = Vec::new();
        let mut atomtypes = Vec::new();

        for raw_line in text.lines() {
            let line = raw_line.trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with('#') || line.starts_with("//") {
                continue;
            }
            if line.starts_with('=') || line.starts_with('-') {
                continue;
            }

            let mut tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            if tokens.last() == Some(&"&") {
                tokens.pop();
            }

            match tokens[0] {
                "WILDATOM" => {
                    if tokens.len() >= 3 {
                        let name = tokens[1].to_string();
                        let elements = tokens[2..].iter().map(|s| (*s).to_string()).collect();
                        wildatoms.push(WildAtom { name, elements });
                    }
                }
                "ATD" => {
                    if tokens.len() < 2 {
                        continue;
                    }

                    let name = tokens[1].to_string();
                    let residue = parse_str_field(tokens.get(2).copied());
                    let atomic_number = tokens.get(3).and_then(|s| parse_u8_field(s));

                    // todo: Placeholder, given Element is missing items.
                    let element = atomic_number
                        .map(|n| Element::from_atomic_number(n).unwrap_or(Element::Zinc));

                    let attached_atoms = tokens.get(4).and_then(|s| parse_u8_field(s));
                    let attached_h = tokens.get(5).and_then(|s| parse_u8_field(s));
                    let ew_count = tokens.get(6).and_then(|s| parse_u8_field(s));
                    let atomic_property = parse_str_field(tokens.get(7).copied());
                    let env = parse_str_field(tokens.get(8).copied());

                    atomtypes.push(AtomTypeDef {
                        name,
                        residue,
                        element,
                        attached_atoms,
                        attached_h,
                        electron_withdrawal_count: ew_count,
                        atomic_property,
                        chem_env: env,
                    });
                }
                _ => {}
            }
        }

        Ok(Self {
            wildatoms,
            atomtypes,
        })
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let data_str = fs::read_to_string(path)?;
        Self::new(&data_str)
    }
}
