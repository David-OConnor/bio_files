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
    pub name: String,            // f2: e.g. "c3", "ca"
    pub residue: Option<String>, // f3: e.g. "AA", or None for "*"
    // todo: Use Element instead of atomic number?
    // pub atomic_number: Option<u8>,     // f4
    pub element: Option<Element>,        // f4
    pub attached_atoms: Option<u8>,      // f5
    pub attached_h: Option<u8>,          // f6
    pub ew_count: Option<u8>,            // f7
    pub atomic_property: Option<String>, // f8: e.g. "[RG3]" or "[sb,db,AR2]"
    pub env: Option<String>,             // f9: e.g. "(C3(C3))"
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
                    let element = match atomic_number {
                        // todo: Placeholder, given Element is missing items.
                        Some(n) => Some(Element::from_atomic_number(n).unwrap_or(Element::Zinc)),
                        None => None,
                    };

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
                        ew_count,
                        atomic_property,
                        env,
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
