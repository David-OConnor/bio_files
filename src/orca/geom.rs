//! https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html

use crate::orca::make_inp_block;

#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum Convergence {
    #[default]
    Normal,
    Loose,
    Tight,
}

impl Convergence {
    pub fn keyword(self) -> String {
        match self {
            Self::Normal => "normal",
            Self::Loose => "loose",
            Self::Tight => "tight",
        }
        .to_string()
    }
}

/// https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html#id8
#[derive(Clone, Debug)]
pub struct Geom {
    // todo: Fill this out a/r
    pub max_iter: u16,
    pub convergence: Convergence,
}

impl Geom {
    pub fn make_inp(&self) -> String {
        let mut contents = vec![("Convergence", self.convergence.keyword())];

        make_inp_block("geom", &contents)
    }
}
