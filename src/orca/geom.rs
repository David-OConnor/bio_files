//! [Geometry Optimizations](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html)

use crate::orca::make_inp_block;

/// [Geometry Optimization Thresholds](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html#geometry-optimization-thresholds)
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

/// [Geometry Optimizations, Table 4.4](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html#id8)
#[derive(Clone, Debug)]
pub struct Geom {
    // todo: Fill this out a/r
    pub max_iter: u16,
    pub convergence: Convergence,
    pub in_hess: Option<String>,
    pub print_internal_hess: bool,
}

impl Geom {
    pub fn make_inp(&self) -> String {
        let mut contents = vec![("Convergence", self.convergence.keyword())];

        let mut keywords = Vec::new();
        if let Some(in_hess_name) = &self.in_hess {
            contents.push(("inhessname", in_hess_name.to_string()));

            keywords.push("inhess");
            keywords.push("read");
        }

        make_inp_block("geom", &contents, &keywords)
    }
}
