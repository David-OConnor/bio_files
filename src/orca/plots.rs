//! [Orbital and Density Plots](https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/plots.html)

use crate::orca::make_inp_block;

#[derive(Clone, Debug)]
pub struct Plots {}

impl Plots {
    pub fn make_inp(&self) -> String {
        let contents = vec![];

        make_inp_block("plots", &contents, &[])
    }
}
