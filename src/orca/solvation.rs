//! Implicit and exlicit solvation
//! [Implicit Solvation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html)

use super::make_inp_block;

#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum SolvatorClusterMode {
    #[default]
    None,
    Stochastic,
}

impl SolvatorClusterMode {
    pub fn keyword(self) -> String {
        match self {
            Self::None => String::new(),
            Self::Stochastic => "stochastic".to_owned(),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ImplicitSolvationSurfaceType {
    // todo: This is a stub
    VdwGaussian,
}

impl ImplicitSolvationSurfaceType {
    pub fn keyword(self) -> String {
        match self {
            Self::VdwGaussian => "vdw_gaussian",
        }
        .to_owned()
    }
}

/// [Implicit Solvation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html)
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ImplicitSolvationModel {
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html#the-conductor-like-polarizable-continuum-model-c-pcm
    // todo: Can also be used as a keyword like `CPCM(Water)`?
    Cpcm,
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html#the-smd-solvation-model
    Smd,
    OpenCosmo,
    Alpb,
}

/// [Implicit Solvation: Keyword List](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html#complete-keyword-list-for-the-cpcm-block)
#[derive(Clone, Debug)]
pub struct SolvatorImplicit {
    pub model: ImplicitSolvationModel,
    pub solvent: Solvent,
    pub surface_type: Option<ImplicitSolvationSurfaceType>,
    // todo: GEneralize these keywords?
    // https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html#complete-keyword-list-for-the-cpcm-block
    pub epsilon: Option<f32>,
    pub rsolv: Option<f32>,
    pub draco: bool,
    pub soln: Option<f32>,
    pub soln25: Option<f32>,
    // todo: More fields A/R
}

impl SolvatorImplicit {
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html#complete-keyword-list-for-the-cpcm-block
    pub fn make_inp(&self) -> String {
        // todo: DRY
        let mut contents: Vec<(&str, String)> = vec![("solvent", self.solvent.keyword())];

        if let Some(st) = &self.surface_type {
            contents.push(("surface_type", st.keyword()));
        }

        if let Some(v) = self.epsilon {
            contents.push(("epsilon", format!("{v:.6}")));
        }

        if let Some(v) = self.rsolv {
            contents.push(("rsolv", format!("{v:.6}")));
        }

        if self.draco {
            contents.push(("draco", "true".to_owned()));
        }

        if let Some(v) = self.soln {
            contents.push(("soln", format!("{v:.6}")));
        }

        if let Some(v) = self.soln25 {
            contents.push(("soln25", format!("{v:.6}")));
        }

        match self.model {
            ImplicitSolvationModel::Cpcm => make_inp_block("cpcm", &contents, &[]),
            ImplicitSolvationModel::Smd => {
                contents.push(("smd", "true".to_owned()));

                make_inp_block("cpcm", &contents, &[])
            }
            _ => unimplemented!(),
        }
    }
}

impl ImplicitSolvationModel {
    pub fn keyword(self) -> String {
        match self {
            // Note: It appears the ORCA API treats SMD as a subset of CPCM. (?)
            Self::Cpcm => "cpcm",
            Self::Smd => "smd",
            Self::OpenCosmo => "cosmors",
            Self::Alpb => "alpb",
        }
        .to_owned()
    }
}

/// [Implicit Solvation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html)
#[derive(Clone, Debug)]
// todo: Consider renaming SolvatorExplicit?
pub struct Solvator {
    pub solvent: Solvent,
    pub num_mols: u16,
    pub cluster_mode: SolvatorClusterMode,
    pub droplet: bool,
}

impl Solvator {
    pub fn make_inp(&self) -> String {
        // todo: Add the Solvent too
        let mut contents = vec![
            ("nsolv", self.num_mols.to_string()),
            ("clustermode", self.cluster_mode.keyword()),
        ];

        if self.droplet {
            contents.push(("droplet", "true".to_string()));
        }

        make_inp_block("solvator", &contents, &[])
    }
}

/// [Implicit Solvation, Table 2.56](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html#id42)
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Solvent {
    // todo: Fill out and cite the source
    Water,
    Ethanol,
    Methanol,
    Phenol,
    Amonia,
}

impl Solvent {
    pub fn keyword(self) -> String {
        match self {
            Self::Water => "water",
            Self::Ethanol => "ethanol",
            Self::Methanol => "methanol",
            Self::Phenol => "phenol",
            Self::Amonia => "amonia",
        }
        .to_owned()
    }
}
