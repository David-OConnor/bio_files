use super::make_inp_block;

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum ScfConvergenceTolerance {
    /// Between medium and strong
    #[default]
    None,
    Sloppy,
    Loose,
    Medium,
    Strong,
    Tight,
    VeryTight,
    Extreme,
}

impl ScfConvergenceTolerance {
    pub fn keyword(self) -> String {
        match self {
            Self::None => "",
            Self::Sloppy => "Sloppy",
            Self::Loose => "Loose",
            Self::Medium => "Medium",
            Self::Strong => "Strong",
            Self::Tight => "Tight",
            Self::VeryTight => "VeryTight",
            Self::Extreme => "Extreme",
        }
        .to_string()
    }
}

#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum ScfMode {
    /// Direct is the ORCA default if not specified.
    #[default]
    Direct,
    Conventional,
}

impl ScfMode {
    pub fn keyword(self) -> String {
        match self {
            Self::Direct => "Direct",
            Self::Conventional => "Conventional",
        }
        .to_string()
    }
}

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ScfGuess {
    HCore,
    Hueckel,
    PAtom,
    PModel,
    MORead,
}

impl ScfGuess {
    pub fn keyword(self) -> String {
        match self {
            Self::HCore => "HCore",
            Self::Hueckel => "Hueckel",
            Self::PAtom => "PAtom",
            Self::PModel => "PModel",
            Self::MORead => "MORead",
        }
        .to_owned()
    }
}

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ScfGuessMode {
    FMatrix,
    CMatrix,
}

impl ScfGuessMode {
    pub fn keyword(self) -> String {
        match self {
            Self::FMatrix => "FMatrix",
            Self::CMatrix => "CMatrix",
        }
        .to_owned()
    }
}

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html
/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/integralhandling.html
#[derive(Clone, Debug)]
pub struct Scf {
    // todo: Damping, level shifting etc. Lots more features to implement
    pub convergence_tolerance: ScfConvergenceTolerance,
    pub mode: ScfMode,
    pub thresh: Option<f32>,
    pub t_cut: Option<f32>,
    pub direct_reset_freq: Option<u16>,
    pub max_disk: Option<u16>,
    pub max_int_mem: Option<u16>,
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/finEfield.html#dipolar-electric-fields
    pub e_field: Option<[f32; 3]>,
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/finEfield.html#quadrupolar-electric-fields
    pub q_field: Option<[f32; 6]>, // todo: Allow custom values. See https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html Table 2.9
    pub guess: Option<ScfGuess>,
    pub guess_mode: Option<ScfGuessMode>,
}

impl Scf {
    pub fn make_inp(&self) -> String {
        let mut contents = vec![
            ("Convergence", self.convergence_tolerance.keyword()),
            ("SCFMode", self.mode.keyword()),
        ];

        if let Some(v) = self.thresh {
            contents.push(("Thresh", format!("{v:.6}")));
        }
        if let Some(v) = self.t_cut {
            contents.push(("TCut", format!("{v:.6}")));
        }
        if let Some(v) = self.thresh {
            contents.push(("Thresh", format!("{v}")));
        }
        if let Some(v) = self.max_disk {
            contents.push(("MaxDisk", format!("{v}")));
        }
        if let Some(v) = self.max_int_mem {
            contents.push(("MaxIntMem", format!("{v}")));
        }
        if let Some([x, y, z]) = self.e_field {
            contents.push(("EField", format!("{x:.6} {y:.6} {z:.6}")))
        }
        if let Some([xx, yy, zz, xy, xz, yz]) = self.q_field {
            contents.push((
                "QField",
                format!("{xx:.6} {yy:.6} {zz:.6} {xy:.6} {xz:.6} {yz:.6}"),
            ))
        }
        if let Some(v) = self.guess {
            contents.push(("Guess", v.keyword()));
        }
        if let Some(v) = self.guess_mode {
            contents.push(("GuessMode", v.keyword()));
        }

        make_inp_block("scf", &contents)
    }
}
