use crate::orca::make_inp_block;

/// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum Method {
    #[default]
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#hartree-fock-hf
    HartreeFock,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#density-functional-theory-dft
    /// todo: Multiple DFT functionals.
    Dft,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#mp2-perturbation-theory
    Mp2Perturbation,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#spin-component-scaled-mp2-scs-mp2
    SpinComponentScaledMp2,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#orbital-optimized-mp2-oo-mp2
    OrbitalOptimzedMp2,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#regularized-mp2
    RegularlizedMp2,
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#double-hybrid-dft-dh-dft
    DoubleHybridDft,
    TripleHybridDft, // todo: QC this one.
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#coupled-cluster-cc
    CoupledCluster,
    Xtb,
    // todo: Support other semiemperical methods; XTB2 is just one.
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#semiempirical-methods-sqm
    SemiEmpericalMethods,

    // --- Gradient corrected functions (GGA) -----
    // [Docs, Table 3.2](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#gradient-corrected-functionals-meta-ggas)
    BP86,
    BLYP,
    OLYP,
    GLYP,
    XLYP,
    PW91,
    MPWPW,
    MPWLYP,
    PBE,
    RPBE,
    REVPBE,
    RPW86PBE,
    PWP,
    B97_3C,
    B97M_V,
    V97M_D3BJ,
    B97M_D4,
    SCANFUNC,
    RSCAN,
    R2SCAN,
    TPSS,
    REVTPSS,
    R2SCAN_3C,
    // --- End GGAs
    // --- Global hybrid functions. [Docs Table 3.4](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#global-hybrid-functionals)
    B1LYP,
    B3LYP,
    B3LYP_G,
    O3LYP,
    X3LYP,
    B1P86,
    B3PW91,
    PW1PW,
    MPW1PW,
    MPW1LYP,
    PBE0,
    REVPBE0,
    REVPBE38,
    BHANDHLYP,
    M06,
    M062X,
    PW6B95,
    TPSSH,
    TPSS0,
    R2SCANH,
    R2SCAN0,
    R2SCAN50,
    PBEH_3C,
    B3LYP_3C,
    // --- End HFX
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/fod.html
    FractionalOccupationDensity,
    None, // todo: I believe this is right?
}

impl Method {
    /// Prefixed with an !, starts the .inp file.
    pub fn keyword(self) -> String {
        match self {
            Self::HartreeFock => "HF",
            Self::Dft => "B3LYP",
            Self::Mp2Perturbation => "RI-MP2", // Or "DLPNO-MP2"
            Self::SpinComponentScaledMp2 => "RI-SCS-MP2",
            Self::OrbitalOptimzedMp2 => "OO-RI-MP2",
            Self::RegularlizedMp2 => "RI-MP2",
            Self::DoubleHybridDft => "B2PLYP",
            Self::TripleHybridDft => "B3LYP-D3(BJ)", // todo: QC this one
            Self::CoupledCluster => "DLPNO-CCSD(T)",
            Self::Xtb => "XTB", // todo?
            Self::SemiEmpericalMethods => "XTB2",
            // GGA
            Self::BP86 => "BP86",
            Self::BLYP => "BLYP",
            Self::OLYP => "OLYP",
            Self::GLYP => "GLYP",
            Self::XLYP => "XLYP",
            Self::PW91 => "PW91",
            Self::MPWPW => "MPWPW",
            Self::MPWLYP => "MPWLYP",
            Self::PBE => "PBE",
            Self::RPBE => "RPBE",
            Self::REVPBE => "REVPBE",
            Self::RPW86PBE => "RPW86PBE",
            Self::PWP => "PWP",
            Self::B97_3C => "B97-3C",
            Self::B97M_V => "B97M-V",
            Self::V97M_D3BJ => "V97M-D3BJ",
            Self::B97M_D4 => "B97M-D4",
            Self::SCANFUNC => "SCANFUNC",
            Self::RSCAN => "RSCAN",
            Self::R2SCAN => "R2SCAN",
            Self::TPSS => "TPSS",
            Self::REVTPSS => "REVTPSS",
            Self::R2SCAN_3C => "R2SCAN-3C",
            // end GGA
            // Start Global Hybrid functionals
            Self::B1LYP => "B1LYP",
            Self::B3LYP => "B3LYP",
            Self::B3LYP_G => "B3LYP-G",
            Self::O3LYP => "O3LYP",
            Self::X3LYP => "X3LYP",
            Self::B1P86 => "B1P86",
            Self::B3PW91 => "B3PW91",
            Self::PW1PW => "PW1PW",
            Self::MPW1PW => "MPW1PW",
            Self::MPW1LYP => "MPW1LYP",
            Self::PBE0 => "PBE0",
            Self::REVPBE0 => "REVPBE0",
            Self::REVPBE38 => "REVPBE38",
            Self::BHANDHLYP => "BHANDHLYP",
            Self::M06 => "M06",
            Self::M062X => "M062X",
            Self::PW6B95 => "PW6B95",
            Self::TPSSH => "TPSSH",
            Self::TPSS0 => "TPSS0",
            Self::R2SCANH => "R2SCANH",
            Self::R2SCAN0 => "R2SCAN0",
            Self::R2SCAN50 => "R2SCAN50",
            Self::PBEH_3C => "PBEH-3C",
            Self::B3LYP_3C => "B3LYP-3C",
            // End Global hybrid functionals
            Self::FractionalOccupationDensity => "FOD",
            Self::None => "",
        }
        .to_string()
    }
}

// todo: Rename etc as required, to distinguish method section
// todo items from first-line ones.

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum FrozenCore {
    FcElectrons,
    FcEwin,
    FcNone,
    // n(u16), // todo: QC this
}

impl FrozenCore {
    pub fn keyword(self) -> String {
        match self {
            Self::FcElectrons => "FC_ELECTRONS",
            Self::FcEwin => "FC_EWIN",
            Self::FcNone => "FC_NONE",
        }
        .to_string()
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
// todo: A/R
pub enum Functional {
    Bp86,
}

impl Functional {
    pub fn keyword(self) -> String {
        match self {
            Self::Bp86 => "BP86",
        }
        .to_string()
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
// todo: A/R
pub enum Correlation {
    C_LYP,
}

impl Correlation {
    pub fn keyword(self) -> String {
        match self {
            Self::C_LYP => "C_LYP",
        }
        .to_string()
    }
}

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html
#[derive(Clone, Debug)]
pub struct MethodSection {
    // todo: How does `funcitonal` differ from `method` key? For example, I see
    // todo BP86 as both in docs examples.
    pub functional: Option<Functional>,
    pub correlation: Option<Correlation>,
    pub switch_to_soscf: Option<bool>,
    pub frozen_core: Option<FrozenCore>,
    // pub new_n_core: Option<()>
    pub check_frozen_core: Option<bool>,
    pub correct_frozen_core: Option<bool>,
}

impl MethodSection {
    pub fn make_inp(&self) -> String {
        let mut contents = Vec::new();

        if let Some(f) = self.functional {
            contents.push(("functional", f.keyword()));
        }

        if let Some(c) = self.correlation {
            contents.push(("correlation", c.keyword()));
        }

        if let Some(v) = self.switch_to_soscf {
            contents.push(("SwitchToSOSCF", format!("{v:?}")));
        }

        if let Some(fc) = self.frozen_core {
            contents.push(("FrozenCore", fc.keyword()));
        }

        if let Some(fc) = self.check_frozen_core {
            contents.push(("CheckFrozenCore", format!("{fc:?}")));
        }

        if let Some(fc) = self.correct_frozen_core {
            contents.push(("CorrectFrozenCore", format!("{fc:?}")));
        }

        make_inp_block("method", &contents)
    }
}
