//! Experimental methods: HF, DFT, etc.

use crate::orca::make_inp_block;

/// [Single point energies](https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html)
/// [Hartree Fock Theory](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/hartreefock.html)
/// [Density Functional Theory](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html)
///  Note that this currently sets the first line of the ORCA input,
/// not the %method% block.
#[derive(Clone, Copy, PartialEq, Debug, Default)]
// todo: Overlap between this and Functional
pub enum Method {
    /// https://www.faccts.de/docs/orca/6.0/tutorials/prop/single_point.html#hartree-fock-hf
    HartreeFock,
    /// https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/3cmethods.html?q=r2scan&n=0#hf-3c
    Hf_3c,
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
    B97_3c,
    wB97x_3c,
    B97M_V,
    V97M_D3BJ,
    B97M,
    SCANFUNC,
    RSCAN,
    R2SCAN,
    TPSS,
    REVTPSS,
    #[default]
    r2SCAN_3c,
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
    r2SCANH,
    r2SCAN0,
    r2SCAN50,
    PBEh_3c,
    B3LYP_3c,
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
            Self::Hf_3c => "HF-3c",
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
            Self::B97_3c => "B97-3c",
            Self::wB97x_3c => "wB97X-3c",
            Self::B97M_V => "B97M-V",
            Self::V97M_D3BJ => "V97M-D3BJ",
            Self::B97M => "B97M",
            Self::SCANFUNC => "SCANFUNC",
            Self::RSCAN => "RSCAN",
            Self::R2SCAN => "R2SCAN",
            Self::TPSS => "TPSS",
            Self::REVTPSS => "REVTPSS",
            Self::r2SCAN_3c => "r2SCAN-3c",
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
            Self::r2SCANH => "r2SCANH",
            Self::r2SCAN0 => "r2SCAN0",
            Self::r2SCAN50 => "r2SCAN50",
            Self::PBEh_3c => "PBEh-3c",
            Self::B3LYP_3c => "B3LYP-3c",
            // End Global hybrid functionals
            Self::FractionalOccupationDensity => "FOD",
            Self::None => "",
        }
        .to_string()
    }

    /// [Manual 3.6: Composite methods) These methods don't use a separate basis set; they have small tailored
    /// ones of their own.
    pub fn is_composite(self) -> bool {
        matches!(
            self,
            Self::Hf_3c
                | Self::B97_3c
                | Self::r2SCAN_3c
                | Self::PBEh_3c
                | Self::B3LYP_3c
                | Self::wB97x_3c
        )
    }
}

/// https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/dispersioncorrections.html#dispersion-corrections
/// https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/dispersioncorrections.html#id41
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum DispersionCorrection {
    DftDopt,
    D3S6,
    D3S8,
    D3A1,
    D3A2,
    D3RS6,
    D3alpha6,
    D4S6,
    D4S8,
    D4A1,
    D4A2,
    D4S9,
    DFTDScaleC6,
}

impl DispersionCorrection {
    // pub fn keyword(self) -> String {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::DftDopt => "DFTDOPT",
            Self::D3S6 => "D3S6",
            Self::D3S8 => "D3S8",
            Self::D3A1 => "D3A1",
            Self::D3A2 => "D3A2",
            Self::D3RS6 => "D3RS6",
            Self::D3alpha6 => "D3alpha6",
            Self::D4S6 => "D4S6",
            Self::D4S8 => "D4S8",
            Self::D4A1 => "D4A1",
            Self::D4A2 => "D4A2",
            Self::D4S9 => "D4S9",
            Self::DFTDScaleC6 => "DFTDScaleC6",
        }
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
    // Pbe,
    // Bp86,
}

impl Functional {
    pub fn keyword(self) -> String {
        match self {
            // Pbe => "PBE",
            // Self::Bp86 => "BP86",
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

/// [General structure of the Input File](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html)
#[derive(Clone, Debug)]
pub struct MethodSection {
    // todo: How does `functional` differ from `method` key? For example, I see
    // todo BP86 as both in docs examples.
    // pub functional: Option<Functional>,
    pub correlation: Option<Correlation>,
    pub switch_to_soscf: Option<bool>,
    pub frozen_core: Option<FrozenCore>,
    // pub new_n_core: Option<()>
    pub check_frozen_core: Option<bool>,
    pub correct_frozen_core: Option<bool>,
    pub dispersion_correction: Option<(DispersionCorrection, f32)>,
}

impl MethodSection {
    pub fn make_inp(&self) -> String {
        let mut contents = Vec::new();

        // if let Some(f) = self.functional {
        //     contents.push(("functional", f.keyword()));
        // }

        if let Some(c) = self.correlation {
            contents.push(("correlation", c.keyword()));
        }

        if let Some(v) = self.switch_to_soscf {
            contents.push(("SwitchToSOSCF", format!("{v:?}")));
        }

        if let Some(v) = self.frozen_core {
            contents.push(("FrozenCore", v.keyword()));
        }

        if let Some(v) = self.check_frozen_core {
            contents.push(("CheckFrozenCore", format!("{v:?}")));
        }

        if let Some(v) = self.correct_frozen_core {
            contents.push(("CorrectFrozenCore", format!("{v:?}")));
        }

        if let Some((kw, v)) = self.dispersion_correction {
            contents.push((kw.keyword(), format!("{v:?}")));
        }

        make_inp_block("method", &contents, &[])
    }
}
