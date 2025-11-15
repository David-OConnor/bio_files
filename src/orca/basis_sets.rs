//! A separate file, as these are quite lengthy!


use std::fmt::{Display, Formatter};
use BasisSet::*;

/// https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html
/// See Table 2.12, 2.13 and further
/// todo: Do we want to split this up into wrapped enums by category for organization,
/// todo: or let downstream applications handle that?
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum BasisSet {
    #[default]
    None,
    // --- Pople start
    Sto3G,
    B3_21G,
    B3_21GSP,
    B4_22GSP,
    B6_31G,
    B6_31GStar,
    M6_31G,
    M6_31GStar,
    B6_31GStarStar,
    B6_31G_d,
    B6_31G_d_p,
    B6_31G_2d,
    B6_31G_2d_p,
    B6_31G_2d_2p,
    B6_31G_2df,
    B6_31G_2df_2p,
    B6_31G_2df_2pd,
    B6_31PlusGStar,
    B6_31PlusGStarStar,
    B6_31PlusG_d,
    B6_31PlusG_d_p,
    B6_31PlusG_2d,
    B6_31PlusG_2d_p,
    B6_31PlusG_2d_2p,
    B6_31PlusG_2df,
    B6_31PlusG_2df_2p,
    B6_31PlusG_2df_2pd,
    B6_31PlusPlusGStarStar,
    B6_31PlusPlusG_d_p,
    B6_31PlusPlusG_2d_p,
    B6_31PlusPlusG_2d_2p,
    B6_31PlusPlusG_2df_2p,
    B6_31PlusPlusG_2df_2pd,
    B6_31PlusPlusG_3df_3pd,
    B6_311G,
    B6_311GStar,
    B6_311GStarStar,
    B6_311G_d,
    B6_311G_d_p,
    B6_311G_2d,
    B6_311G_2d_p,
    B6_311G_2d_2p,
    B6_311G_2df,
    B6_311G_2df_2p,
    B6_311G_2df_2pd,
    B6_311G_3df,
    B6_311G_3df_3pd,
    B6_311PlusGStar,
    B6_311PlusGStarStar,
    B6_311PlusG_d,
    B6_311PlusG_d_p,
    B6_311PlusG_2d,
    B6_311PlusG_2d_p,
    B6_311PlusG_2d_2p,
    B6_311PlusG_2df,
    B6_311PlusG_2df_2p,
    B6_311PlusG_2df_2pd,
    B6_311PlusG_3df,
    B6_311PlusG_3df_2p,
    B6_311PlusG_3df_3pd,
    B6_311PlusPlusGStarStar,
    B6_311PlusPlusG_d_p,
    B6_311PlusPlusG_2d_p,
    B6_311PlusPlusG_2d_2p,
    B6_311PlusPlusG_2df_2p,
    B6_311PlusPlusG_2df_2pd,
    B6_311PlusPlusG_3df_3pd,
    // --- Pople end
    // --- Ahlrich start. Table 2.13
    // Ahlrichs valence basis sets (H–Kr)
    Sv,
    SvP,
    Svp,
    Tzv,
    TzvP,
    Tzvp,
    Tzvpp,
    Qzvp,
    Qzvpp,

    // Ahlrichs def- family (H–Lr, def-ECP for Rb–Lr / Fr–Lr)
    DefSvP,
    DefSvp,
    DefTzvp,
    DefTzvpp,
    MaDefTzvp,
    // --- Ahlrich end
    // --- Karlsruhe def2 start
    // Karlsruhe def2- family
    Def2Svp,
    Def2Svp_,
    Def2Tzvp,
    Def2TzvpMinusF,
    Def2Tzvpp,
    Def2Qzvp,
    Def2Qzvpp,

    // Diffuse def2- sets
    Def2Svpd,
    Def2Tzvpd,
    Def2Tzvppd,
    Def2Qzvpd,
    Def2Qzvppd,

    // Minimally augmented ma-def2- sets
    MaDef2Svp,
    MaDef2SvP,
    MaDef2MSvp,
    MaDef2Tzvp,
    MaDef2TzvpMinusF,
    MaDef2Tzvpp,
    MaDef2Qzvpp,
    // --- Karlsruhe def2 end
    // --- Karlsruhe dhf start
    // Karlsruhe dhf- family
    DhfSvp_,
    DhfSvp,
    DhfTzvp,
    DhfTzvpp,
    DhfQzvp,
    DhfQzvpp,

    // Karlsruhe dhf- two-component variants
    DhfSvp2c,
    DhfTzvp2c,
    DhfTzvpp2c,
    DhfQzvp2c,
    DhfQzvpp2c,
    // --- Karlsruhe dhf end
    // todo: Jensen, hydrogen and others.
    // --- Start Correlation-consistent
    // Correlation-consistent cc-pVnZ
    CcPvdz,
    CcPvtz,
    CcPvqz,
    CcPv5z,
    CcPv6z,

    // Augmented aug-cc-pVnZ
    AugCcPvdz,
    AugCcPvtz,
    AugCcPvqz,
    AugCcPv5z,
    AugCcPv6z,

    // Tight-d variants cc-pVn(+d)Z
    CcPvdPlusdZ,
    CcPvtPlusdZ,
    CcPvqPlusdZ,
    CcPv5PlusdZ,

    // Tight-d augmented aug-cc-pVn(+d)Z
    AugCcPvdPlusdZ,
    AugCcPvtPlusdZ,
    AugCcPvqPlusdZ,
    AugCcPv5PlusdZ,
    AugCcPv6PlusdZ,

    // Partially augmented Truhlar sets
    AprCcPvQPlusdZ,
    MayCcPvTPlusdZ,
    MayCcPvQPlusdZ,
    JunCcPvDPlusdZ,
    JunCcPvTPlusdZ,
    JunCcPvQPlusdZ,
    JulCcPvDPlusdZ,
    JulCcPvTPlusdZ,
    JulCcPvQPlusdZ,
    MaugCcPvDPlusdZ,
    MaugCcPvTPlusdZ,
    MaugCcPvQPlusdZ,

    // Core-valence cc-pCVnZ
    CcPcvdz,
    CcPcvtz,
    CcPcvqz,
    CcPcv5z,
    CcPcv6z,

    // Augmented core-valence aug-cc-pCVnZ
    AugCcPcvdz,
    AugCcPcvtz,
    AugCcPcvqz,
    AugCcPcv5z,
    AugCcPcv6z,

    // Weighted core-valence cc-pwCVnZ
    CcPwCvdz,
    CcPwCvtz,
    CcPwCvqz,
    CcPwCv5z,

    // Augmented weighted core-valence aug-cc-pwCVnZ
    AugCcPwCvdz,
    AugCcPwCvtz,
    AugCcPwCvqz,
    AugCcPwCv5z,

    // Pseudopotential cc-pVnZ-PP
    CcPvdzPp,
    CcPvtzPp,
    CcPvqzPp,
    CcPv5zPp,

    // Augmented pseudopotential aug-cc-pVnZ-PP
    AugCcPvdzPp,
    AugCcPvtzPp,
    AugCcPvqzPp,
    AugCcPv5zPp,

    // Core-valence pseudopotential cc-pCVnZ-PP
    CcPcvdzPp,
    CcPcvtzPp,
    CcPcvqzPp,
    CcPcv5zPp,

    // Augmented core-valence pseudopotential aug-cc-pCVnZ-PP
    AugCcPcvdzPp,
    AugCcPcvtzPp,
    AugCcPcvqzPp,
    AugCcPcv5zPp,

    // Weighted core-valence pseudopotential cc-pwCVnZ-PP
    CcPwCvdzPp,
    CcPwCvtzPp,
    CcPwCvqzPp,
    CcPwCv5zPp,

    // Augmented weighted core-valence pseudopotential aug-cc-pwCVnZ-PP
    AugCcPwCvdzPp,
    AugCcPwCvtzPp,
    AugCcPwCvqzPp,
    AugCcPwCv5zPp,

    // Hyperfine-optimized
    AugCcPvtzJ,

    // W4 theory haV sets
    HaVTPlusdZ,
    HaVQPlusdZ,
    HaV5PlusdZ,

    // --- End Correlation-consistent
}

impl BasisSet {
    /// Prefixed with an !, starts the .inp file.
    pub fn keyword(self) -> String {
        match self {
            None => "",
            // --- Pople start
            Sto3G => "STO-3G",
            B3_21G => "3-21G",
            B3_21GSP => "3-21GSP",
            B4_22GSP => "4-22GSP",
            //
            B6_31G => "6-31G",
            B6_31GStar => "6-31G*",
            M6_31G => "m6-31G",
            M6_31GStar => "m6-31G*",
            B6_31GStarStar => "6-31G**",
            B6_31G_d => "6-31G(d)",
            B6_31G_d_p => "6-31G(d,p)",
            B6_31G_2d => "6-31G(2d)",
            B6_31G_2d_p => "6-31G(2d,p)",
            B6_31G_2d_2p => "6-31G(2d,2p)",
            B6_31G_2df => "6-31G(2df)",
            B6_31G_2df_2p => "6-31G(2df,2p)",
            B6_31G_2df_2pd => "6-31G(2df,2pd)",
            //
            B6_31PlusGStar => "6-31+G*",
            B6_31PlusGStarStar => "6-31+G**",
            B6_31PlusG_d => "6-31+G(d)",
            B6_31PlusG_d_p => "6-31+G(d,p)",
            B6_31PlusG_2d => "6-31+G(2d)",
            B6_31PlusG_2d_p => "6-31+G(2d,p)",
            B6_31PlusG_2d_2p => "6-31+G(2d,2p)",
            B6_31PlusG_2df => "6-31+G(2df)",
            B6_31PlusG_2df_2p => "6-31+G(2df,2p)",
            B6_31PlusG_2df_2pd => "6-31+G(2df,2pd)",
            //
            B6_31PlusPlusGStarStar => "6-31++G**",
            B6_31PlusPlusG_d_p => "6-31++G(d,p)",
            B6_31PlusPlusG_2d_p => "6-31++G(2d,p)",
            B6_31PlusPlusG_2d_2p => "6-31++G(2d,2p)",
            B6_31PlusPlusG_2df_2p => "6-31++G(2df,2p)",
            B6_31PlusPlusG_2df_2pd => "6-31++G(2df,2pd)",
            B6_31PlusPlusG_3df_3pd => "6-31++G(3df,3pd)",
            //
            B6_311G => "6-311G",
            B6_311GStar => "6-311G*",
            B6_311GStarStar => "6-311G**",
            B6_311G_d => "6-311G(d)",
            B6_311G_d_p => "6-311G(d,p)",
            B6_311G_2d => "6-311G(2d)",
            B6_311G_2d_p => "6-311G(2d,p)",
            B6_311G_2d_2p => "6-311G(2d,2p)",
            B6_311G_2df => "6-311G(2df)",
            B6_311G_2df_2p => "6-311G(2df,2p)",
            B6_311G_2df_2pd => "6-311G(2df,2pd)",
            B6_311G_3df => "6-311G(3df)",
            B6_311G_3df_3pd => "6-311G(3df,3pd)",
            //
            B6_311PlusGStar => "6-311+G*",
            B6_311PlusGStarStar => "6-311+G**",
            B6_311PlusG_d => "6-311+G(d)",
            B6_311PlusG_d_p => "6-311+G(d,p)",
            B6_311PlusG_2d => "6-311+G(2d)",
            B6_311PlusG_2d_p => "6-311+G(2d,p)",
            B6_311PlusG_2d_2p => "6-311+G(2d,2p)",
            B6_311PlusG_2df => "6-311+G(2df)",
            B6_311PlusG_2df_2p => "6-311+G(2df,2p)",
            B6_311PlusG_2df_2pd => "6-311+G(2df,2pd)",
            B6_311PlusG_3df => "6-311+G(3df)",
            B6_311PlusG_3df_2p => "6-311+G(3df,2p)",
            B6_311PlusG_3df_3pd => "6-311+G(3df,3pd)",
            //
            B6_311PlusPlusGStarStar => "6-311++G**",
            B6_311PlusPlusG_d_p => "6-311++G(d,p)",
            B6_311PlusPlusG_2d_p => "6-311++G(2d,p)",
            B6_311PlusPlusG_2d_2p => "6-311++G(2d,2p)",
            B6_311PlusPlusG_2df_2p => "6-311++G(2df,2p)",
            B6_311PlusPlusG_2df_2pd => "6-311++G(2df,2pd)",
            B6_311PlusPlusG_3df_3pd => "6-311++G(3df,3pd)",
            // --- Pople end
            // --- Ahlrich start
            // Ahlrichs valence basis sets
            Sv => "SV",
            SvP => "SV(P)",
            Svp => "SVP",
            Tzv => "TZV",
            TzvP => "TZV(P)",
            Tzvp => "TZVP",
            Tzvpp => "TZVPP",
            Qzvp => "QZVP",
            Qzvpp => "QZVPP",

            // Ahlrichs def- family
            DefSvP => "def-SV(P)",
            DefSvp => "def-SVP",
            DefTzvp => "def-TZVP",
            DefTzvpp => "def-TZVPP",
            MaDefTzvp => "ma-def-TZVP",
            // --- Ahlrich end
            // --- Karlsruhe def2 start
            // Karlsruhe def2- family
            Def2Svp => "DEF2-SVP",
            Def2Svp_ => "def2-SV(P)",
            Def2Tzvp => "def2-TZVP",
            Def2TzvpMinusF => "def2-TZVP(-f)",
            Def2Tzvpp => "def2-TZVPP",
            Def2Qzvp => "def2-QZVP",
            Def2Qzvpp => "def2-QZVPP",

            // Diffuse def2- sets
            Def2Svpd => "def2-SVPD",
            Def2Tzvpd => "def2-TZVPD",
            Def2Tzvppd => "def2-TZVPPD",
            Def2Qzvpd => "def2-QZVPD",
            Def2Qzvppd => "def2-QZVPPD",

            // Minimally augmented ma-def2- sets
            MaDef2Svp => "ma-def2-SVP",
            MaDef2SvP => "ma-def2-SV(P)",
            MaDef2MSvp => "ma-def2-mSVP",
            MaDef2Tzvp => "ma-def2-TZVP",
            MaDef2TzvpMinusF => "ma-def2-TZVP(-f)",
            MaDef2Tzvpp => "ma-def2-TZVPP",
            MaDef2Qzvpp => "ma-def2-QZVPP",
            // --- Karlsruhe def2 end
            // --- Karlsruhe dhf start
            // Karlsruhe dhf-
            DhfSvp_ => "dhf-SV(P)",
            DhfSvp => "dhf-SVP",
            DhfTzvp => "dhf-TZVP",
            DhfTzvpp => "dhf-TZVPP",
            DhfQzvp => "dhf-QZVP",
            DhfQzvpp => "dhf-QZVPP",

            // Karlsruhe dhf- two-component
            DhfSvp2c => "dhf-SVP-2c",
            DhfTzvp2c => "dhf-TZVP-2c",
            DhfTzvpp2c => "dhf-TZVPP-2c",
            DhfQzvp2c => "dhf-QZVP-2c",
            DhfQzvpp2c => "dhf-QZVPP-2c",
            // --- Karlsruhe dhf end
            // --- Start Correlation-consistent

            // Correlation-consistent cc-pVnZ
            CcPvdz => "cc-pVDZ",
            CcPvtz => "cc-pVTZ",
            CcPvqz => "cc-pVQZ",
            CcPv5z => "cc-pV5Z",
            CcPv6z => "cc-pV6Z",

            // Augmented aug-cc-pVnZ
            AugCcPvdz => "aug-cc-pVDZ",
            AugCcPvtz => "aug-cc-pVTZ",
            AugCcPvqz => "aug-cc-pVQZ",
            AugCcPv5z => "aug-cc-pV5Z",
            AugCcPv6z => "aug-cc-pV6Z",

            // Tight-d variants
            CcPvdPlusdZ => "cc-pVD(+d)Z",
            CcPvtPlusdZ => "cc-pVT(+d)Z",
            CcPvqPlusdZ => "cc-pVQ(+d)Z",
            CcPv5PlusdZ => "cc-pV5(+d)Z",

            // Tight-d augmented
            AugCcPvdPlusdZ => "aug-cc-pVD(+d)Z",
            AugCcPvtPlusdZ => "aug-cc-pVT(+d)Z",
            AugCcPvqPlusdZ => "aug-cc-pVQ(+d)Z",
            AugCcPv5PlusdZ => "aug-cc-pV5(+d)Z",
            AugCcPv6PlusdZ => "aug-cc-pV6(+d)Z",

            // Partially augmented Truhlar sets
            AprCcPvQPlusdZ => "apr-cc-pV(Q+d)Z",
            MayCcPvTPlusdZ => "may-cc-pV(T+d)Z",
            MayCcPvQPlusdZ => "may-cc-pV(Q+d)Z",
            JunCcPvDPlusdZ => "jun-cc-pV(D+d)Z",
            JunCcPvTPlusdZ => "jun-cc-pV(T+d)Z",
            JunCcPvQPlusdZ => "jun-cc-pV(Q+d)Z",
            JulCcPvDPlusdZ => "jul-cc-pV(D+d)Z",
            JulCcPvTPlusdZ => "jul-cc-pV(T+d)Z",
            JulCcPvQPlusdZ => "jul-cc-pV(Q+d)Z",
            MaugCcPvDPlusdZ => "maug-cc-pV(D+d)Z",
            MaugCcPvTPlusdZ => "maug-cc-pV(T+d)Z",
            MaugCcPvQPlusdZ => "maug-cc-pV(Q+d)Z",

            // Core-valence cc-pCVnZ
            CcPcvdz => "cc-pCVDZ",
            CcPcvtz => "cc-pCVTZ",
            CcPcvqz => "cc-pCVQZ",
            CcPcv5z => "cc-pCV5Z",
            CcPcv6z => "cc-pCV6Z",

            // Augmented core-valence
            AugCcPcvdz => "aug-cc-pCVDZ",
            AugCcPcvtz => "aug-cc-pCVTZ",
            AugCcPcvqz => "aug-cc-pCVQZ",
            AugCcPcv5z => "aug-cc-pCV5Z",
            AugCcPcv6z => "aug-cc-pCV6Z",

            // Weighted core-valence cc-pwCVnZ
            CcPwCvdz => "cc-pwCVDZ",
            CcPwCvtz => "cc-pwCVTZ",
            CcPwCvqz => "cc-pwCVQZ",
            CcPwCv5z => "cc-pwCV5Z",

            // Augmented weighted core-valence
            AugCcPwCvdz => "aug-cc-pwCVDZ",
            AugCcPwCvtz => "aug-cc-pwCVTZ",
            AugCcPwCvqz => "aug-cc-pwCVQZ",
            AugCcPwCv5z => "aug-cc-pwCV5Z",

            // Pseudopotential cc-pVnZ-PP
            CcPvdzPp => "cc-pVDZ-PP",
            CcPvtzPp => "cc-pVTZ-PP",
            CcPvqzPp => "cc-pVQZ-PP",
            CcPv5zPp => "cc-pV5Z-PP",

            // Augmented pseudopotential
            AugCcPvdzPp => "aug-cc-pVDZ-PP",
            AugCcPvtzPp => "aug-cc-pVTZ-PP",
            AugCcPvqzPp => "aug-cc-pVQZ-PP",
            AugCcPv5zPp => "aug-cc-pV5Z-PP",

            // Core-valence pseudopotential cc-pCVnZ-PP
            CcPcvdzPp => "cc-pCVDZ-PP",
            CcPcvtzPp => "cc-pCVTZ-PP",
            CcPcvqzPp => "cc-pCVQZ-PP",
            CcPcv5zPp => "cc-pCV5Z-PP",

            // Augmented core-valence pseudopotential
            AugCcPcvdzPp => "aug-cc-pCVDZ-PP",
            AugCcPcvtzPp => "aug-cc-pCVTZ-PP",
            AugCcPcvqzPp => "aug-cc-pCVQZ-PP",
            AugCcPcv5zPp => "aug-cc-pCV5Z-PP",

            // Weighted core-valence pseudopotential cc-pwCVnZ-PP
            CcPwCvdzPp => "cc-pwCVDZ-PP",
            CcPwCvtzPp => "cc-pwCVTZ-PP",
            CcPwCvqzPp => "cc-pwCVQZ-PP",
            CcPwCv5zPp => "cc-pwCV5Z-PP",

            // Augmented weighted core-valence pseudopotential
            AugCcPwCvdzPp => "aug-cc-pwCVDZ-PP",
            AugCcPwCvtzPp => "aug-cc-pwCVTZ-PP",
            AugCcPwCvqzPp => "aug-cc-pwCVQZ-PP",
            AugCcPwCv5zPp => "aug-cc-pwCV5Z-PP",

            // Hyperfine-optimized
            AugCcPvtzJ => "aug-cc-pVTZ-J",

            // W4 theory haV sets
            HaVTPlusdZ => "haV(T+d)Z",
            HaVQPlusdZ => "haV(Q+d)Z",
            HaV5PlusdZ => "haV(5+d)Z",
            // --- End Correlation-consistent

        }
            .to_string()
    }
}

/// May be useful to organize basis sets.
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub enum BasisSetCategory {
    #[default]
    Pople,
    Ahlrich,
    KarlseruhDef2,
    KarlseruhDhs,
    CorrelationConsistent,
}

impl BasisSetCategory {
    pub fn get_sets(self) -> Vec<BasisSet> {
        match self {
            Self::Pople => vec![
                None,
                // --- Pople start
                Sto3G,
                B3_21G,
                B3_21GSP,
                B4_22GSP,
                B6_31G,
                B6_31GStar,
                M6_31G,
                M6_31GStar,
                B6_31GStarStar,
                B6_31G_d,
                B6_31G_d_p,
                B6_31G_2d,
                B6_31G_2d_p,
                B6_31G_2d_2p,
                B6_31G_2df,
                B6_31G_2df_2p,
                B6_31G_2df_2pd,
                B6_31PlusGStar,
                B6_31PlusGStarStar,
                B6_31PlusG_d,
                B6_31PlusG_d_p,
                B6_31PlusG_2d,
                B6_31PlusG_2d_p,
                B6_31PlusG_2d_2p,
                B6_31PlusG_2df,
                B6_31PlusG_2df_2p,
                B6_31PlusG_2df_2pd,
                B6_31PlusPlusGStarStar,
                B6_31PlusPlusG_d_p,
                B6_31PlusPlusG_2d_p,
                B6_31PlusPlusG_2d_2p,
                B6_31PlusPlusG_2df_2p,
                B6_31PlusPlusG_2df_2pd,
                B6_31PlusPlusG_3df_3pd,
                B6_311G,
                B6_311GStar,
                B6_311GStarStar,
                B6_311G_d,
                B6_311G_d_p,
                B6_311G_2d,
                B6_311G_2d_p,
                B6_311G_2d_2p,
                B6_311G_2df,
                B6_311G_2df_2p,
                B6_311G_2df_2pd,
                B6_311G_3df,
                B6_311G_3df_3pd,
                B6_311PlusGStar,
                B6_311PlusGStarStar,
                B6_311PlusG_d,
                B6_311PlusG_d_p,
                B6_311PlusG_2d,
                B6_311PlusG_2d_p,
                B6_311PlusG_2d_2p,
                B6_311PlusG_2df,
                B6_311PlusG_2df_2p,
                B6_311PlusG_2df_2pd,
                B6_311PlusG_3df,
                B6_311PlusG_3df_2p,
                B6_311PlusG_3df_3pd,
                B6_311PlusPlusGStarStar,
                B6_311PlusPlusG_d_p,
                B6_311PlusPlusG_2d_p,
                B6_311PlusPlusG_2d_2p,
                B6_311PlusPlusG_2df_2p,
                B6_311PlusPlusG_2df_2pd,
                B6_311PlusPlusG_3df_3pd,
            ],
            Self::Ahlrich => vec![
                DefSvP,
                DefSvp,
                DefTzvp,
                DefTzvpp,
                MaDefTzvp,
                Def2Svp,
                Def2Svp_,
                Def2Tzvp,
                Def2TzvpMinusF,
                Def2Tzvpp,
                Def2Qzvp,
                Def2Qzvpp,
                Def2Svpd,
                Def2Tzvpd,
                Def2Tzvppd,
                Def2Qzvpd,
                Def2Qzvppd,
                MaDef2Svp,
                MaDef2SvP,
                MaDef2MSvp,
                MaDef2Tzvp,
                MaDef2TzvpMinusF,
                MaDef2Tzvpp,
                MaDef2Qzvpp,
            ],
            Self::KarlseruhDef2 => vec![
                Def2Svp,
                Def2Svp_,
                Def2Tzvp,
                Def2TzvpMinusF,
                Def2Tzvpp,
                Def2Qzvp,
                Def2Qzvpp,

                // Diffuse def2- sets
                Def2Svpd,
                Def2Tzvpd,
                Def2Tzvppd,
                Def2Qzvpd,
                Def2Qzvppd,

                // Minimally augmented ma-def2- sets
                MaDef2Svp,
                MaDef2SvP,
                MaDef2MSvp,
                MaDef2Tzvp,
                MaDef2TzvpMinusF,
                MaDef2Tzvpp,
                MaDef2Qzvpp,
            ],
            Self::KarlseruhDhf => vec![
                // Karlsruhe dhf- family
                DhfSvp_,
                DhfSvp,
                DhfTzvp,
                DhfTzvpp,
                DhfQzvp,
                DhfQzvpp,

                // Karlsruhe dhf- two-component variants
                DhfSvp2c,
                DhfTzvp2c,
                DhfTzvpp2c,
                DhfQzvp2c,
                DhfQzvpp2c,
            ],
            Self::CorrelationConsistent => vec![
                CcPvdz,
                CcPvtz,
                CcPvqz,
                CcPv5z,
                CcPv6z,

                // Augmented aug-cc-pVnZ
                AugCcPvdz,
                AugCcPvtz,
                AugCcPvqz,
                AugCcPv5z,
                AugCcPv6z,

                // Tight-d variants cc-pVn(+d)Z
                CcPvdPlusdZ,
                CcPvtPlusdZ,
                CcPvqPlusdZ,
                CcPv5PlusdZ,

                // Tight-d augmented aug-cc-pVn(+d)Z
                AugCcPvdPlusdZ,
                AugCcPvtPlusdZ,
                AugCcPvqPlusdZ,
                AugCcPv5PlusdZ,
                AugCcPv6PlusdZ,

                // Partially augmented Truhlar sets
                AprCcPvQPlusdZ,
                MayCcPvTPlusdZ,
                MayCcPvQPlusdZ,
                JunCcPvDPlusdZ,
                JunCcPvTPlusdZ,
                JunCcPvQPlusdZ,
                JulCcPvDPlusdZ,
                JulCcPvTPlusdZ,
                JulCcPvQPlusdZ,
                MaugCcPvDPlusdZ,
                MaugCcPvTPlusdZ,
                MaugCcPvQPlusdZ,

                // Core-valence cc-pCVnZ
                CcPcvdz,
                CcPcvtz,
                CcPcvqz,
                CcPcv5z,
                CcPcv6z,

                // Augmented core-valence aug-cc-pCVnZ
                AugCcPcvdz,
                AugCcPcvtz,
                AugCcPcvqz,
                AugCcPcv5z,
                AugCcPcv6z,

                // Weighted core-valence cc-pwCVnZ
                CcPwCvdz,
                CcPwCvtz,
                CcPwCvqz,
                CcPwCv5z,

                // Augmented weighted core-valence aug-cc-pwCVnZ
                AugCcPwCvdz,
                AugCcPwCvtz,
                AugCcPwCvqz,
                AugCcPwCv5z,

                // Pseudopotential cc-pVnZ-PP
                CcPvdzPp,
                CcPvtzPp,
                CcPvqzPp,
                CcPv5zPp,

                // Augmented pseudopotential aug-cc-pVnZ-PP
                AugCcPvdzPp,
                AugCcPvtzPp,
                AugCcPvqzPp,
                AugCcPv5zPp,

                // Core-valence pseudopotential cc-pCVnZ-PP
                CcPcvdzPp,
                CcPcvtzPp,
                CcPcvqzPp,
                CcPcv5zPp,

                // Augmented core-valence pseudopotential aug-cc-pCVnZ-PP
                AugCcPcvdzPp,
                AugCcPcvtzPp,
                AugCcPcvqzPp,
                AugCcPcv5zPp,

                // Weighted core-valence pseudopotential cc-pwCVnZ-PP
                CcPwCvdzPp,
                CcPwCvtzPp,
                CcPwCvqzPp,
                CcPwCv5zPp,

                // Augmented weighted core-valence pseudopotential aug-cc-pwCVnZ-PP
                AugCcPwCvdzPp,
                AugCcPwCvtzPp,
                AugCcPwCvqzPp,
                AugCcPwCv5zPp,

                // Hyperfine-optimized
                AugCcPvtzJ,

                // W4 theory haV sets
                HaVTPlusdZ,
                HaVQPlusdZ,
                HaV5PlusdZ,
            ]
        }
    }
}

impl Display for BasisSetCategory {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let v = match self {
            Self::Pople => "Pople",
            Self::Ahlrich => "Ahlrich",
            Self::KarlseruhDef2 => "KarlseruhDef2",
            Self::KarlseruhDhs => "karlseruhDhs",
            Self::CorrelationConsistent => "CorrelationConsistent",
        };
        write!(f, "{v}")
    }
}