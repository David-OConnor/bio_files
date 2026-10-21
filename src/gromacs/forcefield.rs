//! Force-field conventions for self-contained GROMACS topologies. Parameters
//! supplied by the caller must belong to the selected family.

use std::fmt::Write;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ForceField {
    #[default]
    Amber,
    Charmm36,
    Opls,
}

/// Rigid water and neutralizing Na/Cl parameters, in nm, kJ/mol, amu and e.
/// A zero O-M distance represents a three-site model.
#[derive(Clone, Debug)]
pub struct WaterModel {
    pub mass_o: f64,
    pub mass_h: f64,
    pub oh: f64,
    pub om: f64,
    pub angle: f64,
    pub q_h: f64,
    pub sigma_o: f64,
    pub epsilon_o: f64,
    pub sigma_h: f64,
    pub epsilon_h: f64,
    /// (mass, sigma, epsilon).
    pub sodium: (f64, f64, f64),
    pub chloride: (f64, f64, f64),
}

#[derive(Clone, Debug)]
pub struct TopologyConfig {
    pub force_field: ForceField,
    /// Override the family's default water and ions.
    pub water_model: Option<WaterModel>,
}

impl Default for TopologyConfig {
    fn default() -> Self {
        Self {
            force_field: ForceField::Amber,
            water_model: None,
        }
    }
}

impl TopologyConfig {
    pub fn for_family(force_field: ForceField) -> Self {
        Self {
            force_field,
            water_model: None,
        }
    }

    /// GROMACS combination rule, fudgeLJ, fudgeQQ. The exact Amber ratio
    /// matters for agreement with native calculations.
    pub fn defaults(&self) -> (u8, f64, f64) {
        match self.force_field {
            ForceField::Amber => (2, 0.5, 1. / 1.2),
            ForceField::Charmm36 => (2, 1., 1.),
            ForceField::Opls => (3, 0.5, 0.5),
        }
    }

    pub fn water(&self) -> WaterModel {
        if let Some(w) = &self.water_model {
            return w.clone();
        }
        let sigma_factor = 2. / 2_f64.powf(1. / 6.);
        match self.force_field {
            ForceField::Amber => WaterModel {
                mass_o: 16.,
                mass_h: 1.008,
                oh: 0.087243313,
                om: 0.015939833,
                angle: 103.6_f64.to_radians(),
                q_h: 0.6791,
                sigma_o: 1.777167268 * sigma_factor / 10.,
                epsilon_o: 0.212800813 * 4.184,
                sigma_h: 0.,
                epsilon_h: 0.,
                sodium: (22.99, 0.2439, 0.1065 * 4.184),
                chloride: (35.45, 0.4478, 0.0073 * 4.184),
            },
            ForceField::Charmm36 => WaterModel {
                mass_o: 15.9994,
                mass_h: 1.008,
                oh: 0.09572,
                om: 0.,
                angle: 104.52_f64.to_radians(),
                q_h: 0.417,
                sigma_o: 1.7682 * sigma_factor / 10.,
                epsilon_o: 0.1521 * 4.184,
                sigma_h: 0.2245 * sigma_factor / 10.,
                epsilon_h: 0.046 * 4.184,
                sodium: (22.98977, 1.41075 * sigma_factor / 10., 0.0469 * 4.184),
                chloride: (35.45, 2.27 * sigma_factor / 10., 0.150 * 4.184),
            },
            ForceField::Opls => WaterModel {
                mass_o: 15.9994,
                mass_h: 1.008,
                oh: 0.09572,
                om: 0.015,
                angle: 104.52_f64.to_radians(),
                q_h: 0.52,
                sigma_o: 0.315365,
                epsilon_o: 0.64852,
                sigma_h: 0.,
                epsilon_h: 0.,
                sodium: (22.98977, 0.333045, 0.011598),
                chloride: (35.453, 0.441724, 0.492833),
            },
        }
    }

    pub fn defaults_text(&self) -> String {
        let (rule, lj, qq) = self.defaults();
        format!(
            "[ defaults ]\n; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n1 {rule} yes {lj:.12} {qq:.12}\n\n"
        )
    }

    /// Atom types for the water/ion blocks. Types are private to this exporter.
    pub fn solvent_atomtypes(&self, include_water: bool) -> String {
        let w = self.water();
        let mut s = String::new();
        let types = [
            ("BF_OW", 8, w.mass_o, w.sigma_o, w.epsilon_o),
            ("BF_HW", 1, w.mass_h, w.sigma_h, w.epsilon_h),
            ("BF_MW", 0, 0., 0., 0.),
            ("BF_NA", 11, w.sodium.0, w.sodium.1, w.sodium.2),
            ("BF_CL", 17, w.chloride.0, w.chloride.1, w.chloride.2),
        ];
        for (i, (name, z, mass, sigma, epsilon)) in types.into_iter().enumerate() {
            if i < 3 && (!include_water || (i == 2 && w.om == 0.)) {
                continue;
            }
            writeln!(s, "{name} {z} {mass:.8} 0 A {sigma:.12e} {epsilon:.12e}").unwrap();
        }
        s
    }

    pub fn solvent_moleculetypes(&self, include_water: bool) -> String {
        let w = self.water();
        let mut s = String::new();
        for (res, kind, q, m) in [
            ("NA", "BF_NA", 1., w.sodium.0),
            ("CL", "BF_CL", -1., w.chloride.0),
        ] {
            writeln!(
                s,
                "\n[ moleculetype ]\n{res} 1\n[ atoms ]\n1 {kind} 1 {res} {res} 1 {q} {m}"
            )
            .unwrap();
        }
        if !include_water {
            return s;
        }
        writeln!(s, "\n[ moleculetype ]\nSOL 2\n[ atoms ]").unwrap();
        writeln!(
            s,
            "1 BF_OW 1 SOL OW 1 {:.12} {:.8}",
            if w.om == 0. { -2. * w.q_h } else { 0. },
            w.mass_o
        )
        .unwrap();
        writeln!(s, "2 BF_HW 1 SOL HW1 1 {:.12} {:.8}", w.q_h, w.mass_h).unwrap();
        writeln!(s, "3 BF_HW 1 SOL HW2 1 {:.12} {:.8}", w.q_h, w.mass_h).unwrap();
        if w.om != 0. {
            writeln!(s, "4 BF_MW 1 SOL MW 1 {:.12} 0", -2. * w.q_h).unwrap();
        }
        writeln!(
            s,
            "[ settles ]\n1 1 {:.12} {:.12}",
            w.oh,
            2. * w.oh * (w.angle / 2.).sin()
        )
        .unwrap();
        let n = if w.om == 0. { 3 } else { 4 };
        s.push_str("[ exclusions ]\n");
        for i in 1..=n {
            write!(s, "{i}").unwrap();
            for j in 1..=n {
                if i != j {
                    write!(s, " {j}").unwrap();
                }
            }
            s.push('\n');
        }
        if w.om != 0. {
            let a = w.om / (2. * w.oh * (w.angle / 2.).cos());
            writeln!(s, "[ virtual_sites3 ]\n4 1 2 3 1 {a:.12} {a:.12}").unwrap();
        }
        s.push('\n');
        s
    }
}
