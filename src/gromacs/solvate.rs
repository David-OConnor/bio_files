//! [Docs](https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html)
//! [More docs](https://manual.gromacs.org/current/how-to/topology.html#water-solvation)

use std::{fmt::Write as _, io, path::Path};

use lin_alg::f32::Vec3;

use crate::gromacs::{
    CUSTOM_SOLVENT_GRO, GMX_BUILTIN_TIP4P_WATER, MoleculeInput,
    gro::{AtomGro, Gro},
    run_gmx,
};

pub(in crate::gromacs) const CUSTOM_REGIONS_SOLVENT_GRO: &str = "water_custom_regions.gro";
const CUSTOM_REGIONS_FULL_BOX_GRO: &str = "water_custom_regions_full_box.gro";

/// Water model to use for explicit solvation via `gmx solvate`.
/// [solvate docs](https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html#gmx-solvate)
///
/// Selecting a model causes the pipeline to:
/// 1. Include the model's atom-type and molecule-type sections in the topology.
/// 2. Run `gmx solvate` to fill the simulation box with water before `grompp`.
#[derive(Clone, Debug, Default)]
pub enum Solvent {
    #[default]
    /// Fill the entire sim box with rigid water molecules, at a realistic density.
    ///
    /// OPC 4-point rigid water (recommended with Amber ff19SB).
    /// Requires GROMACS 2022+ which ships `tip4p.gro` and `opc.itp` in its data directory.
    /// We use Gromac's TIP4P water molecule, which is suitable for some
    /// 4-point rigid water models including OPC.
    WaterOpc,
    /// Fill sub-regions of the full sim box with rigid water molecules at a realistic density.
    ///
    /// Bounds use the GROMACS box coordinate frame and units: nm in `[0, box_nm]`.
    /// Each region must have positive dimensions and fit inside the full simulation box.
    WaterOpcCustomRegions(Vec<(Vec3, Vec3)>), // (bounds low, bounds high)
    Custom(CustomSolventTemplate),
}

impl Solvent {
    /// Pre-equilibrated water-box filename written by `save()` and passed to `gmx solvate -cs`.
    pub(in crate::gromacs) fn gro_filename(&self) -> &'static str {
        match self {
            // 4-site model is fine for OPC; we just use a different force field.
            Self::WaterOpc => GMX_BUILTIN_TIP4P_WATER,
            Self::WaterOpcCustomRegions(_) => GMX_BUILTIN_TIP4P_WATER,
            Self::Custom(_) => CUSTOM_SOLVENT_GRO,
        }
    }

    pub(in crate::gromacs) fn prepare_gro(
        &self,
        dir: &Path,
        box_nm: Option<(f64, f64, f64)>,
    ) -> io::Result<Option<&'static str>> {
        match self {
            Self::WaterOpcCustomRegions(regions) => {
                make_water_opc_custom_regions_gro(dir, box_nm, regions)
            }
            _ => Ok(Some(self.gro_filename())),
        }
    }
}

/// Ask GROMACS to generate a full periodic OPC-compatible water box, retain only
/// molecules wholly inside one requested region, then use that as the solvent
/// template for the final `gmx solvate -cp ...` pass.
///
/// Starting from one full-cell template preserves the periodic packing phase
/// across region boundaries. The final `gmx solvate` invocation remains
/// responsible for removing waters that clash with the solute and updating the
/// topology molecule count.
fn make_water_opc_custom_regions_gro(
    dir: &Path,
    box_nm: Option<(f64, f64, f64)>,
    regions: &[(Vec3, Vec3)],
) -> io::Result<Option<&'static str>> {
    if regions.is_empty() {
        return Ok(None);
    }

    let box_vec = validate_regions(box_nm, regions)?;
    let x = box_vec.x.to_string();
    let y = box_vec.y.to_string();
    let z = box_vec.z.to_string();

    run_gmx(
        dir,
        &[
            "solvate",
            "-cs",
            GMX_BUILTIN_TIP4P_WATER,
            "-box",
            &x,
            &y,
            &z,
            "-o",
            CUSTOM_REGIONS_FULL_BOX_GRO,
        ],
    )?;

    let full_box = Gro::load(&dir.join(CUSTOM_REGIONS_FULL_BOX_GRO))?;
    let regional_box = retain_waters_in_regions(full_box, regions);
    regional_box.save(&dir.join(CUSTOM_REGIONS_SOLVENT_GRO))?;

    Ok(Some(CUSTOM_REGIONS_SOLVENT_GRO))
}

fn validate_regions(
    box_nm: Option<(f64, f64, f64)>,
    regions: &[(Vec3, Vec3)],
) -> io::Result<lin_alg::f64::Vec3> {
    let Some((x, y, z)) = box_nm else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "WaterOpcCustomRegions requires box_nm.",
        ));
    };

    let box_vec = lin_alg::f64::Vec3::new(x, y, z);
    if !is_positive_finite(box_vec) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Simulation box must have positive finite dimensions; got {box_vec:?}."),
        ));
    }

    for (i, &(low, high)) in regions.iter().enumerate() {
        if !is_finite(low)
            || !is_finite(high)
            || low.x < 0.0
            || low.y < 0.0
            || low.z < 0.0
            || high.x <= low.x
            || high.y <= low.y
            || high.z <= low.z
            || high.x as f64 > box_vec.x
            || high.y as f64 > box_vec.y
            || high.z as f64 > box_vec.z
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Water OPC custom region {i} must have positive finite dimensions and fit \
                     inside box {box_vec:?}; got low={low:?}, high={high:?}."
                ),
            ));
        }
    }

    Ok(box_vec)
}

fn retain_waters_in_regions(full_box: Gro, regions: &[(Vec3, Vec3)]) -> Gro {
    let mut atoms = Vec::new();
    let mut start = 0;

    while start < full_box.atoms.len() {
        let first = &full_box.atoms[start];
        let mut end = start + 1;

        while end < full_box.atoms.len()
            && full_box.atoms[end].mol_id == first.mol_id
            && full_box.atoms[end].mol_name == first.mol_name
        {
            end += 1;
        }

        let water = &full_box.atoms[start..end];
        if regions
            .iter()
            .any(|&(low, high)| water.iter().all(|atom| region_contains(low, high, atom)))
        {
            atoms.extend_from_slice(water);
        }

        start = end;
    }

    Gro {
        atoms,
        head_text: "OPC water custom regions generated by Bio Files".to_owned(),
        box_vec: full_box.box_vec,
    }
}

fn is_positive_finite(v: lin_alg::f64::Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite() && v.x > 0.0 && v.y > 0.0 && v.z > 0.0
}

fn is_finite(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

fn region_contains(low: Vec3, high: Vec3, atom: &AtomGro) -> bool {
    atom.posit.x >= low.x as f64
        && atom.posit.y >= low.y as f64
        && atom.posit.z >= low.z as f64
        && atom.posit.x <= high.x as f64
        && atom.posit.y <= high.y as f64
        && atom.posit.z <= high.z as f64
}

/// C+P from `dynamics`.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterInitTemplate {
    // velocity is o velocity, instead of 3 separate velocities
    pub o_posits: Vec<Vec3>,
    pub h0_posits: Vec<Vec3>,
    pub h1_posits: Vec<Vec3>,
    pub o_velocities: Vec<Vec3>,
    pub h0_velocities: Vec<Vec3>,
    pub h1_velocities: Vec<Vec3>,
    // todo: Cache these, or infer?
    /// This must correspond to the positions. Cached.
    pub bounds: (Vec3, Vec3),
}

impl WaterInitTemplate {
    /// Create `.gro` file text for use as the `gmx solvate -cs` template.
    ///
    /// Writes 4 atoms per molecule (OW, HW1, HW2, MW) in GROMACS `.gro` format.
    /// The MW virtual-site position is computed from the OW/HW1/HW2 positions
    /// using the OPC `[virtual_sites3]` coefficients (a = b = 0.147803).
    ///
    /// Positions are converted Å → nm. Velocities (Å/ps → nm/ps) are included
    /// so that the template is also usable as a pre-equilibrated starting frame.
    pub fn to_gro(&self) -> String {
        // OPC virtual-site coefficients: r_MW = r_OW + a*(r_HW1-r_OW) + a*(r_HW2-r_OW)
        const VS_A: f32 = 0.147803;

        let n_mols = self.o_posits.len();
        let n_atoms = n_mols * 4; // OW + HW1 + HW2 + MW

        let mut s = String::from("OPC water template — generated by Bio Files\n");
        let _ = writeln!(s, "{}", n_atoms);

        let mut atom_serial = 1usize;

        for i in 0..n_mols {
            let res = (i + 1) % 100_000;

            let ow = self.o_posits[i];
            let hw1 = self.h0_posits[i];
            let hw2 = self.h1_posits[i];
            let mw = ow + (hw1 - ow) * VS_A + (hw2 - ow) * VS_A;

            let v_ow = self.o_velocities[i];
            let v_hw1 = self.h0_velocities[i];
            let v_hw2 = self.h1_velocities[i];

            for (name, pos, vel) in [
                ("OW", ow, v_ow),
                ("HW1", hw1, v_hw1),
                ("HW2", hw2, v_hw2),
                ("MW", mw, Vec3::new_zero()), // virtual site — no independent velocity
            ] {
                let _ = writeln!(
                    s,
                    "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}",
                    res,
                    "SOL",
                    name,
                    atom_serial % 100_000,
                    pos.x / 10.0, // Å → nm
                    pos.y / 10.0,
                    pos.z / 10.0,
                    vel.x / 10.0, // Å/ps → nm/ps
                    vel.y / 10.0,
                    vel.z / 10.0,
                );
                atom_serial += 1;
            }
        }

        // Box vector line (Å → nm)
        let lo = self.bounds.0;
        let hi = self.bounds.1;
        let _ = writeln!(
            s,
            "{:>10.5}{:>10.5}{:>10.5}",
            (hi.x - lo.x).abs() / 10.0,
            (hi.y - lo.y).abs() / 10.0,
            (hi.z - lo.z).abs() / 10.0,
        );

        s
    }
}

/// Generic custom solvent-box template passed to `gmx solvate -cs`.
///
/// `gro_text` is written verbatim to `solvent.gro`. Any residue names present in that
/// coordinate file that are not built into the selected force field must also be described
/// in `topology_molecules` so the generated topology remains self-contained after
/// `gmx solvate -p` appends the molecule counts.
#[derive(Clone, Debug)]
pub struct CustomSolventTemplate {
    pub gro_text: String,
    pub topology_molecules: Vec<MoleculeInput>,
    /// Include the Amber OPC water molecule type (`SOL`) in the topology.
    pub include_opc_water: bool,
}

impl CustomSolventTemplate {
    pub fn from_water(template: WaterInitTemplate) -> Self {
        Self {
            gro_text: template.to_gro(),
            topology_molecules: Vec::new(),
            include_opc_water: false,
        }
    }
}
