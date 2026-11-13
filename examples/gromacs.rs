use std::{io, path::Path};

use bio_apis::pubchem::StructureSearchNamespace;
use bio_files::{
    Sdf,
    gromacs::{GromacsInput, MdpParams, MoleculeInput, OutputControl, mdp::Integrator},
    md_params::ForceFieldParams,
};

fn main() -> io::Result<()> {
    // The default implementation for each field will use GROMACS defaults; override
    // ones we wish to change.
    let mdp = MdpParams {
        integrator: Integrator::MdVv,
        nsteps: 1_000,
        dt: 0.001,
        output_control: OutputControl {
            nstxout: Some(100),
            nstenergy: Some(500),
            ..Default::default()
        },
        ..Default::default()
    };

    let molecules = {
        let mol = Sdf::load_pubchem(StructureSearchNamespace::Cid, "2244")?;

        // todo: Demonstrate how to load mol specific params, e.g. with antechamber

        vec![MoleculeInput {
            name: mol.ident.clone(),
            atoms: mol.atoms.clone(),
            bonds: mol.bonds.clone(),
            // ff_params: Some(ff_params_mol_specific),
            ff_params: None,
            count: 1,
            copy_atom_posits: None,
        }]
    };

    let params = ForceFieldParams::load_dat(Path::new("gaff2.dat"))?;

    let input = GromacsInput {
        mdp,
        molecules,
        box_nm: Some((3., 3., 3.)),
        ff_global: Some(params),
        ..Default::default()
    };

    // Writes to a new `md_out/gromacs_{timestamp}/` folder.
    let output = input.run()?;

    println!("Gromacs output log: \n\n{}", output.log_text);

    for frame in &output.trajectory {
        println!("Frame time: {:.3}", frame.time);
        println!("Energy: {:?}", frame.energy);

        // Do something with `frame.atom_posits` and `frame.atom_velocities` as required.
    }

    Ok(())
}
