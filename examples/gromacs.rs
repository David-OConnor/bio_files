use bio_files::{
    ForceFieldParams,
    gromacs::{GromacsInput, MdpParams, MoleculeInput, OutputControl},
};

fn main() {
    let mols = ();

    // The default implementation for each field will use GROMACS defaults; override
    // ones we wish to change.
    let mdp = MdpParams {
        integrator: Integrator::MdVv,
        nsteps: 1_000,
        dt: 0.001,
        output_control: OutputControl {
            txout: Some(100),
            nstenergy: Some(500),
            ..Default::default()
        }..Default::default(),
    };

    let molecules = {
        let mol = Sdf::load_pubchem(2244).unwrap();

        // todo: Demonstrate how to load mol specific params, e.g. with antechamber

        vec![MoleculeInput {
            name: &mol.ident,
            atoms: mol.atoms.clone(),
            bonds: mol.bonds.clone(),
            // ff_params: Some(ff_params_mol_specific),
            ff_params: None,
            count: 0,
        }]
    };

    let p = Path::new("gaff2.dat");
    let params = ForceFieldParams::load_dat(p).unwrap();

    let input = GromacsInput {
        mdp,
        molecules,
        box_nm: Some((30., 30., 30.)),
        ff_global: Some(params),
        water_model,
    };

    let output = input.run()?;

    println!("Gromacs output log: \n\n{}", output.log_text);

    for frame in &output.trajectory {
        println!("Frame time: {:?.3}", frame.time);

        println!("Energy: {:?}", frame.energy);

        // Do something with `frame.atom_posits` and `frame.atom_velocities` as required.
    }
}
