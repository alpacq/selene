use nalgebra::dvector;

use crate::{
    linearize::LinearizationProblemBuilder,
    model::{Transport, dynamicmodel::fixedwing3dof::FixedWing3DoF},
    trim::TrimProblemBuilder,
};

pub fn transport_3dof_linearization_example() -> Result<(), Box<dyn std::error::Error>> {
    // initially we need to trim aircraft in given state
    let setpoints = dvector![51.816, 0.0, 0.0]; // setpoints: [vt, altitude, gamma]
    let init_params = dvector![0.1, -10.0, 0.1]; // initial params: [throttle, elevator, alpha]
    let trim_problem = TrimProblemBuilder::new()
        .for_system(Transport::new())
        .with_model(FixedWing3DoF)
        .with_setpoints(setpoints)
        .with_initial_params(init_params)
        .build();
    let (x, u, _cost) = trim_problem.trim()?;

    // then we can linearize around the trimmed state
    let lin_problem = LinearizationProblemBuilder::new()
        .for_system(Transport::new())
        .with_model(FixedWing3DoF)
        .with_trimmed_input_and_state(x, u)
        .build();

    let a_matrix = lin_problem.jacobian_a();
    let b_matrix = lin_problem.jacobian_b();

    println!("A matrix:\n{:?}", a_matrix);
    println!("B matrix:\n{:?}", b_matrix);

    Ok(())
}
