use nalgebra::dvector;

use crate::{
    linearize::LinearizationProblemBuilder,
    model::{F16, dynamicmodel::fixedwing6dof::FixedWing6DoF},
    trim::TrimProblemBuilder,
};

pub fn f16_6dof_linearization_example_pullup() -> Result<(), Box<dyn std::error::Error>> {
    let pullup_setpoints = dvector![152.4, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0];
    let pullup_init_params = dvector![0.2, 1.0, 0.02, 0.0, 0.0, 0.0];
    let pullup_trim_problem = TrimProblemBuilder::new()
        .for_system(F16::new())
        .with_model(FixedWing6DoF)
        .with_setpoints(pullup_setpoints)
        .with_initial_params(pullup_init_params)
        .build();
    let (pullup_x, pullup_u, _cost) = pullup_trim_problem.trim()?;

    let pullup_lin_problem = LinearizationProblemBuilder::new()
        .for_system(F16::new())
        .with_model(FixedWing6DoF)
        .with_trimmed_input_and_state(pullup_x, pullup_u)
        .build();

    eprintln!("Wings-level pull-up:");
    eprintln!("A matrix:\n{:.9}", pullup_lin_problem.jacobian_a());
    eprintln!("B matrix:\n{:.9}", pullup_lin_problem.jacobian_b());

    Ok(())
}

pub fn f16_6dof_linearization_example_turn() -> Result<(), Box<dyn std::error::Error>> {
    let turn_setpoints = dvector![152.4, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 1.0];
    let turn_init_params = dvector![0.85, -6.0, 0.24, 0.1, -0.4, 0.0];
    let turn_trim_problem = TrimProblemBuilder::new()
        .for_system(F16::new())
        .with_model(FixedWing6DoF)
        .with_setpoints(turn_setpoints)
        .with_initial_params(turn_init_params)
        .build();
    let (turn_x, turn_u, _cost) = turn_trim_problem.trim()?;

    let turn_lin_problem = LinearizationProblemBuilder::new()
        .for_system(F16::new())
        .with_model(FixedWing6DoF)
        .with_trimmed_input_and_state(turn_x, turn_u)
        .build();

    eprintln!("Coordinated turn:");
    eprintln!("A matrix:\n{:.9}", turn_lin_problem.jacobian_a());
    eprintln!("B matrix:\n{:.9}", turn_lin_problem.jacobian_b());

    Ok(())
}
