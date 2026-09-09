use nalgebra::dvector;

use crate::{
    linearize::linearize,
    model::{F16, dynamicmodel::fixedwing6dof::FixedWing6DoF},
    trim::create_trim_problem_and_trim,
};

/// Example 3.7-2 from "Aircraft Control and Simulation - 3rd edition" by Brian L. Stevens, Frank L. Lewis and Eric N. Johnson.
pub fn f16_6dof_linearization_example_pullup() -> Result<(), Box<dyn std::error::Error>> {
    let pullup_setpoints = dvector![152.4, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0];
    let pullup_init_params = dvector![0.2, 1.0, 0.02, 0.0, 0.0, 0.0];
    let (pullup_x, pullup_u, _cost) = create_trim_problem_and_trim(
        F16::new(),
        FixedWing6DoF,
        pullup_setpoints,
        pullup_init_params,
    )?;

    let linearized = linearize(F16::new(), FixedWing6DoF, pullup_x, pullup_u)?;

    eprintln!("Wings-level pull-up:");
    eprintln!("A matrix:\n{:.9}", linearized.a());
    eprintln!("B matrix:\n{:.9}", linearized.b());

    Ok(())
}

/// Example 3.7-3 from "Aircraft Control and Simulation - 3rd edition" by Brian L. Stevens, Frank L. Lewis and Eric N. Johnson.
pub fn f16_6dof_linearization_example_turn() -> Result<(), Box<dyn std::error::Error>> {
    let turn_setpoints = dvector![152.4, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 1.0];
    let turn_init_params = dvector![0.85, -6.0, 0.24, 0.1, -0.4, 0.0];
    let (turn_x, turn_u, _cost) =
        create_trim_problem_and_trim(F16::new(), FixedWing6DoF, turn_setpoints, turn_init_params)?;

    let linearized = linearize(F16::new(), FixedWing6DoF, turn_x, turn_u)?;

    eprintln!("Coordinated turn:");
    eprintln!("A matrix:\n{:.9}", linearized.a());
    eprintln!("B matrix:\n{:.9}", linearized.b());

    Ok(())
}
