use nalgebra::dvector;

use crate::{
    linearize::linearize,
    model::{Transport, dynamicmodel::fixedwing3dof::FixedWing3DoF},
    trim::create_trim_problem_and_trim,
};

/// Example 3.7-1 from "Aircraft Control and Simulation - 3rd edition" by Brian L. Stevens, Frank L. Lewis and Eric N. Johnson.
pub fn transport_3dof_linearization_example() -> Result<(), Box<dyn std::error::Error>> {
    // initially we need to trim aircraft in given state
    let (x, u, _cost) = create_trim_problem_and_trim(
        Transport::new(),
        FixedWing3DoF,
        dvector![60.96, 0.0, 15.0], // setpoints: [vt, altitude, gamma]
        dvector![1.0, 10.0, 14.0],  // initial params: [throttle, elevator, alpha]
    )?;

    // then we can linearize around the trimmed state
    let linearized = linearize(Transport::new(), FixedWing3DoF, x, u)?;

    eprintln!("A matrix:\n{:.9}", linearized.a());
    eprintln!("B matrix:\n{:.9}", linearized.b());

    Ok(())
}
