use nalgebra::dvector;

use crate::{
    linearize::linearize,
    math::modal_analysis::{get_6dof_lateral, get_6dof_longitudal},
    model::{F16, dynamicmodel::fixedwing6dof::FixedWing6DoF},
    trim::create_trim_problem_and_trim,
};

/// Example from section 3.8 from "Aircraft Control and Simulation - 3rd edition" by Brian L. Stevens, Frank L. Lewis and Eric N. Johnson.
pub fn f16_6dof_dynamic_behavior_example() -> Result<(), Box<dyn std::error::Error>> {
    let setpoints = dvector![
        153.01, // vt [m/s] — 502 ft/s
        0.0,    // altitude [m]
        0.0,    // gamma [deg]
        0.0,    // roll rate [rad/s]
        0.0,    // pitch rate [rad/s]
        0.0,    // turn rate [rad/s]
        0.0,    // phi [rad]
        0.0,    // coordinated turn flag
        0.30,   // Xcg
    ];
    let init_params = dvector![
        0.15,  // throttle
        -2.0,  // elevator
        0.039, // alpha
        0.0,   // aileron
        0.0,   // rudder
        0.0,   // beta
    ];

    let (x, u, _cost) =
        create_trim_problem_and_trim(F16::new(), FixedWing6DoF, setpoints, init_params)?;

    let linearized = linearize(F16::new(), FixedWing6DoF, x, u)?;

    let longitudal = get_6dof_longitudal(linearized.a());
    let lateral = get_6dof_lateral(linearized.a());

    eprintln!("Longitudal matrix:\n{:.9}", longitudal);
    eprintln!("Lateral matrix:\n{:.9}", lateral);

    Ok(())
}
