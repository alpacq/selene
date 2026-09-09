use nalgebra::{DVector, dvector};

use crate::{
    math::{SizedVector, input_functions::doublet, timestep::TimeStep},
    model::{
        Transport,
        dynamicmodel::fixedwing3dof::{FixedWing3DoF, FixedWing3DoFInput, FixedWing3DoFStates},
    },
    plots::state_variables_plot,
    sim::simulator::create_sim_and_run,
    trim::create_trim_problem_and_trim,
};

/// Example of a throttle doublet input applied to trimmed transport aircraft 3DoF-model.
/// Example 3.6-4 from "Aircraft Control and Simulation - 3rd edition" by Brian L. Stevens, Frank L. Lewis and Eric N. Johnson.
pub fn transport_3dof_throttle_doublet_example() -> Result<(), Box<dyn std::error::Error>> {
    let (x, u, _cost) = create_trim_problem_and_trim(
        Transport::new(),
        FixedWing3DoF,
        dvector![76.2, 0.0, 0.0],  // setpoints: [vt, altitude, gamma]
        dvector![0.1, -10.0, 0.1], // initial params: [throttle, elevator, alpha]
    )?;

    let throttle_doublet = |params: &DVector<f64>, time: f64| -> FixedWing3DoFInput {
        FixedWing3DoFInput::from_vector(dvector![
            doublet(params[0], 0.0, 3.0, 0.1, time),
            params[1],
            params[2],
            params[3]
        ])
    };

    let simulator = create_sim_and_run(
        Transport::new(),
        FixedWing3DoF,
        x,
        u,
        Some(throttle_doublet),
        60.0,
        TimeStep::new(0.001),
    );

    state_variables_plot(
        vec![
            FixedWing3DoFStates::Alpha as usize,
            FixedWing3DoFStates::Theta as usize,
        ],
        simulator.output,
        vec!["alpha(t)".into(), "theta(t)".into()],
        "time (s)".into(),
        "angle (rad)".into(),
        "Throttle doublet response".into(),
    )
}
