use nalgebra::{DVector, dvector};

use crate::{
    math::{SizedVector, input_functions::doublet, timestep::TimeStep},
    model::{
        RAD_TO_DEG, Transport,
        dynamicmodel::fixedwing3dof::{FixedWing3DoF, FixedWing3DoFInput, FixedWing3DoFStates},
    },
    plots::state_variables_plot,
    sim::simulator::create_sim_and_run,
    trim::create_trim_problem_and_trim,
};

/// Example of an elevator doublet input applied to trimmed transport aircraft 3DoF-model.
pub fn transport_3dof_elevator_doublet_example() -> Result<(), Box<dyn std::error::Error>> {
    let (x, u, _cost) = create_trim_problem_and_trim(
        Transport::new(),
        FixedWing3DoF,
        dvector![76.2, 0.0, 0.0],  // setpoints: [vt, altitude, gamma]
        dvector![0.1, -10.0, 0.1], // initial params: [throttle, elevator, alpha]
    )?;

    let elevator_doublet = |params: &DVector<f64>, time: f64| -> FixedWing3DoFInput {
        FixedWing3DoFInput::from_vector(dvector![
            params[0],
            doublet(params[1], 1.0, 0.5, 2.0 / RAD_TO_DEG, time),
            params[2],
            params[3]
        ])
    };

    let simulator = create_sim_and_run(
        Transport::new(),
        FixedWing3DoF,
        x,
        u,
        Some(elevator_doublet),
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
        "Elevator doublet response".into(),
    )
}
