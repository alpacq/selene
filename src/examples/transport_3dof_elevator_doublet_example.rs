use nalgebra::{DVector, dvector};

use crate::{
    math::{SizedVector, input_functions::doublet, timestep::TimeStep},
    model::{
        RAD_TO_DEG, Transport,
        dynamicmodel::fixedwing3dof::{FixedWing3DoF, FixedWing3DoFInput, FixedWing3DoFStates},
    },
    plots::state_variables_plot,
    sim::simulator::SimulatorBuilder,
    trim::TrimProblemBuilder,
};

pub fn transport_3dof_elevator_doublet_example() -> Result<(), Box<dyn std::error::Error>> {
    let setpoints = dvector![76.2, 0.0, 0.0]; // setpoints: [vt, altitude, gamma]
    let init_params = dvector![0.1, -10.0, 0.1]; // initial params: [throttle, elevator, alpha]

    let system = Transport::new();
    let model = FixedWing3DoF;

    let problem = TrimProblemBuilder::new()
        .for_system(system)
        .with_model(model)
        .with_setpoints(setpoints)
        .with_initial_params(init_params)
        .build();
    let (x, u, _cost) = problem.trim()?;

    let system = Transport::new();
    let model = FixedWing3DoF;

    let mut simulator = SimulatorBuilder::new()
        .for_system(system)
        .with_model(model)
        .with_state(x)
        .build();

    let elevator_doublet = |params: &DVector<f64>, time: f64| -> FixedWing3DoFInput {
        FixedWing3DoFInput::from_vector(dvector![
            params[0],
            doublet(params[1], 1.0, 0.5, 2.0 / RAD_TO_DEG, time),
            params[2],
            params[3]
        ])
    };

    simulator.run(u, Some(elevator_doublet), 60.0, TimeStep::new(0.001));

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
