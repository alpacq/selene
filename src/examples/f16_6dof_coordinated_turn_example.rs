use nalgebra::dvector;

use crate::{
    math::timestep::TimeStep,
    model::{
        F16,
        dynamicmodel::fixedwing6dof::{FixedWing6DoF, FixedWing6DoFStates},
    },
    plots::state_variable_of_state_variable_plot,
    sim::simulator::SimulatorBuilder,
    trim::TrimProblemBuilder,
};

/// Example of a coordinated turn simulation for trimmed F-16 aircraft 6DoF-model.
pub fn f16_6dof_coordinated_turn_example() -> Result<(), Box<dyn std::error::Error>> {
    let setpoints = dvector![
        152.1, // vt [m/s]
        0.0,   // altitude [m]
        0.0,   // gamma [deg]
        0.0,   // roll rate [rad/s]
        0.0,   // pitch rate [rad/s]
        0.3,   // turn rate [rad/s]
        0.0,   // phi [rad] — calculated by solver for coordinated turn
        1.0,   // setpoint for coordinated turn
    ];
    let init_params = dvector![
        0.85, // throttle
        -6.0, // elevator
        0.24, // alpha
        0.1,  // aileron
        -0.4, // rudder
        0.0,  // beta
    ];

    let system = F16::new();
    let model = FixedWing6DoF;

    let problem = TrimProblemBuilder::new()
        .for_system(system)
        .with_model(model)
        .with_setpoints(setpoints)
        .with_initial_params(init_params)
        .build();
    let (x, u, _cost) = problem.trim()?;

    let system = F16::new();
    let model = FixedWing6DoF;

    let mut simulator = SimulatorBuilder::new()
        .for_system(system)
        .with_model(model)
        .with_state(x)
        .build();

    simulator.run(u, None, 180.0, TimeStep::new(0.001));

    state_variable_of_state_variable_plot(
        simulator.output,
        FixedWing6DoFStates::PosE as usize,
        FixedWing6DoFStates::PosN as usize,
        "Position E".into(),
        "Position N".into(),
        "F-16 Coordinated Turn".into(),
    )
}
