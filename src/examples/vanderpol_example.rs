use nalgebra::dvector;

use crate::{
    math::timestep::TimeStep,
    model::{
        VanDerPol,
        dynamicmodel::state2::{State2, State2Input, State2State},
    },
    plots::phase_portrait,
    sim::simulator::SimulatorBuilder,
};

/// Example of simulating the dynamic system - Van der Pol Oscillator
pub fn vanderpol_example() -> Result<(), Box<dyn std::error::Error>> {
    let mut simulator = SimulatorBuilder::new()
        .for_system(VanDerPol {})
        .with_model(State2 {})
        .with_state(State2State::new(dvector![0.1, 0.1]))
        .build();

    simulator.run(
        State2Input::new(dvector![0.8]),
        None,
        60.0,
        TimeStep::new(0.001),
    );

    phase_portrait(simulator.output, "Van der Pol Oscillator".into())
}
