use nalgebra::dvector;

use crate::{
    math::timestep::TimeStep,
    model::{
        VanDerPol,
        dynamicmodel::state2::{State2, State2Input, State2State},
    },
    plots::phase_portrait,
    sim::simulator::create_sim_and_run,
};

/// Example of simulating the dynamic system - Van der Pol Oscillator
/// Example 3.4-1 from "Aircraft Control and Simulation - 3rd edition" by Brian L. Stevens, Frank L. Lewis and Eric N. Johnson.
pub fn vanderpol_example() -> Result<(), Box<dyn std::error::Error>> {
    let simulator = create_sim_and_run(
        VanDerPol {},
        State2 {},
        State2State::new(dvector![0.1, 0.1]),
        State2Input::new(dvector![0.8]),
        None,
        60.0,
        TimeStep::new(0.001),
    );

    phase_portrait(simulator.output, "Van der Pol Oscillator".into())
}
