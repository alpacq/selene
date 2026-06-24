pub mod math;
pub mod model;
pub mod plots;
pub mod sim;

use nalgebra::dvector;

use crate::{
    math::timestep::TimeStep,
    model::{
        VanDerPol,
        dynamicmodel::state2::{State2, State2Input, State2State},
    },
    plots::phase_portrait,
    sim::{output::SimOutput, simulator::Simulator},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut simulator = Simulator::<VanDerPol, State2> {
        system: VanDerPol {},
        model: State2 {},
        time: 0.0,
        state: State2State::new(dvector![0.1, 0.1]),
        output: SimOutput::default(),
    };
    let initial_state = simulator.state.clone();
    simulator.run(
        initial_state,
        State2Input::new(dvector![0.8]),
        60.0,
        TimeStep::new(0.001),
    );

    phase_portrait(simulator.output, "Van der Pol Oscillator".into())
}
