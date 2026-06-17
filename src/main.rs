use math::{input::Input, state::State, timestep::TimeStep};
use model::VanDerPol;
use nalgebra::dvector;
use plots::phase_portrait;
use sim::run;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output = run(
        &VanDerPol {},
        State::new(dvector![0.1, 0.1]),
        Input::new(dvector![0.8]),
        60.0,
        &TimeStep::new(0.001),
    );

    phase_portrait(output, "Van der Pol Oscillator".into())
}
