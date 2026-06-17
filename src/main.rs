use math::{input::Input, state::State, timestep::TimeStep};
use nalgebra::dvector;
use plots::phase_portait;
use sim::run;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output = run(
        &van_der_pol,
        State::new(2, dvector![0.1, 0.1]),
        Input::new(1, dvector![0.8]),
        60.0,
        &TimeStep::new(0.001),
    )?;

    phase_portait(output, "Van der Pol Oscillator".into())
}

pub fn van_der_pol(x: &State, u: &Input) -> State {
    let u = u.input_vector[0];
    let x1 = x.state_vector[0];
    let x2 = x.state_vector[1];

    State::new(2, dvector![x2, u * (1.0 - x1 * x1) * x2 - x1])
}
