//! Sim is a crate when flight simulation takes place.

use math::{input::Input, rk4::rk4, state::State, timestep::TimeStep};
use nalgebra::{DVector, dvector};
use std::io::Error;

/// Runs the simulation for the given time step and returns the final state.
///
/// # Arguments
///
/// * `state_equations` - The state equations to use for the simulation.
/// * `initial_state` - The initial state of the simulation.
/// * `initial_input` - The initial input to the simulation.
/// * `time` - The total time to simulate.
/// * `dt` - The time step to use for the simulation.
///
/// # Returns
///
/// Returns `Ok(Vec<DVector<f64>>)` if the simulation completes successfully, or `Err` if an error occurs.
/// Returned vector represents the simulation output, where first element of each DVector is the time,
/// and the remaining elements are the state variables.
pub fn run<F>(
    state_equations: &F,
    initial_state: State,
    initial_input: Input,
    duration: f64,
    dt: &TimeStep,
) -> Result<Vec<DVector<f64>>, Error>
where
    F: Fn(&State, &Input) -> State,
{
    let number_of_steps = (duration / dt.seconds()) as usize;
    let mut time = 0.0_f64;
    let mut state = initial_state;
    let mut input = initial_input;

    let mut output = Vec::<DVector<f64>>::with_capacity(number_of_steps);

    for _ in 0..number_of_steps {
        state = rk4(state_equations, &state, &input, dt);
        output.push(dvector![time, state[0], state[1]]);
        time += dt.seconds();
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    fn van_der_pol(x: &State, u: &Input) -> State {
        let u = u[0];
        let x1 = x[0];
        let x2 = x[1];

        dvector![x2, u * (1.0 - x1 * x1) * x2 - x1]
    }

    #[test]
    fn it_works() {
        let result = run(
            &van_der_pol,
            dvector![0.1, 0.1],
            dvector![0.8],
            1.0,
            &TimeStep::new(0.001),
        );
        assert!(result.is_ok());
    }
}
