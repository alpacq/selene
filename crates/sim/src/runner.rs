use crate::output::SimOutput;
use math::{input::Input, rk4::rk4, state::State, timestep::TimeStep};
use model::Model;

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
pub fn run(
    model: &dyn Model,
    initial_state: State,
    initial_input: Input,
    duration: f64,
    dt: &TimeStep,
) -> SimOutput {
    let system_rank = initial_state.size();
    let number_of_steps = (duration / dt.seconds()) as usize;
    let mut time = 0.0_f64;
    let mut state = initial_state;
    let input = initial_input;

    let mut output = SimOutput::with_capacity(number_of_steps);

    for _ in 0..number_of_steps {
        state = rk4(
            |state, input| model.state_equations(state, input),
            &state,
            &input,
            dt,
        );
        let mut current_output = Vec::<f64>::with_capacity(system_rank);
        output.time.push(time);
        for state_variable in state.state_vector.iter() {
            current_output.push(*state_variable);
        }
        output.output_vector.push(current_output);
        time += dt.seconds();
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::VanDerPol;
    use nalgebra::dvector;

    #[test]
    fn it_works() {
        let number_of_steps = (1.0 / 0.001) as usize;
        let result = run(
            &VanDerPol {},
            State::new(dvector![0.1, 0.1]),
            Input::new(dvector![0.8]),
            1.0,
            &TimeStep::new(0.001),
        );
        assert_eq!(result.len(), number_of_steps)
    }
}
