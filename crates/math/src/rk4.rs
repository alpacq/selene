//! Runge-Kutta 4th order integration algorithm

use crate::{input::Input, state::State, timestep::TimeStep};

/// Integrates the state equations using the Runge-Kutta 4th order method.
///
/// # Arguments
///
/// * `state_equations` - A function that computes the state equations.
///
/// * `x` - The  state vector.
/// * `u` - The input vector.
/// * `dt` - The time step.
///
/// # Returns
///
/// Returns the state vector in the next time step.
pub fn rk4<F>(state_equations: F, x: &State, u: &Input, dt: &TimeStep) -> State
where
    F: Fn(&State, &Input) -> State,
{
    let dt = dt.seconds();
    let dt_half = dt / 2.0;

    let k1 = state_equations(x, u);

    let x2 = &x.state_vector + &k1.state_vector * dt_half;
    let k2 = state_equations(&State::new(x.size, x2), u);

    let x3 = &x.state_vector + &k2.state_vector * dt_half;
    let k3 = state_equations(&State::new(x.size, x3), u);

    let x4 = &x.state_vector + &k3.state_vector * dt;
    let k4 = state_equations(&State::new(x.size, x4), u);

    let result_vector = &x.state_vector
        + &(k1.state_vector + &k2.state_vector * 2.0 + &k3.state_vector * 2.0 + &k4.state_vector)
            * (dt / 6.0);

    State::new(x.size, result_vector)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn it_works() {
        let state = State::new(1, dvector![0.0]);
        let input = Input::new(1, dvector![0.0]);
        let state_equation = |x: &State, _u: &Input| State::new(1, dvector![x.state_vector[0]]);
        let result = rk4(state_equation, &state, &input, &TimeStep::new(0.001));
        assert_eq!(result, State::new(1, dvector![0.0]));
    }
}
