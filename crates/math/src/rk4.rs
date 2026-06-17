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
    let k2 = state_equations(&(x + &k1 * dt_half), u);
    let k3 = state_equations(&(x + &k2 * dt_half), u);
    let k4 = state_equations(&(x + &k3 * dt), u);

    x + &(k1 + &k2 * 2.0 + &k3 * 2.0 + &k4) * (dt / 6.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn it_works() {
        let state = dvector![0.0];
        let input = dvector![0.0];
        let state_equation = |x: &State, _u: &Input| dvector![x[0]];
        let result = rk4(state_equation, &state, &input, &TimeStep::new(0.001));
        assert_eq!(result, dvector![0.0]);
    }
}
