//! Runge-Kutta 4th order integration algorithm

use crate::math::{IntegrableState, timestep::TimeStep};

/// Integrates the state equations using the Runge-Kutta 4th order method.
///
/// # Arguments
///
/// * `state_equations` - A function that computes the state equations.
/// * `system` - The system to simulate.
/// * `x` - The  state vector.
/// * `u` - The input vector.
/// * `dt` - The time step.
///
/// # Returns
///
/// Returns the state vector in the next time step.
pub fn rk4<F, S, St, In>(state_equations: F, system: &S, x: &St, u: &In, dt: &TimeStep) -> St
where
    St: IntegrableState,
    F: Fn(&S, &St, &In) -> St,
{
    let dt = dt.seconds();
    let dt_half = dt / 2.0;

    let k1 = state_equations(system, x, u);

    let x2 = St::from_vector(x.vector() + k1.vector() * dt_half);
    let k2 = state_equations(system, &x2, u);

    let x3 = St::from_vector(x.vector() + k2.vector() * dt_half);
    let k3 = state_equations(system, &x3, u);

    let x4 = St::from_vector(x.vector() + k3.vector() * dt);
    let k4 = state_equations(system, &x4, u);

    let result_vector = x.vector()
        + (k1.vector() + k2.vector() * 2.0 + k3.vector() * 2.0 + k4.vector()) * (dt / 6.0);

    St::from_vector(result_vector)
}

#[cfg(test)]
mod tests {
    use crate::math::{input::Input, state::State, timestep::TimeStep};
    use crate::model::VanDerPol;

    use super::*;
    use nalgebra::dvector;

    #[test]
    fn it_works() {
        let state = State::new(dvector![0.0]);
        let input = Input::new(dvector![0.0]);
        let vanderpol = VanDerPol {};
        let state_equation =
            |_system: &VanDerPol, x: &State, _u: &Input| State::new(dvector![x.state_vector[0]]);
        let result = rk4(
            state_equation,
            &vanderpol,
            &state,
            &input,
            &TimeStep::new(0.001),
        );
        assert_eq!(result, State::new(dvector![0.0]));
    }
}
