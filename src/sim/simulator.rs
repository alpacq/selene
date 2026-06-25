use crate::math::{SizedVector, timestep::TimeStep};
use crate::model::DynamicModel;

use crate::sim::output::SimOutput;
use crate::sim::rk4::rk4;

/// The simulator for a dynamic model.
///
/// # Example
///
/// ```
/// use selene::sim::Simulator;
/// use selene::model::dynamicmodel::State2;
/// use selene::model::VanDerPol;
///
/// let mut simulator = Simulator::new(VanDerPol::default(), State2);
/// ```
pub struct Simulator<S, DM>
where
    DM: DynamicModel<S>,
{
    /// The system being simulated.
    pub system: S,
    /// The model being simulated.
    pub model: DM,
    /// The current time of the simulation.
    pub time: f64,
    /// The current state of the simulation.
    pub state: DM::State,
    /// The output of the simulation.
    pub output: SimOutput,
}

impl<S, DM> Simulator<S, DM>
where
    DM: DynamicModel<S>,
{
    /// Steps the simulation forward by one time step using the RK4 method.
    ///
    /// # Arguments
    ///
    /// * `state` - The current state of the simulation.
    /// * `input` - The input to the simulation.
    /// * `dt` - The time step to use for the simulation.
    ///
    /// # Returns
    ///
    /// The state of the simulation after one time step.
    pub fn step(&mut self, state: DM::State, input: &DM::Input, dt: &TimeStep) -> DM::State {
        let state = rk4(
            |system: &S, state: &DM::State, input: &DM::Input| {
                self.model.state_equations(system, state, input)
            },
            &self.system,
            &state,
            input,
            dt,
        );
        state
    }

    /// Runs the simulation for the given duration and time step and saves the output.
    ///
    /// # Arguments
    ///
    /// * `initial_state` - The initial state of the simulation.
    /// * `initial_input` - The initial input to the simulation.
    /// * `time` - The total time to simulate.
    /// * `dt` - The time step to use for the simulation.
    pub fn run(
        &mut self,
        initial_state: DM::State,
        initial_input: DM::Input,
        duration: f64,
        dt: TimeStep,
    ) {
        let system_rank = self.state.size();
        let number_of_steps = (duration / dt.seconds()) as usize;
        let mut time = 0.0_f64;
        let mut state = initial_state;
        let input = initial_input;

        self.output = SimOutput::with_capacity(number_of_steps);

        for _ in 0..number_of_steps {
            state = self.step(state, &input, &dt);
            let mut current_output = Vec::<f64>::with_capacity(system_rank);
            self.output.time.push(time);
            for state_variable in state.vector().iter() {
                current_output.push(*state_variable);
            }
            self.output.output_vector.push(current_output);
            time += dt.seconds();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::timestep::TimeStep;
    use crate::model::VanDerPol;
    use crate::model::dynamicmodel::state2::{State2, State2Input, State2State};
    use nalgebra::dvector;

    #[test]
    fn it_works() {
        let mut simulator = Simulator::<VanDerPol, State2> {
            system: VanDerPol {},
            model: State2 {},
            time: 0.0,
            state: State2State::new(dvector![0.1, 0.1]),
            output: SimOutput::default(),
        };
        let number_of_steps = (1.0 / 0.001) as usize;
        let input = State2Input::new(dvector![0.8]);
        let initial_state = simulator.state.clone();
        simulator.run(initial_state, input, 1.0, TimeStep::new(0.001));
        assert_eq!(simulator.output.len(), number_of_steps)
    }
}
