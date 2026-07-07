use nalgebra::DVector;

use crate::math::{SizedVector, timestep::TimeStep};
use crate::model::DynamicModel;

use crate::sim::output::SimOutput;
use crate::sim::rk4::rk4;

/// The simulator for a dynamic model.
///
/// Use SimulatorBuilder to create a simulator.
pub struct Simulator<Sys, M>
where
    M: DynamicModel<Sys>,
{
    /// The system being simulated.
    system: Sys,
    /// The model being simulated.
    model: M,
    /// The current time of the simulation.
    pub time: f64,
    /// The current state of the simulation.
    pub state: M::State,
    /// The output of the simulation.
    pub output: SimOutput,
}

impl<Sys, M> Simulator<Sys, M>
where
    M: DynamicModel<Sys>,
{
    /// Returns a builder for creating a simulator from specific
    /// system and dynamic model.
    ///
    /// # Example
    ///
    /// ```
    /// use selene::sim::Simulator;
    /// use selene::model::dynamicmodel::State2;
    /// use selene::model::VanDerPol;
    ///
    /// let mut simulator = Simulator::builder()
    ///     .for_system(VanDerPol {})
    ///     .with_model(State2 {})
    ///     .with_state(State2State::new(dvector![0.1, 0.1]))
    ///     .build();
    /// ```
    pub fn builder() -> SimulatorBuilder {
        SimulatorBuilder::new()
    }
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
    pub fn step(&mut self, state: M::State, input: &M::Input, dt: &TimeStep) -> M::State {
        let state = rk4(
            |system: &Sys, state: &M::State, input: &M::Input| {
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
    /// * `initial_input` - The initial input to the simulation.
    /// * `input_fn` - An optional function to compute the input at each time step.
    ///     If `None`, `initial_input` is used for all time steps.
    /// * `duration` - The total time to simulate.
    /// * `dt` - The time step to use for the simulation.
    pub fn run(
        &mut self,
        initial_input: M::Input,
        input_fn: Option<fn(&DVector<f64>, f64) -> M::Input>,
        duration: f64,
        dt: TimeStep,
    ) {
        let system_rank = self.state.size();
        let number_of_steps = (duration / dt.seconds()) as usize;
        let mut time = 0.0_f64;
        let mut state = M::State::from_vector(self.state.vector().clone());
        let input_params = initial_input.vector().clone();
        let mut input = initial_input;

        self.output = SimOutput::with_capacity(number_of_steps);

        for _ in 0..number_of_steps {
            if let Some(f) = input_fn {
                input = f(&input_params, time);
            }
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

/// First step in building the simulator process
/// Empty builder
pub struct SimulatorBuilder;

/// Second step in building the simulator process
/// Requires a system to be specified
pub struct SimulatorBuilderWithSystem<Sys> {
    system: Sys,
}

/// Third step in building the simulator process
/// Requires a dynamic model to be specified
pub struct SimulatorBuilderWithModel<Sys, M: DynamicModel<Sys>> {
    system: Sys,
    model: M,
}

/// Fourth step in building the simulator process
/// Requires an initial state to be specified
pub struct SimulatorBuilderWithInitialState<Sys, M: DynamicModel<Sys>> {
    system: Sys,
    model: M,
    time: f64,
    state: M::State,
    output: SimOutput,
}

impl SimulatorBuilder {
    /// Returns a new empty builder for creating a simulator.
    pub fn new() -> Self {
        SimulatorBuilder
    }

    /// Specifies the system to use for the simulator.
    pub fn for_system<Sys>(self, system: Sys) -> SimulatorBuilderWithSystem<Sys> {
        SimulatorBuilderWithSystem { system }
    }
}

impl<Sys> SimulatorBuilderWithSystem<Sys> {
    /// Specifies the model to use for the simulator.
    pub fn with_model<M: DynamicModel<Sys>>(self, model: M) -> SimulatorBuilderWithModel<Sys, M> {
        SimulatorBuilderWithModel {
            system: self.system,
            model,
        }
    }
}

impl<Sys, M: DynamicModel<Sys>> SimulatorBuilderWithModel<Sys, M> {
    /// Specifies the initial state to use for the simulator.
    pub fn with_state(self, state: M::State) -> SimulatorBuilderWithInitialState<Sys, M> {
        SimulatorBuilderWithInitialState {
            system: self.system,
            model: self.model,
            time: 0.0,
            state,
            output: SimOutput::default(),
        }
    }
}

impl<Sys, M: DynamicModel<Sys>> SimulatorBuilderWithInitialState<Sys, M> {
    /// Specifies the time to use for the simulator.
    pub fn time(mut self, time: f64) -> Self {
        self.time = time;
        self
    }

    /// Returns the simulator with given parameters
    pub fn build(self) -> Simulator<Sys, M> {
        Simulator {
            system: self.system,
            model: self.model,
            time: self.time,
            state: self.state,
            output: self.output,
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
        let mut simulator = SimulatorBuilder::new()
            .for_system(VanDerPol {})
            .with_model(State2 {})
            .with_state(State2State::new(dvector![0.1, 0.1]))
            .build();
        let number_of_steps = (1.0 / 0.001) as usize;
        let input = State2Input::new(dvector![0.8]);
        simulator.run(input, None, 1.0, TimeStep::new(0.001));
        assert_eq!(simulator.output.len(), number_of_steps)
    }
}
