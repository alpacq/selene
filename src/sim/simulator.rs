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
    pub fn step(&self, state: M::State, input: &M::Input, dt: &TimeStep) -> M::State {
        rk4(
            |system: &Sys, state: &M::State, input: &M::Input| {
                self.model.state_equations(system, state, input).0
            },
            &self.system,
            &state,
            input,
            dt,
        )
    }

    /// Appends one sample of the time history to the simulation output.
    ///
    /// The model output `y = g(x, u)` is evaluated at the sampled state and
    /// input, so it is consistent with the recorded state.
    ///
    /// # Arguments
    ///
    /// * `time` - The time of the sample.
    /// * `state` - The state at that time.
    /// * `input` - The input applied at that time.
    fn record(&mut self, time: f64, state: &M::State, input: &M::Input) {
        let output = self.model.state_equations(&self.system, state, input).1;
        self.output.time.push(time);
        self.output
            .state_vector
            .push(state.vector().iter().copied().collect());
        self.output
            .output_vector
            .push(output.vector().iter().copied().collect());
    }

    /// Runs the simulation for the given duration and time step and saves the output.
    ///
    /// The run starts at the simulator's current time and state, and both are
    /// updated to the final values when it completes. The initial condition is
    /// recorded as the first sample, so `n` steps produce `n + 1` samples.
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
        let number_of_steps = (duration / dt.seconds()) as usize;
        let mut time = self.time;
        let mut state = M::State::from_vector(self.state.vector().clone());
        let input_params = initial_input.vector().clone();
        let mut input = initial_input;

        self.output = SimOutput::with_capacity(number_of_steps + 1);

        if let Some(f) = input_fn {
            input = f(&input_params, time);
        }
        self.record(time, &state, &input);

        for _ in 0..number_of_steps {
            // The input is held constant over the step, then re-evaluated at
            // the new time so that the recorded sample is consistent.
            state = self.step(state, &input, &dt);
            time += dt.seconds();
            if let Some(f) = input_fn {
                input = f(&input_params, time);
            }
            self.record(time, &state, &input);
        }

        self.state = state;
        self.time = time;
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

    use crate::math::test_utils::assert_approx;

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
        // The initial condition is recorded too, hence the extra sample.
        assert_eq!(simulator.output.len(), number_of_steps + 1)
    }

    /// Every recorded vector must describe the same instant, and the first
    /// sample must be the initial condition at the starting time.
    #[test]
    fn output_starts_at_the_initial_condition() {
        let mut simulator = SimulatorBuilder::new()
            .for_system(VanDerPol {})
            .with_model(State2 {})
            .with_state(State2State::new(dvector![0.1, 0.2]))
            .build();

        simulator.run(
            State2Input::new(dvector![0.8]),
            None,
            1.0,
            TimeStep::new(0.001),
        );

        let output = &simulator.output;
        assert_eq!(output.time.len(), output.state_vector.len());
        assert_eq!(output.time.len(), output.output_vector.len());
        assert_approx(output.time[0], 0.0, 1e-12, "t0");
        assert_approx(output.state_vector[0][0], 0.1, 1e-12, "x1(0)");
        assert_approx(output.state_vector[0][1], 0.2, 1e-12, "x2(0)");
        // State2 defines y = x1.
        assert_approx(output.output_vector[0][0], 0.1, 1e-12, "y(0)");
    }

    /// The last recorded sample must sit at the end of the requested interval
    /// and match the simulator's own final state.
    #[test]
    fn run_leaves_the_simulator_at_the_final_state() {
        let mut simulator = SimulatorBuilder::new()
            .for_system(VanDerPol {})
            .with_model(State2 {})
            .with_state(State2State::new(dvector![0.1, 0.1]))
            .build();

        simulator.run(
            State2Input::new(dvector![0.8]),
            None,
            1.0,
            TimeStep::new(0.001),
        );

        let last = simulator.output.state_vector.last().expect("no samples");
        assert_approx(simulator.time, 1.0, 1e-9, "final time");
        assert_approx(*simulator.output.time.last().unwrap(), 1.0, 1e-9, "t_end");
        assert_approx(simulator.state.x1(), last[0], 1e-12, "final x1");
        assert_approx(simulator.state.x2(), last[1], 1e-12, "final x2");
    }
}
