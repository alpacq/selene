use math::{input::Input, state::State};

pub trait DynamicModel {
    /// Computes the step of model's state equations given the current state and input.
    ///
    /// # Arguments
    ///
    /// * `x` - The current state of the model.
    /// * `u` - The input to the model.
    ///
    /// # Returns
    ///
    /// The state of the model after one step.
    fn state_equations(&self, x: &State, u: &Input) -> State;
}
