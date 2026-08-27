use crate::math::SizedVector;

/// Defines how a dynamic system is modeled with state equations
///
/// The state equations are used to compute the state of the system after one step.
pub trait DynamicModel<Sys> {
    type State: SizedVector;
    type Input: SizedVector;
    type Output: SizedVector;
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
    fn state_equations(
        &self,
        system: &Sys,
        x: &Self::State,
        u: &Self::Input,
    ) -> (Self::State, Self::Output);

    /// Returns the rank of the system.
    ///
    /// The rank is the size of the state vector.
    fn system_rank(&self) -> usize;

    /// Returns the number of inputs to the system.
    fn num_inputs(&self) -> usize;

    /// Returns the number of outputs from the system.
    fn num_outputs(&self) -> usize;
}
