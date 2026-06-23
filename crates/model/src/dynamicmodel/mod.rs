pub mod fixedwing3dof;

/// Defines how a dynamic system is modeled with state equations
///
/// The state equations are used to compute the state of the system after one step.
pub trait DynamicModel<System> {
    type State;
    type Input;
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
    fn state_equations(&self, system: &System, x: &Self::State, u: &Self::Input) -> Self::State;
}
