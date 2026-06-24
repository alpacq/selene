use crate::math::{IntegrableState, SizedVector};

pub mod fixedwing3dof;
pub mod fixedwing6dof;
pub mod state2;

/// Defines how a dynamic system is modeled with state equations
///
/// The state equations are used to compute the state of the system after one step.
pub trait DynamicModel<System> {
    type State: IntegrableState;
    type Input: SizedVector;
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

    /// Returns the rank of the system.
    ///
    /// The rank is the size of the state vector.
    fn system_rank(&self) -> usize;
}
