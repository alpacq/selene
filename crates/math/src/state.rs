//! Represents the state of a system.

use nalgebra::DVector;

#[derive(Debug, Clone, PartialEq)]
pub struct State {
    pub state_vector: DVector<f64>,
}

impl State {
    /// Creates a new `State` with the given state vector.
    ///
    /// # Arguments
    ///
    /// * `state_vector` - The state vector.
    pub fn new(state_vector: DVector<f64>) -> Self {
        Self { state_vector }
    }

    /// Returns the size of the state vector
    ///
    /// # Returns
    ///
    /// The size of the state vector.
    pub fn size(&self) -> usize {
        self.state_vector.len()
    }
}
