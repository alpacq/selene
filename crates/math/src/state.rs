//! Represents the state of a system.

use nalgebra::DVector;

#[derive(Debug, Clone, PartialEq)]
pub struct State {
    pub size: usize,
    pub state_vector: DVector<f64>,
}

impl State {
    /// Creates a new `State` with the given size and state vector.
    ///
    /// # Arguments
    ///
    /// * `size` - The size of the state vector.
    /// * `state_vector` - The state vector.
    ///
    /// # Panics
    ///
    /// Panics if `size` does not match `state_vector` length.
    pub fn new(size: usize, state_vector: DVector<f64>) -> Self {
        if size != state_vector.len() {
            panic!("size must match state_vector length");
        }
        Self { size, state_vector }
    }
}
