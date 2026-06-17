//! Input struct for the simulation.

use nalgebra::DVector;

#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    pub input_vector: DVector<f64>,
}

impl Input {
    /// Creates a new `Input` with the given input vector.
    ///
    /// # Arguments
    ///
    /// * `input_vector` - The input vector.
    pub fn new(input_vector: DVector<f64>) -> Self {
        Self { input_vector }
    }

    /// Returns the size of the input vector
    ///
    /// # Returns
    ///
    /// The size of the input vector.
    pub fn size(&self) -> usize {
        self.input_vector.len()
    }
}
