//! Input struct for the simulation.

use nalgebra::DVector;

#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    pub size: usize,
    pub input_vector: DVector<f64>,
}

impl Input {
    /// Creates a new `Input` with the given size and input vector.
    ///
    /// # Arguments
    ///
    /// * `size` - The size of the input vector.
    /// * `input_vector` - The input vector.
    ///
    /// # Panics
    ///
    /// Panics if `size` does not match `input_vector` length.
    pub fn new(size: usize, input_vector: DVector<f64>) -> Self {
        if size != input_vector.len() {
            panic!("size must match input_vector length");
        }
        Self { size, input_vector }
    }
}
