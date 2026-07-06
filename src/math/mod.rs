//! Necessary low level math abstractions and traits

use nalgebra::DVector;

pub mod input_functions;
pub mod interpolation;
pub mod timestep;

/// A trait for sized vectors that are used for states and inputs
pub trait SizedVector {
    /// Returns the size of the vector
    fn size(&self) -> usize;

    /// Returns a reference to the underlying [`DVector`]
    fn vector(&self) -> &DVector<f64>;

    /// Creates a [`SizedVector`] from a [`DVector`]
    fn from_vector(vector: DVector<f64>) -> Self;
}

#[cfg(test)]
pub mod test_utils;
