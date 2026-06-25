//! Necessary low level math abstractions and traits

use nalgebra::DVector;

pub mod timestep;

/// A trait for sized vectors that are used for states and inputs
pub trait SizedVector {
    /// Returns the size of the vector
    fn size(&self) -> usize;

    /// Returns a reference to the underlying [`DVector`]
    fn vector(&self) -> &DVector<f64>;
}

/// A trait for integrable states that can be created from [`DVector`]
pub trait IntegrableState: SizedVector {
    /// Creates an [`IntegrableState`] from a [`DVector`]
    fn from_vector(vector: DVector<f64>) -> Self;
}
