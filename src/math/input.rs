//! Input struct for the simulation.

use nalgebra::DVector;

use crate::math::SizedVector;

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
}

impl SizedVector for Input {
    /// Returns the size of the input vector
    ///
    /// # Returns
    ///
    /// The size of the input vector.
    fn size(&self) -> usize {
        self.input_vector.len()
    }

    /// Returns the input vector.
    ///
    /// # Returns
    ///
    /// The input vector.
    fn vector(&self) -> DVector<f64> {
        self.input_vector.clone()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn four_element_vector_returns_size_four() {
        let input = Input::new(dvector![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(input.size(), 4);
    }

    #[test]
    fn can_access_input_vector_elements() {
        let input = Input::new(dvector![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(input.input_vector[0], 1.0);
        assert_eq!(input.input_vector[1], 2.0);
        assert_eq!(input.input_vector[2], 3.0);
        assert_eq!(input.input_vector[3], 4.0);
    }
}
