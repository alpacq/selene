//! Represents the state of a system.

use nalgebra::DVector;

use crate::math::{IntegrableState, SizedVector};

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
}

impl SizedVector for State {
    /// Returns the size of the state vector
    ///
    /// # Returns
    ///
    /// The size of the state vector.
    fn size(&self) -> usize {
        self.state_vector.len()
    }

    /// Returns the state vector.
    ///
    /// # Returns
    ///
    /// The state vector.
    fn vector(&self) -> DVector<f64> {
        self.state_vector.clone()
    }
}

impl IntegrableState for State {
    fn from_vector(vector: DVector<f64>) -> Self {
        Self::new(vector)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn four_element_vector_returns_size_four() {
        let state = State::new(dvector![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(state.size(), 4);
    }

    #[test]
    fn can_access_state_vector_elements() {
        let state = State::new(dvector![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(state.state_vector[0], 1.0);
        assert_eq!(state.state_vector[1], 2.0);
        assert_eq!(state.state_vector[2], 3.0);
        assert_eq!(state.state_vector[3], 4.0);
    }
}
