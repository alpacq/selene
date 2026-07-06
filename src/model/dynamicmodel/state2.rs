use crate::{
    math::SizedVector,
    model::{DynamicModel, VanDerPol},
};
use nalgebra::{DVector, dvector};

/// A struct representing the state of a 2-dimensional system
#[derive(Debug, Clone)]
pub struct State2State {
    state_vector: DVector<f64>,
}

impl State2State {
    /// Creates a new `State2State` with the given state vector
    pub fn new(state_vector: DVector<f64>) -> Self {
        Self { state_vector }
    }

    /// Returns the first state variable
    pub fn x1(&self) -> f64 {
        self.state_vector[0]
    }

    /// Returns the second state variable
    pub fn x2(&self) -> f64 {
        self.state_vector[1]
    }
}

impl SizedVector for State2State {
    /// Returns the size of the state vector
    fn size(&self) -> usize {
        self.state_vector.len()
    }

    /// Returns a reference to the underlying [`DVector`]
    fn vector(&self) -> &DVector<f64> {
        &self.state_vector
    }

    // Creates a new `State2State` from the given vector
    fn from_vector(vector: DVector<f64>) -> Self {
        Self::new(vector)
    }
}
/// Input for the `State2` model
pub struct State2Input {
    input_vector: DVector<f64>,
}

impl State2Input {
    /// Creates a new `State2Input` with the given input vector
    pub fn new(input_vector: DVector<f64>) -> Self {
        Self { input_vector }
    }

    /// Returns the input value
    pub fn u(&self) -> f64 {
        self.input_vector[0]
    }
}

impl SizedVector for State2Input {
    /// Returns the size of the input vector
    fn size(&self) -> usize {
        self.input_vector.len()
    }

    /// Returns a reference to the underlying [`DVector`]
    fn vector(&self) -> &DVector<f64> {
        &self.input_vector
    }

    // Creates a new `State2Input` from the given vector
    fn from_vector(vector: DVector<f64>) -> Self {
        Self::new(vector)
    }
}

/// The 2-dimensional dynamic system model
/// with 1-dimensional input
pub struct State2;

impl DynamicModel<VanDerPol> for State2 {
    type State = State2State;
    type Input = State2Input;

    fn state_equations(
        &self,
        _system: &VanDerPol,
        x: &Self::State,
        u: &Self::Input,
    ) -> Self::State {
        State2State::new(dvector![
            x.x2(),
            u.u() * (1.0 - x.x1() * x.x1()) * x.x2() - x.x1()
        ])
    }

    fn system_rank(&self) -> usize {
        2
    }
}
