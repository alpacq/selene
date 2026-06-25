use crate::{
    math::{IntegrableState, SizedVector},
    model::{DynamicModel, VanDerPol},
};
use nalgebra::{DVector, dvector};

#[derive(Debug, Clone)]
pub struct State2State {
    state_vector: DVector<f64>,
}

impl State2State {
    pub fn new(state_vector: DVector<f64>) -> Self {
        Self { state_vector }
    }

    pub fn x1(&self) -> f64 {
        self.state_vector[0]
    }

    pub fn x2(&self) -> f64 {
        self.state_vector[1]
    }
}

impl SizedVector for State2State {
    fn size(&self) -> usize {
        self.state_vector.len()
    }

    fn vector(&self) -> &DVector<f64> {
        &self.state_vector
    }
}

impl IntegrableState for State2State {
    fn from_vector(vector: DVector<f64>) -> Self {
        Self::new(vector)
    }
}

pub struct State2Input {
    input_vector: DVector<f64>,
}

impl State2Input {
    pub fn new(input_vector: DVector<f64>) -> Self {
        Self { input_vector }
    }

    pub fn u(&self) -> f64 {
        self.input_vector[0]
    }
}

impl SizedVector for State2Input {
    fn size(&self) -> usize {
        self.input_vector.len()
    }

    fn vector(&self) -> &DVector<f64> {
        &self.input_vector
    }
}

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
