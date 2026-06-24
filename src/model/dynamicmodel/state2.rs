use crate::{
    math::{IntegrableState, SizedVector, input::Input, state::State},
    model::{DynamicModel, VanDerPol},
};
use nalgebra::{DVector, dvector};

#[derive(Clone)]
pub struct State2State(State);

impl State2State {
    pub fn new(state_vector: DVector<f64>) -> Self {
        Self(State::new(state_vector))
    }

    pub fn x1(&self) -> f64 {
        self.0.state_vector[0]
    }

    pub fn x2(&self) -> f64 {
        self.0.state_vector[1]
    }
}

impl SizedVector for State2State {
    fn size(&self) -> usize {
        self.0.size()
    }

    fn vector(&self) -> DVector<f64> {
        self.0.vector()
    }
}

impl IntegrableState for State2State {
    fn from_vector(vector: DVector<f64>) -> Self {
        Self::new(vector)
    }
}

pub struct State2Input(Input);

impl State2Input {
    pub fn new(input_vector: DVector<f64>) -> Self {
        Self(Input::new(input_vector))
    }

    pub fn u(&self) -> f64 {
        self.0.input_vector[0]
    }
}

impl SizedVector for State2Input {
    fn size(&self) -> usize {
        self.0.size()
    }

    fn vector(&self) -> DVector<f64> {
        self.0.vector()
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
