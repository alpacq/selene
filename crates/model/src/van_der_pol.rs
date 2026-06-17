use crate::model::Model;
use math::{input::Input, state::State};
use nalgebra::dvector;

pub struct VanDerPol {}

impl Model for VanDerPol {
    fn state_equations(&self, x: &State, u: &Input) -> State {
        let u = u.input_vector[0];
        let x1 = x.state_vector[0];
        let x2 = x.state_vector[1];

        State::new(dvector![x2, u * (1.0 - x1 * x1) * x2 - x1])
    }
}
