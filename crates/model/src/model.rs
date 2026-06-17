use math::{input::Input, state::State};

pub trait Model {
    fn state_equations(&self, x: &State, u: &Input) -> State;
}
