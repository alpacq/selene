use nalgebra::DVector;

#[derive(Debug, Clone, PartialEq)]
pub struct State {
    pub size: usize,
    pub state_vector: DVector<f64>,
}

impl State {
    pub fn new(size: usize, state_vector: DVector<f64>) -> Self {
        if size != state_vector.len() {
            panic!("size must match state_vector length");
        }
        Self { size, state_vector }
    }
}
