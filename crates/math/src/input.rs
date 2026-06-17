use nalgebra::DVector;

#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    pub size: usize,
    pub input_vector: DVector<f64>,
}

impl Input {
    pub fn new(size: usize, input_vector: DVector<f64>) -> Self {
        if size != input_vector.len() {
            panic!("size must match input_vector length");
        }
        Self { size, input_vector }
    }
}
