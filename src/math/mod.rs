use nalgebra::DVector;

pub mod timestep;

pub trait SizedVector {
    fn size(&self) -> usize;

    fn vector(&self) -> &DVector<f64>;
}

pub trait IntegrableState: SizedVector {
    fn from_vector(vector: DVector<f64>) -> Self;
}
