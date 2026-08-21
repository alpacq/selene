use crate::examples::transport_3dof_linearization_example::transport_3dof_linearization_example;

pub mod error;
pub mod examples;
pub mod linearize;
pub mod math;
pub mod model;
pub mod plots;
pub mod sim;
pub mod trim;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    transport_3dof_linearization_example()
}
