use crate::examples::f16_6dof_coordinated_turn_example::f16_6dof_coordinated_turn_example;

pub mod error;
pub mod examples;
pub mod linearize;
pub mod math;
pub mod model;
pub mod plots;
pub mod sim;
pub mod trim;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    f16_6dof_coordinated_turn_example()
}
