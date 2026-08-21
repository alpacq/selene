use crate::examples::f16_6dof_linearization_example::{
    f16_6dof_linearization_example_pullup, f16_6dof_linearization_example_turn,
};

pub mod error;
pub mod examples;
pub mod linearize;
pub mod math;
pub mod model;
pub mod plots;
pub mod sim;
pub mod trim;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    f16_6dof_linearization_example_pullup()?;
    f16_6dof_linearization_example_turn()?;
    Ok(())
}
