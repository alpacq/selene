use crate::examples::f16_6dof_dynamic_behavior_example::f16_6dof_dynamic_behavior_example;

pub mod error;
pub mod examples;
pub mod linearize;
pub mod math;
pub mod model;
pub mod plots;
pub mod sim;
pub mod trim;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    f16_6dof_dynamic_behavior_example()?;
    Ok(())
}
