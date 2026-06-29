pub mod engineparams;
pub mod f100pw220;
pub mod staticthrust;

pub trait Engine {
    /// Converts throttle input to set power output.
    ///
    /// # Arguments
    ///
    /// * `throttle` - The throttle input, between 0 and 1.
    ///
    /// # Returns
    ///
    /// The power setpoint corresponding to the throttle input.
    fn throttle_to_power(&self, throttle: f64) -> f64;

    /// Computes the power dynamics of the engine.
    ///
    /// # Arguments
    ///
    /// * `power` - The current power output of the engine.
    /// * `set_power` - The power setpoint to achieve.
    ///
    /// # Returns
    ///
    /// The power output after one step of the dynamics.
    fn power_dynamics(&self, power: f64, set_power: f64) -> f64;

    /// Computes the thrust produced by the engine.
    ///
    /// # Arguments
    ///
    /// * `power` - The current power output of the engine.
    /// * `altitude` - The current altitude.
    /// * `mach` - The current Mach number.
    ///
    /// # Returns
    ///
    /// The thrust produced by the engine.
    fn thrust(&self, power: f64, altitude: f64, mach: f64) -> f64;

    /// Reverse time constant of the engine.
    ///
    /// # Arguments
    ///
    /// * `delta_power` - difference between set and current power output of the engine.
    ///
    /// # Returns
    /// The reversed time constant of the engine.
    fn tau_inverse(&self, delta_power: f64) -> f64;

    /// Engine angular momentum.
    ///
    /// # Returns
    /// Engine angular momentum.
    fn hx(&self) -> f64;
}
