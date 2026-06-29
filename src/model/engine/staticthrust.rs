use crate::model::{
    atmosphere::tas_from_mach,
    engine::{Engine, engineparams::EngineParams},
};

/// A struct representing a static thrust engine
/// Simple, constant thrust turboprop engine model
pub struct StaticThrust {
    pub params: EngineParams,
}

impl StaticThrust {
    /// Creates a new static thrust engine with default parameters
    pub fn new() -> Self {
        Self {
            params: EngineParams {
                hx: 0.0, //unused,
                tstat: 266893.0,
                dtdv: -554.6,
            },
        }
    }
}

impl Engine for StaticThrust {
    /// Converts throttle to power
    fn throttle_to_power(&self, throttle: f64) -> f64 {
        throttle
    }

    /// Power dynamics model (constant power output)
    fn power_dynamics(&self, _power: f64, set_power: f64) -> f64 {
        set_power
    }

    /// Thrust model (constant thrust)
    fn thrust(&self, power: f64, altitude: f64, mach: f64) -> f64 {
        (self.params.tstat + tas_from_mach(mach, altitude) * self.params.dtdv) * power
    }

    /// Inverse of tau model (constant time constant)
    fn tau_inverse(&self, _delta_power: f64) -> f64 {
        1.0
    }

    /// Engine angular momentum
    fn hx(&self) -> f64 {
        0.0
    }
}
