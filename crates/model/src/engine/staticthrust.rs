use crate::{
    atmosphere::tas_from_mach,
    engine::{Engine, engineparams::EngineParams},
};

pub struct StaticThrust {
    pub params: EngineParams,
}

impl StaticThrust {
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
    fn throttle_to_power(&self, throttle: f64) -> f64 {
        throttle
    }

    fn power_dynamics(&self, _power: f64, set_power: f64) -> f64 {
        set_power
    }

    fn thrust(&self, power: f64, altitude: f64, mach: f64) -> f64 {
        (self.params.tstat + tas_from_mach(mach, altitude) * self.params.dtdv) * power
    }

    fn tau_inverse(&self, _delta_power: f64) -> f64 {
        1.0
    }
}
