use crate::model::{
    aircraft::Aircraft, airframeparams::AirframeParams, engine::staticthrust::StaticThrust,
};

/// Simple model of a medium size transport fixed-wing aircraft
/// powered by two turboprop engines
pub type Transport = Aircraft<StaticThrust>;

impl Transport {
    /// Creates a new `Transport` model with default parameters
    pub fn new() -> Self {
        Self {
            airframe: AirframeParams {
                s: 201.6,
                b: 0.0, //unused
                cbar: 5.334,
                mass: 72970.0,
                ixx: 0.0, //unused
                iyy: 5.559e6,
                izz: 0.0, //unused
                ixz: 0.0, //unused
                ze: 0.6096,
                cdcls: 0.042,
                cla: 0.085,
                cma: -0.022,
                cmde: -0.016,
                cmq: -16.0,
                cmadot: -6.0,
                cladot: 0.0,
            },
            engine: StaticThrust::new(),
        }
    }
}
