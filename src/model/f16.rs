use crate::model::{
    aerodynamics::f16aero::F16Aero, aircraft::Aircraft, airframeparams::AirframeParams,
    engine::f100pw220::F100PW220,
};

/// A struct representing an F-16 aircraft
/// based on the Stevens & Lewis book's
/// F-16 model parameters
pub type F16 = Aircraft<F16Aero, F100PW220>;

impl F16 {
    /// Creates a new F-16 aircraft with default parameters
    pub fn new() -> Self {
        Self {
            airframe: AirframeParams {
                s: 27.87,
                b: 9.144,
                cbar: 3.45,
                mass: 11340.5,
                ixx: 12874.9,
                iyy: 75674.3,
                izz: 85552.4,
                ixz: 1331.4,
                ze: 0.6096,
                cdcls: 0.042,
                cla: 0.085,
                cma: -0.022,
                cmde: -0.016,
                cmadot: -6.0,
                cladot: 0.0,
            },
            aerodynamics: F16Aero {},
            engine: F100PW220::new(),
        }
    }
}
