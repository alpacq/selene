use crate::model::{
    aircraft::Aircraft, airframeparams::AirframeParams, engine::f100pw220::F100PW220,
};

pub type F16 = Aircraft<F100PW220>;

impl F16 {
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
                cmq: -16.0,
                cmadot: -6.0,
                cladot: 0.0,
            },
            engine: F100PW220::new(),
        }
    }
}
