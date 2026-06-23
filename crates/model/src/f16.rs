use crate::{
    DynamicModel,
    aircraft::Aircraft,
    airframeparams::AirframeParams,
    atmosphere::{air_pressure, mach},
    engine::f100pw220::F100PW220,
};
use math::{input::Input, state::State};
use nalgebra::dvector;

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

impl DynamicModel for F16 {
    fn state_equations(&self, x: &State, u: &Input) -> State {
        // input vector: throttle, elevator, aileron, rudder
        let throttle = u.input_vector[0]; // throttle position [0..1]
        let elevator = u.input_vector[1]; // elevator position [-1..1]
        let aileron = u.input_vector[2]; // aileron position [-1..1]
        let rudder = u.input_vector[3]; // rudder position [-1..1]

        // state vector: vt, alpha, beta, phi, theta, psi, p, q, r, position N (unused), position E (unused), altitude, power
        let vt = x.state_vector[0]; // TAS [m/s]
        let alpha = x.state_vector[1]; // angle of attack [rad]
        let beta = x.state_vector[2]; // sideslip angle [rad]
        let phi = x.state_vector[3]; // roll angle [rad]
        let theta = x.state_vector[4]; // pitch angle [rad]
        let psi = x.state_vector[5]; // yaw angle [rad]
        let p = x.state_vector[6]; // roll rate [rad/s]
        let q = x.state_vector[7]; // pitch rate [rad/s]
        let r = x.state_vector[8]; // yaw rate [rad/s]
        let altitude = x.state_vector[9]; // altitude [m]
        let power = x.state_vector[10]; // power [W]

        // additional temporary variables
        let xpq = self.airframe.ixz * (self.airframe.ixx - self.airframe.iyy + self.airframe.izz);
        let det = self.airframe.ixx * self.airframe.izz - (self.airframe.ixz * self.airframe.ixz);
        let xqr = self.airframe.izz * (self.airframe.izz - self.airframe.iyy)
            + (self.airframe.ixz * self.airframe.ixz);
        let zpq = (self.airframe.ixx - self.airframe.iyy) * self.airframe.ixx
            + (self.airframe.ixz * self.airframe.ixz);
        let ypr = self.airframe.izz - self.airframe.ixx;

        let pressure = air_pressure(vt, altitude);
        let mach = mach(vt, altitude);

        State::new(dvector![0.0])
    }
}
