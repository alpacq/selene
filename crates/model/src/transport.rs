use crate::{
    GD, RTOD,
    aircraft::Aircraft,
    airframeparams::AirframeParams,
    atmosphere::{dynamic_pressure, mach},
    dynamicmodel::DynamicModel,
    engine::{Engine, staticthrust::StaticThrust},
};
use math::{input::Input, state::State};
use nalgebra::dvector;

/// Simple longitudinal model of a medium size transport fixed-wing aircraft
/// powered by two turboprop engines
/// only 3-DoF !!! - translation and pitching motion in the vertical plane
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

impl DynamicModel for Transport {
    fn state_equations(&self, x: &State, u: &Input) -> State {
        // input vector: throttle, elevator, x_cg, landing_gear
        let throttle = if u.input_vector[0] >= 0.0 {
            u.input_vector[0]
        } else {
            0.0
        }; // throttle position [0..1]
        let elevator = u.input_vector[1]; // elevator position [-1..1]
        let x_cg = u.input_vector[2]; // x-axis position of the center of gravity [m]
        let landing_gear = u.input_vector[3]; // landing gear position [0, 1]

        // state vector: vt, alpha, theta, q, h
        let vt = x.state_vector[0]; // TAS [m/s]
        let alpha = x.state_vector[1]; // angle of attack [rad]
        let theta = x.state_vector[2]; // pitch angle [deg]
        let q = x.state_vector[3]; // pitch rate [deg/s]
        let altitude = x.state_vector[4]; // altitude [m]

        let alpha_deg = RTOD * alpha;
        let pressure = dynamic_pressure(vt, altitude);
        let mach = mach(vt, altitude);

        let qs = pressure * self.airframe.s; // static pressure times wing area, reference force for scaling aerodynamic forces

        let gamma = theta - alpha; // flight path angle = pitch - angle of attack
        let sin_gamma = gamma.sin();
        let cos_gamma = gamma.cos();

        let cl0 = if landing_gear == 0.0 { 0.2 } else { 1.0 };
        let cd0 = if landing_gear == 0.0 { 0.016 } else { 0.08 };
        let cm0 = if landing_gear == 0.0 { 0.05 } else { -0.2 };
        let dcdg = if landing_gear == 0.0 { 0.0 } else { 0.02 };
        let dcmg = if landing_gear == 0.0 { 0.0 } else { -0.05 };

        let thrust = self.engine.thrust(throttle, altitude, mach);
        let cl = cl0 + self.airframe.cla * alpha_deg;
        let cm = dcmg
            + cm0
            + self.airframe.cma * alpha_deg
            + self.airframe.cmde * elevator
            + cl * (x_cg - 0.25);
        let cd = dcdg + cd0 + self.airframe.cdcls * cl * cl;

        let x0_derivative = (thrust * alpha.cos() - qs * cd) / self.airframe.mass - GD * sin_gamma; // vt'
        let x1_derivative = (-thrust * alpha.sin() - qs * cl
            + self.airframe.mass * (vt * q + GD * cos_gamma))
            / (self.airframe.mass * vt + qs * self.airframe.cladot); // alpha'
        let x2_derivative = q; // theta'

        let pitch_damping = 0.5
            * self.airframe.cbar
            * (self.airframe.cmq * q + self.airframe.cmadot * x1_derivative)
            / vt; // pitch damping

        let x3_derivative = (qs * self.airframe.cbar * (cm + pitch_damping)
            + thrust * self.airframe.ze)
            / self.airframe.iyy; // q'
        let x4_derivative = vt * sin_gamma; // h'
        let x5_derivative = vt * cos_gamma; // x'

        State::new(dvector![
            x0_derivative,
            x1_derivative,
            x2_derivative,
            x3_derivative,
            x4_derivative,
            x5_derivative,
        ])
    }
}
