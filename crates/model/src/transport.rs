use crate::{GD, RTOD, atmosphere::air_data, fixedwingparams::AircraftParams, model::Model};
use math::{input::Input, state::State};
use nalgebra::dvector;

/// Simple longitudinal model of a medium size transport fixed-wing aircraft
/// powered by two turboprop engines
/// only 3-DoF !!! - translation and pitching motion in the vertical plane
pub struct Transport {
    pub params: AircraftParams,
}

impl Transport {
    /// Creates a new `Transport` model with default parameters
    pub fn new() -> Self {
        Self {
            params: AircraftParams {
                s: 201.6,
                cbar: 5.334,
                mass: 72970.0,
                iyy: 5.559e6,
                tstat: 266893.0,
                dtdv: -554.6,
                ze: 0.6096,
                cdcls: 0.042,
                cla: 0.085,
                cma: -0.022,
                cmde: -0.016,
                cmq: -16.0,
                cmadot: -6.0,
                cladot: 0.0,
            },
        }
    }
}

impl Model for Transport {
    fn state_equations(&self, x: &State, u: &Input) -> State {
        let throttle = if u.input_vector[0] >= 0.0 {
            u.input_vector[0]
        } else {
            0.0
        }; // throttle position [0..1]
        let elevator = u.input_vector[1]; // elevator position [-1..1]
        let x_cg = u.input_vector[2]; // x-axis position of the center of gravity [m]
        let landing_gear = u.input_vector[3]; // landing gear position [0, 1]

        let vt = x.state_vector[0]; // TAS [m/s]
        let alpha = RTOD * x.state_vector[1]; // angle of attack [deg]
        let theta = x.state_vector[2]; // pitch angle [deg]
        let q = x.state_vector[3]; // pitch rate [deg/s]
        let h = x.state_vector[4]; // altitude [m]

        let (pressure, _) = air_data(vt, h);

        let qs = pressure * self.params.s; // static pressure times wing area, reference force for scaling aerodynamic forces

        let sin_alpha = alpha.sin();
        let cos_alpha = alpha.cos();

        let gamma = theta - x.state_vector[1]; // flight path angle = pitch - angle of attack
        let sin_gamma = gamma.sin();
        let cos_gamma = gamma.cos();

        let cl0 = if landing_gear == 0.0 { 0.2 } else { 1.0 };
        let cd0 = if landing_gear == 0.0 { 0.016 } else { 0.08 };
        let cm0 = if landing_gear == 0.0 { 0.05 } else { -0.2 };
        let dcdg = if landing_gear == 0.0 { 0.0 } else { 0.02 };
        let dcmg = if landing_gear == 0.0 { 0.0 } else { -0.05 };

        let thrust = (self.params.tstat + vt * self.params.dtdv) * throttle;
        let cl = cl0 + self.params.cla * alpha;
        let cm =
            dcmg + cm0 + self.params.cma * alpha + self.params.cmde * elevator + cl * (x_cg - 0.25);
        let cd = dcdg + cd0 + self.params.cdcls * cl * cl;

        let x0_derivative = (thrust * cos_alpha - qs * cd) / self.params.mass - GD * sin_gamma;
        let x1_derivative = (-thrust * sin_alpha - qs * cl
            + self.params.mass * (vt * q + GD * cos_gamma))
            / (self.params.mass * vt + qs * self.params.cladot);
        let x2_derivative = q;

        let d = 0.5 * self.params.cbar * (self.params.cmq * q + self.params.cmadot * x1_derivative)
            / vt; // pitch damping

        let x3_derivative =
            (qs * self.params.cbar * (cm + d) + thrust * self.params.ze) / self.params.iyy;
        let x4_derivative = vt * sin_gamma;
        let x5_derivative = vt * cos_gamma;

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
