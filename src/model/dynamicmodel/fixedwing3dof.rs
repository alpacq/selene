use crate::{
    math::{IntegrableState, SizedVector},
    model::{
        DynamicModel, GRAVITY, RAD_TO_DEG,
        aerodynamics::Aerodynamics,
        aircraft::Aircraft,
        atmosphere::{dynamic_pressure, mach},
        engine::Engine,
    },
};
use nalgebra::{DVector, dvector};

/// A struct representing the state of a fixed-wing aircraft in 3D space
/// for simple longitudinal dynamics
pub struct FixedWing3DoFState {
    state_vector: DVector<f64>,
}

impl FixedWing3DoFState {
    /// Creates a new FixedWing3DoFState with the given state vector
    pub fn new(state_vector: DVector<f64>) -> Self {
        Self { state_vector }
    }

    /// TAS [m/s]
    pub fn vt(&self) -> f64 {
        self.state_vector[0]
    }

    /// angle of attack [rad]
    pub fn alpha(&self) -> f64 {
        self.state_vector[1]
    }

    /// pitch angle [deg]
    pub fn theta(&self) -> f64 {
        self.state_vector[2]
    }

    /// pitch rate [deg/s]
    pub fn q(&self) -> f64 {
        self.state_vector[3]
    }

    /// altitude [m]
    pub fn altitude(&self) -> f64 {
        self.state_vector[4]
    }
}

impl SizedVector for FixedWing3DoFState {
    /// Returns the size of the state vector
    ///
    /// # Returns
    ///
    /// The size of the state vector.
    fn size(&self) -> usize {
        self.state_vector.len()
    }

    /// Returns the state vector.
    ///
    /// # Returns
    ///
    /// The state vector.
    fn vector(&self) -> &DVector<f64> {
        &self.state_vector
    }
}

impl IntegrableState for FixedWing3DoFState {
    /// Creates a new FixedWing3DoFState from the given vector
    fn from_vector(vector: DVector<f64>) -> Self {
        FixedWing3DoFState::new(vector)
    }
}

/// A struct representing the input to a fixed-wing aircraft in 3D space
/// for simple longitudinal dynamics
pub struct FixedWing3DoFInput {
    input_vector: DVector<f64>,
}

impl FixedWing3DoFInput {
    /// throttle position [0..1]
    pub fn throttle(&self) -> f64 {
        self.input_vector[0].clamp(0.0, 1.0)
    }

    /// elevator position [-1..1]
    pub fn elevator(&self) -> f64 {
        self.input_vector[1].clamp(-1.0, 1.0)
    }

    /// x-axis position of the center of gravity [m]
    pub fn x_cg(&self) -> f64 {
        self.input_vector[2]
    }

    /// landing gear position [0..1]
    pub fn landing_gear_position(&self) -> f64 {
        self.input_vector[3]
    }
}

impl SizedVector for FixedWing3DoFInput {
    /// Returns the size of the input vector
    ///
    /// # Returns
    ///
    /// The size of the input vector.
    fn size(&self) -> usize {
        self.input_vector.len()
    }

    /// Returns the input vector.
    ///
    /// # Returns
    ///
    /// The input vector.
    fn vector(&self) -> &DVector<f64> {
        &self.input_vector
    }
}

/// Simple longitudinal 3 degrees of freedom model
/// of a fixed-wing aircraft
/// translation and pitching motion in the vertical plane
pub struct FixedWing3DOF;

impl<A: Aerodynamics, E: Engine> DynamicModel<Aircraft<A, E>> for FixedWing3DOF {
    type State = FixedWing3DoFState;
    type Input = FixedWing3DoFInput;

    fn state_equations(
        &self,
        system: &Aircraft<A, E>,
        x: &Self::State,
        u: &Self::Input,
    ) -> Self::State {
        let alpha_deg = RAD_TO_DEG * x.alpha();
        let pressure = dynamic_pressure(x.vt(), x.altitude());
        let mach = mach(x.vt(), x.altitude());

        let qs = pressure * system.airframe.s; // static pressure times wing area, reference force for scaling aerodynamic forces

        let gamma = x.theta() - x.alpha(); // flight path angle = pitch - angle of attack
        let sin_gamma = gamma.sin();
        let cos_gamma = gamma.cos();

        let cl0 = if u.landing_gear_position() == 0.0 {
            0.2
        } else {
            1.0
        };
        let cd0 = if u.landing_gear_position() == 0.0 {
            0.016
        } else {
            0.08
        };
        let cm0 = if u.landing_gear_position() == 0.0 {
            0.05
        } else {
            -0.2
        };
        let dcdg = if u.landing_gear_position() == 0.0 {
            0.0
        } else {
            0.02
        };
        let dcmg = if u.landing_gear_position() == 0.0 {
            0.0
        } else {
            -0.05
        };

        let thrust = system.engine.thrust(u.throttle(), x.altitude(), mach);
        let cl = cl0 + system.airframe.cla * alpha_deg;
        let cm = dcmg
            + cm0
            + system.airframe.cma * alpha_deg
            + system.airframe.cmde * u.elevator()
            + cl * (u.x_cg() - 0.25);
        let cd = dcdg + cd0 + system.airframe.cdcls * cl * cl;

        let x0_derivative =
            (thrust * x.alpha().cos() - qs * cd) / system.airframe.mass - GRAVITY * sin_gamma; // vt'
        let x1_derivative = (-thrust * x.alpha().sin() - qs * cl
            + system.airframe.mass * (x.vt() * x.q() + GRAVITY * cos_gamma))
            / (system.airframe.mass * x.vt() + qs * system.airframe.cladot); // alpha'
        let x2_derivative = x.q(); // theta'

        let pitch_damping = 0.5
            * system.airframe.cbar
            * (system.aerodynamics.cmq(x.alpha()) * x.q() + system.airframe.cmadot * x1_derivative)
            / x.vt(); // pitch damping

        let x3_derivative = (qs * system.airframe.cbar * (cm + pitch_damping)
            + thrust * system.airframe.ze)
            / system.airframe.iyy; // q'
        let x4_derivative = x.vt() * sin_gamma; // h'
        let x5_derivative = x.vt() * cos_gamma; // x'

        FixedWing3DoFState::new(dvector![
            x0_derivative,
            x1_derivative,
            x2_derivative,
            x3_derivative,
            x4_derivative,
            x5_derivative,
        ])
    }

    fn system_rank(&self) -> usize {
        5
    }
}
