use crate::{
    error::TrimError,
    math::SizedVector,
    model::{
        DynamicModel, GRAVITY, RAD_TO_DEG,
        aerodynamics::Aerodynamics,
        aircraft::Aircraft,
        atmosphere::{dynamic_pressure, mach},
        engine::Engine,
    },
    trim::TrimTarget,
};
use nalgebra::{DVector, dvector};
use num_traits::Pow;

pub enum FixedWing3DoFStates {
    Vt,
    Alpha,
    Theta,
    Q,
    Altitude,
}

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
        self.state_vector[FixedWing3DoFStates::Vt as usize]
    }

    /// angle of attack [rad]
    pub fn alpha(&self) -> f64 {
        self.state_vector[FixedWing3DoFStates::Alpha as usize]
    }

    /// pitch angle [deg]
    pub fn theta(&self) -> f64 {
        self.state_vector[FixedWing3DoFStates::Theta as usize]
    }

    /// pitch rate [deg/s]
    pub fn q(&self) -> f64 {
        self.state_vector[FixedWing3DoFStates::Q as usize]
    }

    /// altitude [m]
    pub fn altitude(&self) -> f64 {
        self.state_vector[FixedWing3DoFStates::Altitude as usize]
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

    /// elevator deflection [deg], limited to +-35 deg of travel
    pub fn elevator(&self) -> f64 {
        self.input_vector[1].clamp(-35.0, 35.0)
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

    /// Creates a new FixedWing3DoFInput from the given vector
    fn from_vector(vector: DVector<f64>) -> Self {
        FixedWing3DoFInput {
            input_vector: vector,
        }
    }
}

/// Simple longitudinal 3 degrees of freedom model
/// of a fixed-wing aircraft
/// translation and pitching motion in the vertical plane
pub struct FixedWing3DoF;

impl<A: Aerodynamics, E: Engine> DynamicModel<Aircraft<A, E>> for FixedWing3DoF {
    type State = FixedWing3DoFState;
    type Input = FixedWing3DoFInput;

    /// State equations for the 3-DoF simple longitudinal model
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

        let (cl0, cd0, cm0, dcdg, dcmg) = if u.landing_gear_position() == 0.0 {
            (0.2, 0.016, 0.05, 0.0, 0.0)
        } else {
            (1.0, 0.08, -0.2, 0.02, -0.05)
        };

        let thrust = system.engine.thrust(u.throttle(), x.altitude(), mach);
        let cl = cl0 + system.aerodynamics.cz(x.alpha(), 0.0, u.elevator()) * alpha_deg;
        let cm = dcmg
            + cm0
            + system.aerodynamics.cm(x.alpha(), 0.0) * alpha_deg
            + system.aerodynamics.cm(0.0, u.elevator()) * u.elevator()
            + cl * (u.x_cg() - 0.25);
        let cd = dcdg + cd0 + system.aerodynamics.cx(x.alpha(), u.elevator()) * cl * cl;

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

    /// Returns the rank of the system (number of state variables)
    /// For this 3-DoF simple longitudinal model, the rank is 5.
    fn system_rank(&self) -> usize {
        5
    }
}

impl<A: Aerodynamics, E: Engine> TrimTarget<Aircraft<A, E>> for FixedWing3DoF {
    /// Sets up the state and input for the given setpoints and parameters.
    ///
    /// # Arguments
    ///
    /// * `system` - The aircraft system to trim.
    /// * `setpoints` - Must be a `DVector` of length 3 containing the setpoints for the trim problem:
    ///     * `setpoints[0]` - The desired velocity setpoint [m/s].
    ///     * `setpoints[1]` - The desired altitude setpoint [m].
    ///     * `setpoints[2]` - The desired gamma angle setpoint [deg].
    /// * `params` - The parameters for the trim problem.
    ///     * `params[0]` - The throttle [0.0 - 1.0].
    ///     * `params[1]` - The elevator [deg].
    ///     * `params[2]` - The alpha angle [rad].
    ///
    /// # Returns
    ///
    /// A tuple containing the state and input for the trim problem.
    fn setup(
        &self,
        _system: &Aircraft<A, E>,
        setpoints: &DVector<f64>,
        params: &DVector<f64>,
    ) -> Result<(FixedWing3DoFState, FixedWing3DoFInput), TrimError> {
        if setpoints.len() != 3 {
            return Err(TrimError::SetpointsError(
                "setpoints must have length 3".to_string(),
            ));
        }
        if params.len() != 3 {
            return Err(TrimError::ParamsError(
                "params must have length 3".to_string(),
            ));
        }

        let set_vt = setpoints[0];
        let set_altitude = setpoints[1];
        let set_gamma = setpoints[2] / RAD_TO_DEG;
        // Clean configuration (gear + flaps up)
        let u = FixedWing3DoFInput {
            input_vector: dvector![params[0], params[1], 0.25, 0.0],
        };
        let x = FixedWing3DoFState::new(dvector![
            set_vt,
            params[2],
            params[2] + set_gamma,
            0.0,
            set_altitude,
            0.0
        ]);
        Ok((x, u))
    }

    /// Cost for 3-DoF fixed-wing model is calculated as
    /// `vt^2 + 100 * alpha^2 + 10 * q^2`.
    fn cost(&self, x_dot: &FixedWing3DoFState) -> f64 {
        x_dot.vt().pow(2.0) + 100.0 * x_dot.alpha().pow(2.0) + 10.0 * x_dot.q().pow(2.0)
    }
}
