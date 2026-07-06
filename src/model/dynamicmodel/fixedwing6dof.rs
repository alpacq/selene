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

pub enum FixedWing6DoFStates {
    Vt,
    Alpha,
    Beta,
    Phi,
    Theta,
    Psi,
    P,
    Q,
    R,
    PosN,
    PosE,
    Altitude,
    Power,
}

/// A struct representing the state of a fixed-wing aircraft in 3D space
/// for 6-DoF full dynamic model
pub struct FixedWing6DoFState {
    state_vector: DVector<f64>,
}

impl FixedWing6DoFState {
    /// Creates a new FixedWing6DoFState with the given state vector
    pub fn new(state_vector: DVector<f64>) -> Self {
        Self { state_vector }
    }

    /// TAS [m/s]
    pub fn vt(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Vt as usize]
    }

    /// angle of attack [rad]
    pub fn alpha(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Alpha as usize]
    }

    /// sideslip angle [rad]
    pub fn beta(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Beta as usize]
    }

    /// roll angle [rad]
    pub fn phi(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Phi as usize]
    }

    /// pitch angle [rad]
    pub fn theta(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Theta as usize]
    }

    /// yaw angle [rad]
    pub fn psi(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Psi as usize]
    }

    /// roll rate [rad/s]
    pub fn p(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::P as usize]
    }

    /// pitch rate [rad/s]
    pub fn q(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Q as usize]
    }

    /// yaw rate [rad/s]
    pub fn r(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::R as usize]
    }

    /// north position [m]
    pub fn pos_n(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::PosN as usize]
    }

    /// east position [m]
    pub fn pos_e(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::PosE as usize]
    }

    /// altitude [m]
    pub fn altitude(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Altitude as usize]
    }

    /// power [0..100]
    pub fn power(&self) -> f64 {
        self.state_vector[FixedWing6DoFStates::Power as usize]
    }
}

impl SizedVector for FixedWing6DoFState {
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

    /// Creates a new FixedWing6DoFState from the given vector
    fn from_vector(vector: DVector<f64>) -> Self {
        FixedWing6DoFState::new(vector)
    }
}

/// A struct representing the input to a fixed-wing aircraft in 3D space
/// for 6-DoF full dynamic model
pub struct FixedWing6DoFInput {
    input_vector: DVector<f64>,
}

impl FixedWing6DoFInput {
    /// throttle position [0..1]
    pub fn throttle(&self) -> f64 {
        self.input_vector[0].clamp(0.0, 1.0)
    }

    /// elevator position [deg]
    pub fn elevator(&self) -> f64 {
        self.input_vector[1]
    }

    /// aileron position [deg]
    pub fn aileron(&self) -> f64 {
        self.input_vector[2]
    }

    /// rudder position [deg]
    pub fn rudder(&self) -> f64 {
        self.input_vector[3]
    }

    /// center of gravity on x-axis
    pub fn x_cg(&self) -> f64 {
        self.input_vector[4]
    }
}

impl SizedVector for FixedWing6DoFInput {
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

    /// Creates a new FixedWing6DoFInput from the given vector
    fn from_vector(vector: DVector<f64>) -> Self {
        FixedWing6DoFInput {
            input_vector: vector,
        }
    }
}

/// A struct representing a 6-DoF fixed-wing aircraft dynamic model
pub struct FixedWing6DoF;

impl<A: Aerodynamics, E: Engine> DynamicModel<Aircraft<A, E>> for FixedWing6DoF {
    type State = FixedWing6DoFState;
    type Input = FixedWing6DoFInput;

    /// State equations for the 6-DoF fixed-wing aircraft model.
    /// Full non-linear model
    fn state_equations(
        &self,
        system: &Aircraft<A, E>,
        x: &Self::State,
        u: &Self::Input,
    ) -> Self::State {
        // rigid body helper variables
        let xpq =
            system.airframe.ixz * (system.airframe.ixx - system.airframe.iyy + system.airframe.izz);
        let det =
            system.airframe.ixx * system.airframe.izz - (system.airframe.ixz * system.airframe.ixz);
        let xqr = system.airframe.izz * (system.airframe.izz - system.airframe.iyy)
            + (system.airframe.ixz * system.airframe.ixz);
        let zpq = (system.airframe.ixx - system.airframe.iyy) * system.airframe.ixx
            + (system.airframe.ixz * system.airframe.ixz);
        let ypr = system.airframe.izz - system.airframe.ixx;

        let alpha_deg = x.alpha() * RAD_TO_DEG;
        let beta_deg = x.beta() * RAD_TO_DEG;

        // values from atmosphere parameters

        let pressure = dynamic_pressure(x.vt(), x.altitude());
        let mach = mach(x.vt(), x.altitude());

        // engine modeling

        let set_power = system.engine.throttle_to_power(u.throttle());
        let power_dot = system.engine.power_dynamics(x.power(), set_power);
        let thrust = system.engine.thrust(x.power(), x.altitude(), mach);

        // control and state helper variables

        let aileron_deg = u.aileron() / 20.0;
        let rudder_deg = u.rudder() / 30.0;

        let tvt = 0.5 / x.vt();
        let b2v = system.airframe.b * tvt;
        let cq = system.airframe.cbar * x.q() * tvt;

        // aerodynamic properties with damping derivatives

        let cxt = system.aerodynamics.cx(alpha_deg, u.elevator())
            + cq * system.aerodynamics.cxq(alpha_deg);
        let cyt = system.aerodynamics.cy(beta_deg, u.aileron(), u.rudder())
            + b2v
                * (system.aerodynamics.cyr(alpha_deg) * x.r()
                    + system.aerodynamics.cyp(alpha_deg) * x.p());
        let czt = system.aerodynamics.cz(alpha_deg, beta_deg, u.elevator())
            + cq * system.aerodynamics.czq(alpha_deg);
        let clt = system.aerodynamics.cl(alpha_deg, beta_deg)
            + system.aerodynamics.dlda(alpha_deg, beta_deg) * aileron_deg
            + system.aerodynamics.dldr(alpha_deg, beta_deg) * rudder_deg
            + b2v
                * (system.aerodynamics.clr(alpha_deg) * x.r()
                    + system.aerodynamics.clp(alpha_deg) * x.p());
        let cmt = system.aerodynamics.cm(alpha_deg, u.elevator())
            + cq * system.aerodynamics.cmq(alpha_deg)
            + czt * (system.airframe.xcg - u.x_cg());
        let cnt = system.aerodynamics.cn(alpha_deg, beta_deg)
            + system.aerodynamics.dnda(alpha_deg, beta_deg) * aileron_deg
            + system.aerodynamics.dndr(alpha_deg, beta_deg) * rudder_deg
            + b2v
                * (system.aerodynamics.cnr(alpha_deg) * x.r()
                    + system.aerodynamics.cnp(alpha_deg) * x.p())
            - cyt * (system.airframe.xcg - u.x_cg()) * system.airframe.cbar / system.airframe.b;

        // helper variables for state equations

        let uu = x.vt() * x.alpha().cos() * x.beta().cos();
        let vv = x.vt() * x.beta().sin();
        let ww = x.vt() * x.alpha().sin() * x.beta().cos();
        let qs = pressure * system.airframe.s;
        let qsb = qs * system.airframe.b;
        let rmqs = qs / system.airframe.mass;

        let ay = rmqs * cyt;
        let az = rmqs * czt;

        // force equations

        let u_dot = x.r() * vv - x.q() * ww - GRAVITY * x.theta().sin()
            + (qs * cxt + thrust) / system.airframe.mass;
        let v_dot = x.p() * ww - x.r() * uu + GRAVITY * x.theta().cos() * x.phi().sin() + ay;
        let w_dot = x.q() * uu - x.p() * vv + GRAVITY * x.theta().cos() * x.phi().cos() + az;
        let dum = uu * uu + ww * ww;

        let vt_dot = (uu * u_dot + vv * v_dot + ww * w_dot) / x.vt();
        let alpha_dot = (uu * w_dot - ww * u_dot) / dum;
        let beta_dot = (x.vt() * v_dot - vv * vt_dot) * (x.beta().cos() / dum);

        // kinematic equations

        let phi_dot = x.p()
            + (x.theta().sin() / x.theta().cos()) * (x.q() * x.phi().sin() + x.r() * x.phi().cos());
        let theta_dot = x.q() * x.phi().cos() - x.r() * x.phi().sin();
        let psi_dot = (x.q() * x.phi().sin() + x.r() * x.phi().cos()) / x.theta().cos();

        // moments equations

        let roll = qsb * clt;
        let pitch = qs * system.airframe.cbar * cmt;
        let yaw = qsb * cnt;

        let p_dot = (xpq * x.p() * x.q() - xqr * x.q() * x.r()
            + system.airframe.izz * roll
            + system.airframe.ixz * (yaw + x.q() * system.engine.hx()))
            / det;
        let q_dot = (ypr * x.p() * x.r() - system.airframe.ixz * (x.p() * x.p() - x.r() * x.r())
            + pitch
            - x.r() * system.engine.hx())
            / system.airframe.iyy;
        let r_dot = (zpq * x.p() * x.q() - xpq * x.q() * x.r()
            + system.airframe.ixz * roll
            + system.airframe.ixx * (yaw + x.q() * system.engine.hx()))
            / det;

        // navigation equations

        let s1 = x.theta().cos() * x.psi().cos();
        let s2 = x.theta().cos() * x.psi().sin();
        let s3 = x.phi().sin() * x.psi().cos() * x.theta().sin() - x.phi().cos() * x.psi().sin();
        let s4 = x.phi().sin() * x.psi().sin() * x.theta().sin() + x.phi().cos() * x.psi().cos();
        let s5 = x.phi().sin() * x.theta().cos();
        let s6 = x.phi().cos() * x.theta().sin() * x.psi().cos() + x.phi().sin() * x.psi().sin();
        let s7 = x.phi().cos() * x.theta().sin() * x.psi().sin() - x.phi().sin() * x.psi().cos();
        let s8 = x.phi().cos() * x.theta().cos();

        let posn_dot = uu * s1 + vv * s3 + ww * s6;
        let pose_dot = uu * s2 + vv * s4 + ww * s7;
        let altitude_dot = uu * x.theta().sin() - vv * s5 - ww * s8;

        // TODO: telemetry output

        let _an = -az / GRAVITY;
        let _alat = ay / GRAVITY;

        FixedWing6DoFState::new(dvector![
            vt_dot,
            alpha_dot,
            beta_dot,
            phi_dot,
            theta_dot,
            psi_dot,
            p_dot,
            q_dot,
            r_dot,
            posn_dot,
            pose_dot,
            altitude_dot,
            power_dot
        ])
    }

    /// Returns the rank of the system (number of state variables)
    /// For this 6-DoF model, the rank is 13.
    fn system_rank(&self) -> usize {
        13
    }
}

impl<A: Aerodynamics, E: Engine> TrimTarget<Aircraft<A, E>> for FixedWing6DoF {
    /// Sets up the state and input for the given setpoints and parameters.
    ///
    /// # Arguments
    ///
    /// * `system` - The aircraft system to trim.
    /// * `setpoints` - Must be a `DVector` of length 8 containing the setpoints for the trim problem:
    ///     * `setpoints[0]` - The desired velocity setpoint [m/s].
    ///     * `setpoints[1]` - The desired altitude setpoint [m].
    ///     * `setpoints[2]` - The desired gamma angle setpoint [deg].
    ///     * `setpoints[3]` - The desired roll rate setpoint [rad/s].
    ///     * `setpoints[4]` - The desired pitch rate setpoint [rad/s].
    ///     * `setpoints[5]` - The desired turn rate setpoint [rad/s].
    ///     * `setpoints[6]` - The desired phi angle setpoint [rad].
    ///     * `setpoints[7]` - The coordinated turn 'boolean' (1.0 for coordinated turn).
    /// * `params` - The parameters for the trim problem.
    ///     * `params[0]` - The throttle [0.0 - 1.0].
    ///     * `params[1]` - The elevator [deg].
    ///     * `params[2]` - The alpha angle [rad]
    ///     * `params[3]` - The aileron [deg].
    ///     * `params[4]` - The rudder [deg].
    ///     * `params[5]` - The beta angle [rad].
    ///
    /// # Returns
    ///
    /// A tuple containing the state and input for the trim problem.
    fn setup(
        &self,
        system: &Aircraft<A, E>,
        setpoints: &DVector<f64>,
        params: &DVector<f64>,
    ) -> Result<(FixedWing6DoFState, FixedWing6DoFInput), TrimError> {
        if setpoints.len() < 8 {
            return Err(TrimError::SetpointsError(
                "setpoints must have length 8 or more".to_string(),
            ));
        }
        if params.len() != 6 {
            return Err(TrimError::ParamsError(
                "params must have length 6".to_string(),
            ));
        }

        // setpoints for steady state flight
        let set_vt = setpoints[0];
        let set_altitude = setpoints[1];
        let set_gamma = setpoints[2] / RAD_TO_DEG;
        let set_roll_rate = setpoints[3];
        let set_pitch_rate = setpoints[4];
        let set_turn_rate = setpoints[5];
        let phi = setpoints[6];

        // control inputs vector [throttle, elevator, aileron, rudder]
        let input_vector = if setpoints.len() > 8 {
            dvector![params[0], params[1], params[3], params[4], setpoints[8]]
        } else {
            dvector![params[0], params[1], params[3], params[4], 0.35]
        };
        let u = FixedWing6DoFInput { input_vector };

        // coordinated turn / skidding turn / non-turning flight
        let (beta, phi, theta, p, q, r) = if setpoints[7] == 1.0 {
            let alpha = params[2];
            let beta = params[5];
            let g = set_turn_rate * set_vt / GRAVITY;

            // turn-coordination constraint -> bank angle phi
            let a = 1.0 - g * alpha.tan() * beta.sin();
            let b = set_gamma.sin() / beta.cos();
            let c = 1.0 + g * g * beta.cos() * beta.cos();
            let phi = ((g
                * beta.cos()
                * ((a - b * b)
                    + b * alpha.tan()
                        * (c * (1.0 - b * b) + g * g * beta.sin() * beta.sin()).sqrt()))
                / (alpha.cos() * (a * a - b * b * (1.0 + c * alpha.tan() * alpha.tan()))))
            .atan();

            // rate-of-climb constraint -> pitch angle theta
            let a = alpha.cos() * beta.cos();
            let b = phi.sin() * beta.sin() + phi.cos() * alpha.sin() * beta.cos();
            let theta = ((a * b
                + set_gamma.sin() * (a * a - set_gamma.sin() * set_gamma.sin() + b * b).sqrt())
                / (a * a - set_gamma.sin() * set_gamma.sin()))
            .atan();

            // steady-turn body-axis rates
            let p = -set_turn_rate * theta.sin();
            let q = set_turn_rate * phi.sin() * theta.cos();
            let r = set_turn_rate * phi.cos() * theta.cos();
            (beta, phi, theta, p, q, r)
        } else if set_turn_rate != 0.0 {
            let phi = 0.0_f64;
            let sgcb = set_gamma.sin() / params[5].cos();
            let theta = params[2] + (sgcb / (1.0 - sgcb * sgcb).sqrt()).atan();
            let p = 0.0;
            let q = set_turn_rate * phi.sin() * theta.cos();
            let r = set_turn_rate * phi.cos() * theta.cos();
            (params[5], phi, theta, p, q, r)
        } else {
            let d = if phi != 0.0 { -params[2] } else { params[2] };
            let theta = if set_gamma.sin() != 0.0 {
                let sgcb = set_gamma.sin() / params[5].cos();
                d + (sgcb / (1.0 - sgcb * sgcb).sqrt()).atan()
            } else {
                d
            };
            (params[5], phi, theta, set_roll_rate, set_pitch_rate, 0.0)
        };

        let x = FixedWing6DoFState::new(dvector![
            set_vt,
            params[2],
            beta,
            phi,
            theta,
            0.0,
            p,
            q,
            r,
            0.0,
            0.0,
            set_altitude,
            system.engine.throttle_to_power(params[0])
        ]);
        Ok((x, u))
    }

    /// Cost for 6-DoF fixed-wing model is calculated as
    /// `vt^2 + 100 * (alpha^2 + beta^2) + 10 * (P^2 + Q^2 + R^2)`.
    fn cost(&self, x_dot: &FixedWing6DoFState) -> f64 {
        x_dot.vt().pow(2.0)
            + 100.0 * (x_dot.alpha().pow(2.0) + x_dot.beta().pow(2.0))
            + 10.0 * (x_dot.p().pow(2.0) + x_dot.q().pow(2.0) + x_dot.r().pow(2.0))
    }
}

#[cfg(test)]
mod tests {
    use crate::model::F16;

    use super::*;
    use crate::math::test_utils::assert_approx;
    use nalgebra::dvector;

    fn prepare_model() -> FixedWing6DoFState {
        let x = FixedWing6DoFState::new(dvector![
            152.4, 0.5, -0.2, -1.0, 1.0, -1.0, 0.7, -0.8, 0.9, 304.8, 274.32, 3048.0, 90.0
        ]);
        let u = FixedWing6DoFInput {
            input_vector: dvector![0.9, 20.0, -15.0, -20.0, 0.4],
        };

        let system = FixedWing6DoF {};

        system.state_equations(&F16::new(), &x, &u)
    }

    const EPSILON: f64 = 1e-3;

    #[test]
    fn textbook_model_test_vt() {
        assert_approx(prepare_model().vt(), -22.933428, EPSILON, "vt_dot");
    }

    #[test]
    fn textbook_model_test_alpha() {
        assert_approx(prepare_model().alpha(), -0.8813417, EPSILON, "alpha_dot");
    }

    #[test]
    fn textbook_model_test_beta() {
        assert_approx(prepare_model().beta(), -0.4760048, EPSILON, "beta_dot");
    }

    #[test]
    fn textbook_model_test_phi() {
        assert_approx(prepare_model().phi(), 2.505734, EPSILON, "phi_dot");
    }

    #[test]
    fn textbook_model_test_theta() {
        assert_approx(prepare_model().theta(), 0.3250820, EPSILON, "theta_dot");
    }

    #[test]
    fn textbook_model_test_psi() {
        assert_approx(prepare_model().psi(), 2.145926, EPSILON, "psi_dot");
    }

    #[test]
    fn textbook_model_test_p() {
        assert_approx(prepare_model().p(), 12.62395, EPSILON, "p_dot");
    }

    #[test]
    fn textbook_model_test_q() {
        assert_approx(prepare_model().q(), 0.9648011, EPSILON, "q_dot");
    }

    #[test]
    fn textbook_model_test_r() {
        assert_approx(prepare_model().r(), 0.5809014, EPSILON, "r_dot");
    }

    #[test]
    fn textbook_model_test_altitude() {
        assert_approx(prepare_model().altitude(), 75.629, EPSILON, "altitude_dot");
    }

    #[test]
    fn textbook_model_test_power() {
        assert_approx(prepare_model().power(), -58.68999, EPSILON, "power_dot");
    }

    /// Trimmed level flight at sea level (nominal cg) must match Stevens &
    /// Lewis Table 3.6-2. Speeds are converted from ft/s to m/s; throttle is
    /// dimensionless and AOA / elevator are in degrees.
    ///
    /// The lowest tabulated speeds (130-170 ft/s) drive AOA toward ~45 deg,
    /// the edge of the aero envelope where the book notes trimming becomes
    /// ill-conditioned, so they are excluded here.
    #[test]
    fn trimmed_level_flight_matches_textbook_table_3_6_2() {
        use crate::trim::TrimProblemBuilder;

        const FT: f64 = 0.3048;
        // (speed [ft/s], throttle [-], AOA [deg], elevator [deg])
        let rows = [
            (200.0, 0.287, 19.70, 0.723),
            (300.0, 0.122, 8.49, -0.591),
            (350.0, 0.107, 5.87, -0.539),
            (500.0, 0.137, 2.14, -0.756),
            (700.0, 0.282, 0.382, -0.900),
            (800.0, 0.378, -0.045, -0.943),
        ];

        for (v_fts, throttle, aoa, elevator) in rows {
            let vt = v_fts * FT;
            let setpoints = dvector![vt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
            // Seed near the tabulated AOA so the optimizer locks onto the
            // level-flight trim rather than a spurious local minimum.
            let init = dvector![0.3, 0.0, aoa / RAD_TO_DEG, 0.0, 0.0, 0.0];
            let problem = TrimProblemBuilder::new()
                .for_system(F16::new())
                .with_model(FixedWing6DoF)
                .with_setpoints(setpoints)
                .with_initial_params(init)
                .build();
            let (x, u, cost) = problem.trim().expect("trim returned an error");

            assert!(cost < 1e-6, "vt={v_fts} ft/s did not converge: cost={cost}");
            assert_approx(
                u.throttle(),
                throttle,
                5e-3,
                format!("throttle at {} m/s", vt).as_str(),
            );
            assert_approx(
                x.alpha() * RAD_TO_DEG,
                aoa,
                5e-2,
                format!("AOA at {} m/s", vt).as_str(),
            );
            assert_approx(
                u.elevator(),
                elevator,
                5e-2,
                format!("elevator at {} m/s", vt).as_str(),
            );
        }
    }
}
