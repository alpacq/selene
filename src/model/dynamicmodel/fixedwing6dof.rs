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
        self.state_vector[0]
    }

    /// angle of attack [rad]
    pub fn alpha(&self) -> f64 {
        self.state_vector[1]
    }

    /// sideslip angle [rad]
    pub fn beta(&self) -> f64 {
        self.state_vector[2]
    }

    /// roll angle [rad]
    pub fn phi(&self) -> f64 {
        self.state_vector[3]
    }

    /// pitch angle [rad]
    pub fn theta(&self) -> f64 {
        self.state_vector[4]
    }

    /// yaw angle [rad]
    pub fn psi(&self) -> f64 {
        self.state_vector[5]
    }

    /// roll rate [rad/s]
    pub fn p(&self) -> f64 {
        self.state_vector[6]
    }

    /// pitch rate [rad/s]
    pub fn q(&self) -> f64 {
        self.state_vector[7]
    }

    /// yaw rate [rad/s]
    pub fn r(&self) -> f64 {
        self.state_vector[8]
    }

    /// north position [m]
    pub fn pos_n(&self) -> f64 {
        self.state_vector[9]
    }

    /// east position [m]
    pub fn pos_e(&self) -> f64 {
        self.state_vector[10]
    }

    /// altitude [m]
    pub fn altitude(&self) -> f64 {
        self.state_vector[11]
    }

    /// power [0..100]
    pub fn power(&self) -> f64 {
        self.state_vector[12]
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
}

impl IntegrableState for FixedWing6DoFState {
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
        self.input_vector[0]
    }

    /// elevator position [-1..1]
    pub fn elevator(&self) -> f64 {
        self.input_vector[1]
    }

    /// aileron position [-1..1]
    pub fn aileron(&self) -> f64 {
        self.input_vector[2]
    }

    /// rudder position [-1..1]
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
}

/// A struct representing a 6-DoF fixed-wing aircraft dynamic model
pub struct FixedWing6DoF;

impl<A: Aerodynamics, E: Engine> DynamicModel<Aircraft<A, E>> for FixedWing6DoF {
    type State = FixedWing6DoFState;
    type Input = FixedWing6DoFInput;

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

    fn system_rank(&self) -> usize {
        13
    }
}

#[cfg(test)]
mod tests {
    use crate::model::F16;

    use super::*;
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

    fn assert_approx(actual: f64, expected: f64, name: &str) {
        assert!(
            (actual - expected).abs() < EPSILON,
            "{name}: expected {expected}, got {actual}"
        );
    }

    #[test]
    fn textbook_model_test_vt() {
        assert_approx(prepare_model().vt(), -25.74150, "vt_dot");
    }

    #[test]
    fn textbook_model_test_alpha() {
        assert_approx(prepare_model().alpha(), -0.8708620, "alpha_dot");
    }

    #[test]
    fn textbook_model_test_beta() {
        assert_approx(prepare_model().beta(), -0.4797399, "beta_dot");
    }

    #[test]
    fn textbook_model_test_phi() {
        assert_approx(prepare_model().phi(), 2.505734, "phi_dot");
    }

    #[test]
    fn textbook_model_test_theta() {
        assert_approx(prepare_model().theta(), 0.3250820, "theta_dot");
    }

    #[test]
    fn textbook_model_test_psi() {
        assert_approx(prepare_model().psi(), 2.145926, "psi_dot");
    }

    #[test]
    fn textbook_model_test_p() {
        assert_approx(prepare_model().p(), 12.62395, "p_dot");
    }

    #[test]
    fn textbook_model_test_q() {
        assert_approx(prepare_model().q(), 0.9648011, "q_dot");
    }

    #[test]
    fn textbook_model_test_r() {
        assert_approx(prepare_model().r(), 0.5809014, "r_dot");
    }

    #[test]
    fn textbook_model_test_altitude() {
        assert_approx(prepare_model().altitude(), 75.629, "altitude_dot");
    }

    #[test]
    fn textbook_model_test_power() {
        assert_approx(prepare_model().power(), -58.68999, "power_dot");
    }
}
