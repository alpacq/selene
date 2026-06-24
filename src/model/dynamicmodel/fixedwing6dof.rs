use crate::{
    math::{IntegrableState, SizedVector, input::Input, state::State},
    model::{
        DynamicModel,
        aircraft::Aircraft,
        atmosphere::{dynamic_pressure, mach},
        engine::Engine,
    },
};
use nalgebra::{DVector, dvector};

pub struct FixedWing6DoFState(State);

impl FixedWing6DoFState {
    pub fn new(state_vector: DVector<f64>) -> Self {
        Self(State::new(state_vector))
    }

    /// TAS [m/s]
    pub fn vt(&self) -> f64 {
        self.0.state_vector[0]
    }

    /// angle of attack [rad]
    pub fn alpha(&self) -> f64 {
        self.0.state_vector[1]
    }

    /// sideslip angle [rad]
    pub fn beta(&self) -> f64 {
        self.0.state_vector[2]
    }

    /// roll angle [rad]
    pub fn phi(&self) -> f64 {
        self.0.state_vector[3]
    }

    /// pitch angle [rad]
    pub fn theta(&self) -> f64 {
        self.0.state_vector[4]
    }

    /// yaw angle [rad]
    pub fn psi(&self) -> f64 {
        self.0.state_vector[5]
    }

    /// roll rate [rad/s]
    pub fn p(&self) -> f64 {
        self.0.state_vector[6]
    }

    /// pitch rate [rad/s]
    pub fn q(&self) -> f64 {
        self.0.state_vector[7]
    }

    /// yaw rate [rad/s]
    pub fn r(&self) -> f64 {
        self.0.state_vector[8]
    }

    /// north position [m]
    pub fn pos_n(&self) -> f64 {
        self.0.state_vector[9]
    }

    /// east position [m]
    pub fn pos_e(&self) -> f64 {
        self.0.state_vector[10]
    }

    /// altitude [m]
    pub fn altitude(&self) -> f64 {
        self.0.state_vector[11]
    }

    /// power [0..100]
    pub fn power(&self) -> f64 {
        self.0.state_vector[12]
    }
}

impl SizedVector for FixedWing6DoFState {
    /// Returns the size of the state vector
    ///
    /// # Returns
    ///
    /// The size of the state vector.
    fn size(&self) -> usize {
        self.0.size()
    }

    /// Returns the state vector.
    ///
    /// # Returns
    ///
    /// The state vector.
    fn vector(&self) -> DVector<f64> {
        self.0.vector()
    }
}

impl IntegrableState for FixedWing6DoFState {
    fn from_vector(vector: DVector<f64>) -> Self {
        FixedWing6DoFState::new(vector)
    }
}

pub struct FixedWing6DoFInput(Input);

impl FixedWing6DoFInput {
    /// throttle position [0..1]
    pub fn throttle(&self) -> f64 {
        self.0.input_vector[0]
    }

    /// elevator position [-1..1]
    pub fn elevator(&self) -> f64 {
        self.0.input_vector[1]
    }

    /// aileron position [-1..1]
    pub fn aileron(&self) -> f64 {
        self.0.input_vector[2]
    }

    /// rudder position [-1..1]
    pub fn rudder(&self) -> f64 {
        self.0.input_vector[3]
    }
}

impl SizedVector for FixedWing6DoFInput {
    /// Returns the size of the input vector
    ///
    /// # Returns
    ///
    /// The size of the input vector.
    fn size(&self) -> usize {
        self.0.size()
    }

    /// Returns the input vector.
    ///
    /// # Returns
    ///
    /// The input vector.
    fn vector(&self) -> DVector<f64> {
        self.0.vector()
    }
}

pub struct FixedWing6DoF;

impl<E: Engine> DynamicModel<Aircraft<E>> for FixedWing6DoF {
    type State = FixedWing6DoFState;
    type Input = FixedWing6DoFInput;

    fn state_equations(
        &self,
        system: &Aircraft<E>,
        x: &Self::State,
        _u: &Self::Input,
    ) -> Self::State {
        // temporary helper variables
        let xpq =
            system.airframe.ixz * (system.airframe.ixx - system.airframe.iyy + system.airframe.izz);
        let det =
            system.airframe.ixx * system.airframe.izz - (system.airframe.ixz * system.airframe.ixz);
        let xqr = system.airframe.izz * (system.airframe.izz - system.airframe.iyy)
            + (system.airframe.ixz * system.airframe.ixz);
        let zpq = (system.airframe.ixx - system.airframe.iyy) * system.airframe.ixx
            + (system.airframe.ixz * system.airframe.ixz);
        let ypr = system.airframe.izz - system.airframe.ixx;

        let pressure = dynamic_pressure(x.vt(), x.altitude());
        let mach = mach(x.vt(), x.altitude());

        FixedWing6DoFState::new(dvector![0.0])
    }

    fn system_rank(&self) -> usize {
        13
    }
}
