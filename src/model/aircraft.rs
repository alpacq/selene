use crate::model::{aerodynamics::Aerodynamics, airframeparams::AirframeParams, engine::Engine};

/// A struct representing an aircraft with an airframe, aerodynamic properties and engine
pub struct Aircraft<A: Aerodynamics, E: Engine> {
    /// The airframe parameters of the aircraft
    pub airframe: AirframeParams,
    /// The aerodynamic properties of the aircraft
    pub aerodynamics: A,
    /// The engine of the aircraft
    pub engine: E,
}
