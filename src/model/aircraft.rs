use crate::model::{airframeparams::AirframeParams, engine::Engine};

/// A struct representing an aircraft with an airframe and engine
pub struct Aircraft<E: Engine> {
    /// The airframe parameters of the aircraft
    pub airframe: AirframeParams,
    /// The engine of the aircraft
    pub engine: E,
}
