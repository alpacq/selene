use crate::{airframeparams::AirframeParams, engine::Engine};

pub struct Aircraft<E: Engine> {
    pub airframe: AirframeParams,
    pub engine: E,
}
