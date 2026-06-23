/// Parameters of the engine
pub struct EngineParams {
    /// engine angular momentum [kg * m^2 / s]
    pub hx: f64,
    /// Static thrust, N
    pub tstat: f64,
    /// Thrust derivative w.r.t. velocity, N/(m/s)
    pub dtdv: f64,
}
