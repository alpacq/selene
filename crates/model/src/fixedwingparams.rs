/// Parameters for a fixed-wing aircraft model
pub struct AircraftParams {
    /// Reference wing area, m²
    pub s: f64,
    /// Mean aerodynamic chord, m
    pub cbar: f64,
    /// Mass, kg
    pub mass: f64,
    /// Pitch-axis moment of inertia, kg·m²
    pub iyy: f64,
    /// Static thrust, N
    pub tstat: f64,
    /// Thrust derivative w.r.t. velocity, N/(m/s)
    pub dtdv: f64,
    /// Thrust line offset from CG, m
    pub ze: f64,
    /// Induced drag coefficient
    pub cdcls: f64,
    /// dCL/dalpha, 1/deg
    pub cla: f64,
    /// dCM/dalpha, 1/deg
    pub cma: f64,
    /// dCM/delevator, 1/deg
    pub cmde: f64,
    /// dCM/dq, 1/rad
    pub cmq: f64,
    /// dCM/dalpha_dot, 1/rad
    pub cmadot: f64,
    /// dCL/dalpha_dot, 1/rad
    pub cladot: f64,
}
