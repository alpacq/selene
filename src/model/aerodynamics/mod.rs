pub mod f16aero;
pub mod transportaero;

/// Trait for aerodynamics dynamic properties
pub trait Aerodynamics {
    /// damping derivative CXq - CX derivative with respect to pitch rate q
    fn cxq(&self, alpha: f64) -> f64;

    /// damping derivative CYr - CY derivative with respect to yaw rate r
    fn cyr(&self, alpha: f64) -> f64;

    /// damping derivative CYp - CY derivative with respect to roll rate p
    fn cyp(&self, alpha: f64) -> f64;

    /// damping derivative CZq - CZ derivative with respect to pitch rate q
    fn czq(&self, alpha: f64) -> f64;

    /// damping derivative Clr - Cl(roll) derivative with respect to yaw rate r
    fn clr(&self, alpha: f64) -> f64;

    /// damping derivative Clp - Cl(roll) derivative with respect to roll rate p
    fn clp(&self, alpha: f64) -> f64;

    /// damping derivative Cmq - Cm(pitch) derivative with respect to pitch rate q
    fn cmq(&self, alpha: f64) -> f64;

    /// damping derivative Cnr - Cn(yaw) derivative with respect to yaw rate r
    fn cnr(&self, alpha: f64) -> f64;

    /// damping derivative Cnp - Cn(yaw) derivative with respect to roll rate p
    fn cnp(&self, alpha: f64) -> f64;

    /// x-axis aerodynamic force coefficient
    fn cx(&self, alpha: f64, elevator: f64) -> f64;

    /// y-axis aerodynamic force (sideforce) coefficient
    fn cy(&self, beta: f64, aileron: f64, rudder: f64) -> f64;

    /// z-axis aerodynamic force coefficient
    fn cz(&self, alpha: f64, beta: f64, elevator: f64) -> f64;

    /// pitching moment coefficient
    fn cm(&self, alpha: f64, elevator: f64) -> f64;

    /// rolling moment coefficient
    fn cl(&self, alpha: f64, beta: f64) -> f64;

    /// yawing moment coefficient
    fn cn(&self, alpha: f64, beta: f64) -> f64;

    /// rolling moment due to ailerons
    fn dlda(&self, alpha: f64, beta: f64) -> f64;

    /// rolling moment due to rudder
    fn dldr(&self, alpha: f64, beta: f64) -> f64;

    /// yawing moment due to ailerons
    fn dnda(&self, alpha: f64, beta: f64) -> f64;

    /// yawing moment due to rudder
    fn dndr(&self, alpha: f64, beta: f64) -> f64;
}
