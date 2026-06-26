use crate::math::SizedVector;

/// Trait for aerodynamics dynamic properties
pub trait Aerodynamics {
    /// Returns the damping derivatives
    /// in form of a 9-element vector
    /// [0] - CXq - CX derivative with respect to pitch rate q
    /// [1] - CYr - CY derivative with respect to yaw rate r
    /// [2] - CYp - CY derivative with respect to roll rate p
    /// [3] - CZq - CZ derivative with respect to pitch rate q
    /// [4] - Clr - Cl(roll) derivative with respect to yaw rate r
    /// [5] - Clp - Cl(roll) derivative with respect to roll rate p
    /// [6] - Cmq - Cm(pitch) derivative with respect to pitch rate q
    /// [7] - Cnr - Cn(yaw) derivative with respect to yaw rate r
    /// [8] - Cnp - Cn(yaw) derivative with respect to roll rate p
    fn damping_derivatives(&self) -> impl SizedVector;

    /// x-axis aerodynamic force coefficient
    fn cx(&self) -> f64;

    /// y-axis aerodynamic force (sideforce) coefficient
    fn cy(&self) -> f64;

    /// z-axis aerodynamic force coefficient
    fn cz(&self) -> f64;

    /// pitching moment coefficient
    fn cm(&self) -> f64;

    /// rolling moment coefficient
    fn cl(&self) -> f64;

    /// yawing moment coefficient
    fn cn(&self) -> f64;

    /// rolling moment due to ailerons
    fn dlda(&self) -> f64;

    /// rolling moment due to rudder
    fn dldr(&self) -> f64;

    /// yawing moment due to ailerons
    fn dnda(&self) -> f64;

    /// yawing moment due to rudder
    fn dndr(&self) -> f64;
}
