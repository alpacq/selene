use crate::model::aerodynamics::Aerodynamics;

pub struct TransportAero;

impl Aerodynamics for TransportAero {
    /// damping derivative CXq - CX derivative with respect to pitch rate q
    fn cxq(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// damping derivative CYr - CY derivative with respect to yaw rate r
    fn cyr(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// damping derivative CYp - CY derivative with respect to roll rate p
    fn cyp(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// damping derivative CZq - CZ derivative with respect to pitch rate q
    fn czq(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// damping derivative Clr - Cl(roll) derivative with respect to yaw rate r
    fn clr(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// damping derivative Clp - Cl(roll) derivative with respect to roll rate p
    fn clp(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// damping derivative Cmq - Cm(pitch) derivative with respect to pitch rate q
    fn cmq(&self, _alpha: f64) -> f64 {
        -16.0 // cmq
    }

    /// damping derivative Cnr - Cn(yaw) derivative with respect to yaw rate r
    fn cnr(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// damping derivative Cnp - Cn(yaw) derivative with respect to roll rate p
    fn cnp(&self, _alpha: f64) -> f64 {
        0.0
    }

    /// x-axis aerodynamic force coefficient
    fn cx(&self, _alpha: f64, _elevator: f64) -> f64 {
        0.042 // cdcls
    }

    /// y-axis aerodynamic force (sideforce) coefficient
    fn cy(&self, _beta: f64, _aileron: f64, _rudder: f64) -> f64 {
        0.0
    }

    /// z-axis aerodynamic force coefficient
    fn cz(&self, _alpha: f64, _beta: f64, _elevator: f64) -> f64 {
        0.085 //cla
    }

    /// pitching moment coefficient
    fn cm(&self, alpha: f64, elevator: f64) -> f64 {
        if alpha != 0.0 {
            -0.022 // cma
        } else if elevator != 0.0 {
            -0.016 //cmde
        } else {
            0.0
        }
    }

    /// rolling moment coefficient
    fn cl(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    /// yawing moment coefficient
    fn cn(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    /// rolling moment due to ailerons
    fn dlda(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    /// rolling moment due to rudder
    fn dldr(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    /// yawing moment due to ailerons
    fn dnda(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    /// yawing moment due to rudder
    fn dndr(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }
}
