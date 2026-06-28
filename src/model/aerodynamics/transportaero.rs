use crate::model::aerodynamics::Aerodynamics;

pub struct TransportAero;

impl Aerodynamics for TransportAero {
    fn cxq(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn cyr(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn cyp(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn czq(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn clr(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn clp(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn cmq(&self, _alpha: f64) -> f64 {
        -16.0
    }

    fn cnr(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn cnp(&self, _alpha: f64) -> f64 {
        0.0
    }

    fn cx(&self, _alpha: f64, _elevator: f64) -> f64 {
        0.0
    }

    fn cy(&self, _beta: f64, _aileron: f64, _rudder: f64) -> f64 {
        0.0
    }

    fn cz(&self, _alpha: f64, _beta: f64, _elevator: f64) -> f64 {
        0.0
    }

    fn cm(&self, _alpha: f64, _elevator: f64) -> f64 {
        0.0
    }

    fn cl(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    fn cn(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    fn dlda(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    fn dldr(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    fn dnda(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }

    fn dndr(&self, _alpha: f64, _beta: f64) -> f64 {
        0.0
    }
}
