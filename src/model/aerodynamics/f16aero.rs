use num_traits::pow::Pow;

use crate::{
    math::interpolation::{lut_interpolation_1d, lut_interpolation_2d},
    model::aerodynamics::Aerodynamics,
};

/// ------------------------------------------------------------
/// dynamic damping derivatives
/// Index: alpha 12 points, -10 to 45 deg
/// ------------------------------------------------------------
/// CXq — x-force due to pitch rate
const CXQ_LUT: [f64; 12] = [
    -0.267, -0.110, 0.308, 1.34, 2.08, 2.91, 2.76, 2.05, 1.50, 1.49, 1.83, 1.21,
];
/// CYr — y-force due to yaw rate
pub const CYR_LUT: [f64; 12] = [
    0.882, 0.852, 0.876, 0.958, 0.962, 0.974, 0.819, 0.483, 0.590, 1.21, -0.493, -1.04,
];
/// CYp — y-force due to roll rate
pub const CYP_LUT: [f64; 12] = [
    -0.108, -0.108, -0.188, 0.110, 0.258, 0.226, 0.344, 0.362, 0.611, 0.529, 0.298, -2.27,
];
/// CZq — z-force due to pitch rate
pub const CZQ_LUT: [f64; 12] = [
    -8.80, -25.8, -28.9, -31.4, -31.2, -30.7, -27.7, -28.2, -29.0, -29.8, -38.3, -35.3,
];
/// Clr — rolling moment due to yaw rate
pub const CLR_LUT: [f64; 12] = [
    -0.126, -0.026, 0.063, 0.113, 0.208, 0.230, 0.319, 0.437, 0.680, 0.100, 0.447, -0.330,
];
/// Clp — rolling moment due to roll rate
pub const CLP_LUT: [f64; 12] = [
    -0.360, -0.359, -0.443, -0.420, -0.383, -0.375, -0.329, -0.294, -0.230, -0.210, -0.120, -0.100,
];
/// Cmq — pitching moment due to pitch rate
pub const CMQ_LUT: [f64; 12] = [
    -7.21, -0.540, -5.23, -5.26, -6.11, -6.64, -5.69, -6.00, -6.20, -6.40, -6.60, -6.00,
];
/// Cnr — yawing moment due to yaw rate
pub const CNR_LUT: [f64; 12] = [
    -0.380, -0.363, -0.378, -0.386, -0.370, -0.453, -0.550, -0.582, -0.595, -0.637, -1.02, -0.840,
];
/// Cnp — yawing moment due to roll rate
pub const CNP_LUT: [f64; 12] = [
    0.061, 0.052, 0.052, -0.012, -0.013, -0.024, 0.050, 0.150, 0.130, 0.158, 0.240, 0.150,
];

/// ------------------------------------------------------------
/// CX — longitudinal force coefficient
///  Rows: alpha 12 points, -10 to 45 deg
/// Cols: elevator 5 points, -24 to 24 deg, step 12 deg
/// ------------------------------------------------------------
pub const CX: [&[f64]; 12] = [
    // EL=-24   EL=-12   EL=0    EL=12   EL=24
    &[-0.099, -0.048, -0.022, -0.040, -0.083], // alpha = -10 deg
    &[-0.081, -0.038, -0.020, -0.038, -0.073], // alpha =  -5 deg
    &[-0.081, -0.040, -0.021, -0.039, -0.076], // alpha =   0 deg
    &[-0.063, -0.021, -0.004, -0.025, -0.072], // alpha =   5 deg
    &[-0.025, 0.016, 0.032, 0.006, -0.046],    // alpha =  10 deg
    &[0.044, 0.083, 0.094, 0.062, 0.012],      // alpha =  15 deg
    &[0.097, 0.127, 0.128, 0.087, 0.024],      // alpha =  20 deg
    &[0.113, 0.137, 0.130, 0.085, 0.025],      // alpha =  25 deg
    &[0.145, 0.162, 0.154, 0.100, 0.043],      // alpha =  30 deg
    &[0.167, 0.177, 0.161, 0.110, 0.053],      // alpha =  35 deg
    &[0.174, 0.179, 0.155, 0.104, 0.047],      // alpha =  40 deg
    &[0.166, 0.167, 0.138, 0.091, 0.040],      // alpha =  45 deg
];

/// ------------------------------------------------------------
/// CZ — normal force coefficient (lift axis)
/// Rows: alpha 12 points, -10 to 45 deg
/// ------------------------------------------------------------
pub const CZ: [f64; 12] = [
    0.770, 0.241, -0.100, -0.416, -0.731, -1.053, -1.366, -1.646, -1.917, -2.120, -2.248, -2.229,
];

/// ------------------------------------------------------------
/// CM — pitching moment coefficient
/// Rows: alpha 12 points, -10 to 45 deg
/// Cols: elevator 5 points, -24 to 24 deg, step 12 deg
/// ------------------------------------------------------------
pub const CM: [&[f64]; 12] = [
    // EL=-24   EL=-12   EL=0    EL=12   EL=24
    &[0.205, 0.081, -0.046, -0.174, -0.259], // alpha = -10 deg
    &[0.168, 0.077, -0.020, -0.145, -0.202], // alpha =  -5 deg
    &[0.186, 0.107, -0.009, -0.121, -0.184], // alpha =   0 deg
    &[0.196, 0.110, -0.005, -0.127, -0.193], // alpha =   5 deg
    &[0.213, 0.110, -0.006, -0.129, -0.199], // alpha =  10 deg
    &[0.251, 0.141, 0.010, -0.102, -0.150],  // alpha =  15 deg
    &[0.245, 0.127, 0.006, -0.097, -0.160],  // alpha =  20 deg
    &[0.238, 0.119, -0.001, -0.113, -0.167], // alpha =  25 deg
    &[0.252, 0.133, 0.014, -0.087, -0.104],  // alpha =  30 deg
    &[0.231, 0.108, 0.000, -0.084, -0.076],  // alpha =  35 deg
    &[0.198, 0.081, -0.013, -0.069, -0.041], // alpha =  40 deg
    &[0.192, 0.093, 0.032, -0.006, -0.005],  // alpha =  45 deg
];

/// ------------------------------------------------------------
/// CL — rolling moment coefficient
/// Rows: alpha 12 points, -10 to 45 deg
/// Cols: |beta| 7 points, 0 to 30 deg, step 5 deg
/// ------------------------------------------------------------
pub const CL: [&[f64]; 12] = [
    // |beta|=0  5      10      15      20      25      30 deg
    &[0.000, -0.001, -0.003, -0.001, 0.000, 0.007, 0.009], // alpha = -10 deg
    &[0.000, -0.004, -0.009, -0.010, -0.010, -0.010, -0.011], // alpha =  -5 deg
    &[0.000, -0.008, -0.017, -0.020, -0.022, -0.023, -0.023], // alpha =   0 deg
    &[0.000, -0.012, -0.024, -0.030, -0.034, -0.034, -0.037], // alpha =   5 deg
    &[0.000, -0.016, -0.030, -0.039, -0.047, -0.049, -0.050], // alpha =  10 deg
    &[0.000, -0.019, -0.034, -0.044, -0.046, -0.046, -0.047], // alpha =  15 deg
    &[0.000, -0.020, -0.040, -0.050, -0.059, -0.068, -0.074], // alpha =  20 deg
    &[0.000, -0.020, -0.037, -0.049, -0.061, -0.071, -0.079], // alpha =  25 deg
    &[0.000, -0.015, -0.016, -0.023, -0.033, -0.060, -0.091], // alpha =  30 deg
    &[0.000, -0.008, -0.002, -0.006, -0.036, -0.058, -0.076], // alpha =  35 deg
    &[0.000, -0.013, -0.010, -0.014, -0.035, -0.062, -0.077], // alpha =  40 deg
    &[0.000, -0.015, -0.019, -0.027, -0.035, -0.059, -0.076], // alpha =  45 deg
];

/// ------------------------------------------------------------
/// CN — yawing moment coefficient
/// Rows: alpha 12 points, -10 to 45 deg
/// Cols: |beta| 7 points, 0 to 30 deg, step 5 deg
/// ------------------------------------------------------------
pub const CN: [&[f64]; 12] = [
    // |beta|=0   5      10      15      20      25      30 deg
    &[0.000, 0.018, 0.038, 0.056, 0.064, 0.074, 0.079], // alpha = -10 deg
    &[0.000, 0.019, 0.042, 0.057, 0.077, 0.086, 0.090], // alpha =  -5 deg
    &[0.000, 0.018, 0.042, 0.059, 0.076, 0.093, 0.106], // alpha =   0 deg
    &[0.000, 0.019, 0.042, 0.058, 0.074, 0.089, 0.106], // alpha =   5 deg
    &[0.000, 0.019, 0.043, 0.058, 0.073, 0.080, 0.096], // alpha =  10 deg
    &[0.000, 0.018, 0.039, 0.053, 0.057, 0.062, 0.080], // alpha =  15 deg
    &[0.000, 0.013, 0.030, 0.032, 0.029, 0.049, 0.068], // alpha =  20 deg
    &[0.000, 0.007, 0.017, 0.012, 0.007, 0.022, 0.030], // alpha =  25 deg
    &[0.000, 0.004, 0.004, 0.002, 0.012, 0.028, 0.064], // alpha =  30 deg
    &[0.000, -0.014, -0.035, -0.046, -0.034, -0.012, 0.015], // alpha =  35 deg
    &[0.000, -0.017, -0.047, -0.071, -0.065, -0.002, 0.011], // alpha =  40 deg
    &[0.000, -0.033, -0.057, -0.073, -0.041, -0.013, -0.001], // alpha =  45 deg
];

/// ------------------------------------------------------------
/// DLDA — rolling moment due to ailerons
/// Rows: alpha 12 points, -10 to 45 deg
/// Cols: beta 7 points, -30 to 30 deg, step 10 deg
/// ------------------------------------------------------------
pub const DLDA: [&[f64]; 12] = [
    // beta=-30   -20    -10      0     10     20     30 deg
    &[-0.041, -0.041, -0.042, -0.040, -0.043, -0.044, -0.043], // alpha = -10 deg
    &[-0.052, -0.053, -0.053, -0.052, -0.049, -0.048, -0.049], // alpha =  -5 deg
    &[-0.053, -0.053, -0.052, -0.051, -0.048, -0.048, -0.047], // alpha =   0 deg
    &[-0.056, -0.053, -0.051, -0.052, -0.049, -0.047, -0.045], // alpha =   5 deg
    &[-0.050, -0.050, -0.049, -0.048, -0.043, -0.042, -0.042], // alpha =  10 deg
    &[-0.056, -0.051, -0.049, -0.048, -0.042, -0.041, -0.037], // alpha =  15 deg
    &[-0.082, -0.066, -0.043, -0.042, -0.042, -0.020, -0.003], // alpha =  20 deg
    &[-0.059, -0.043, -0.035, -0.037, -0.036, -0.028, -0.013], // alpha =  25 deg
    &[-0.042, -0.038, -0.026, -0.031, -0.025, -0.013, -0.010], // alpha =  30 deg
    &[-0.038, -0.027, -0.016, -0.026, -0.021, -0.014, -0.003], // alpha =  35 deg
    &[-0.027, -0.023, -0.018, -0.017, -0.016, -0.011, -0.007], // alpha =  40 deg
    &[-0.017, -0.016, -0.014, -0.012, -0.011, -0.010, -0.008], // alpha =  45 deg
];

/// ------------------------------------------------------------
/// DLDR — rolling moment due to rudder
/// Rows: alpha 12 points, -10 to 45 deg
/// Cols: beta 7 points, -30 to 30 deg, step 10 deg
/// ------------------------------------------------------------
pub const DLDR: [&[f64]; 12] = [
    // beta=-30   -20    -10      0     10     20     30 deg
    &[0.005, 0.007, 0.013, 0.018, 0.015, 0.021, 0.023], // alpha = -10 deg
    &[0.017, 0.016, 0.013, 0.015, 0.014, 0.011, 0.010], // alpha =  -5 deg
    &[0.014, 0.014, 0.011, 0.015, 0.013, 0.010, 0.011], // alpha =   0 deg
    &[0.010, 0.014, 0.012, 0.014, 0.013, 0.011, 0.011], // alpha =   5 deg
    &[-0.005, 0.013, 0.011, 0.014, 0.012, 0.010, 0.011], // alpha =  10 deg
    &[0.009, 0.009, 0.009, 0.014, 0.011, 0.009, 0.010], // alpha =  15 deg
    &[0.019, 0.012, 0.008, 0.014, 0.011, 0.008, 0.008], // alpha =  20 deg
    &[0.005, 0.005, 0.005, 0.015, 0.010, 0.010, 0.010], // alpha =  25 deg
    &[-0.000, 0.000, -0.002, 0.013, 0.008, 0.006, 0.006], // alpha =  30 deg
    &[-0.005, 0.004, 0.005, 0.011, 0.008, 0.005, 0.014], // alpha =  35 deg
    &[-0.011, 0.009, 0.003, 0.006, 0.007, 0.000, 0.020], // alpha =  40 deg
    &[0.008, 0.007, 0.005, 0.001, 0.003, 0.001, 0.000], // alpha =  45 deg
];

/// ------------------------------------------------------------
/// DNDA — yawing moment due to ailerons
/// Rows: alpha 12 points, -10 to 45 deg
/// Cols: beta 7 points, -30 to 30 deg, step 10 deg
/// ------------------------------------------------------------
pub const DNDA: [&[f64]; 12] = [
    // beta=-30   -20    -10      0     10     20     30 deg
    &[0.001, 0.002, -0.006, -0.011, -0.015, -0.024, -0.022], // alpha = -10 deg
    &[-0.027, -0.014, -0.008, -0.011, -0.015, -0.010, 0.002], // alpha =  -5 deg
    &[-0.017, -0.016, -0.006, -0.010, -0.014, -0.004, -0.003], // alpha =   0 deg
    &[-0.013, -0.016, -0.006, -0.009, -0.012, -0.002, -0.005], // alpha =   5 deg
    &[-0.012, -0.014, -0.005, -0.008, -0.011, -0.001, -0.003], // alpha =  10 deg
    &[-0.016, -0.019, -0.008, -0.006, -0.008, 0.003, -0.001], // alpha =  15 deg
    &[0.001, -0.021, -0.005, 0.000, -0.002, 0.014, -0.009],  // alpha =  20 deg
    &[0.017, 0.002, 0.007, 0.004, 0.002, 0.006, -0.009],     // alpha =  25 deg
    &[0.011, 0.012, 0.004, 0.007, 0.006, -0.001, -0.001],    // alpha =  30 deg
    &[0.017, 0.015, 0.007, 0.010, 0.012, 0.004, 0.003],      // alpha =  35 deg
    &[0.008, 0.015, 0.006, 0.004, 0.011, 0.004, -0.002],     // alpha =  40 deg
    &[0.016, 0.011, 0.006, 0.010, 0.011, 0.006, 0.001],      // alpha =  45 deg
];

/// ------------------------------------------------------------
/// DNDR — yawing moment due to rudder
/// Rows: alpha 12 points, -10 to 45 deg
/// Cols: beta 7 points, -30 to 30 deg, step 10 deg
/// ------------------------------------------------------------
pub const DNDR: [&[f64]; 12] = [
    // beta=-30   -20    -10      0     10     20     30 deg
    &[-0.018, -0.028, -0.037, -0.048, -0.043, -0.052, -0.062], // alpha = -10 deg
    &[-0.052, -0.051, -0.041, -0.045, -0.044, -0.034, -0.034], // alpha =  -5 deg
    &[-0.052, -0.043, -0.038, -0.045, -0.041, -0.036, -0.027], // alpha =   0 deg
    &[-0.052, -0.046, -0.040, -0.045, -0.041, -0.036, -0.028], // alpha =   5 deg
    &[-0.054, -0.045, -0.040, -0.044, -0.040, -0.035, -0.027], // alpha =  10 deg
    &[-0.049, -0.049, -0.038, -0.045, -0.038, -0.028, -0.027], // alpha =  15 deg
    &[-0.059, -0.057, -0.037, -0.047, -0.034, -0.024, -0.023], // alpha =  20 deg
    &[-0.051, -0.052, -0.030, -0.048, -0.035, -0.023, -0.023], // alpha =  25 deg
    &[-0.030, -0.030, -0.027, -0.049, -0.035, -0.020, -0.019], // alpha =  30 deg
    &[-0.037, -0.033, -0.024, -0.045, -0.029, -0.016, -0.009], // alpha =  35 deg
    &[-0.026, -0.030, -0.019, -0.033, -0.022, -0.010, -0.025], // alpha =  40 deg
    &[-0.013, -0.008, -0.013, -0.016, -0.009, -0.014, -0.010], // alpha =  45 deg
];

pub struct F16Aero;

impl Aerodynamics for F16Aero {
    /// damping derivative CXq - CX derivative with respect to pitch rate q
    fn cxq(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CXQ_LUT, 2, 0.2)
    }

    /// damping derivative CYr - CY derivative with respect to yaw rate r
    fn cyr(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CYR_LUT, 2, 0.2)
    }

    /// damping derivative CYP - CY derivative with respect to roll rate p
    fn cyp(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CYP_LUT, 2, 0.2)
    }

    /// damping derivative CZq - CZ derivative with respect to pitch rate q
    fn czq(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CZQ_LUT, 2, 0.2)
    }

    /// damping derivative Clr - Cl(roll) derivative with respect to yaw rate r
    fn clr(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CLR_LUT, 2, 0.2)
    }

    /// damping derivative Clp - Cl(roll) derivative with respect to roll rate p
    fn clp(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CLP_LUT, 2, 0.2)
    }

    /// damping derivative Cmq - Cm(pitch) derivative with respect to pitch rate q
    fn cmq(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CMQ_LUT, 2, 0.2)
    }

    /// damping derivative Cnr - Cn(yaw) derivative with respect to yaw rate r
    fn cnr(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CNR_LUT, 2, 0.2)
    }

    /// damping derivative Cnp - Cn(yaw) derivative with respect to roll rate p
    fn cnp(&self, alpha: f64) -> f64 {
        lut_interpolation_1d(alpha, &CNP_LUT, 2, 0.2)
    }

    /// x-axis aerodynamic force coefficient
    fn cx(&self, alpha: f64, elevator: f64) -> f64 {
        lut_interpolation_2d(alpha, elevator, 2, 2, 0.2, 1.0 / 12.0, &CX)
    }

    /// y-axis aerodynamic force (sideforce) coefficient
    fn cy(&self, beta: f64, aileron: f64, rudder: f64) -> f64 {
        beta * -0.02 + (aileron / 20.0) * 0.021 + (rudder / 30.0) * 0.086
    }

    /// z-axis aerodynamic force coefficient
    fn cz(&self, alpha: f64, beta: f64, elevator: f64) -> f64 {
        let s = lut_interpolation_1d(alpha, &CZ, 2, 0.2);
        s * (1.0 - (beta / 57.3).pow(2.0)) - (elevator / 25.0) * 0.19
    }

    /// pitching moment coefficient
    fn cm(&self, alpha: f64, elevator: f64) -> f64 {
        lut_interpolation_2d(alpha, elevator, 2, 2, 0.2, 1.0 / 12.0, &CM)
    }

    /// rolling moment coefficient
    fn cl(&self, alpha: f64, beta: f64) -> f64 {
        let s = lut_interpolation_2d(alpha, beta.abs(), 2, 0, 0.2, 0.2, &CL);
        s * beta.signum()
    }

    /// yawing moment coefficient
    fn cn(&self, alpha: f64, beta: f64) -> f64 {
        let s = lut_interpolation_2d(alpha, beta.abs(), 2, 0, 0.2, 0.2, &CN);
        s * beta.signum()
    }

    /// rolling moment due to ailerons
    fn dlda(&self, alpha: f64, beta: f64) -> f64 {
        lut_interpolation_2d(alpha, beta, 2, 3, 0.2, 0.1, &DLDA)
    }

    /// rolling moment due to rudder
    fn dldr(&self, alpha: f64, beta: f64) -> f64 {
        lut_interpolation_2d(alpha, beta, 2, 3, 0.2, 0.1, &DLDR)
    }

    /// yawing moment due to ailerons
    fn dnda(&self, alpha: f64, beta: f64) -> f64 {
        lut_interpolation_2d(alpha, beta, 2, 3, 0.2, 0.1, &DNDA)
    }

    /// yawing moment due to rudder
    fn dndr(&self, alpha: f64, beta: f64) -> f64 {
        lut_interpolation_2d(alpha, beta, 2, 3, 0.2, 0.1, &DNDR)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::test_utils::assert_approx;

    const EPSILON: f64 = 1e-6;

    fn aero() -> F16Aero {
        F16Aero
    }

    // ----------------------------------------------------------------
    // CXQ — x-force due to pitch rate
    // ----------------------------------------------------------------

    #[test]
    fn cxq_node_at_0_deg() {
        // alpha=0 deg is a grid node, no interpolation
        assert_approx(aero().cxq(0.0), 0.308, EPSILON, "cxq node 0deg");
    }

    #[test]
    fn cxq_node_at_25_deg() {
        assert_approx(aero().cxq(25.0), 2.05, EPSILON, "cxq node 25deg");
    }

    #[test]
    fn cxq_midpoint_7_5_deg() {
        // midpoint between 5 deg (1.34) and 10 deg (2.08) => 1.71
        assert_approx(aero().cxq(7.5), 1.71, EPSILON, "cxq mid 7.5deg");
    }

    #[test]
    fn cxq_clamp_above_45_deg() {
        // alpha=50 clamps to 45 deg => 0.59
        assert_approx(aero().cxq(50.0), 0.590, EPSILON, "cxq clamp 50deg");
    }

    // ----------------------------------------------------------------
    // CYR — y-force due to yaw rate
    // ----------------------------------------------------------------

    #[test]
    fn cyr_node_at_0_deg() {
        assert_approx(aero().cyr(0.0), 0.876, EPSILON, "cyr node 0deg");
    }

    #[test]
    fn cyr_node_at_25_deg() {
        assert_approx(aero().cyr(25.0), 0.483, EPSILON, "cyr node 25deg");
    }

    #[test]
    fn cyr_clamp_above_45_deg() {
        assert_approx(aero().cyr(50.0), -1.587, EPSILON, "cyr clamp 50deg");
    }

    // ----------------------------------------------------------------
    // CYP — y-force due to roll rate
    // ----------------------------------------------------------------

    #[test]
    fn cyp_node_at_0_deg() {
        assert_approx(aero().cyp(0.0), -0.188, EPSILON, "cyp node 0deg");
    }

    #[test]
    fn cyp_node_at_25_deg() {
        assert_approx(aero().cyp(25.0), 0.362, EPSILON, "cyp node 25deg");
    }

    #[test]
    fn cyp_midpoint_7_5_deg() {
        // midpoint between 5 deg (0.110) and 10 deg (0.258) => 0.184
        assert_approx(aero().cyp(7.5), 0.184, EPSILON, "cyp mid 7.5deg");
    }

    // ----------------------------------------------------------------
    // CZQ — z-force due to pitch rate
    // ----------------------------------------------------------------

    #[test]
    fn czq_node_at_0_deg() {
        assert_approx(aero().czq(0.0), -28.9, EPSILON, "czq node 0deg");
    }

    #[test]
    fn czq_node_at_25_deg() {
        assert_approx(aero().czq(25.0), -28.2, EPSILON, "czq node 25deg");
    }

    #[test]
    fn czq_midpoint_7_5_deg() {
        // midpoint between 5 deg (-31.4) and 10 deg (-31.2) => -31.3
        assert_approx(aero().czq(7.5), -31.3, EPSILON, "czq mid 7.5deg");
    }

    // ----------------------------------------------------------------
    // CLR — rolling moment due to yaw rate
    // ----------------------------------------------------------------

    #[test]
    fn clr_node_at_0_deg() {
        assert_approx(aero().clr(0.0), 0.063, EPSILON, "clr node 0deg");
    }

    #[test]
    fn clr_node_at_25_deg() {
        assert_approx(aero().clr(25.0), 0.437, EPSILON, "clr node 25deg");
    }

    #[test]
    fn clr_midpoint_7_5_deg() {
        // midpoint between 5 deg (0.113) and 10 deg (0.208) => 0.1605
        assert_approx(aero().clr(7.5), 0.1605, EPSILON, "clr mid 7.5deg");
    }

    // ----------------------------------------------------------------
    // CLP — rolling moment due to roll rate
    // ----------------------------------------------------------------

    #[test]
    fn clp_node_at_0_deg() {
        assert_approx(aero().clp(0.0), -0.443, EPSILON, "clp node 0deg");
    }

    #[test]
    fn clp_node_at_25_deg() {
        assert_approx(aero().clp(25.0), -0.294, EPSILON, "clp node 25deg");
    }

    #[test]
    fn clp_midpoint_7_5_deg() {
        // midpoint between 5 deg (-0.420) and 10 deg (-0.383) => -0.4015
        assert_approx(aero().clp(7.5), -0.4015, EPSILON, "clp mid 7.5deg");
    }

    // ----------------------------------------------------------------
    // CMQ — pitching moment due to pitch rate
    // ----------------------------------------------------------------

    #[test]
    fn cmq_node_at_0_deg() {
        assert_approx(aero().cmq(0.0), -5.23, EPSILON, "cmq node 0deg");
    }

    #[test]
    fn cmq_node_at_25_deg() {
        assert_approx(aero().cmq(25.0), -6.0, EPSILON, "cmq node 25deg");
    }

    #[test]
    fn cmq_midpoint_7_5_deg() {
        // midpoint between 5 deg (-5.26) and 10 deg (-6.11) => -5.685
        assert_approx(aero().cmq(7.5), -5.685, EPSILON, "cmq mid 7.5deg");
    }

    // ----------------------------------------------------------------
    // CNR — yawing moment due to yaw rate
    // ----------------------------------------------------------------

    #[test]
    fn cnr_node_at_0_deg() {
        assert_approx(aero().cnr(0.0), -0.378, EPSILON, "cnr node 0deg");
    }

    #[test]
    fn cnr_node_at_25_deg() {
        assert_approx(aero().cnr(25.0), -0.582, EPSILON, "cnr node 25deg");
    }

    #[test]
    fn cnr_clamp_above_45_deg() {
        assert_approx(aero().cnr(50.0), -0.660, EPSILON, "cnr clamp 50deg");
    }

    // ----------------------------------------------------------------
    // CNP — yawing moment due to roll rate
    // ----------------------------------------------------------------

    #[test]
    fn cnp_node_at_0_deg() {
        assert_approx(aero().cnp(0.0), 0.052, EPSILON, "cnp node 0deg");
    }

    #[test]
    fn cnp_node_at_25_deg() {
        assert_approx(aero().cnp(25.0), 0.150, EPSILON, "cnp node 25deg");
    }

    #[test]
    fn cnp_midpoint_7_5_deg() {
        // midpoint between 5 deg (-0.012) and 10 deg (-0.013) => -0.0125
        assert_approx(aero().cnp(7.5), -0.0125, EPSILON, "cnp mid 7.5deg");
    }

    // ----------------------------------------------------------------
    // CX — longitudinal force coefficient
    // ----------------------------------------------------------------

    #[test]
    fn cx_node_alpha0_el0() {
        // both axes on grid nodes => exact table value
        assert_approx(aero().cx(0.0, 0.0), -0.021, EPSILON, "cx node alpha=0 el=0");
    }

    #[test]
    fn cx_node_alpha20_el12() {
        assert_approx(
            aero().cx(20.0, 12.0),
            0.087,
            EPSILON,
            "cx node alpha=20 el=12",
        );
    }

    #[test]
    fn cx_midpoint_alpha7_5_el6() {
        // alpha midpoint between 5 and 10 deg, el midpoint between 0 and 12 deg
        // => 0.00225
        assert_approx(
            aero().cx(7.5, 6.0),
            0.00225,
            EPSILON,
            "cx mid alpha=7.5 el=6",
        );
    }

    #[test]
    fn cx_clamp_alpha50_el30() {
        // alpha clamps to 45 deg, el clamps to 24 deg => 0.0105
        assert_approx(
            aero().cx(50.0, 30.0),
            0.0105,
            EPSILON,
            "cx clamp alpha=50 el=30",
        );
    }

    // ----------------------------------------------------------------
    // CY — sideforce (analytic formula, no LUT)
    // ----------------------------------------------------------------

    #[test]
    fn cy_all_zero() {
        assert_approx(aero().cy(0.0, 0.0, 0.0), 0.0, EPSILON, "cy all zero");
    }

    #[test]
    fn cy_positive_inputs() {
        // beta=5, ail=10, rdr=15 => -0.1 + 0.0105 + 0.043 = -0.0465
        assert_approx(
            aero().cy(5.0, 10.0, 15.0),
            -0.0465,
            EPSILON,
            "cy positive inputs",
        );
    }

    #[test]
    fn cy_mixed_sign_inputs() {
        // beta=-10, ail=-20, rdr=30 => 0.2 - 0.021 + 0.086 = 0.265
        assert_approx(
            aero().cy(-10.0, -20.0, 30.0),
            0.265,
            EPSILON,
            "cy mixed sign inputs",
        );
    }

    // ----------------------------------------------------------------
    // CZ — normal force coefficient
    // ----------------------------------------------------------------

    #[test]
    fn cz_node_alpha0_beta0_el0() {
        // base value from LUT at alpha=0 is -0.100, no beta/el correction
        assert_approx(
            aero().cz(0.0, 0.0, 0.0),
            -0.100,
            EPSILON,
            "cz node alpha=0 beta=0 el=0",
        );
    }

    #[test]
    fn cz_with_beta_and_elevator() {
        // alpha=10: base=-0.731, beta=5: factor=(1-(5/57.3)^2)=0.992394, el=12
        // => -0.731 * 0.992394 - (12/25)*0.19 = -0.816634
        assert_approx(
            aero().cz(10.0, 5.0, 12.0),
            -0.816634,
            EPSILON,
            "cz alpha=10 beta=5 el=12",
        );
    }

    #[test]
    fn cz_clamp_alpha50() {
        // alpha=50 clamps to 45 deg => base=-2.21, no beta/el correction
        assert_approx(
            aero().cz(50.0, 0.0, 0.0),
            -2.21,
            EPSILON,
            "cz clamp alpha=50",
        );
    }

    // ----------------------------------------------------------------
    // CM — pitching moment coefficient
    // ----------------------------------------------------------------

    #[test]
    fn cm_node_alpha0_el0() {
        assert_approx(aero().cm(0.0, 0.0), -0.009, EPSILON, "cm node alpha=0 el=0");
    }

    #[test]
    fn cm_node_alpha20_el12() {
        assert_approx(
            aero().cm(20.0, 12.0),
            -0.097,
            EPSILON,
            "cm node alpha=20 el=12",
        );
    }

    #[test]
    fn cm_midpoint_alpha7_5_el6() {
        // => -0.06675
        assert_approx(
            aero().cm(7.5, 6.0),
            -0.06675,
            EPSILON,
            "cm mid alpha=7.5 el=6",
        );
    }

    // ----------------------------------------------------------------
    // CL — rolling moment coefficient
    // NOTE: cl() adds beta.signum() — verify this is intentional,
    // as the original Fortran adds signum only to the interpolated
    // aerodynamic value, not to a unit offset.
    // ----------------------------------------------------------------

    #[test]
    fn cl_zero_beta_returns_zero_aero_plus_signum() {
        // beta=0 => signum=0, table value at alpha=0/beta=0 is 0.0 => total 0.0
        assert_approx(aero().cl(0.0, 0.0), 0.0, EPSILON, "cl beta=0");
    }

    #[test]
    fn cl_positive_beta_node() {
        // alpha=20, beta=10 => table s=-0.040, sign=+1 => total=-0.040
        assert_approx(
            aero().cl(20.0, 10.0),
            -0.040,
            EPSILON,
            "cl alpha=20 beta=10",
        );
    }

    #[test]
    fn cl_negative_beta_node() {
        // alpha=10, beta=-10 => table s=-0.030, sign=-1 => total=+0.030
        assert_approx(
            aero().cl(10.0, -10.0),
            0.030,
            EPSILON,
            "cl alpha=10 beta=-10",
        );
    }

    // ----------------------------------------------------------------
    // CN — yawing moment coefficient
    // ----------------------------------------------------------------

    #[test]
    fn cn_zero_beta() {
        // beta=0 => signum=0, table value at alpha=0/|beta|=0 is 0.0
        assert_approx(aero().cn(0.0, 0.0), 0.0, EPSILON, "cn beta=0");
    }

    #[test]
    fn cn_positive_beta_node() {
        // alpha=15, beta=15 => table s=0.053, sign=+1 => total=0.053
        assert_approx(aero().cn(15.0, 15.0), 0.053, EPSILON, "cn alpha=15 beta=15");
    }

    #[test]
    fn cn_negative_beta_node() {
        // alpha=0, beta=-10 => table s=0.042, sign=-1 => total=-0.042
        assert_approx(
            aero().cn(0.0, -10.0),
            -0.042,
            EPSILON,
            "cn alpha=0 beta=-10",
        );
    }

    // ----------------------------------------------------------------
    // DLDA — rolling moment due to ailerons
    // ----------------------------------------------------------------

    #[test]
    fn dlda_node_alpha0_beta0() {
        assert_approx(
            aero().dlda(0.0, 0.0),
            -0.051,
            EPSILON,
            "dlda node alpha=0 beta=0",
        );
    }

    #[test]
    fn dlda_node_alpha20_beta10() {
        assert_approx(
            aero().dlda(20.0, 10.0),
            -0.042,
            EPSILON,
            "dlda node alpha=20 beta=10",
        );
    }

    #[test]
    fn dlda_node_alpha0_neg_beta() {
        // negative beta grid node
        assert_approx(
            aero().dlda(0.0, -10.0),
            -0.052,
            EPSILON,
            "dlda node alpha=0 beta=-10",
        );
    }

    #[test]
    fn dlda_midpoint_alpha7_5_beta5() {
        // => -0.048
        assert_approx(
            aero().dlda(7.5, 5.0),
            -0.048,
            EPSILON,
            "dlda mid alpha=7.5 beta=5",
        );
    }

    // ----------------------------------------------------------------
    // DLDR — rolling moment due to rudder
    // ----------------------------------------------------------------

    #[test]
    fn dldr_node_alpha0_beta0() {
        assert_approx(
            aero().dldr(0.0, 0.0),
            0.015,
            EPSILON,
            "dldr node alpha=0 beta=0",
        );
    }

    #[test]
    fn dldr_node_alpha10_neg_beta20() {
        assert_approx(
            aero().dldr(10.0, -20.0),
            0.013,
            EPSILON,
            "dldr node alpha=10 beta=-20",
        );
    }

    #[test]
    fn dldr_midpoint_alpha12_5_beta15() {
        // alpha midpoint between 10 and 15 deg, beta midpoint between 10 and 20 => 0.0105
        assert_approx(
            aero().dldr(12.5, 15.0),
            0.0105,
            EPSILON,
            "dldr mid alpha=12.5 beta=15",
        );
    }

    // ----------------------------------------------------------------
    // DNDA — yawing moment due to ailerons
    // ----------------------------------------------------------------

    #[test]
    fn dnda_node_alpha0_beta0() {
        assert_approx(
            aero().dnda(0.0, 0.0),
            -0.010,
            EPSILON,
            "dnda node alpha=0 beta=0",
        );
    }

    #[test]
    fn dnda_node_alpha30_beta20() {
        assert_approx(
            aero().dnda(30.0, 20.0),
            -0.001,
            EPSILON,
            "dnda node alpha=30 beta=20",
        );
    }

    #[test]
    fn dnda_midpoint_alpha17_5_neg_beta10() {
        // alpha midpoint between 15 and 20, beta=-10 node => -0.0065
        assert_approx(
            aero().dnda(17.5, -10.0),
            -0.0065,
            EPSILON,
            "dnda mid alpha=17.5 beta=-10",
        );
    }

    // ----------------------------------------------------------------
    // DNDR — yawing moment due to rudder
    // ----------------------------------------------------------------

    #[test]
    fn dndr_node_alpha0_beta0() {
        assert_approx(
            aero().dndr(0.0, 0.0),
            -0.045,
            EPSILON,
            "dndr node alpha=0 beta=0",
        );
    }

    #[test]
    fn dndr_node_alpha25_neg_beta30() {
        assert_approx(
            aero().dndr(25.0, -30.0),
            -0.051,
            EPSILON,
            "dndr node alpha=25 beta=-30",
        );
    }

    #[test]
    fn dndr_midpoint_alpha22_5_beta10() {
        // alpha midpoint between 20 and 25 deg, beta=10 node => -0.0345
        assert_approx(
            aero().dndr(22.5, 10.0),
            -0.0345,
            EPSILON,
            "dndr mid alpha=22.5 beta=10",
        );
    }
}
