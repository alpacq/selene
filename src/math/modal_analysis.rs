use std::f64::consts::{LN_2, PI};

use nalgebra::{Complex, DMatrix, DVector};

const IMAG_EPSILON: f64 = 1e-6;
const SHORT_PERIOD_PHUGOID_SPLIT_S: f64 = 15.0;
const ROLL_SUBSIDENCE_MAX_TAU_S: f64 = 1.0;
const SPIRAL_MIN_TAU_S: f64 = 20.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    Longitudinal,
    Lateral,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModeKind {
    ShortPeriod,
    Phugoid,
    DutchRoll,
    RollSubsidence,
    Spiral,
    Unclassified,
}

#[derive(Debug, Clone, Copy)]
pub struct ModeAnalysis {
    pub eigenvalue: Complex<f64>,
    pub natural_frequency: f64,
    pub period: Option<f64>,
    pub damping_ratio: Option<f64>,
    pub doubling_time: Option<f64>,
    pub time_constant: Option<f64>,
    pub stable: bool,
    pub kind: ModeKind,
}

fn get_a_submatrix(a: &DMatrix<f64>, indices: &[usize]) -> DMatrix<f64> {
    a.select_rows(indices).select_columns(indices)
}

pub fn get_6dof_longitudal(a: &DMatrix<f64>) -> DMatrix<f64> {
    get_a_submatrix(a, &[0, 1, 4, 7])
}

pub fn get_6dof_lateral(a: &DMatrix<f64>) -> DMatrix<f64> {
    get_a_submatrix(a, &[2, 3, 6, 8])
}

fn calculate_mode_period(eigenvalue: Complex<f64>) -> f64 {
    2.0 * PI / eigenvalue.im.abs()
}

fn calculate_natural_frequency(eigenvalue: Complex<f64>) -> f64 {
    (eigenvalue.re * eigenvalue.re + eigenvalue.im * eigenvalue.im).sqrt()
}

fn calculate_damping_ratio(eigenvalue: Complex<f64>) -> f64 {
    -1.0 * eigenvalue.re / calculate_natural_frequency(eigenvalue)
}

fn calculate_time_constant(eigenvalue: Complex<f64>) -> f64 {
    -1.0 / eigenvalue.re
}

fn calculate_amplitude_doubling_time(eigenvalue: Complex<f64>) -> f64 {
    LN_2 / eigenvalue.re
}

fn determine_mode_stability(eigenvalue: Complex<f64>) -> bool {
    eigenvalue.re < 0.0
}

fn classify_complex_mode(axis: Axis, period: f64) -> ModeKind {
    match axis {
        Axis::Longitudinal if period < SHORT_PERIOD_PHUGOID_SPLIT_S => ModeKind::ShortPeriod,
        Axis::Longitudinal => ModeKind::Phugoid,
        Axis::Lateral => ModeKind::DutchRoll,
    }
}

fn classify_real_mode(axis: Axis, stable: bool, time_constant: f64) -> ModeKind {
    if axis != Axis::Lateral || !stable {
        return ModeKind::Unclassified;
    }
    if time_constant < ROLL_SUBSIDENCE_MAX_TAU_S {
        ModeKind::RollSubsidence
    } else if time_constant > SPIRAL_MIN_TAU_S {
        ModeKind::Spiral
    } else {
        ModeKind::Unclassified
    }
}

pub fn perform_modal_analysis(
    axis: Axis,
    eigenvalues: &DVector<Complex<f64>>,
) -> Vec<ModeAnalysis> {
    let mut modes = Vec::new();

    for eigenvalue in eigenvalues.iter().copied() {
        if eigenvalue.im.abs() < IMAG_EPSILON || eigenvalue.im < 0.0 {
            continue;
        }

        let stable = determine_mode_stability(eigenvalue);
        let period = calculate_mode_period(eigenvalue);
        let (damping_ratio, doubling_time) = if stable {
            (Some(calculate_damping_ratio(eigenvalue)), None)
        } else {
            (None, Some(calculate_amplitude_doubling_time(eigenvalue)))
        };

        modes.push(ModeAnalysis {
            eigenvalue,
            natural_frequency: calculate_natural_frequency(eigenvalue),
            period: Some(period),
            damping_ratio,
            doubling_time,
            time_constant: None,
            stable,
            kind: classify_complex_mode(axis, period),
        });
    }

    for eigenvalue in eigenvalues.iter().copied() {
        if eigenvalue.im.abs() >= IMAG_EPSILON {
            continue;
        }

        let stable = determine_mode_stability(eigenvalue);
        let time_constant = calculate_time_constant(eigenvalue);

        modes.push(ModeAnalysis {
            eigenvalue,
            natural_frequency: calculate_natural_frequency(eigenvalue),
            period: None,
            damping_ratio: None,
            doubling_time: None,
            time_constant: Some(time_constant),
            stable,
            kind: classify_real_mode(axis, stable, time_constant),
        });
    }

    modes
}

#[cfg(test)]
mod tests {
    use nalgebra::dvector;

    use crate::{
        linearize::linearize,
        math::test_utils::assert_approx,
        model::{F16, dynamicmodel::fixedwing6dof::FixedWing6DoF},
        trim::create_trim_problem_and_trim,
    };

    use super::*;

    /// Stevens & Lewis tabulate the F-16 Jacobian in English units (ft,
    /// slug), while this model works in SI (m, kg). `Vt` is the only
    /// longitudinal state carrying a length dimension — `Alpha`, `Theta`
    /// and `Q` are angles/angular rates and are unit-system independent.
    /// So any entry in the `Vt` row or column of the book's table must be
    /// rescaled before comparing against this model:
    /// - row `Vt`, other column: book value is an acceleration in ft/s²,
    ///   multiply by `FT_TO_M` to get m/s².
    /// - other row, column `Vt`: book value is "per ft/s", multiply by
    ///   `1 / FT_TO_M` to get "per m/s".
    /// - `[Vt, Vt]`, or any entry involving only angle-like states: the
    ///   length unit cancels out, so no conversion is needed.
    const FT_TO_M: f64 = 0.3048;

    /// Trims and linearizes the F-16 at 502 ft/s (≈153.01 m/s), straight and
    /// level, cg at 0.3c̄ — the flight condition of Stevens & Lewis Examples
    /// 3.8-1 (longitudinal A matrix) and 3.8-2 (lateral-directional A
    /// matrix).
    #[test]
    fn check_matrices() -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![
            153.01, // vt [m/s]
            0.0,    // altitude [m]
            0.0,    // gamma [deg]
            0.0,    // roll rate [rad/s]
            0.0,    // pitch rate [rad/s]
            0.0,    // turn rate [rad/s]
            0.0,    // phi [rad]
            0.0,    // coordinated turn flag
            0.30,   // Xcg
        ];
        let init_params = dvector![
            0.15,  // throttle
            -2.0,  // elevator
            0.039, // alpha
            0.0,   // aileron
            0.0,   // rudder
            0.0,   // beta
        ];

        let (x, u, _cost) =
            create_trim_problem_and_trim(F16::new(), FixedWing6DoF, setpoints, init_params)?;

        let linearized = linearize(F16::new(), FixedWing6DoF, x, u)?;

        let longitudal = get_6dof_longitudal(linearized.a());
        let lateral = get_6dof_lateral(linearized.a());

        // Example 3.8-1, longitudinal A matrix (states: Vt, Alpha, Theta, Q).
        assert_approx(longitudal[(0, 0)], -2.0244e-02, 1e-4, "A_long[Vt, Vt]");
        assert_approx(
            longitudal[(0, 1)],
            7.8763 * FT_TO_M,
            1e-2,
            "A_long[Vt, Alpha]",
        );
        assert_approx(
            longitudal[(0, 2)],
            -3.2170e+01 * FT_TO_M,
            2e-3,
            "A_long[Vt, Theta]",
        );
        assert_approx(
            longitudal[(0, 3)],
            -6.5020e-01 * FT_TO_M,
            1e-3,
            "A_long[Vt, Q]",
        );

        assert_approx(
            longitudal[(1, 0)],
            -2.5372e-04 / FT_TO_M,
            1e-5,
            "A_long[Alpha, Vt]",
        );
        assert_approx(longitudal[(1, 1)], -1.0190, 1e-3, "A_long[Alpha, Alpha]");
        assert_approx(longitudal[(1, 2)], 0.0, 1e-6, "A_long[Alpha, Theta]");
        assert_approx(longitudal[(1, 3)], 9.0484e-01, 1e-3, "A_long[Alpha, Q]");

        assert_approx(longitudal[(2, 0)], 0.0, 1e-6, "A_long[Theta, Vt]");
        assert_approx(longitudal[(2, 1)], 0.0, 1e-6, "A_long[Theta, Alpha]");
        assert_approx(longitudal[(2, 2)], 0.0, 1e-6, "A_long[Theta, Theta]");
        assert_approx(longitudal[(2, 3)], 1.0, 1e-6, "A_long[Theta, Q]");

        // Book value is 7.9472E-11 (ft-based) — numerical noise around zero.
        assert_approx(longitudal[(3, 0)], 0.0, 1e-6, "A_long[Q, Vt]");
        assert_approx(longitudal[(3, 1)], -2.4982, 1e-3, "A_long[Q, Alpha]");
        assert_approx(longitudal[(3, 2)], 0.0, 1e-6, "A_long[Q, Theta]");
        assert_approx(longitudal[(3, 3)], -1.3861, 1e-3, "A_long[Q, Q]");

        // Example 3.8-2, lateral-directional A matrix (states: Beta, Phi, P,
        // R). None of these states carry a length dimension, so the book
        // values compare directly — no unit conversion needed.
        assert_approx(lateral[(0, 0)], -3.2200e-01, 1e-3, "A_lat[Beta, Beta]");
        assert_approx(lateral[(0, 1)], 6.4032e-02, 1e-3, "A_lat[Beta, Phi]");
        assert_approx(lateral[(0, 2)], 3.8904e-02, 1e-3, "A_lat[Beta, P]");
        assert_approx(lateral[(0, 3)], -9.9156e-01, 1e-3, "A_lat[Beta, R]");

        assert_approx(lateral[(1, 0)], 0.0, 1e-6, "A_lat[Phi, Beta]");
        assert_approx(lateral[(1, 1)], 0.0, 1e-6, "A_lat[Phi, Phi]");
        assert_approx(lateral[(1, 2)], 1.0, 1e-6, "A_lat[Phi, P]");
        assert_approx(lateral[(1, 3)], 3.9385e-02, 1e-3, "A_lat[Phi, R]");

        assert_approx(lateral[(2, 0)], -3.0919e+01, 5e-2, "A_lat[P, Beta]");
        assert_approx(lateral[(2, 1)], 0.0, 1e-6, "A_lat[P, Phi]");
        assert_approx(lateral[(2, 2)], -3.6730, 1e-2, "A_lat[P, P]");
        assert_approx(lateral[(2, 3)], 6.7425e-01, 1e-2, "A_lat[P, R]");

        assert_approx(lateral[(3, 0)], 9.4724, 5e-2, "A_lat[R, Beta]");
        assert_approx(lateral[(3, 1)], 0.0, 1e-6, "A_lat[R, Phi]");
        assert_approx(lateral[(3, 2)], -2.6358e-02, 1e-3, "A_lat[R, P]");
        assert_approx(lateral[(3, 3)], -4.9849e-01, 1e-2, "A_lat[R, R]");

        Ok(())
    }
}
