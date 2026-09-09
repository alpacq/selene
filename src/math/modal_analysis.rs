use nalgebra::DMatrix;

/// Extracts the square submatrix formed by the given state indices, keeping
/// only the rows *and* columns that belong to the reduced state set.
///
/// This is valid because a finite-difference partial derivative `∂ẋᵢ/∂xⱼ` is
/// evaluated independently of the other state dimensions, so slicing the
/// full Jacobian gives exactly the same numbers as building a dedicated
/// reduced-order model would.
fn get_a_submatrix(a: &DMatrix<f64>, indices: &[usize]) -> DMatrix<f64> {
    a.select_rows(indices).select_columns(indices)
}

/// Longitudinal states: `Vt, Alpha, Theta, Q` (indices 0, 1, 4, 7 in
/// [`FixedWing6DoFStates`](crate::model::dynamicmodel::fixedwing6dof::FixedWing6DoFStates)).
pub fn get_6dof_longitudal(a: &DMatrix<f64>) -> DMatrix<f64> {
    get_a_submatrix(a, &[0, 1, 4, 7])
}

/// Lateral-directional states: `Beta, Phi, P, R` (indices 2, 3, 6, 8 in
/// [`FixedWing6DoFStates`](crate::model::dynamicmodel::fixedwing6dof::FixedWing6DoFStates)).
pub fn get_6dof_lateral(a: &DMatrix<f64>) -> DMatrix<f64> {
    get_a_submatrix(a, &[2, 3, 6, 8])
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
