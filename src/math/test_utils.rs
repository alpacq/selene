#[cfg(test)]
pub fn assert_approx(actual: f64, expected: f64, epsilon: f64, label: &str) {
    assert!(
        (actual - expected).abs() < epsilon,
        "{label}: expected {expected}, got {actual} (diff={diff:.2e}, epsilon={epsilon:.2e})",
        diff = (actual - expected).abs()
    );
}
