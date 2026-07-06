/// Returns a doublet input function that alternates between two values over a given duration.
///
/// # Arguments
///
/// * `initial_value` - The initial value of the doublet.
/// * `start_time` - The time at which the doublet starts.
/// * `pulse_duration` - The duration of each pulse.
/// * `pulse_value` - The value of each pulse.
/// * `time` - The current time.
pub fn doublet(
    initial_value: f64,
    start_time: f64,
    pulse_duration: f64,
    pulse_value: f64,
    time: f64,
) -> f64 {
    if time >= start_time && time < start_time + pulse_duration {
        initial_value + pulse_value
    } else if time >= start_time + pulse_duration && time < start_time + 2.0 * pulse_duration {
        initial_value - pulse_value
    } else {
        initial_value
    }
}
