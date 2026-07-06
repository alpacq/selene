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
