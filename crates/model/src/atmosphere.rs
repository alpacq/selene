use crate::{GAMMA, GC, R0};

/// Returns the pressure and mach number for given true air speed and altitude
///
/// # Arguments
///
/// * `vt` - true air speed [m/s]
/// * `altitude` - altitude [m]
///
/// # Returns
///
/// * `(pressure, mach)` - pressure [Pa] and mach number [-]
pub fn air_data(vt: f64, altitude: f64) -> (f64, f64) {
    let temperature_drop_factor = 1.0 - 2.306e-5 * altitude; // temperature drop factor
    let temperature = if altitude >= 10668.0 {
        216.5
    } else {
        288.15 * temperature_drop_factor
    };

    let rho = R0 * temperature_drop_factor.powf(4.14);
    let mach = vt / ((GAMMA * GC * temperature).sqrt());
    let pressure = 0.5 * rho * vt * vt;
    (pressure, mach)
}
