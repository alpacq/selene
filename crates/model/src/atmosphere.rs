use crate::{GAMMA, GC, R0};

/// Returns the pressure for given true air speed and altitude
///
/// # Arguments
///
/// * `vt` - true air speed [m/s]
/// * `altitude` - altitude [m]
///
/// # Returns
///
/// * `pressure` - pressure [Pa] and mach number [-]
pub fn air_pressure(vt: f64, altitude: f64) -> f64 {
    let temperature_drop_factor = 1.0 - 2.306e-5 * altitude; // temperature drop factor

    let rho = R0 * temperature_drop_factor.powf(4.14);

    0.5 * rho * vt * vt
}

/// Returns the mach number for given true air speed and altitude
///
/// # Arguments
///
/// * `vt` - true air speed [m/s]
/// * `altitude` - altitude [m]
///
/// # Returns
///
/// * `mach` - mach number [-]
pub fn mach(vt: f64, altitude: f64) -> f64 {
    let temperature = temperature_at_altitude(altitude);

    vt / ((GAMMA * GC * temperature).sqrt())
}

/// Returns the true air speed for given mach number and altitude
///
/// # Arguments
///
/// * `mach` - mach number [-]
/// * `altitude` - altitude [m]
///
/// # Returns
///
/// * `tas` - true air speed [m/s]
pub fn tas_from_mach(mach: f64, altitude: f64) -> f64 {
    let temperature = temperature_at_altitude(altitude);

    mach * ((GAMMA * GC * temperature).sqrt())
}

fn temperature_at_altitude(altitude: f64) -> f64 {
    let temperature_drop_factor = 1.0 - 2.306e-5 * altitude; // temperature drop factor
    if altitude >= 10668.0 {
        216.5
    } else {
        288.15 * temperature_drop_factor
    }
}
