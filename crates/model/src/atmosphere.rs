use crate::{GAMMA, GC, R0};

/// Returns the dynamic pressure for given true air speed and altitude
///
/// # Arguments
///
/// * `vt` - true air speed [m/s]
/// * `altitude` - altitude [m]
///
/// # Returns
///
/// * `pressure` - pressure [Pa] and mach number [-]
pub fn dynamic_pressure(vt: f64, altitude: f64) -> f64 {
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

/// Returns the temperature at a given altitude
///
/// # Arguments
///
/// * `altitude` - altitude [m]
///
/// # Returns
///
/// * `temperature` - temperature [K]
fn temperature_at_altitude(altitude: f64) -> f64 {
    let temperature_drop_factor = 1.0 - 2.306e-5 * altitude; // temperature drop factor
    if altitude >= 10668.0 {
        216.5
    } else {
        288.15 * temperature_drop_factor
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temperature_at_mount_everest() {
        assert_eq!(format!("{:.2}", temperature_at_altitude(8848.0)), "229.36");
    }

    #[test]
    fn temperature_at_airliner_cruising_altitude() {
        assert_eq!(temperature_at_altitude(11000.0), 216.5);
    }

    #[test]
    fn temperature_at_sea_level() {
        assert_eq!(temperature_at_altitude(0.0), 288.15);
    }

    #[test]
    fn speed_of_sound() {
        assert_eq!(format!("{:.2}", tas_from_mach(1.0, 0.0)), "340.26");
    }

    #[test]
    fn mach_number() {
        assert_eq!(mach(340.2626485525556, 0.0), 1.0);
    }

    #[test]
    fn dynamic_pressure_at_mount_everest() {
        assert_eq!(
            format!("{:.2}", dynamic_pressure(300.0, 8848.0)),
            "21431.37"
        );
    }

    #[test]
    fn dynamic_pressure_at_sea_level() {
        assert_eq!(
            format!("{:.2}", dynamic_pressure(340.2626485525556, 0.0)),
            "70914.44"
        );
    }
}
