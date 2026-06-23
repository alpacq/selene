use crate::engine::{Engine, engineparams::EngineParams};

const IDLE_THRUST_LUT: [[f64; 6]; 6] = [
    [4715.1, 2825.6, 3914.4, 5071.0, 6672.3, 8273.7],
    [2824.6, 1890.5, 3069.3, 4492.7, 5915.5, 7561.9],
    [266.9, 111.2, 1534.6, 3358.4, 5025.9, 6782.6],
    [-4537.6, -3158.0, -1334.5, 1556.9, 4047.5, 6049.6],
    [-12010.2, -8451.6, -5782.7, -1098.5, 2668.9, 4893.0],
    [-16013.6, -6227.5, -2646.7, -1521.3, -889.6, 3113.7],
];

const MIL_THRUST_LUT: [[f64; 6]; 6] = [
    [56403.4, 40701.2, 27579.0, 17570.5, 10898.1, 6227.5],
    [56403.4, 40701.2, 28082.0, 17970.9, 10987.0, 6227.5],
    [56092.8, 41407.5, 29401.3, 19083.1, 11565.4, 6938.2],
    [56226.1, 43770.5, 31537.8, 20728.8, 12632.9, 7384.0],
    [55118.4, 45262.8, 34473.7, 23665.3, 14456.7, 8585.1],
    [51955.4, 43803.9, 35808.2, 27134.2, 16903.2, 10275.5],
];

const MAX_THRUST_LUT: [[f64; 6]; 6] = [
    [88964.0, 66723.3, 48040.8, 31137.5, 17792.9, 11120.6],
    [95280.8, 69836.9, 49924.0, 32572.5, 19727.5, 11565.4],
    [100973.9, 74986.9, 54491.5, 36263.5, 22241.1, 12610.9],
    [107824.5, 84096.6, 61207.9, 41304.5, 25354.9, 14298.9],
    [115944.5, 93730.7, 71069.7, 49447.3, 30514.8, 17570.5],
    [128458.9, 103736.7, 81402.7, 59972.4, 38435.1, 22492.7],
];

pub struct F100PW220 {
    pub params: EngineParams,
}

impl F100PW220 {
    pub fn new() -> Self {
        Self {
            params: EngineParams {
                hx: 216.9,
                tstat: 0.0, //unused
                dtdv: 0.0,  //unused
            },
        }
    }
}

impl Engine for F100PW220 {
    /// For F-100-PW-220 there there is linear characteristic but with two different slope angles
    /// depending on wheter standard or military power is used
    /// (IDLE/MIL is up to 0.77 throttle value, and above it is MIL/MAX)
    /// MAX is with afterburner
    fn throttle_to_power(&self, throttle: f64) -> f64 {
        if throttle <= 0.77 {
            throttle * 64.94
        } else {
            throttle * 217.38 - 117.38
        }
    }

    /// For F-100-PW-220 it is first-order system
    /// Dynamics is different for different ranges of throttle
    /// and for the direction of power change
    fn power_dynamics(&self, power: f64, set_power: f64) -> f64 {
        let (t, actual_set_power) = if set_power >= 50.0 {
            if power >= 50.0 {
                (5.0, set_power)
            } else {
                (self.tau_inverse(set_power - power), 60.0)
            }
        } else {
            if power >= 50.0 {
                (5.0, 40.0)
            } else {
                (self.tau_inverse(set_power - power), set_power)
            }
        };
        t * (actual_set_power - power)
    }

    /// For F-100-PW-220 linear interpolation in power ranges
    /// 0-50.0 IDLE/MIL
    /// 50.0-100.0 MIL/MAX
    /// using thrust LUTs defined as constants for this engine
    fn thrust(&self, power: f64, altitude: f64, mach: f64) -> f64 {
        let altitude = altitude / 3048.0;
        let (altitude_index, altitude_remainder): (usize, f64) = if altitude >= 5.0 {
            (4, altitude - 4.0)
        } else {
            (altitude.floor() as usize, altitude - altitude.floor())
        };
        let altitude_complement = 1.0 - altitude_remainder;
        let mach = 5.0 * mach;
        let (mach_index, mach_remainder): (usize, f64) = if mach >= 5.0 {
            (4, mach - 4.0)
        } else {
            (mach.floor() as usize, mach - mach.floor())
        };

        let military_thrust_s = MIL_THRUST_LUT[altitude_index][mach_index] * altitude_complement
            + MIL_THRUST_LUT[altitude_index + 1][mach_index] * altitude_remainder;
        let military_thrust_t = MIL_THRUST_LUT[altitude_index][mach_index + 1]
            * altitude_complement
            + MIL_THRUST_LUT[altitude_index + 1][mach_index + 1] * altitude_remainder;
        let military_thrust =
            military_thrust_s + (military_thrust_t - military_thrust_s) * mach_remainder;

        if power <= 50.0 {
            let idle_thrust_s = IDLE_THRUST_LUT[altitude_index][mach_index] * altitude_complement
                + IDLE_THRUST_LUT[altitude_index + 1][mach_index] * altitude_remainder;
            let idle_thrust_t = IDLE_THRUST_LUT[altitude_index][mach_index + 1]
                * altitude_complement
                + IDLE_THRUST_LUT[altitude_index + 1][mach_index + 1] * altitude_remainder;
            let idle_thrust = idle_thrust_s + (idle_thrust_t - idle_thrust_s) * mach_remainder;
            idle_thrust * (military_thrust - idle_thrust) * power * 0.02
        } else {
            let max_thrust_s = MAX_THRUST_LUT[altitude_index][mach_index] * altitude_complement
                + MAX_THRUST_LUT[altitude_index + 1][mach_index] * altitude_remainder;
            let max_thrust_t = MAX_THRUST_LUT[altitude_index][mach_index + 1] * altitude_complement
                + MAX_THRUST_LUT[altitude_index + 1][mach_index + 1] * altitude_remainder;
            let max_thrust = max_thrust_s + (max_thrust_t - max_thrust_s) * mach_remainder;
            military_thrust * (max_thrust - military_thrust) * (power - 50.0) * 0.02
        }
    }

    /// For F-100-PW-220 it reacts faster for little differences
    /// between set and current power
    /// and slower for larger throttle changes
    fn tau_inverse(&self, delta_power: f64) -> f64 {
        if delta_power <= 25.0 {
            1.0
        } else if delta_power >= 50.0 {
            0.1
        } else {
            1.9 - 0.036 * delta_power
        }
    }
}
