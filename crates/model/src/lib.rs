mod atmosphere;
pub mod fixedwingparams;
pub mod model;
pub mod transport;
pub mod van_der_pol;

pub use model::Model;
pub use transport::Transport;
pub use van_der_pol::VanDerPol;

/// Sea-level air density [kg/m^3]
const R0: f64 = 1225.0;

/// air adiabatic constant
const GAMMA: f64 = 1.4;

/// air gas constant // J / (kg * K)
const GC: f64 = 287.0;

/// Gravitational acceleration constant (m / s^2)
const GD: f64 = 9.80665;

/// Radians to degrees conversion factor
const RTOD: f64 = 57.29578;
