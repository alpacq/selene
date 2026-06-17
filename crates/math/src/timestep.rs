//! TimeStep is a struct that holds the time step for the simulation.

pub struct TimeStep {
    dt: f64,
}

impl TimeStep {
    /// Creates a new TimeStep with the given time step in seconds.
    ///
    /// # Arguments
    ///
    /// * `dt` - The time step in seconds.
    ///
    /// # Panics
    ///
    /// Panics if `dt` is not positive.
    pub fn new(dt: f64) -> Self {
        if dt > 0.0 {
            TimeStep { dt }
        } else {
            panic!("dt must be positive.");
        }
    }

    /// Returns the time step in seconds.
    pub fn seconds(&self) -> f64 {
        self.dt
    }
}
