//! Output data from a simulation.

/// Time history recorded by the simulator.
///
/// The three vectors are always the same length and share the same sample
/// index: `time[k]`, `state_vector[k]` and `output_vector[k]` all describe the
/// system at the same instant. The first sample is the initial condition, so a
/// run of `n` steps produces `n + 1` samples.
#[derive(Default)]
pub struct SimOutput {
    /// Time vector [t0, t1, ..., tn]
    pub time: Vec<f64>,
    /// State vector [x0, x1, ..., xn] where x0 = [x00, x01, ..., x0n] etc
    pub state_vector: Vec<Vec<f64>>,
    /// Output vector [y0, y1, ..., yn] where y0 = [y00, y01, ..., y0m] etc
    pub output_vector: Vec<Vec<f64>>,
}

impl SimOutput {
    /// Returns the time history of the state variable at the given index.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the state variable to return.
    pub fn state_variable_at(&self, index: usize) -> Vec<f64> {
        self.state_vector.iter().map(|x| x[index]).collect()
    }

    /// Returns the time history of the output variable at the given index.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the output variable to return.
    pub fn output_variable_at(&self, index: usize) -> Vec<f64> {
        self.output_vector.iter().map(|y| y[index]).collect()
    }

    /// Returns the number of recorded samples.
    ///
    /// # Returns
    ///
    /// The number of recorded samples.
    pub fn len(&self) -> usize {
        self.time.len()
    }

    /// Returns whether the output is empty.
    ///
    /// # Returns
    ///
    /// `true` if the output is empty, `false` otherwise.
    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    /// Creates a new empty `SimOutput` with the given capacity.
    ///
    /// # Arguments
    ///
    /// * `capacity` - The capacity of the output vectors.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            time: Vec::with_capacity(capacity),
            state_vector: Vec::with_capacity(capacity),
            output_vector: Vec::with_capacity(capacity),
        }
    }
}
