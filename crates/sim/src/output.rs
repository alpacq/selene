//! Output data from a simulation.

pub struct SimOutput {
    /// Time vector [t0, t1, ..., tn]
    pub time: Vec<f64>,
    /// Output vector [y0, y1, ..., yn]
    pub output_vector: Vec<Vec<f64>>,
}

impl SimOutput {
    /// Returns the output variable at the given index.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the output variable to return.
    pub fn output_variable(&self, index: usize) -> Vec<f64> {
        self.output_vector.iter().map(|y| y[index]).collect()
    }

    /// Returns the number of time steps.
    pub fn len(&self) -> usize {
        self.time.len()
    }

    /// Returns whether the output is empty.
    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    /// Creates a new empty `SimOutput` with the given capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            time: Vec::with_capacity(capacity),
            output_vector: Vec::with_capacity(capacity),
        }
    }
}
