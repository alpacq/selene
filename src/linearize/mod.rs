//! This module provides linearization features for dynamic systems around selected equilibrium points.

use nalgebra::DMatrix;

use crate::error::LinearizationError;
use crate::math::SizedVector;
use crate::model::DynamicModel;

const TOLERANCE_MIN: f64 = 3.3e-5;
const TOLERANCE_OK: f64 = 8.1e-4;

/// Represents a linearization problem for a dynamic system around an equilibrium point.
pub struct LinearizationProblem<Sys, M>
where
    M: DynamicModel<Sys>,
{
    system: Sys,
    model: M,
    x_trimmed: M::State,
    u_trimmed: M::Input,
}

impl<Sys, M: DynamicModel<Sys>> LinearizationProblem<Sys, M> {
    pub fn fdx(&self, i: usize, j: usize, delta: f64) -> f64 {
        let mut x_minus = self.x_trimmed.vector().clone();
        let mut x_plus = self.x_trimmed.vector().clone();
        x_minus[j] -= delta;
        x_plus[j] += delta;
        let xd_minus = self.model.state_equations(
            &self.system,
            &M::State::from_vector(x_minus),
            &self.u_trimmed,
        );
        let xd_plus = self.model.state_equations(
            &self.system,
            &M::State::from_vector(x_plus),
            &self.u_trimmed,
        );
        (xd_plus.vector()[i] - xd_minus.vector()[i]) / (2.0 * delta)
    }

    pub fn fdu(&self, i: usize, j: usize, delta: f64) -> f64 {
        let mut u_minus = self.u_trimmed.vector().clone();
        let mut u_plus = self.u_trimmed.vector().clone();
        u_minus[j] -= delta;
        u_plus[j] += delta;
        let xd_minus = self.model.state_equations(
            &self.system,
            &self.x_trimmed,
            &M::Input::from_vector(u_minus),
        );
        let xd_plus = self.model.state_equations(
            &self.system,
            &self.x_trimmed,
            &M::Input::from_vector(u_plus),
        );
        (xd_plus.vector()[i] - xd_minus.vector()[i]) / (2.0 * delta)
    }

    fn adaptive_derivative(
        &self,
        partial_fn: impl Fn(usize, usize, f64) -> f64,
        i: usize,
        j: usize,
        v_j: f64,
    ) -> Result<f64, LinearizationError> {
        let initial_delta = (0.01 * v_j).abs();
        let mut delta = if initial_delta > 0.1 {
            initial_delta
        } else {
            0.1
        };
        let mut tolerance = 0.1;
        let (mut a0, mut a1, mut a2) = (0.0, 0.0, 0.0);
        let mut best = 0.0;
        let mut best_tolerance = f64::INFINITY;

        for k in 0..18 {
            a0 = partial_fn(i, j, delta);

            if k >= 2 {
                if a0 == a1 && a1 == a2 {
                    return Ok(a1);
                } else {
                    let d0 = (a0 - a1).abs(); // change between last 2
                    let d1 = (a1 - a2).abs(); // change between previous 2
                    let b0 = if a0 > a1 { a1 } else { a0 }; //scale to norm d0
                    let b1 = if a1 > a2 { a2 } else { a1 }; //scale to norm d1

                    if b0 > 0.0 && d0 <= tolerance * b0 && d1 <= tolerance * b1 {
                        best = a1;
                        best_tolerance = tolerance;
                        tolerance *= 0.2;
                        if tolerance <= TOLERANCE_MIN {
                            return Ok(a1);
                        }
                    }
                }
            }
            a2 = a1;
            a1 = a0;
            delta *= 0.6;
        }
        if best_tolerance <= TOLERANCE_OK {
            Ok(best)
        } else {
            return Err(LinearizationError::ConvergenceError(
                "Failed to converge".to_string(),
            ));
        }
    }

    pub fn jacobian_A(&self) -> DMatrix<f64> {
        let n = self.x_trimmed.vector().len();
        let mut result = DMatrix::<f64>::zeros(n, n);

        for j in 0..n {
            let v_j = self.x_trimmed.vector()[j];
            for i in 0..n {
                let partial = |i: usize, j: usize, delta: f64| -> f64 { self.fdx(i, j, delta) };
                result[(i, j)] = self.adaptive_derivative(partial, i, j, v_j).unwrap_or(0.0);
            }
        }

        result
    }

    pub fn jacobian_B(&self) -> DMatrix<f64> {
        let n = self.x_trimmed.vector().len();
        let m = self.u_trimmed.vector().len();
        let mut result = DMatrix::<f64>::zeros(n, m);

        for j in 0..m {
            let v_j = self.u_trimmed.vector()[j];
            for i in 0..n {
                let partial = |i: usize, j: usize, delta: f64| -> f64 { self.fdu(i, j, delta) };
                result[(i, j)] = self.adaptive_derivative(partial, i, j, v_j).unwrap_or(0.0);
            }
        }

        result
    }
}

/// First step in building the linearization problem process
/// Empty builder
pub struct LinearizationProblemBuilder;

/// Second step in building the linearization problem process
/// Requires a system to be specified
pub struct LinearizationProblemBuilderWithSystem<Sys> {
    system: Sys,
}

/// Third step in building the linearization problem process
/// Requires a model to be specified
pub struct LinearizationProblemBuilderWithModel<Sys, M: DynamicModel<Sys>> {
    system: Sys,
    model: M,
}

/// Fourth step in building the linearization problem process
/// Requires a trimmed state and input to be specified
pub struct LinearizationProblemBuilderWithTrimmedInputAndState<Sys, M: DynamicModel<Sys>> {
    system: Sys,
    model: M,
    x_trimmed: M::State,
    u_trimmed: M::Input,
}

impl LinearizationProblemBuilder {
    /// Returns a new empty builder for creating a linearization problem.
    pub fn new() -> Self {
        LinearizationProblemBuilder
    }

    /// Specifies the system to use for the linearization problem.
    pub fn for_system<Sys>(self, system: Sys) -> LinearizationProblemBuilderWithSystem<Sys> {
        LinearizationProblemBuilderWithSystem { system }
    }
}

impl<Sys> LinearizationProblemBuilderWithSystem<Sys> {
    /// Specifies the model to use for the linearization problem.
    pub fn with_model<M: DynamicModel<Sys>>(
        self,
        model: M,
    ) -> LinearizationProblemBuilderWithModel<Sys, M> {
        LinearizationProblemBuilderWithModel {
            system: self.system,
            model,
        }
    }
}

impl<Sys, M: DynamicModel<Sys>> LinearizationProblemBuilderWithModel<Sys, M> {
    /// Specifies the trimmed state and input to use for the linearization problem.
    pub fn with_trimmed_input_and_state(
        self,
        x_trimmed: M::State,
        u_trimmed: M::Input,
    ) -> LinearizationProblemBuilderWithTrimmedInputAndState<Sys, M> {
        LinearizationProblemBuilderWithTrimmedInputAndState {
            system: self.system,
            model: self.model,
            x_trimmed,
            u_trimmed,
        }
    }
}

impl<Sys, M: DynamicModel<Sys>> LinearizationProblemBuilderWithTrimmedInputAndState<Sys, M> {
    /// Returns the linearization problem built from this builder.
    pub fn build(self) -> LinearizationProblem<Sys, M> {
        LinearizationProblem {
            system: self.system,
            model: self.model,
            x_trimmed: self.x_trimmed,
            u_trimmed: self.u_trimmed,
        }
    }
}
