//! This module provides linearization features for dynamic systems around selected equilibrium points.

use nalgebra::DMatrix;

use crate::error::LinearizationError;
use crate::math::SizedVector;
use crate::model::DynamicModel;
use crate::model::dynamicmodel::linearizeddynamicmodel::LinearizedDynamicModel;

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
    /// Computes the partial derivatives of the state vs state (A matrix element)
    ///
    /// # Arguments
    ///
    /// * `i` - the row index of the A matrix element.
    /// * `j` - the column index of the A matrix element.
    /// * `delta` - the finite difference step size.
    ///
    /// # Returns
    ///
    /// The [i,j] element of A matrix.
    fn fdx(&self, i: usize, j: usize, delta: f64) -> f64 {
        let mut x_minus = self.x_trimmed.vector().clone();
        let mut x_plus = self.x_trimmed.vector().clone();
        x_minus[j] -= delta;
        x_plus[j] += delta;
        let xd_minus = self
            .model
            .state_equations(
                &self.system,
                &M::State::from_vector(x_minus),
                &self.u_trimmed,
            )
            .0;
        let xd_plus = self
            .model
            .state_equations(
                &self.system,
                &M::State::from_vector(x_plus),
                &self.u_trimmed,
            )
            .0;
        (xd_plus.vector()[i] - xd_minus.vector()[i]) / (2.0 * delta)
    }

    /// Computes the partial derivatives of the state vs input (B matrix element)
    ///
    /// # Arguments
    ///
    /// * `i` - the row index of the B matrix element.
    /// * `j` - the column index of the B matrix element.
    /// * `delta` - the finite difference step size.
    ///
    /// # Returns
    ///
    /// The [i,j] element of B matrix.
    fn fdu(&self, i: usize, j: usize, delta: f64) -> f64 {
        let mut u_minus = self.u_trimmed.vector().clone();
        let mut u_plus = self.u_trimmed.vector().clone();
        u_minus[j] -= delta;
        u_plus[j] += delta;
        let xd_minus = self
            .model
            .state_equations(
                &self.system,
                &self.x_trimmed,
                &M::Input::from_vector(u_minus),
            )
            .0;
        let xd_plus = self
            .model
            .state_equations(
                &self.system,
                &self.x_trimmed,
                &M::Input::from_vector(u_plus),
            )
            .0;
        (xd_plus.vector()[i] - xd_minus.vector()[i]) / (2.0 * delta)
    }

    /// Computes the partial derivatives of the output vs state (C matrix element)
    ///
    /// # Arguments
    ///
    /// * `i` - the row index of the C matrix element.
    /// * `j` - the column index of the C matrix element.
    /// * `delta` - the finite difference step size
    ///
    /// # Returns
    ///
    /// The [i,j] element of the C matrix.
    fn ydx(&self, i: usize, j: usize, delta: f64) -> f64 {
        let mut x_minus = self.x_trimmed.vector().clone();
        let mut x_plus = self.x_trimmed.vector().clone();
        x_minus[j] -= delta;
        x_plus[j] += delta;
        let yd_minus = self
            .model
            .state_equations(
                &self.system,
                &M::State::from_vector(x_minus),
                &self.u_trimmed,
            )
            .1;
        let yd_plus = self
            .model
            .state_equations(
                &self.system,
                &M::State::from_vector(x_plus),
                &self.u_trimmed,
            )
            .1;
        (yd_plus.vector()[i] - yd_minus.vector()[i]) / (2.0 * delta)
    }

    /// Computes the partial derivatives of the output vs input (D matrix element)
    ///
    /// # Arguments
    ///
    /// * `i` - the row index of the D matrix element.
    /// * `j` - the column index of the D matrix element.
    /// * `delta` - the finite difference step size.
    ///
    /// # Returns
    ///
    /// The [i,j] element of the D matrix.
    fn ydu(&self, i: usize, j: usize, delta: f64) -> f64 {
        let mut u_minus = self.u_trimmed.vector().clone();
        let mut u_plus = self.u_trimmed.vector().clone();
        u_minus[j] -= delta;
        u_plus[j] += delta;
        let yd_minus = self
            .model
            .state_equations(
                &self.system,
                &self.x_trimmed,
                &M::Input::from_vector(u_minus),
            )
            .1;
        let yd_plus = self
            .model
            .state_equations(
                &self.system,
                &self.x_trimmed,
                &M::Input::from_vector(u_plus),
            )
            .1;
        (yd_plus.vector()[i] - yd_minus.vector()[i]) / (2.0 * delta)
    }

    /// Adaptive algorithm that doesn't use constant delta step,
    /// but instead adaptively decreases perturbation until reaches convergence.
    ///
    /// # Arguments
    ///
    /// * `partial_fn` - The function to compute the partial derivative.
    /// * `i` - The row index of the matrix element.
    /// * `j` - The column index of the matrix element.
    /// * `v_j` - The value of the j-th element of th trimmed state or input vector.
    ///
    /// # Returns
    ///
    /// The computed derivative value.
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
        let mut a0;
        let (mut a1, mut a2) = (0.0, 0.0);
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
                    // Scale to the magnitude of the smaller of each pair, regardless of
                    // sign, so negative derivatives converge just as well as positive
                    // ones. Floor the scale at TOLERANCE_MIN so near-zero derivatives
                    // fall back to an absolute (rather than relative) tolerance test.
                    let b0 = a0.abs().min(a1.abs()).max(TOLERANCE_MIN); //scale to norm d0
                    let b1 = a1.abs().min(a2.abs()).max(TOLERANCE_MIN); //scale to norm d1

                    if d0 <= tolerance * b0 && d1 <= tolerance * b1 {
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

    /// Computes the Jacobian matrix A using adaptive finite difference.
    ///
    /// # Returns
    ///
    /// The Jacobian matrix A.
    pub fn jacobian_a(&self) -> DMatrix<f64> {
        let n = self.model.system_rank();
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

    /// Computes the Jacobian matrix B using adaptive finite difference.
    ///
    /// # Returns
    ///
    /// The Jacobian matrix B.
    pub fn jacobian_b(&self) -> DMatrix<f64> {
        let n = self.model.system_rank();
        let m = self.model.num_inputs();
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

    /// Computes the Jacobian matrix C using adaptive finite difference.
    ///
    /// # Returns
    ///
    /// The Jacobian matrix C.
    pub fn jacobian_c(&self) -> DMatrix<f64> {
        let p = self.model.num_outputs();
        let n = self.model.system_rank();
        let mut result = DMatrix::<f64>::zeros(p, n);

        for j in 0..n {
            let v_j = self.x_trimmed.vector()[j];
            for i in 0..p {
                let partial = |i: usize, j: usize, delta: f64| -> f64 { self.ydx(i, j, delta) };
                result[(i, j)] = self.adaptive_derivative(partial, i, j, v_j).unwrap_or(0.0);
            }
        }

        result
    }

    /// Computes the Jacobian matrix D using adaptive finite difference.
    ///
    /// # Returns
    ///
    /// The Jacobian matrix D.
    pub fn jacobian_d(&self) -> DMatrix<f64> {
        let p = self.model.num_outputs();
        let m = self.model.num_inputs();
        let mut result = DMatrix::<f64>::zeros(p, m);

        for j in 0..m {
            let v_j = self.u_trimmed.vector()[j];
            for i in 0..p {
                let partial = |i: usize, j: usize, delta: f64| -> f64 { self.ydu(i, j, delta) };
                result[(i, j)] = self.adaptive_derivative(partial, i, j, v_j).unwrap_or(0.0);
            }
        }

        result
    }

    pub fn to_linearized_model(
        &self,
    ) -> Result<LinearizedDynamicModel<Sys, M>, LinearizationError> {
        let a = self.jacobian_a();
        let b = self.jacobian_b();
        let c = self.jacobian_c();
        let d = self.jacobian_d();

        let n = self.model.system_rank();
        if a.nrows() != n {
            return Err(LinearizationError::DimensionMismatch(format!(
                "A rows {} ≠ model.system_rank() {}",
                a.nrows(),
                n
            )));
        }
        let m = self.model.num_inputs();
        if b.ncols() != m {
            return Err(LinearizationError::DimensionMismatch(format!(
                "B cols {} ≠ model.num_inputs() {}",
                b.ncols(),
                m
            )));
        }
        let p = self.model.num_outputs();
        if c.nrows() != p {
            return Err(LinearizationError::DimensionMismatch(format!(
                "C rows {} ≠ model.num_outputs() {}",
                c.nrows(),
                p
            )));
        }

        LinearizedDynamicModel::new(a, b, c, d)
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

/// Abstraction over the linearization process.
///
/// # Inputs:
/// - `system`: The system to linearize.
/// - `model`: The model to linearize.
/// - `x_trimmed`: The trimmed state vector.
/// - `u_trimmed`: The trimmed input vector.
///
/// # Returns:
/// - `Ok(LinearizedDynamicModel)`: The linearized model.
/// - `Err(LinearizationError)`: An error occurred during linearization.
pub fn linearize<Sys, M: DynamicModel<Sys>>(
    system: Sys,
    model: M,
    x_trimmed: M::State,
    u_trimmed: M::Input,
) -> Result<LinearizedDynamicModel<Sys, M>, LinearizationError> {
    let problem = LinearizationProblemBuilder::new()
        .for_system(system)
        .with_model(model)
        .with_trimmed_input_and_state(x_trimmed, u_trimmed)
        .build();

    problem.to_linearized_model()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::test_utils::assert_approx;
    use crate::model::VanDerPol;
    use crate::model::dynamicmodel::state2::{State2, State2Input, State2State};
    use nalgebra::{dmatrix, dvector};

    fn assert_matrix_approx(actual: &DMatrix<f64>, expected: &DMatrix<f64>, label: &str) {
        assert_eq!(
            actual.shape(),
            expected.shape(),
            "{label}: shape {:?} != {:?}",
            actual.shape(),
            expected.shape()
        );
        for i in 0..expected.nrows() {
            for j in 0..expected.ncols() {
                assert_approx(
                    actual[(i, j)],
                    expected[(i, j)],
                    1e-9,
                    format!("{label}[{i},{j}]").as_str(),
                );
            }
        }
    }

    /// The Van der Pol oscillator ẋ = [x2, u(1 - x1²)x2 - x1], y = x1 has the
    /// analytic Jacobian A = [[0, 1], [-1, u]] at the origin, where B vanishes
    /// because the damping term is multiplied by x2 = 0.
    #[test]
    fn van_der_pol_at_the_origin_matches_the_analytic_jacobian()
    -> Result<(), Box<dyn std::error::Error>> {
        let problem = LinearizationProblemBuilder::new()
            .for_system(VanDerPol {})
            .with_model(State2)
            .with_trimmed_input_and_state(
                State2State::new(dvector![0.0, 0.0]),
                State2Input::new(dvector![0.8]),
            )
            .build();

        let linear = problem.to_linearized_model()?;

        assert_matrix_approx(linear.a(), &dmatrix![0.0, 1.0; -1.0, 0.8], "A");
        assert_matrix_approx(linear.b(), &dmatrix![0.0; 0.0], "B");
        assert_matrix_approx(linear.c(), &dmatrix![1.0, 0.0], "C");
        assert_matrix_approx(linear.d(), &dmatrix![0.0], "D");
        Ok(())
    }

    /// Linearizing a model that is already linear must return the original
    /// matrices, regardless of the point it is linearized about. This exercises
    /// the finite-difference Jacobians against an exactly known answer.
    #[test]
    fn linearizing_a_linear_model_recovers_its_matrices() -> Result<(), Box<dyn std::error::Error>>
    {
        let a = dmatrix![0.0, 1.0; -2.0, -3.0];
        let b = dmatrix![0.0; 1.5];
        let c = dmatrix![1.0, 0.25];
        let d = dmatrix![0.5];

        let linear = LinearizedDynamicModel::<VanDerPol, State2>::new(
            a.clone(),
            b.clone(),
            c.clone(),
            d.clone(),
        )?;

        // Deliberately off-equilibrium: a linear model has the same Jacobians
        // everywhere, so the trim point must not influence the result.
        let problem = LinearizationProblemBuilder::new()
            .for_system(VanDerPol {})
            .with_model(linear)
            .with_trimmed_input_and_state(
                State2State::new(dvector![0.3, -0.2]),
                State2Input::new(dvector![0.7]),
            )
            .build();

        let relinearized = problem.to_linearized_model()?;

        assert_matrix_approx(relinearized.a(), &a, "A");
        assert_matrix_approx(relinearized.b(), &b, "B");
        assert_matrix_approx(relinearized.c(), &c, "C");
        assert_matrix_approx(relinearized.d(), &d, "D");
        Ok(())
    }
}
