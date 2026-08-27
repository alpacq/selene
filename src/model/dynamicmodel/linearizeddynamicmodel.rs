use crate::{error::LinearizationError, math::SizedVector, model::DynamicModel};
use nalgebra::DMatrix;
use std::marker::PhantomData;

pub struct LinearizedDynamicModel<Sys, M: DynamicModel<Sys>> {
    a: DMatrix<f64>,
    b: DMatrix<f64>,
    c: DMatrix<f64>,
    d: DMatrix<f64>,
    _phantom: PhantomData<(Sys, M)>,
}

impl<Sys, M: DynamicModel<Sys>> LinearizedDynamicModel<Sys, M> {
    pub fn new(
        a: DMatrix<f64>,
        b: DMatrix<f64>,
        c: DMatrix<f64>,
        d: DMatrix<f64>,
    ) -> Result<Self, LinearizationError> {
        if !a.is_square() {
            return Err(LinearizationError::DimensionMismatch(
                "A must be square".to_string(),
            ));
        }

        let n = a.nrows();
        let m = b.ncols();
        let p = c.nrows();

        if b.nrows() != n {
            return Err(LinearizationError::DimensionMismatch(format!(
                "B rows {} ≠ A rows {}",
                b.nrows(),
                n
            )));
        }
        if c.ncols() != n {
            return Err(LinearizationError::DimensionMismatch(format!(
                "C cols {} ≠ A rows {}",
                c.ncols(),
                n
            )));
        }
        if d.nrows() != p || d.ncols() != m {
            return Err(LinearizationError::DimensionMismatch(format!(
                "D shape {}×{} ≠ C rows {} × B cols {}",
                d.nrows(),
                d.ncols(),
                p,
                m
            )));
        }

        Ok(Self {
            a,
            b,
            c,
            d,
            _phantom: PhantomData,
        })
    }

    pub fn a(&self) -> &DMatrix<f64> {
        &self.a
    }

    pub fn b(&self) -> &DMatrix<f64> {
        &self.b
    }

    pub fn c(&self) -> &DMatrix<f64> {
        &self.c
    }

    pub fn d(&self) -> &DMatrix<f64> {
        &self.d
    }
}

impl<Sys, M: DynamicModel<Sys>> DynamicModel<Sys> for LinearizedDynamicModel<Sys, M> {
    type State = M::State;
    type Input = M::Input;
    type Output = M::Output;

    fn state_equations(
        &self,
        _system: &Sys,
        x: &Self::State,
        u: &Self::Input,
    ) -> (Self::State, Self::Output) {
        let dx_vec = &self.a * x.vector() + &self.b * u.vector();
        let dx = Self::State::from_vector(dx_vec);

        let y_vec = &self.c * x.vector() + &self.d * u.vector();
        let y = Self::Output::from_vector(y_vec);

        (dx, y)
    }

    fn system_rank(&self) -> usize {
        self.a.nrows()
    }

    fn num_inputs(&self) -> usize {
        self.b.ncols()
    }

    fn num_outputs(&self) -> usize {
        self.c.nrows()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::test_utils::assert_approx;
    use crate::math::timestep::TimeStep;
    use crate::model::VanDerPol;
    use crate::model::dynamicmodel::state2::{State2, State2Input, State2State};
    use crate::sim::simulator::SimulatorBuilder;
    use nalgebra::{dmatrix, dvector};

    /// Any linear model borrows the state / input / output types of the model
    /// it is a linearization of, so a 2-state, 1-input, 1-output alias is
    /// enough to exercise the whole trait implementation.
    type LinearState2 = LinearizedDynamicModel<VanDerPol, State2>;

    /// ẋ = [[0, 1], [-2, -3]]·x + [[0], [1]]·u,  y = [1, 0]·x + 0.5·u
    fn model() -> LinearState2 {
        LinearState2::new(
            dmatrix![0.0, 1.0; -2.0, -3.0],
            dmatrix![0.0; 1.0],
            dmatrix![1.0, 0.0],
            dmatrix![0.5],
        )
        .expect("matrix dimensions are consistent")
    }

    #[test]
    fn state_derivative_is_a_x_plus_b_u() {
        let x = State2State::new(dvector![1.0, 2.0]);
        let u = State2Input::new(dvector![3.0]);

        let (x_dot, _) = model().state_equations(&VanDerPol {}, &x, &u);

        // A·x = [2, -8], B·u = [0, 3]
        assert_approx(x_dot.x1(), 2.0, 1e-12, "x1_dot");
        assert_approx(x_dot.x2(), -5.0, 1e-12, "x2_dot");
    }

    #[test]
    fn output_is_c_x_plus_d_u() {
        let x = State2State::new(dvector![1.0, 2.0]);
        let u = State2Input::new(dvector![3.0]);

        let (_, y) = model().state_equations(&VanDerPol {}, &x, &u);

        // C·x = 1.0, D·u = 1.5
        assert_approx(y.y(), 2.5, 1e-12, "y");
    }

    #[test]
    fn system_is_linear_in_state_and_input() {
        let model = model();
        let system = VanDerPol {};
        let x = State2State::new(dvector![0.4, -1.1]);
        let u = State2Input::new(dvector![0.7]);

        let (x_dot, y) = model.state_equations(&system, &x, &u);
        let (x_dot_scaled, y_scaled) = model.state_equations(
            &system,
            &State2State::new(x.vector() * 2.0),
            &State2Input::new(u.vector() * 2.0),
        );

        assert_approx(x_dot_scaled.x1(), 2.0 * x_dot.x1(), 1e-12, "x1_dot scaling");
        assert_approx(x_dot_scaled.x2(), 2.0 * x_dot.x2(), 1e-12, "x2_dot scaling");
        assert_approx(y_scaled.y(), 2.0 * y.y(), 1e-12, "y scaling");
    }

    #[test]
    fn dimensions_are_derived_from_the_matrices() {
        let model = model();

        assert_eq!(model.system_rank(), 2);
        assert_eq!(model.num_inputs(), 1);
        assert_eq!(model.num_outputs(), 1);
    }

    #[test]
    fn non_square_state_matrix_is_rejected() {
        let result = LinearState2::new(
            dmatrix![0.0, 1.0, 0.0; -2.0, -3.0, 0.0],
            dmatrix![0.0; 1.0],
            dmatrix![1.0, 0.0],
            dmatrix![0.0],
        );

        assert!(matches!(
            result,
            Err(LinearizationError::DimensionMismatch(_))
        ));
    }

    #[test]
    fn input_matrix_with_wrong_row_count_is_rejected() {
        let result = LinearState2::new(
            dmatrix![0.0, 1.0; -2.0, -3.0],
            dmatrix![0.0; 1.0; 0.0],
            dmatrix![1.0, 0.0],
            dmatrix![0.0],
        );

        assert!(matches!(
            result,
            Err(LinearizationError::DimensionMismatch(_))
        ));
    }

    #[test]
    fn output_matrix_with_wrong_column_count_is_rejected() {
        let result = LinearState2::new(
            dmatrix![0.0, 1.0; -2.0, -3.0],
            dmatrix![0.0; 1.0],
            dmatrix![1.0, 0.0, 0.0],
            dmatrix![0.0],
        );

        assert!(matches!(
            result,
            Err(LinearizationError::DimensionMismatch(_))
        ));
    }

    #[test]
    fn feedthrough_matrix_with_wrong_shape_is_rejected() {
        let result = LinearState2::new(
            dmatrix![0.0, 1.0; -2.0, -3.0],
            dmatrix![0.0; 1.0],
            dmatrix![1.0, 0.0],
            dmatrix![0.0, 0.0],
        );

        assert!(matches!(
            result,
            Err(LinearizationError::DimensionMismatch(_))
        ));
    }

    /// A linear model must be a drop-in replacement for the model it was
    /// derived from, so it has to run in the existing simulator. The harmonic
    /// oscillator ẋ = [[0, 1], [-1, 0]]·x has the closed-form solution
    /// x(t) = [cos t, -sin t] for x(0) = [1, 0].
    #[test]
    fn simulated_harmonic_oscillator_matches_analytic_solution() {
        let oscillator = LinearState2::new(
            dmatrix![0.0, 1.0; -1.0, 0.0],
            dmatrix![0.0; 0.0],
            dmatrix![1.0, 0.0],
            dmatrix![0.0],
        )
        .expect("matrix dimensions are consistent");

        let mut simulator = SimulatorBuilder::new()
            .for_system(VanDerPol {})
            .with_model(oscillator)
            .with_state(State2State::new(dvector![1.0, 0.0]))
            .build();

        simulator.run(
            State2Input::new(dvector![0.0]),
            None,
            1.0,
            TimeStep::new(0.001),
        );

        assert_approx(simulator.state.x1(), 1.0_f64.cos(), 1e-9, "x1(1 s)");
        assert_approx(simulator.state.x2(), -1.0_f64.sin(), 1e-9, "x2(1 s)");
    }

    /// With C = [1, 0] and D = 0 the recorded output history must be the first
    /// state variable, which is what separates it from the state history.
    #[test]
    fn recorded_output_history_follows_the_c_matrix() {
        let oscillator = LinearState2::new(
            dmatrix![0.0, 1.0; -1.0, 0.0],
            dmatrix![0.0; 0.0],
            dmatrix![1.0, 0.0],
            dmatrix![0.0],
        )
        .expect("matrix dimensions are consistent");

        let mut simulator = SimulatorBuilder::new()
            .for_system(VanDerPol {})
            .with_model(oscillator)
            .with_state(State2State::new(dvector![1.0, 0.0]))
            .build();

        simulator.run(
            State2Input::new(dvector![0.0]),
            None,
            0.5,
            TimeStep::new(0.001),
        );

        let states = simulator.output.state_variable_at(0);
        let outputs = simulator.output.output_variable_at(0);

        assert_eq!(states.len(), outputs.len());
        for (k, (x1, y)) in states.iter().zip(outputs.iter()).enumerate() {
            assert_approx(*y, *x1, 1e-12, format!("y[{k}] vs x1[{k}]").as_str());
        }
        assert_approx(*outputs.last().unwrap(), 0.5_f64.cos(), 1e-9, "y(0.5 s)");
    }
}
