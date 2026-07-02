use argmin::{
    core::{CostFunction, Error, Executor, State},
    solver::neldermead::NelderMead,
};
use nalgebra::DVector;

use crate::{error::TrimError, model::DynamicModel};

/// Model-specific glue required to trim a [`DynamicModel`].
///
/// `Sys` is the system the model runs against (e.g. an aircraft). The trait
/// reuses the model's [`DynamicModel::State`] / [`DynamicModel::Input`] types,
/// so the generic [`TrimProblem`] never needs to know their concrete shape.
pub trait TrimTarget<Sys>: DynamicModel<Sys> {
    /// Maps the optimizer `params` and the `setpoints` into a concrete state
    /// and input for a single cost evaluation.
    fn setup(
        &self,
        setpoints: &DVector<f64>,
        params: &DVector<f64>,
    ) -> Result<(Self::State, Self::Input), TrimError>;

    /// Scalarizes the state derivative into the value to minimize.
    fn cost(&self, x_dot: &Self::State) -> f64;

    /// Initial Nelder-Mead simplex (`n + 1` vertices) in parameter space.
    fn initial_simplex(&self) -> Vec<DVector<f64>>;
}

/// A generic trim problem binding a system, a model and its setpoints.
///
/// It minimizes the model's [`TrimTarget::objective`] with Nelder-Mead
/// (argmin's equivalent of MATLAB's `fminsearch`).
pub struct TrimProblem<Sys, M>
where
    M: TrimTarget<Sys>,
{
    system: Sys,
    model: M,
    setpoints: DVector<f64>,
}

impl<Sys, M> TrimProblem<Sys, M>
where
    M: TrimTarget<Sys>,
{
    pub fn new(system: Sys, model: M, setpoints: DVector<f64>) -> Self {
        Self {
            system,
            model,
            setpoints,
        }
    }

    /// Runs the trim and returns the trimmed state, input and final cost.
    pub fn trim(self) -> Result<(M::State, M::Input, f64), Error> {
        let solver = NelderMead::new(self.model.initial_simplex()).with_sd_tolerance(1e-8)?;

        let res = Executor::new(self, solver)
            .configure(|state| state.max_iters(1000))
            .run()?;

        let best = res
            .state()
            .get_best_param()
            .ok_or_else(|| Error::msg("trim: no best parameter found"))?
            .clone();
        let cost = res.state().get_best_cost();

        // Rebuild the trimmed state and input from the minimizing vector,
        // reusing the exact same mapping as the cost evaluation.
        let problem = res
            .problem
            .problem
            .expect("problem is returned by executor");
        let (x, u) = problem.model.setup(&problem.setpoints, &best)?;

        Ok((x, u, cost))
    }
}

impl<Sys, M> CostFunction for TrimProblem<Sys, M>
where
    M: TrimTarget<Sys>,
{
    type Param = DVector<f64>;
    type Output = f64;

    fn cost(&self, params: &Self::Param) -> Result<Self::Output, Error> {
        let (x, u) = self.model.setup(&self.setpoints, params)?;
        let x_dot = self.model.state_equations(&self.system, &x, &u);
        Ok(self.model.cost(&x_dot))
    }
}
