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

    /// Returns the cost function value for the given state derivative.
    fn cost(&self, x_dot: &Self::State) -> f64;
}

/// Builds the initial Nelder-Mead simplex from a single seed point using the
/// same heuristic as MATLAB's `fminsearch`: perturb each coordinate by 5%
/// (or a small fixed step for coordinates that are exactly zero).
fn build_simplex(x0: &DVector<f64>) -> Vec<DVector<f64>> {
    const USUAL_DELTA: f64 = 0.05;
    const ZERO_DELTA: f64 = 0.00025;

    let mut simplex = Vec::with_capacity(x0.len() + 1);
    simplex.push(x0.clone());
    for i in 0..x0.len() {
        let mut vertex = x0.clone();
        vertex[i] = if vertex[i] != 0.0 {
            (1.0 + USUAL_DELTA) * vertex[i]
        } else {
            ZERO_DELTA
        };
        simplex.push(vertex);
    }
    simplex
}

/// A generic trim problem binding a system, a model and its setpoints.
///
/// It minimizes the model's [`TrimTarget::cost`] with Nelder-Mead
/// (argmin's equivalent of MATLAB's `fminsearch`).
pub struct TrimProblem<Sys, M>
where
    M: TrimTarget<Sys>,
{
    system: Sys,
    model: M,
    setpoints: DVector<f64>,
    /// Initial guess for the optimizer parameters (the seed `s0`), in the
    /// same layout as [`TrimTarget::setup`]'s `params`.
    initial_params: DVector<f64>,
}

impl<Sys, M> TrimProblem<Sys, M>
where
    M: TrimTarget<Sys>,
{
    /// Creates a new trim problem for the given system, model, setpoints and
    /// initial parameter guess.
    pub fn new(
        system: Sys,
        model: M,
        setpoints: DVector<f64>,
        initial_params: DVector<f64>,
    ) -> Self {
        Self {
            system,
            model,
            setpoints,
            initial_params,
        }
    }

    /// Runs the trim and returns the trimmed state, input and final cost.
    pub fn trim(self) -> Result<(M::State, M::Input, f64), Error> {
        let simplex = build_simplex(&self.initial_params);
        let solver = NelderMead::new(simplex).with_sd_tolerance(1e-8)?;

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

    /// Evaluates the cost function of the trim problem for the given parameters.
    fn cost(&self, params: &Self::Param) -> Result<Self::Output, Error> {
        let (x, u) = self.model.setup(&self.setpoints, params)?;
        let x_dot = self.model.state_equations(&self.system, &x, &u);
        Ok(self.model.cost(&x_dot))
    }
}
