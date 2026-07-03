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
    ///
    /// # Arguments
    ///
    /// * `setpoints` - The setpoints for the trim problem.
    /// * `params` - The optimizer parameters.
    ///
    /// # Returns
    ///
    /// A tuple of the concrete state and input for a single cost evaluation.
    fn setup(
        &self,
        setpoints: &DVector<f64>,
        params: &DVector<f64>,
    ) -> Result<(Self::State, Self::Input), TrimError>;

    /// Returns the cost function value for the given state derivative.
    ///
    /// # Arguments
    ///
    /// * `x_dot` - The state derivative.
    ///
    /// # Returns
    ///
    /// The cost function value.
    fn cost(&self, x_dot: &Self::State) -> f64;
}

/// Builds the initial Nelder-Mead simplex from a single seed point using the
/// same heuristic as MATLAB's `fminsearch`: perturb each coordinate by 5%
/// (or a small fixed step for coordinates that are exactly zero).
///
/// # Arguments
///
/// * `x0` - The seed point for the simplex.
///
/// # Returns
///
/// A vector of `x0.len() + 1` vertices forming the initial simplex.
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
    /// Returns a builder for creating a trim problem from specific
    /// system and dynamic model.
    ///
    /// # Example
    ///
    /// ```
    /// use selene::trim::TrimProblem;
    /// use selene::model::dynamicmodel::State2;
    /// use selene::model::VanDerPol;
    ///
    /// let mut simulator = Simulator::builder()
    ///     .for_system(VanDerPol {})
    ///     .with_model(State2 {})
    ///     .with_setpoints(dvector![])
    ///     .with_initial_params(dvector![0.1, 0.1])
    ///     .build();
    /// ```
    pub fn builder() -> TrimProblemBuilder {
        TrimProblemBuilder::new()
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

/// CostFunction trait from argmin implementation for TrimProblem
/// It is necessary to implement this trait so that the `Executor` can evaluate the cost function
/// and perform the Nelder-Mead optimization.
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

/// First step in building the trim problem process
/// Empty builder
pub struct TrimProblemBuilder;

/// Second step in building the trim problem process
/// Requires a system to be specified
pub struct TrimProblemBuilderWithSystem<Sys> {
    system: Sys,
}

/// Third step in building the trim problem process
/// Requires a dynamic model to be specified
pub struct TrimProblemBuilderWithModel<Sys, M: TrimTarget<Sys>> {
    system: Sys,
    model: M,
}

/// Fourth step in building the trim problem process
/// Requires setpoints to be specified
pub struct TrimProblemBuilderWithSetpoints<Sys, M: TrimTarget<Sys>> {
    system: Sys,
    model: M,
    setpoints: DVector<f64>,
}

/// Fifth step in building the trim problem process
/// Requires initial parameters to be specified
pub struct TrimProblemBuilderWithInitialParams<Sys, M: TrimTarget<Sys>> {
    system: Sys,
    model: M,
    setpoints: DVector<f64>,
    initial_params: DVector<f64>,
}

impl TrimProblemBuilder {
    /// Returns a new empty builder for creating a trim problem.
    pub fn new() -> Self {
        TrimProblemBuilder
    }

    /// Specifies the system to use for the trim problem.
    pub fn for_system<Sys>(self, system: Sys) -> TrimProblemBuilderWithSystem<Sys> {
        TrimProblemBuilderWithSystem { system }
    }
}

impl<Sys> TrimProblemBuilderWithSystem<Sys> {
    /// Specifies the model to use for the trim problem.
    pub fn with_model<M: TrimTarget<Sys>>(self, model: M) -> TrimProblemBuilderWithModel<Sys, M> {
        TrimProblemBuilderWithModel {
            system: self.system,
            model,
        }
    }
}

impl<Sys, M: TrimTarget<Sys>> TrimProblemBuilderWithModel<Sys, M> {
    /// Specifies the setpoints to use for the trim problem.
    pub fn with_setpoints(
        self,
        setpoints: DVector<f64>,
    ) -> TrimProblemBuilderWithSetpoints<Sys, M> {
        TrimProblemBuilderWithSetpoints {
            system: self.system,
            model: self.model,
            setpoints,
        }
    }
}

impl<Sys, M: TrimTarget<Sys>> TrimProblemBuilderWithSetpoints<Sys, M> {
    /// Specifies the initial parameter guess to use for the trim problem.
    pub fn with_initial_params(
        self,
        initial_params: DVector<f64>,
    ) -> TrimProblemBuilderWithInitialParams<Sys, M> {
        TrimProblemBuilderWithInitialParams {
            system: self.system,
            model: self.model,
            setpoints: self.setpoints,
            initial_params,
        }
    }
}

impl<Sys, M: TrimTarget<Sys>> TrimProblemBuilderWithInitialParams<Sys, M> {
    /// Builds the trim problem with the specified initial parameter guess.
    pub fn build(self) -> TrimProblem<Sys, M> {
        TrimProblem {
            system: self.system,
            model: self.model,
            setpoints: self.setpoints,
            initial_params: self.initial_params,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::RAD_TO_DEG;
    use crate::model::Transport;
    use crate::model::dynamicmodel::fixedwing3dof::FixedWing3DoF;
    use nalgebra::dvector;

    #[test]
    fn test_transport_3dof_model_at_0_altitude_and_51_816_velocity()
    -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![51.816, 0.0, 0.0]; // setpoints: [vt, altitude, gamma]
        let init_params = dvector![0.1, -10.0, 0.1]; // initial params: [throttle, elevator, alpha]
        let problem = TrimProblemBuilder::new()
            .for_system(Transport::new())
            .with_model(FixedWing3DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build();
        let (x, u, cost) = problem.trim()?;
        assert!(cost < 1e-6, "trim did not converge: cost = {cost}");
        assert!((u.throttle() - 0.297).abs() < 5e-3);
        assert!((u.elevator() + 25.7).abs() < 5e-2);
        assert!((x.alpha() * RAD_TO_DEG - 22.1).abs() < 5e-2);
        Ok(())
    }

    #[test]
    fn test_transport_3dof_model_at_0_altitude_and_152_4_velocity()
    -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![152.4, 0.0, 0.0]; // setpoints: [vt, altitude, gamma]
        let init_params = dvector![0.1, -10.0, 0.1]; // initial params: [throttle, elevator, alpha]
        let problem = TrimProblemBuilder::new()
            .for_system(Transport::new())
            .with_model(FixedWing3DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build();
        let (x, u, cost) = problem.trim()?;
        assert!(cost < 1e-6, "trim did not converge: cost = {cost}");
        assert!((u.throttle() - 0.293).abs() < 5e-3);
        assert!((u.elevator() - 2.46).abs() < 5e-2);
        assert!((x.alpha() * RAD_TO_DEG - 0.58).abs() < 5e-2);
        Ok(())
    }

    #[test]
    fn test_transport_3dof_model_at_9144_altitude_and_152_4_velocity()
    -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![152.4, 9144.0, 0.0]; // setpoints: [vt, altitude, gamma]
        let init_params = dvector![0.1, -10.0, 0.1]; // initial params: [throttle, elevator, alpha]
        let problem = TrimProblemBuilder::new()
            .for_system(Transport::new())
            .with_model(FixedWing3DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build();
        let (x, u, cost) = problem.trim()?;
        assert!(cost < 1e-6, "trim did not converge: cost = {cost}");
        assert!((u.throttle() - 0.204).abs() < 5e-3);
        assert!((u.elevator() + 4.1).abs() < 5e-2);
        assert!((x.alpha() * RAD_TO_DEG - 5.43).abs() < 5e-2);
        Ok(())
    }
}
