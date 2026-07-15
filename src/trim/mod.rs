//! Trimming dynamics models in steady state

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
    /// * `system` - The system the model runs against.
    /// * `setpoints` - The setpoints for the trim problem.
    /// * `params` - The optimizer parameters.
    ///
    /// # Returns
    ///
    /// A tuple of the concrete state and input for a single cost evaluation.
    fn setup(
        &self,
        system: &Sys,
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
    /// Runs the trim and returns the trimmed state, input and final cost.
    ///
    /// A single Nelder-Mead run frequently terminates on a collapsed simplex
    /// that is close to, but short of, the true minimum. To recover, the search
    /// is restarted from the current best point (rebuilding a fresh simplex)
    /// until it converges or a restart no longer makes progress. This mirrors
    /// the restart loop used by the Stevens & Lewis trimmer.
    pub fn trim(self) -> Result<(M::State, M::Input, f64), Error> {
        /// Maximum number of simplex restarts before giving up.
        const MAX_RESTARTS: usize = 50;
        /// Cost below which the trim is considered converged.
        const COST_TOLERANCE: f64 = 1e-12;

        let mut problem = self;
        let mut seed = problem.initial_params.clone();
        let mut best_param = seed.clone();
        let mut best_cost = f64::INFINITY;

        for _ in 0..MAX_RESTARTS {
            let simplex = build_simplex(&seed);
            let solver = NelderMead::new(simplex).with_sd_tolerance(1e-10)?;

            let res = Executor::new(problem, solver)
                .configure(|state| state.max_iters(1000))
                .run()?;

            let cost = res.state().get_best_cost();
            let param = res
                .state()
                .get_best_param()
                .ok_or_else(|| Error::msg("trim: no best parameter found"))?
                .clone();

            // Reclaim the problem so the next restart (or the final state
            // reconstruction) can reuse it.
            problem = res
                .problem
                .problem
                .expect("problem is returned by executor");

            let improved = cost + f64::EPSILON < best_cost;
            if cost < best_cost {
                best_cost = cost;
                best_param = param;
            }

            // Stop once converged, or when a restart stops making progress.
            if best_cost < COST_TOLERANCE || !improved {
                break;
            }
            seed = best_param.clone();
        }

        // Rebuild the trimmed state and input from the minimizing vector,
        // reusing the exact same mapping as the cost evaluation.
        let (x, u) = problem
            .model
            .setup(&problem.system, &problem.setpoints, &best_param)?;

        Ok((x, u, best_cost))
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
        let (x, u) = self.model.setup(&self.system, &self.setpoints, params)?;
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
    /// Builds the trim problem with the provided properties..
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
    use crate::math::test_utils::assert_approx;
    use crate::model::F16;
    use crate::model::RAD_TO_DEG;
    use crate::model::Transport;
    use crate::model::dynamicmodel::fixedwing3dof::FixedWing3DoF;
    use crate::model::dynamicmodel::fixedwing6dof::FixedWing6DoF;
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
        assert_approx(u.throttle(), 0.297, 5e-3, "throttle");
        assert_approx(u.elevator(), -25.7, 5e-2, "elevator");
        assert_approx(x.alpha() * RAD_TO_DEG, 22.1, 5e-2, "alpha");
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
        assert_approx(u.throttle(), 0.293, 5e-3, "throttle");
        assert_approx(u.elevator(), 2.46, 5e-2, "elevator");
        assert_approx(x.alpha() * RAD_TO_DEG, 0.58, 5e-2, "alpha");
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
        assert_approx(u.throttle(), 0.204, 5e-3, "throttle");
        assert_approx(u.elevator(), -4.1, 5e-2, "elevator");
        assert_approx(x.alpha() * RAD_TO_DEG, 5.43, 5e-2, "alpha");
        Ok(())
    }

    #[test]
    fn test_f16_6dof_model_at_0_altitude_and_152_1_velocity_coordinated_turn_0_15()
    -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![
            152.1, // vt [m/s]
            0.0,   // altitude [m]
            0.0,   // gamma [deg]
            0.0,   // roll rate [rad/s]
            0.0,   // pitch rate [rad/s]
            0.15,  // turn rate [rad/s]
            0.0,   // phi [rad] — calculated by solver for coordinated turn
            1.0,   // setpoint for coordinated turn
        ];
        let init_params = dvector![
            0.2,  // throttle
            1.0,  // elevator
            0.02, // alpha
            1.0,  // aileron
            1.0,  // rudder
            0.02, // beta
        ];
        let problem = TrimProblemBuilder::new()
            .for_system(F16::new())
            .with_model(FixedWing6DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build();
        let (x, u, cost) = problem.trim()?;
        assert!(cost < 1e-6, "trim did not converge: cost = {cost}");
        // Coordinated turn => essentially zero sideslip.
        assert!(x.beta().abs() * RAD_TO_DEG < 1.0);
        // Bank angle is fixed by the turn-coordination kinematics (~67 deg).
        assert_approx(x.phi() * RAD_TO_DEG, 66.9, 1.0, "phi");
        // Steady-state angle of attack and throttle for the turn.
        assert_approx(x.alpha() * RAD_TO_DEG, 7.0, 0.5, "alpha");
        assert_approx(u.throttle(), 0.334, 0.03, "throttle");
        Ok(())
    }

    #[test]
    fn test_f16_6dof_model_at_0_altitude_and_152_1_velocity_coordinated_turn_0_3()
    -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![
            152.1, // vt [m/s]
            0.0,   // altitude [m]
            0.0,   // gamma [deg]
            0.0,   // roll rate [rad/s]
            0.0,   // pitch rate [rad/s]
            0.3,   // turn rate [rad/s]
            0.0,   // phi [rad] — calculated by solver for coordinated turn
            1.0,   // setpoint for coordinated turn
        ];
        let init_params = dvector![
            0.2,  // throttle
            1.0,  // elevator
            0.02, // alpha
            1.0,  // aileron
            1.0,  // rudder
            0.02, // beta
        ];
        let problem = TrimProblemBuilder::new()
            .for_system(F16::new())
            .with_model(FixedWing6DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build();
        let (x, u, cost) = problem.trim()?;
        assert!(cost < 1e-6, "trim did not converge: cost = {cost}");
        // Coordinated turn => essentially zero sideslip.
        assert!(x.beta().abs() * RAD_TO_DEG < 1.0);
        // Bank angle is fixed by the turn-coordination kinematics (~78 deg).
        assert_approx(x.phi() * RAD_TO_DEG, 78.2, 1.0, "phi");
        // Steady-state angle of attack and throttle for the turn.
        assert_approx(x.alpha() * RAD_TO_DEG, 13.79, 0.5, "alpha");
        assert_approx(u.throttle(), 0.835, 0.03, "throttle");
        Ok(())
    }

    #[test]
    fn trim_6dof_f16_nominal_xcg_035() -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![
            152.9, // vt [m/s]
            0.0,   // altitude [m]
            0.0,   // gamma [deg]
            0.0,   // roll rate [rad/s]
            0.0,   // pitch rate [rad/s]
            0.0,   // turn rate [rad/s]
            0.0,   // phi [rad]
            0.0,   // coordinated turn flag (0 = none)
        ];
        // Xcg = 0.35 is the default in F16::new(), so no extra param needed.
        let init_params = dvector![
            0.15,  // throttle — close to expected 0.1385
            -1.0,  // elevator — close to expected -0.7588 deg
            0.037, // alpha — close to expected 0.03691 rad
            0.0,   // aileron
            0.0,   // rudder
            0.0,   // beta
        ];

        let (x, u, cost) = TrimProblemBuilder::new()
            .for_system(F16::new())
            .with_model(FixedWing6DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build()
            .trim()?;

        assert!(cost < 1e-8, "trim did not converge: cost = {cost:.3e}");

        // State
        assert_approx(x.alpha(), 0.03691, 1e-3, "alpha [rad]");
        assert_approx(x.beta(), 0.0, 1e-4, "beta [rad]");
        assert_approx(x.phi(), 0.0, 1e-6, "phi [rad]");
        assert_approx(x.theta(), 0.03691, 1e-3, "theta [rad]");
        assert_approx(x.p(), 0.0, 1e-6, "P [rad/s]");
        assert_approx(x.q(), 0.0, 1e-6, "Q [rad/s]");
        assert_approx(x.r(), 0.0, 1e-6, "R [rad/s]");

        // Controls
        assert_approx(u.throttle(), 0.1385, 1e-3, "throttle");
        assert_approx(u.elevator(), -0.7588, 5e-3, "elevator [deg]");
        // AIL and RDR are ~1e-7 deg at trim (numerical noise) — just check near zero
        assert!(
            u.aileron().abs() < 1e-3,
            "aileron should be ~0, got {}",
            u.aileron()
        );
        assert!(
            u.rudder().abs() < 2e-3,
            "rudder should be ~0, got {}",
            u.rudder()
        );

        Ok(())
    }

    #[test]
    fn trim_6dof_f16_xcg_030() -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![
            152.9, // vt [m/s]
            0.0,   // altitude [m]
            0.0,   // gamma [deg]
            0.0,   // roll rate [rad/s]
            0.0,   // pitch rate [rad/s]
            0.0,   // turn rate [rad/s]
            0.0,   // phi [rad]
            0.0,   // coordinated turn flag
            0.30,  // Xcg — adjust index to match your setpoints layout
        ];
        let init_params = dvector![
            0.15,  // throttle
            -2.0,  // elevator
            0.039, // alpha
            0.0,   // aileron
            0.0,   // rudder
            0.0,   // beta
        ];

        let (x, u, cost) = TrimProblemBuilder::new()
            .for_system(F16::new())
            .with_model(FixedWing6DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build()
            .trim()?;

        assert!(cost < 1e-8, "trim did not converge: cost = {cost:.3e}");

        assert_approx(x.alpha(), 0.03936, 1e-3, "alpha [rad]");
        assert_approx(x.beta(), 0.0, 1e-6, "beta [rad]");
        assert_approx(x.phi(), 0.0, 1e-6, "phi [rad]");
        assert_approx(x.theta(), 0.03936, 1e-3, "theta [rad]");
        assert_approx(x.p(), 0.0, 1e-6, "P [rad/s]");
        assert_approx(x.q(), 0.0, 1e-6, "Q [rad/s]");
        assert_approx(x.r(), 0.0, 1e-6, "R [rad/s]");

        assert_approx(u.throttle(), 0.1485, 1e-3, "throttle");
        assert_approx(u.elevator(), -1.931, 5e-3, "elevator [deg]");
        assert!(
            u.aileron().abs() < 1e-4,
            "aileron should be ~0, got {}",
            u.aileron()
        );
        assert!(
            u.rudder().abs() < 1e-3,
            "rudder should be ~0, got {}",
            u.rudder()
        );

        Ok(())
    }

    #[test]
    fn trim_6dof_f16_xcg_038() -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![
            152.9, // vt [m/s]
            0.0,   // altitude [m]
            0.0,   // gamma [deg]
            0.0,   // roll rate [rad/s]
            0.0,   // pitch rate [rad/s]
            0.0,   // turn rate [rad/s]
            0.0,   // phi [rad]
            0.0,   // coordinated turn flag
            0.38,  // Xcg — adjust index to match your setpoints layout
        ];
        let init_params = dvector![
            0.13,  // throttle
            -0.1,  // elevator
            0.035, // alpha
            0.0,   // aileron
            0.0,   // rudder
            0.0,   // beta
        ];

        let (x, u, cost) = TrimProblemBuilder::new()
            .for_system(F16::new())
            .with_model(FixedWing6DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build()
            .trim()?;

        assert!(cost < 1e-8, "trim did not converge: cost = {cost:.3e}");

        assert_approx(x.alpha(), 0.03544, 1e-4, "alpha [rad]");
        assert_approx(x.beta(), 0.0, 1e-5, "beta [rad]");
        assert_approx(x.phi(), 0.0, 1e-6, "phi [rad]");
        assert_approx(x.theta(), 0.03544, 1e-4, "theta [rad]");
        assert_approx(x.p(), 0.0, 1e-6, "P [rad/s]");
        assert_approx(x.q(), 0.0, 1e-6, "Q [rad/s]");
        assert_approx(x.r(), 0.0, 1e-6, "R [rad/s]");

        assert_approx(u.throttle(), 0.1325, 1e-3, "throttle");
        assert_approx(u.elevator(), -0.05590, 5e-3, "elevator [deg]");
        assert!(
            u.aileron().abs() < 1e-4,
            "aileron should be ~0, got {}",
            u.aileron()
        );
        assert!(
            u.rudder().abs() < 1e-3,
            "rudder should be ~0, got {}",
            u.rudder()
        );

        Ok(())
    }

    #[test]
    fn trim_6dof_f16_xcg_030_coordinated_turn_03_rads() -> Result<(), Box<dyn std::error::Error>> {
        let setpoints = dvector![
            152.9, // vt [m/s]
            0.0,   // altitude [m]
            0.0,   // gamma [deg]
            0.0,   // roll rate — set by CONSTR
            0.0,   // pitch rate — set by CONSTR
            0.3,   // turn rate (psi_dot) [rad/s]
            0.0,   // phi — determined by coordinated turn kinematics
            1.0,   // coordinated turn flag
            0.30,  // Xcg
        ];
        let init_params = dvector![
            0.85, // throttle — close to expected 0.8499
            -6.0, // elevator — close to expected -6.256 deg
            0.25, // alpha — close to expected 0.2485 rad
            0.1,  // aileron — close to expected 0.09891 deg
            -0.4, // rudder — close to expected -0.4218 deg
            0.0,  // beta — coordinated => ~0
        ];

        let (x, u, cost) = TrimProblemBuilder::new()
            .for_system(F16::new())
            .with_model(FixedWing6DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params)
            .build()
            .trim()?;

        assert!(cost < 1e-6, "trim did not converge: cost = {cost:.3e}");

        // State
        assert_approx(x.alpha(), 0.2485, 2e-3, "alpha [rad]");
        assert_approx(x.beta(), 0.0, 1e-3, "beta [rad]"); // coordinated => ~0
        assert_approx(x.phi(), 1.367, 1e-2, "phi [rad]");
        assert_approx(x.theta(), 0.05185, 1e-3, "theta [rad]");
        assert_approx(x.p(), -0.01555, 1e-3, "P [rad/s]");
        assert_approx(x.q(), 0.2934, 2e-3, "Q [rad/s]");
        assert_approx(x.r(), 0.06071, 1e-3, "R [rad/s]");

        // Controls
        assert_approx(u.throttle(), 0.8499, 5e-3, "throttle");
        assert_approx(u.elevator(), -6.256, 5e-2, "elevator [deg]");
        assert_approx(u.aileron(), 0.09891, 5e-3, "aileron [deg]");
        assert_approx(u.rudder(), -0.4218, 5e-3, "rudder [deg]");

        Ok(())
    }
}
