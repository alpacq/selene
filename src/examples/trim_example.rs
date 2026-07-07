use nalgebra::dvector;

use crate::{
    model::{F16, dynamicmodel::fixedwing6dof::FixedWing6DoF},
    plots::yx,
    trim::TrimProblemBuilder,
};

/// Example of a trim problem for the F-16 model.
pub fn trim_example() -> Result<(), Box<dyn std::error::Error>> {
    let set_vts = (42..=244).map(|i| i as f64).collect::<Vec<_>>();
    let mut outs = Vec::new();
    let mut init_params = dvector![
        1.0, // throttle
        0.0, // elevator
        0.1, // alpha
        0.0, // aileron
        0.0, // rudder
        0.0, // beta
    ];
    for set_vt in set_vts {
        let setpoints = dvector![
            set_vt, // vt [m/s]
            0.0,    // altitude [m]
            0.0,    // gamma [deg]
            0.0,    // roll rate [rad/s]
            0.0,    // pitch rate [rad/s]
            0.0,    // turn rate [rad/s]
            0.0,    // phi [rad] — calculated by solver for coordinated turn
            0.0,    // setpoint for coordinated turn
        ];
        let problem = TrimProblemBuilder::new()
            .for_system(F16::new())
            .with_model(FixedWing6DoF)
            .with_setpoints(setpoints)
            .with_initial_params(init_params.clone())
            .build();
        match problem.trim() {
            Ok((x, u, cost)) if cost < 2.5 => {
                init_params = dvector![
                    u.throttle(),
                    u.elevator(),
                    x.alpha(),
                    u.aileron(),
                    u.rudder(),
                    x.beta(),
                ];
                outs.push((set_vt, u.throttle()));
            }
            Ok((_, _, cost)) => {
                eprintln!("No convergence for vt={set_vt}, cost={cost:.3e}, skipping");
            }
            Err(e) => {
                eprintln!("Error for vt={set_vt}: {e}, skipping");
            }
        }
    }

    let (vts, throttles) = outs.into_iter().unzip();

    yx(
        vts,
        throttles,
        "vt".into(),
        "throttle".into(),
        "Trimmed power curve".into(),
    )
}
