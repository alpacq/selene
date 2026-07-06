pub mod error;
pub mod math;
pub mod model;
pub mod plots;
pub mod sim;
pub mod trim;

use nalgebra::{DVector, dvector};

use crate::{
    error::TrimError::{self, ConvergenceError},
    math::{SizedVector, input_functions::doublet, timestep::TimeStep},
    model::{
        F16, RAD_TO_DEG, Transport, VanDerPol,
        dynamicmodel::{
            fixedwing3dof::{FixedWing3DoF, FixedWing3DoFInput, FixedWing3DoFStates},
            fixedwing6dof::FixedWing6DoF,
            state2::{State2, State2Input, State2State},
        },
    },
    plots::{phase_portrait, state_variables_plot, yx},
    sim::simulator::SimulatorBuilder,
    trim::TrimProblemBuilder,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // A: just simulation of a dynamic system
    //
    // let mut simulator = SimulatorBuilder::new()
    //     .for_system(VanDerPol {})
    //     .with_model(State2 {})
    //     .with_state(State2State::new(dvector![0.1, 0.1]))
    //     .build();
    // let initial_state = simulator.state.clone();
    // simulator.run(
    //     initial_state,
    //     State2Input::new(dvector![0.8]),
    //     None,
    //     60.0,
    //     TimeStep::new(0.001),
    // );

    // phase_portrait(simulator.output, "Van der Pol Oscillator".into())

    // B: Trim of an aircraft, plotting throttle vs TAS at sea level
    // let set_vts = (42..=244).map(|i| i as f64).collect::<Vec<_>>();
    // let mut outs = Vec::new();
    // let mut init_params = dvector![
    //     1.0, // throttle
    //     0.0, // elevator
    //     0.1, // alpha
    //     0.0, // aileron
    //     0.0, // rudder
    //     0.0, // beta
    // ];
    // for set_vt in set_vts {
    //     let setpoints = dvector![
    //         set_vt, // vt [m/s]
    //         0.0,    // altitude [m]
    //         0.0,    // gamma [deg]
    //         0.0,    // roll rate [rad/s]
    //         0.0,    // pitch rate [rad/s]
    //         0.0,    // turn rate [rad/s]
    //         0.0,    // phi [rad] — calculated by solver for coordinated turn
    //         0.0,    // setpoint for coordinated turn
    //     ];
    //     let problem = TrimProblemBuilder::new()
    //         .for_system(F16::new())
    //         .with_model(FixedWing6DoF)
    //         .with_setpoints(setpoints)
    //         .with_initial_params(init_params.clone())
    //         .build();
    //     match problem.trim() {
    //         Ok((x, u, cost)) if cost < 2.5 => {
    //             init_params = dvector![
    //                 u.throttle(),
    //                 u.elevator(),
    //                 x.alpha(),
    //                 u.aileron(),
    //                 u.rudder(),
    //                 x.beta(),
    //             ];
    //             outs.push((set_vt, u.throttle()));
    //         }
    //         Ok((_, _, cost)) => {
    //             eprintln!("No convergence for vt={set_vt}, cost={cost:.3e}, skipping");
    //         }
    //         Err(e) => {
    //             eprintln!("Error for vt={set_vt}: {e}, skipping");
    //         }
    //     }
    // }

    // let (vts, throttles) = outs.into_iter().unzip();

    // yx(
    //     vts,
    //     throttles,
    //     "vt".into(),
    //     "throttle".into(),
    //     "Trimmed power curve".into(),
    // )

    let setpoints = dvector![76.2, 0.0, 0.0]; // setpoints: [vt, altitude, gamma]
    let init_params = dvector![0.1, -10.0, 0.1]; // initial params: [throttle, elevator, alpha]

    let system = Transport::new();
    let model = FixedWing3DoF;

    let problem = TrimProblemBuilder::new()
        .for_system(system)
        .with_model(model)
        .with_setpoints(setpoints)
        .with_initial_params(init_params)
        .build();
    let (x, u, _cost) = problem.trim()?;

    let system = Transport::new();
    let model = FixedWing3DoF;

    let mut simulator = SimulatorBuilder::new()
        .for_system(system)
        .with_model(model)
        .with_state(x)
        .build();

    let elevator_doublet = |params: &DVector<f64>, time: f64| -> FixedWing3DoFInput {
        FixedWing3DoFInput::from_vector(dvector![
            params[0],
            doublet(params[1], 1.0, 0.5, 2.0 / RAD_TO_DEG, time),
            params[2],
            params[3]
        ])
    };

    let throttle_doublet = |params: &DVector<f64>, time: f64| -> FixedWing3DoFInput {
        FixedWing3DoFInput::from_vector(dvector![
            doublet(params[0], 0.0, 3.0, 0.1, time),
            params[1],
            params[2],
            params[3]
        ])
    };

    // simulator.run(u, Some(elevator_doublet), 60.0, TimeStep::new(0.001));

    // state_variables_plot(
    //     vec![
    //         FixedWing3DoFStates::Alpha as usize,
    //         FixedWing3DoFStates::Theta as usize,
    //     ],
    //     simulator.output,
    //     vec!["alpha(t)".into(), "theta(t)".into()],
    //     "time (s)".into(),
    //     "angle (rad)".into(),
    //     "Elevator doublet response".into(),
    // )

    simulator.run(u, Some(throttle_doublet), 60.0, TimeStep::new(0.001));

    state_variables_plot(
        vec![
            FixedWing3DoFStates::Alpha as usize,
            FixedWing3DoFStates::Theta as usize,
        ],
        simulator.output,
        vec!["alpha(t)".into(), "theta(t)".into()],
        "time (s)".into(),
        "angle (rad)".into(),
        "Throttle doublet response".into(),
    )
}
