use math::{input::Input, state::State, timestep::TimeStep};
use nalgebra::dvector;
use plotters::prelude::*;
use sim::run;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output = run(
        &van_der_pol,
        dvector![0.1, 0.1],
        dvector![0.8],
        60.0,
        &TimeStep::new(0.001),
    )?;

    let root = BitMapBackend::new("plot.png", (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Van der Pol Oscillator", ("sans-serif", 24).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-3f64..3f64, -4f64..5f64)?;

    chart.configure_mesh().x_desc("x1").y_desc("x2").draw()?;

    chart
        .draw_series(LineSeries::new(output.iter().map(|o| (o[1], o[2])), &BLACK))?
        .label("x2(x1)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

pub fn van_der_pol(x: &State, u: &Input) -> State {
    let u = u[0];
    let x1 = x[0];
    let x2 = x[1];

    dvector![x2, u * (1.0 - x1 * x1) * x2 - x1]
}
