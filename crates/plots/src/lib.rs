//! This module provides plotting utilities for the simulation results.

use plotters::prelude::*;
use sim::output::SimOutput;

/// General helper function to plot a 2D line series.
///
/// # Arguments
///
/// * `x_axis_values` - The values for the x-axis.
/// * `y_axis_values` - The values for the y-axis.
/// * `x_axis_label` - The label for the x-axis.
/// * `y_axis_label` - The label for the y-axis.
/// * `title` - The title of the plot.
///
/// # Returns
///
/// Returns `Ok(())` if the plot is successfully created, or an error if it fails.
fn plot(
    x_axis_values: Vec<f64>,
    y_axis_values: Vec<f64>,
    x_axis_label: String,
    y_axis_label: String,
    title: String,
) -> Result<(), Box<dyn std::error::Error>> {
    let name = format!("{title}.png");
    let label = format!("{y_axis_label}({x_axis_label})");

    let x_min = x_axis_values
        .iter()
        .map(|x| *x)
        .fold(f64::INFINITY, f64::min)
        - 1.0;
    let x_max = x_axis_values
        .iter()
        .map(|x| *x)
        .fold(f64::NEG_INFINITY, f64::max)
        + 1.0;
    let y_min = y_axis_values
        .iter()
        .map(|y| *y)
        .fold(f64::INFINITY, f64::min)
        - 1.0;
    let y_max = y_axis_values
        .iter()
        .map(|y| *y)
        .fold(f64::NEG_INFINITY, f64::max)
        + 1.0;

    let root = BitMapBackend::new(name.as_str(), (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc(x_axis_label)
        .y_desc(y_axis_label)
        .draw()?;

    chart
        .draw_series(LineSeries::new(
            x_axis_values
                .iter()
                .zip(y_axis_values.iter())
                .map(|(x, y)| (*x, *y)),
            &BLACK,
        ))?
        .label(label)
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

/// Plots a phase portrait of the simulation results.
///
/// # Arguments
///
/// * `data` - The simulation results to plot.
/// * `title` - The title of the plot.
///
/// # Returns
///
/// Returns `Ok(())` if the plot is successfully created, or an error if it fails.
pub fn phase_portrait(data: SimOutput, title: String) -> Result<(), Box<dyn std::error::Error>> {
    let x_axis_values = data.output_variable(0);
    let y_axis_values = data.output_variable(1);
    let x_axis_label = "x1".into();
    let y_axis_label = "x2".into();

    plot(
        x_axis_values,
        y_axis_values,
        x_axis_label,
        y_axis_label,
        title,
    )
}
