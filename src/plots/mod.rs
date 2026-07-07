//! This module provides plotting utilities for the simulation results.

use crate::sim::output::SimOutput;
use plotters::prelude::*;

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

    let x_min = x_axis_values.iter().copied().fold(f64::INFINITY, f64::min) - 1.0;
    let x_max = x_axis_values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
        + 1.0;
    let y_min = y_axis_values.iter().copied().fold(f64::INFINITY, f64::min) - 1.0;
    let y_max = y_axis_values
        .iter()
        .copied()
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
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

/// General helper function to plot multiple 2D line series.
///
/// # Arguments
///
/// * `x_axis_values` - The values for the x-axis.
/// * `y_axis_values` - Vector of vectors of values for the y-axis.
/// * `line_labels` - Vector of labels for the x-axis.
/// * `x_axis_label` - The label for the x-axis.
/// * `y_axis_label` - The label for the y-axis.
/// * `title` - The title of the plot.
///
/// # Returns
///
/// Returns `Ok(())` if the plot is successfully created, or an error if it fails.
fn multiple_plot(
    x_axis_values: Vec<f64>,
    y_axis_values: Vec<Vec<f64>>,
    line_labels: Vec<String>,
    x_axis_label: String,
    y_axis_label: String,
    title: String,
) -> Result<(), Box<dyn std::error::Error>> {
    let name = format!("{title}.png");

    let x_data_min = x_axis_values.iter().copied().fold(f64::INFINITY, f64::min);
    let x_data_max = x_axis_values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let y_data_min = y_axis_values
        .iter()
        .flatten()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let y_data_max = y_axis_values
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    // Pad the axes proportionally to the data range so small variations
    // (e.g. a fraction of a degree) stay visible instead of being dwarfed
    // by a fixed +-1.0 margin.
    let x_span = x_data_max - x_data_min;
    let y_span = y_data_max - y_data_min;
    let x_pad = if x_span > 0.0 { x_span * 0.05 } else { 1.0 };
    let y_pad = if y_span > 0.0 { y_span * 0.1 } else { 1.0 };

    let x_min = x_data_min - x_pad;
    let x_max = x_data_max + x_pad;
    let y_min = y_data_min - y_pad;
    let y_max = y_data_max + y_pad;

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
        .x_desc(&x_axis_label)
        .y_desc(&y_axis_label)
        .draw()?;

    let colors = vec![BLACK, RED, BLUE, GREEN];
    let mut i = 0;

    for (y_vector, label) in y_axis_values.iter().zip(line_labels.iter()) {
        let current_color = colors[i];
        chart
            .draw_series(LineSeries::new(
                x_axis_values
                    .iter()
                    .zip(y_vector.iter())
                    .map(|(x, y)| (*x, *y)),
                &current_color,
            ))?
            .label(label)
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], current_color));
        if i == 3 {
            i = 0;
        } else {
            i += 1;
        }
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
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
    let x_axis_values = data.output_variable_at(0);
    let y_axis_values = data.output_variable_at(1);
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

/// Plots a standard y(x) plot of the given x and y.
///
/// # Arguments
///
/// * `x` - The x values vector.
/// * `y` - The y values vector.
/// * `x_label` - The label for the x-axis.
/// * `y_label` - The label for the y-axis.
/// * `title` - The title of the plot.
///
/// # Returns
///
/// Returns `Ok(())` if the plot is successfully created, or an error if it fails.
pub fn yx(
    x: Vec<f64>,
    y: Vec<f64>,
    x_label: String,
    y_label: String,
    title: String,
) -> Result<(), Box<dyn std::error::Error>> {
    plot(x, y, x_label, y_label, title)
}

/// Plots multiple y(x) plots for the given state variables
/// based on simulation output results.
///
/// # Arguments
///
/// * `data` - The simulation results to plot.
/// * `line_labels` - The labels for the lines on the plot.
/// * `x_axis_label` - The label for the x-axis.
/// * `y_axis_label` - The label for the y-axis.
/// * `title` - The title of the plot.
///
/// # Returns
///
/// Returns `Ok(())` if the plot is successfully created, or an error if it fails.
pub fn state_variables_plot(
    variables: Vec<usize>,
    data: SimOutput,
    line_labels: Vec<String>,
    x_axis_label: String,
    y_axis_label: String,
    title: String,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut ys = Vec::new();
    for variable in variables {
        ys.push(data.output_variable_at(variable).clone())
    }
    let x = data.time;

    multiple_plot(x, ys, line_labels, x_axis_label, y_axis_label, title)
}

/// Plots a state variable against another state variable based on simulation output results.
///
/// # Arguments
///
/// * `data` - The simulation results to plot.
/// * `x_variable` - The index of the x-axis variable.
/// * `y_variable` - The index of the y-axis variable.
/// * `x_axis_label` - The label for the x-axis.
/// * `y_axis_label` - The label for the y-axis.
/// * `title` - The title of the plot.
///
/// # Returns
///
/// Returns `Ok(())` if the plot is successfully created, or an error if it fails.
pub fn state_variable_of_state_variable_plot(
    data: SimOutput,
    x_variable: usize,
    y_variable: usize,
    x_axis_label: String,
    y_axis_label: String,
    title: String,
) -> Result<(), Box<dyn std::error::Error>> {
    let x_axis_values = data.output_variable_at(x_variable);
    let y_axis_values = data.output_variable_at(y_variable);

    plot(
        x_axis_values,
        y_axis_values,
        x_axis_label,
        y_axis_label,
        title,
    )
}
