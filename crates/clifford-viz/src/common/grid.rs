//! Grid and axis rendering utilities for 2D plots.
//!
//! Provides functions to draw Cartesian grids and coordinate axes
//! with consistent styling across all visualization demos.

use egui_plot::{Line, PlotPoints};

use super::colors::palette;

/// Create a 2D Cartesian grid centered at the origin.
///
/// # Arguments
/// * `bounds` - Half-width/height of the grid (grid spans from -bounds to +bounds)
/// * `step` - Distance between grid lines
///
/// # Returns
/// A vector of [`Line`] objects representing the grid lines.
/// The axis lines (at x=0 and y=0) are drawn thicker.
#[must_use]
pub fn grid_2d(bounds: f64, step: f64) -> Vec<Line> {
    let mut lines = Vec::new();
    let n = (bounds / step).ceil() as i32;

    for i in -n..=n {
        let v = f64::from(i) * step;
        let width = if i == 0 { 1.5 } else { 0.5 };

        // Vertical line
        lines.push(
            Line::new(PlotPoints::new(vec![[v, -bounds], [v, bounds]]))
                .color(palette::GRID)
                .width(width as f32),
        );
        // Horizontal line
        lines.push(
            Line::new(PlotPoints::new(vec![[-bounds, v], [bounds, v]]))
                .color(palette::GRID)
                .width(width as f32),
        );
    }
    lines
}

/// Create coordinate axes with colored axis lines.
///
/// # Arguments
/// * `bounds` - Half-length of each axis (axes span from -bounds to +bounds)
///
/// # Returns
/// A vector of two [`Line`] objects: the X axis (red) and Y axis (green).
#[must_use]
pub fn axes_2d(bounds: f64) -> Vec<Line> {
    vec![
        Line::new(PlotPoints::new(vec![[-bounds, 0.0], [bounds, 0.0]]))
            .color(palette::X_AXIS)
            .width(2.0)
            .name("x"),
        Line::new(PlotPoints::new(vec![[0.0, -bounds], [0.0, bounds]]))
            .color(palette::Y_AXIS)
            .width(2.0)
            .name("y"),
    ]
}

/// Create a fine grid with major and minor lines.
///
/// # Arguments
/// * `bounds` - Half-width/height of the grid
/// * `major_step` - Distance between major (thick) grid lines
/// * `minor_divisions` - Number of minor lines between each major line
///
/// # Returns
/// Grid lines with major lines drawn thicker than minor lines.
#[must_use]
pub fn grid_2d_with_minor(bounds: f64, major_step: f64, minor_divisions: u32) -> Vec<Line> {
    let mut lines = Vec::new();
    let minor_step = major_step / f64::from(minor_divisions);
    let n = (bounds / minor_step).ceil() as i32;

    for i in -n..=n {
        let v = f64::from(i) * minor_step;
        let is_major = i % (minor_divisions as i32) == 0;
        let is_axis = i == 0;

        let (color, width) = if is_axis {
            (palette::GRID, 2.0)
        } else if is_major {
            (palette::GRID, 1.0)
        } else {
            (super::colors::with_alpha(palette::GRID, 80), 0.5)
        };

        // Vertical line
        lines.push(
            Line::new(PlotPoints::new(vec![[v, -bounds], [v, bounds]]))
                .color(color)
                .width(width),
        );
        // Horizontal line
        lines.push(
            Line::new(PlotPoints::new(vec![[-bounds, v], [bounds, v]]))
                .color(color)
                .width(width),
        );
    }
    lines
}

/// Create a polar grid (concentric circles and radial lines).
///
/// # Arguments
/// * `max_radius` - Maximum radius of the grid
/// * `radial_step` - Distance between concentric circles
/// * `angular_divisions` - Number of radial lines (evenly spaced)
/// * `segments_per_circle` - Number of line segments to approximate each circle
///
/// # Returns
/// Grid lines forming a polar coordinate system.
#[must_use]
pub fn polar_grid(
    max_radius: f64,
    radial_step: f64,
    angular_divisions: u32,
    segments_per_circle: usize,
) -> Vec<Line> {
    use std::f64::consts::PI;

    let mut lines = Vec::new();

    // Concentric circles
    let num_circles = (max_radius / radial_step).ceil() as u32;
    for i in 1..=num_circles {
        let radius = f64::from(i) * radial_step;
        let points: Vec<[f64; 2]> = (0..=segments_per_circle)
            .map(|j| {
                let angle = 2.0 * PI * j as f64 / segments_per_circle as f64;
                [radius * angle.cos(), radius * angle.sin()]
            })
            .collect();
        lines.push(
            Line::new(PlotPoints::new(points))
                .color(palette::GRID)
                .width(0.5),
        );
    }

    // Radial lines
    for i in 0..angular_divisions {
        let angle = 2.0 * PI * f64::from(i) / f64::from(angular_divisions);
        let dx = max_radius * angle.cos();
        let dy = max_radius * angle.sin();
        lines.push(
            Line::new(PlotPoints::new(vec![[0.0, 0.0], [dx, dy]]))
                .color(palette::GRID)
                .width(0.5),
        );
    }

    lines
}

/// Default bounds for most visualizations.
pub const DEFAULT_BOUNDS: f64 = 5.0;

/// Default grid step size.
pub const DEFAULT_STEP: f64 = 1.0;
