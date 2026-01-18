//! Grid and axis rendering utilities for 2D plots.
//!
//! Provides functions to draw Cartesian grids and coordinate axes
//! with consistent styling across all visualization demos.
//!
//! # Visual Hierarchy
//!
//! Grids use a three-level hierarchy to avoid competing with visualization content:
//! 1. **Minor lines** - Very subtle, `GRID_MINOR` color, `HAIRLINE` weight
//! 2. **Major lines** - Slightly visible, `GRID_MAJOR` color, `THIN` weight
//! 3. **Axes** - Clear but desaturated, semantic colors, `NORMAL` weight
//!
//! # Theme Support
//!
//! All grid functions take an `&egui::Context` parameter and automatically
//! adapt colors to the current light/dark theme.

use egui::Context;
use egui_plot::{Line, PlotPoints};

use super::colors::{grid, grid_major, grid_minor, line_weights, x_axis, y_axis};

/// Create a 2D Cartesian grid centered at the origin.
///
/// # Arguments
/// * `ctx` - The egui context (for theme-aware colors)
/// * `bounds` - Half-width/height of the grid (grid spans from -bounds to +bounds)
/// * `step` - Distance between grid lines
///
/// # Returns
/// A vector of [`Line`] objects representing the grid lines.
/// Uses subtle colors that don't compete with visualization content.
#[must_use]
pub fn grid_2d(ctx: &Context, bounds: f64, step: f64) -> Vec<Line> {
    let grid_color = grid(ctx);
    let mut lines = Vec::new();
    let n = (bounds / step).ceil() as i32;

    for i in -n..=n {
        let v = f64::from(i) * step;
        // Skip axis lines (drawn separately with colors)
        if i == 0 {
            continue;
        }

        // Vertical line
        lines.push(
            Line::new(PlotPoints::new(vec![[v, -bounds], [v, bounds]]))
                .color(grid_color)
                .width(line_weights::HAIRLINE),
        );
        // Horizontal line
        lines.push(
            Line::new(PlotPoints::new(vec![[-bounds, v], [bounds, v]]))
                .color(grid_color)
                .width(line_weights::HAIRLINE),
        );
    }
    lines
}

/// Create coordinate axes with colored axis lines.
///
/// # Arguments
/// * `ctx` - The egui context (for theme-aware colors)
/// * `bounds` - Half-length of each axis (axes span from -bounds to +bounds)
///
/// # Returns
/// A vector of two [`Line`] objects: the X axis (muted red) and Y axis (muted green).
/// Uses desaturated colors for a professional, non-fatiguing appearance.
#[must_use]
pub fn axes_2d(ctx: &Context, bounds: f64) -> Vec<Line> {
    vec![
        Line::new(PlotPoints::new(vec![[-bounds, 0.0], [bounds, 0.0]]))
            .color(x_axis(ctx))
            .width(line_weights::NORMAL)
            .name("x"),
        Line::new(PlotPoints::new(vec![[0.0, -bounds], [0.0, bounds]]))
            .color(y_axis(ctx))
            .width(line_weights::NORMAL)
            .name("y"),
    ]
}

/// Create a fine grid with major and minor lines.
///
/// This creates a clear visual hierarchy:
/// - Minor lines are very subtle (every `major_step / minor_divisions`)
/// - Major lines are slightly more visible (every `major_step`)
/// - Axis lines are skipped (use `axes_2d` separately)
///
/// # Arguments
/// * `ctx` - The egui context (for theme-aware colors)
/// * `bounds` - Half-width/height of the grid
/// * `major_step` - Distance between major (more visible) grid lines
/// * `minor_divisions` - Number of minor lines between each major line
///
/// # Returns
/// Grid lines with major lines drawn more prominently than minor lines.
#[must_use]
pub fn grid_2d_with_minor(
    ctx: &Context,
    bounds: f64,
    major_step: f64,
    minor_divisions: u32,
) -> Vec<Line> {
    let major_color = grid_major(ctx);
    let minor_color = grid_minor(ctx);
    let mut lines = Vec::new();
    let minor_step = major_step / f64::from(minor_divisions);
    let n = (bounds / minor_step).ceil() as i32;

    for i in -n..=n {
        let v = f64::from(i) * minor_step;
        let is_major = i % (minor_divisions as i32) == 0;
        let is_axis = i == 0;

        // Skip axis lines (drawn separately)
        if is_axis {
            continue;
        }

        let (color, width) = if is_major {
            (major_color, line_weights::THIN)
        } else {
            (minor_color, line_weights::HAIRLINE)
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
/// * `ctx` - The egui context (for theme-aware colors)
/// * `max_radius` - Maximum radius of the grid
/// * `radial_step` - Distance between concentric circles
/// * `angular_divisions` - Number of radial lines (evenly spaced)
/// * `segments_per_circle` - Number of line segments to approximate each circle
///
/// # Returns
/// Grid lines forming a polar coordinate system.
#[must_use]
pub fn polar_grid(
    ctx: &Context,
    max_radius: f64,
    radial_step: f64,
    angular_divisions: u32,
    segments_per_circle: usize,
) -> Vec<Line> {
    use std::f64::consts::PI;

    let grid_color = grid(ctx);
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
                .color(grid_color)
                .width(line_weights::HAIRLINE),
        );
    }

    // Radial lines
    for i in 0..angular_divisions {
        let angle = 2.0 * PI * f64::from(i) / f64::from(angular_divisions);
        let dx = max_radius * angle.cos();
        let dy = max_radius * angle.sin();
        lines.push(
            Line::new(PlotPoints::new(vec![[0.0, 0.0], [dx, dy]]))
                .color(grid_color)
                .width(line_weights::HAIRLINE),
        );
    }

    lines
}

/// Create an origin marker (small crosshair at the origin).
///
/// This provides a subtle visual anchor at the coordinate origin.
///
/// # Arguments
/// * `ctx` - The egui context (for theme-aware colors)
/// * `size` - Half-length of the crosshair arms
///
/// # Returns
/// Two short line segments forming a crosshair.
#[must_use]
pub fn origin_marker(ctx: &Context, size: f64) -> Vec<Line> {
    let marker_color = grid_major(ctx);
    vec![
        Line::new(PlotPoints::new(vec![[-size, 0.0], [size, 0.0]]))
            .color(marker_color)
            .width(line_weights::THIN),
        Line::new(PlotPoints::new(vec![[0.0, -size], [0.0, size]]))
            .color(marker_color)
            .width(line_weights::THIN),
    ]
}

/// Default bounds for most visualizations.
pub const DEFAULT_BOUNDS: f64 = 5.0;

/// Default grid step size.
pub const DEFAULT_STEP: f64 = 1.0;
