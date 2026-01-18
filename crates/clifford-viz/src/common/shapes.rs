//! 2D shape primitives for geometric visualizations.
//!
//! Provides functions to draw common geometric shapes using `egui_plot`:
//! - Points (markers)
//! - Arrows (vectors)
//! - Lines (finite segments and infinite lines)
//! - Circles and arcs
//!
//! All shapes return [`egui_plot`] primitives that can be added to a plot.
//!
//! # Line Weights
//!
//! Shapes use semantic line weights from [`line_weights`]:
//! - `THICK` (2.0) for primary elements like arrows and shapes
//! - `NORMAL` (1.5) for standard outlines
//! - `THIN` (1.0) for supporting elements

use std::f64::consts::PI;

use egui::Color32;
use egui_plot::{Line, PlotPoints, Points};

use super::colors::line_weights;

/// Draw a point marker at the specified coordinates.
///
/// # Arguments
/// * `x`, `y` - Coordinates of the point
/// * `color` - Fill color of the marker
/// * `radius` - Visual radius of the marker in plot coordinates
#[must_use]
pub fn point_marker(x: f64, y: f64, color: Color32, radius: f32) -> Points {
    Points::new(vec![[x, y]])
        .color(color)
        .radius(radius)
        .filled(true)
}

/// Draw a point marker with a label.
#[must_use]
pub fn labeled_point(x: f64, y: f64, color: Color32, radius: f32, name: &str) -> Points {
    Points::new(vec![[x, y]])
        .color(color)
        .radius(radius)
        .filled(true)
        .name(name)
}

/// Draw an arrow from a starting point in a given direction.
///
/// # Arguments
/// * `ox`, `oy` - Origin (starting point) of the arrow
/// * `dx`, `dy` - Direction vector (arrow points from origin to origin + direction)
/// * `color` - Color of the arrow
///
/// # Returns
/// A vector of [`Line`] objects forming the arrow shaft and head.
#[must_use]
pub fn arrow_2d(ox: f64, oy: f64, dx: f64, dy: f64, color: Color32) -> Vec<Line> {
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return vec![];
    }

    let head_size = len * 0.15;
    let angle = dy.atan2(dx);
    let head_angle = 0.4; // radians (~23 degrees)

    // Main shaft
    let shaft = Line::new(PlotPoints::new(vec![[ox, oy], [ox + dx, oy + dy]]))
        .color(color)
        .width(line_weights::THICK);

    // Arrowhead wings
    let tip_x = ox + dx;
    let tip_y = oy + dy;

    let h1 = Line::new(PlotPoints::new(vec![
        [tip_x, tip_y],
        [
            tip_x - head_size * (angle + head_angle).cos(),
            tip_y - head_size * (angle + head_angle).sin(),
        ],
    ]))
    .color(color)
    .width(line_weights::THICK);

    let h2 = Line::new(PlotPoints::new(vec![
        [tip_x, tip_y],
        [
            tip_x - head_size * (angle - head_angle).cos(),
            tip_y - head_size * (angle - head_angle).sin(),
        ],
    ]))
    .color(color)
    .width(line_weights::THICK);

    vec![shaft, h1, h2]
}

/// Draw a labeled arrow (vector) with the given name.
#[must_use]
pub fn labeled_arrow(ox: f64, oy: f64, dx: f64, dy: f64, color: Color32, name: &str) -> Vec<Line> {
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return vec![];
    }

    let head_size = len * 0.15;
    let angle = dy.atan2(dx);
    let head_angle = 0.4;

    // Main shaft with name
    let shaft = Line::new(PlotPoints::new(vec![[ox, oy], [ox + dx, oy + dy]]))
        .color(color)
        .width(line_weights::THICK)
        .name(name);

    // Arrowhead wings
    let tip_x = ox + dx;
    let tip_y = oy + dy;

    let h1 = Line::new(PlotPoints::new(vec![
        [tip_x, tip_y],
        [
            tip_x - head_size * (angle + head_angle).cos(),
            tip_y - head_size * (angle + head_angle).sin(),
        ],
    ]))
    .color(color)
    .width(line_weights::THICK);

    let h2 = Line::new(PlotPoints::new(vec![
        [tip_x, tip_y],
        [
            tip_x - head_size * (angle - head_angle).cos(),
            tip_y - head_size * (angle - head_angle).sin(),
        ],
    ]))
    .color(color)
    .width(line_weights::THICK);

    vec![shaft, h1, h2]
}

/// Draw a circle centered at the given point.
///
/// # Arguments
/// * `cx`, `cy` - Center coordinates
/// * `radius` - Radius of the circle
/// * `color` - Line color
/// * `segments` - Number of line segments to approximate the circle
#[must_use]
pub fn circle_2d(cx: f64, cy: f64, radius: f64, color: Color32, segments: usize) -> Line {
    let points: Vec<[f64; 2]> = (0..=segments)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / segments as f64;
            [cx + radius * angle.cos(), cy + radius * angle.sin()]
        })
        .collect();
    Line::new(PlotPoints::new(points))
        .color(color)
        .width(line_weights::NORMAL)
}

/// Draw a filled circle (disk) as a polygon.
#[must_use]
pub fn filled_circle(cx: f64, cy: f64, radius: f64, color: Color32, segments: usize) -> Line {
    // For egui_plot, we use a closed Line with fill
    let points: Vec<[f64; 2]> = (0..=segments)
        .map(|i| {
            let angle = 2.0 * PI * i as f64 / segments as f64;
            [cx + radius * angle.cos(), cy + radius * angle.sin()]
        })
        .collect();
    Line::new(PlotPoints::new(points))
        .color(color)
        .width(line_weights::THIN)
        .fill(0.0) // Fill to y=0 baseline - this is a limitation of egui_plot
}

/// Draw a line segment between two points.
#[must_use]
pub fn line_segment(x1: f64, y1: f64, x2: f64, y2: f64, color: Color32) -> Line {
    Line::new(PlotPoints::new(vec![[x1, y1], [x2, y2]]))
        .color(color)
        .width(line_weights::THICK)
}

/// Draw a labeled line segment.
#[must_use]
pub fn labeled_line_segment(
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    color: Color32,
    name: &str,
) -> Line {
    line_segment(x1, y1, x2, y2, color).name(name)
}

/// Clip a line segment to a square boundary centered at the origin.
///
/// Uses the Liang-Barsky algorithm.
///
/// Returns `Some((x1, y1, x2, y2))` if any part of the line is visible,
/// or `None` if the line is completely outside the bounds.
#[must_use]
fn clip_line_to_bounds(
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    bounds: f64,
) -> Option<(f64, f64, f64, f64)> {
    let dx = x2 - x1;
    let dy = y2 - y1;

    // Parameters for the four boundaries: -bounds <= x,y <= bounds
    let p = [-dx, dx, -dy, dy];
    let q = [x1 + bounds, bounds - x1, y1 + bounds, bounds - y1];

    let mut t0 = 0.0f64;
    let mut t1 = 1.0f64;

    for i in 0..4 {
        if p[i].abs() < 1e-10 {
            // Line is parallel to this boundary
            if q[i] < 0.0 {
                return None; // Line is outside and parallel
            }
        } else {
            let t = q[i] / p[i];
            if p[i] < 0.0 {
                // Entering boundary
                t0 = t0.max(t);
            } else {
                // Leaving boundary
                t1 = t1.min(t);
            }
        }
    }

    if t0 > t1 {
        return None; // Line is completely outside
    }

    // Compute clipped endpoints
    let clipped_x1 = x1 + t0 * dx;
    let clipped_y1 = y1 + t0 * dy;
    let clipped_x2 = x1 + t1 * dx;
    let clipped_y2 = y1 + t1 * dy;

    Some((clipped_x1, clipped_y1, clipped_x2, clipped_y2))
}

/// Draw an "infinite" line given a point and direction, clipped to bounds.
///
/// The line is properly clipped to the viewport boundaries to prevent
/// automatic viewport resizing when lines extend beyond the visible area.
///
/// # Arguments
/// * `px`, `py` - A point on the line
/// * `dx`, `dy` - Direction vector of the line
/// * `bounds` - The viewport extent (line is clipped to [-bounds, bounds] on both axes)
/// * `color` - Line color
#[must_use]
pub fn infinite_line_2d(px: f64, py: f64, dx: f64, dy: f64, bounds: f64, color: Color32) -> Line {
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        // Degenerate case: return a point
        return Line::new(PlotPoints::new(vec![[px, py], [px, py]]))
            .color(color)
            .width(line_weights::THICK);
    }

    // Normalize direction and compute extended line endpoints
    let nx = dx / len;
    let ny = dy / len;
    let t_max = bounds * 3.0; // Extend far enough to cover any viewport position

    let x1 = px - t_max * nx;
    let y1 = py - t_max * ny;
    let x2 = px + t_max * nx;
    let y2 = py + t_max * ny;

    // Clip to viewport bounds
    if let Some((cx1, cy1, cx2, cy2)) = clip_line_to_bounds(x1, y1, x2, y2, bounds) {
        Line::new(PlotPoints::new(vec![[cx1, cy1], [cx2, cy2]]))
            .color(color)
            .width(line_weights::THICK)
    } else {
        // Line is completely outside viewport - return empty line
        Line::new(PlotPoints::new(vec![]))
            .color(color)
            .width(line_weights::THICK)
    }
}

/// Draw a line from its homogeneous coordinates (ax + by + c = 0).
///
/// This is useful for visualizing projective geometry where lines
/// are represented as [a, b, c] vectors.
#[must_use]
pub fn line_from_homogeneous(a: f64, b: f64, c: f64, bounds: f64, color: Color32) -> Line {
    let norm = (a * a + b * b).sqrt();
    if norm < 1e-10 {
        // Degenerate line
        return Line::new(PlotPoints::new(vec![[0.0, 0.0], [0.0, 0.0]]))
            .color(color)
            .width(line_weights::THICK);
    }

    // Direction is perpendicular to (a, b)
    let dx = -b / norm;
    let dy = a / norm;

    // Find a point on the line
    let px = -a * c / (a * a + b * b);
    let py = -b * c / (a * a + b * b);

    infinite_line_2d(px, py, dx, dy, bounds, color)
}

/// Draw an arc (portion of a circle).
///
/// # Arguments
/// * `cx`, `cy` - Center of the arc
/// * `radius` - Radius of the arc
/// * `start_angle` - Starting angle in radians (0 = positive x direction)
/// * `end_angle` - Ending angle in radians
/// * `color` - Line color
#[must_use]
pub fn arc_2d(
    cx: f64,
    cy: f64,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    color: Color32,
) -> Line {
    let arc_length = (end_angle - start_angle).abs();
    let segments = ((arc_length * 20.0).ceil() as usize).max(2);

    let points: Vec<[f64; 2]> = (0..=segments)
        .map(|i| {
            let t = i as f64 / segments as f64;
            let angle = start_angle + t * (end_angle - start_angle);
            [cx + radius * angle.cos(), cy + radius * angle.sin()]
        })
        .collect();

    Line::new(PlotPoints::new(points))
        .color(color)
        .width(line_weights::THICK)
}

/// Draw an arc with an arrow at the end, useful for showing rotation direction.
#[must_use]
pub fn arc_with_arrow(
    cx: f64,
    cy: f64,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    color: Color32,
) -> Vec<Line> {
    let arc = arc_2d(cx, cy, radius, start_angle, end_angle, color);

    // Add arrowhead at end
    let head_size = radius * 0.2;
    let end_x = cx + radius * end_angle.cos();
    let end_y = cy + radius * end_angle.sin();

    // Tangent direction at end (perpendicular to radius, in direction of increasing angle)
    let sign = if end_angle > start_angle { 1.0 } else { -1.0 };
    let tangent_angle = end_angle + sign * PI / 2.0;

    let head_angle = 0.4;
    let h1 = Line::new(PlotPoints::new(vec![
        [end_x, end_y],
        [
            end_x - head_size * (tangent_angle + head_angle).cos(),
            end_y - head_size * (tangent_angle + head_angle).sin(),
        ],
    ]))
    .color(color)
    .width(line_weights::THICK);

    let h2 = Line::new(PlotPoints::new(vec![
        [end_x, end_y],
        [
            end_x - head_size * (tangent_angle - head_angle).cos(),
            end_y - head_size * (tangent_angle - head_angle).sin(),
        ],
    ]))
    .color(color)
    .width(line_weights::THICK);

    vec![arc, h1, h2]
}

/// Draw a bivector as an oriented circular arc at the origin.
///
/// In 2D, a bivector represents an oriented area. This function
/// visualizes it as a partial circle with an arrow indicating orientation.
///
/// # Arguments
/// * `magnitude` - The magnitude determines the arc length (larger = more arc)
/// * `color` - Color of the arc
#[must_use]
pub fn bivector_2d(magnitude: f64, color: Color32) -> Vec<Line> {
    let radius = 0.5;
    let arc_extent = magnitude.clamp(-PI, PI); // Clamp to one full turn
    arc_with_arrow(0.0, 0.0, radius, 0.0, arc_extent, color)
}

/// Draw a polygon from a list of vertices.
#[must_use]
pub fn polygon(vertices: &[[f64; 2]], color: Color32) -> Line {
    let mut points = vertices.to_vec();
    if !points.is_empty() {
        points.push(points[0]); // Close the polygon
    }
    Line::new(PlotPoints::new(points))
        .color(color)
        .width(line_weights::THICK)
}

/// Draw a regular polygon centered at the given point.
#[must_use]
pub fn regular_polygon(
    cx: f64,
    cy: f64,
    radius: f64,
    sides: usize,
    rotation: f64,
    color: Color32,
) -> Line {
    let points: Vec<[f64; 2]> = (0..=sides)
        .map(|i| {
            let angle = rotation + 2.0 * PI * i as f64 / sides as f64;
            [cx + radius * angle.cos(), cy + radius * angle.sin()]
        })
        .collect();
    Line::new(PlotPoints::new(points))
        .color(color)
        .width(line_weights::THICK)
}
