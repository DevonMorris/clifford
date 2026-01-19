//! 3D shape primitives for visualization.
//!
//! This module provides functions to render 3D shapes by projecting them
//! through a [`Camera3D`] to 2D coordinates for egui_plot.
//!
//! All functions return `Vec<Line>` or `Line` that can be added to a plot.
//!
//! # Example
//!
//! ```ignore
//! use clifford_viz::common::camera3d::Camera3D;
//! use clifford_viz::common::shapes3d::*;
//!
//! let camera = Camera3D::default();
//!
//! Plot::new("3d_view")
//!     .data_aspect(1.0)
//!     .show(ui, |plot_ui| {
//!         // Draw coordinate axes
//!         for line in coordinate_axes(&camera, 2.0) {
//!             plot_ui.line(line);
//!         }
//!
//!         // Draw a wireframe box
//!         for line in wireframe_box(&camera, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], Color32::WHITE) {
//!             plot_ui.line(line);
//!         }
//!     });
//! ```

use clifford::specialized::euclidean::dim3::Vector;
use egui::Color32;
use egui_plot::{Line, PlotPoints};
use std::f32::consts::PI;

use super::camera3d::Camera3D;

// =============================================================================
// Coordinate Axes
// =============================================================================

/// Draw RGB coordinate axes at the origin.
///
/// - X axis: Red
/// - Y axis: Green
/// - Z axis: Blue
///
/// Returns 6 lines (3 axis lines + 3 arrowheads).
#[must_use]
pub fn coordinate_axes(camera: &Camera3D, length: f32) -> Vec<Line> {
    let mut lines = Vec::new();

    // X axis (red)
    lines.extend(arrow_3d(
        camera,
        [0.0, 0.0, 0.0],
        [length, 0.0, 0.0],
        Color32::from_rgb(244, 67, 54),
    ));

    // Y axis (green)
    lines.extend(arrow_3d(
        camera,
        [0.0, 0.0, 0.0],
        [0.0, length, 0.0],
        Color32::from_rgb(76, 175, 80),
    ));

    // Z axis (blue)
    lines.extend(arrow_3d(
        camera,
        [0.0, 0.0, 0.0],
        [0.0, 0.0, length],
        Color32::from_rgb(33, 150, 243),
    ));

    lines
}

/// Draw coordinate axes with labels.
///
/// Returns lines and label positions for each axis.
#[must_use]
pub fn coordinate_axes_labeled(camera: &Camera3D, length: f32) -> (Vec<Line>, Vec<([f64; 2], &'static str)>) {
    let lines = coordinate_axes(camera, length);

    let labels = vec![
        (camera.project([length * 1.1, 0.0, 0.0]), "X"),
        (camera.project([0.0, length * 1.1, 0.0]), "Y"),
        (camera.project([0.0, 0.0, length * 1.1]), "Z"),
    ];

    (lines, labels)
}

// =============================================================================
// Basic Shapes
// =============================================================================

/// Draw a line segment in 3D.
#[must_use]
pub fn line_3d(camera: &Camera3D, start: [f32; 3], end: [f32; 3], color: Color32) -> Line {
    let p1 = camera.project(start);
    let p2 = camera.project(end);

    Line::new(PlotPoints::new(vec![p1, p2]))
        .color(color)
        .width(1.5)
}

/// Draw an arrow in 3D (line with arrowhead).
///
/// Returns multiple lines: the shaft and arrowhead.
#[must_use]
pub fn arrow_3d(camera: &Camera3D, start: [f32; 3], end: [f32; 3], color: Color32) -> Vec<Line> {
    let mut lines = Vec::new();

    // Main shaft
    lines.push(line_3d(camera, start, end, color));

    // Arrowhead using clifford Vector
    let start_v = Vector::new(start[0], start[1], start[2]);
    let end_v = Vector::new(end[0], end[1], end[2]);
    let dir = end_v - start_v;
    let len = dir.norm();

    if len > 1e-6 {
        let head_size = len * 0.1;
        let dir_norm = dir.normalized();

        // Find perpendicular vectors using GA operations
        let perp1 = perpendicular_vector(dir_norm);
        let perp2 = dir_norm.cross(perp1);

        // Arrowhead points
        let head_base = end_v - dir_norm * head_size;

        for angle in [0.0, PI * 2.0 / 3.0, PI * 4.0 / 3.0] {
            let c = angle.cos() * head_size * 0.4;
            let s = angle.sin() * head_size * 0.4;
            let point = head_base + perp1 * c + perp2 * s;
            lines.push(line_3d(camera, end, [point.x(), point.y(), point.z()], color));
        }
    }

    lines
}

/// Draw a point marker in 3D (rendered as a small circle).
#[must_use]
pub fn point_3d(camera: &Camera3D, position: [f32; 3], radius: f32, color: Color32) -> Vec<Line> {
    // Project the point and draw a small circle around it
    let center = camera.project(position);

    // Check if point is visible
    if center[0].is_infinite() || center[1].is_infinite() {
        return Vec::new();
    }

    // Draw a circle in screen space (not a true 3D sphere projection)
    // Scale radius by distance for approximate perspective
    let (_, depth) = camera.project_with_depth(position);
    let screen_radius = (radius / depth * 50.0) as f64;
    let screen_radius = screen_radius.clamp(2.0, 20.0);

    let segments = 16;
    let points: Vec<[f64; 2]> = (0..=segments)
        .map(|i| {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / segments as f64;
            [
                center[0] + screen_radius * 0.02 * angle.cos(),
                center[1] + screen_radius * 0.02 * angle.sin(),
            ]
        })
        .collect();

    vec![Line::new(PlotPoints::new(points)).color(color).width(2.0)]
}

// =============================================================================
// Wireframe Shapes
// =============================================================================

/// Draw a wireframe box (12 edges).
///
/// # Arguments
/// * `center` - Center of the box
/// * `size` - Half-extents in each dimension [half_x, half_y, half_z]
/// * `color` - Line color
#[must_use]
pub fn wireframe_box(
    camera: &Camera3D,
    center: [f32; 3],
    size: [f32; 3],
    color: Color32,
) -> Vec<Line> {
    let [cx, cy, cz] = center;
    let [sx, sy, sz] = size;

    // 8 vertices of the box
    let vertices = [
        [cx - sx, cy - sy, cz - sz], // 0: back-bottom-left
        [cx + sx, cy - sy, cz - sz], // 1: back-bottom-right
        [cx + sx, cy + sy, cz - sz], // 2: back-top-right
        [cx - sx, cy + sy, cz - sz], // 3: back-top-left
        [cx - sx, cy - sy, cz + sz], // 4: front-bottom-left
        [cx + sx, cy - sy, cz + sz], // 5: front-bottom-right
        [cx + sx, cy + sy, cz + sz], // 6: front-top-right
        [cx - sx, cy + sy, cz + sz], // 7: front-top-left
    ];

    // 12 edges
    let edges = [
        // Back face
        (0, 1), (1, 2), (2, 3), (3, 0),
        // Front face
        (4, 5), (5, 6), (6, 7), (7, 4),
        // Connecting edges
        (0, 4), (1, 5), (2, 6), (3, 7),
    ];

    edges
        .iter()
        .map(|&(a, b)| line_3d(camera, vertices[a], vertices[b], color))
        .collect()
}

/// Draw a wireframe box with specified vertices (8 points).
///
/// Useful when the box has been transformed by a rotor/motor.
#[must_use]
pub fn wireframe_box_vertices(camera: &Camera3D, vertices: &[[f32; 3]; 8], color: Color32) -> Vec<Line> {
    let edges = [
        // Back face
        (0, 1), (1, 2), (2, 3), (3, 0),
        // Front face
        (4, 5), (5, 6), (6, 7), (7, 4),
        // Connecting edges
        (0, 4), (1, 5), (2, 6), (3, 7),
    ];

    edges
        .iter()
        .map(|&(a, b)| line_3d(camera, vertices[a], vertices[b], color))
        .collect()
}

/// Generate the 8 vertices of a unit cube centered at origin.
///
/// Useful as input to transformations.
#[must_use]
pub fn unit_cube_vertices() -> [[f32; 3]; 8] {
    [
        [-0.5, -0.5, -0.5],
        [ 0.5, -0.5, -0.5],
        [ 0.5,  0.5, -0.5],
        [-0.5,  0.5, -0.5],
        [-0.5, -0.5,  0.5],
        [ 0.5, -0.5,  0.5],
        [ 0.5,  0.5,  0.5],
        [-0.5,  0.5,  0.5],
    ]
}

/// Draw a wireframe sphere (latitude/longitude lines).
///
/// # Arguments
/// * `center` - Center of the sphere
/// * `radius` - Sphere radius
/// * `color` - Line color
/// * `segments` - Number of latitude/longitude lines (more = smoother)
#[must_use]
pub fn wireframe_sphere(
    camera: &Camera3D,
    center: [f32; 3],
    radius: f32,
    color: Color32,
    segments: usize,
) -> Vec<Line> {
    let mut lines = Vec::new();
    let segments = segments.max(4);

    // Longitude lines (vertical great circles)
    for i in 0..segments {
        let angle = 2.0 * PI * i as f32 / segments as f32;
        let cos_a = angle.cos();
        let sin_a = angle.sin();

        let points: Vec<[f64; 2]> = (0..=segments * 2)
            .map(|j| {
                let phi = PI * j as f32 / (segments * 2) as f32 - PI / 2.0;
                let cos_phi = phi.cos();
                let sin_phi = phi.sin();
                let point = [
                    center[0] + radius * cos_phi * cos_a,
                    center[1] + radius * sin_phi,
                    center[2] + radius * cos_phi * sin_a,
                ];
                camera.project(point)
            })
            .filter(|p| p[0].is_finite() && p[1].is_finite())
            .collect();

        if points.len() >= 2 {
            lines.push(Line::new(PlotPoints::new(points)).color(color).width(1.0));
        }
    }

    // Latitude lines (horizontal circles)
    for i in 1..segments {
        let phi = PI * i as f32 / segments as f32 - PI / 2.0;
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let r = radius * cos_phi;
        let y = center[1] + radius * sin_phi;

        let points: Vec<[f64; 2]> = (0..=segments * 2)
            .map(|j| {
                let angle = 2.0 * PI * j as f32 / (segments * 2) as f32;
                let point = [
                    center[0] + r * angle.cos(),
                    y,
                    center[2] + r * angle.sin(),
                ];
                camera.project(point)
            })
            .filter(|p| p[0].is_finite() && p[1].is_finite())
            .collect();

        if points.len() >= 2 {
            lines.push(Line::new(PlotPoints::new(points)).color(color).width(1.0));
        }
    }

    lines
}

/// Draw a circle in 3D space.
///
/// # Arguments
/// * `center` - Center of the circle
/// * `normal` - Normal vector to the circle's plane
/// * `radius` - Circle radius
/// * `color` - Line color
/// * `segments` - Number of segments (more = smoother)
#[must_use]
pub fn circle_3d(
    camera: &Camera3D,
    center: [f32; 3],
    normal: [f32; 3],
    radius: f32,
    color: Color32,
    segments: usize,
) -> Line {
    let normal_v = Vector::new(normal[0], normal[1], normal[2]).normalized();
    let center_v = Vector::new(center[0], center[1], center[2]);

    // Find two perpendicular vectors in the circle's plane using GA
    let u = perpendicular_vector(normal_v);
    let v = normal_v.cross(u);

    let points: Vec<[f64; 2]> = (0..=segments)
        .map(|i| {
            let angle = 2.0 * PI * i as f32 / segments as f32;
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            let point = center_v + u * (radius * cos_a) + v * (radius * sin_a);
            camera.project([point.x(), point.y(), point.z()])
        })
        .collect();

    Line::new(PlotPoints::new(points)).color(color).width(1.5)
}

/// Draw a plane (as a bounded quad with optional grid).
///
/// # Arguments
/// * `center` - Center of the plane
/// * `normal` - Normal vector
/// * `size` - Half-extent of the visible plane
/// * `color` - Line color
/// * `grid_lines` - Number of grid lines (0 for just the outline)
#[must_use]
pub fn plane_3d(
    camera: &Camera3D,
    center: [f32; 3],
    normal: [f32; 3],
    size: f32,
    color: Color32,
    grid_lines: usize,
) -> Vec<Line> {
    let mut lines = Vec::new();
    let normal_v = Vector::new(normal[0], normal[1], normal[2]).normalized();
    let center_v = Vector::new(center[0], center[1], center[2]);

    let u = perpendicular_vector(normal_v);
    let v = normal_v.cross(u);

    // Corner vertices
    let corners: [Vector<f32>; 4] = [
        center_v - u * size - v * size,
        center_v + u * size - v * size,
        center_v + u * size + v * size,
        center_v - u * size + v * size,
    ];

    // Outline
    for i in 0..4 {
        let a = corners[i];
        let b = corners[(i + 1) % 4];
        lines.push(line_3d(
            camera,
            [a.x(), a.y(), a.z()],
            [b.x(), b.y(), b.z()],
            color,
        ));
    }

    // Grid lines
    if grid_lines > 0 {
        let faded = Color32::from_rgba_unmultiplied(
            color.r(),
            color.g(),
            color.b(),
            color.a() / 2,
        );

        for i in 1..=grid_lines {
            let t = i as f32 / (grid_lines + 1) as f32;
            let t = -size + 2.0 * size * t;

            // Lines parallel to u
            let start = center_v - u * size + v * t;
            let end = center_v + u * size + v * t;
            lines.push(line_3d(
                camera,
                [start.x(), start.y(), start.z()],
                [end.x(), end.y(), end.z()],
                faded,
            ));

            // Lines parallel to v
            let start = center_v + u * t - v * size;
            let end = center_v + u * t + v * size;
            lines.push(line_3d(
                camera,
                [start.x(), start.y(), start.z()],
                [end.x(), end.y(), end.z()],
                faded,
            ));
        }
    }

    lines
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Find a vector perpendicular to the given unit vector using GA operations.
fn perpendicular_vector(v: Vector<f32>) -> Vector<f32> {
    // Choose the axis least aligned with v
    let abs_x = v.x().abs();
    let abs_y = v.y().abs();
    let abs_z = v.z().abs();

    let axis = if abs_x <= abs_y && abs_x <= abs_z {
        Vector::unit_x()
    } else if abs_y <= abs_z {
        Vector::unit_y()
    } else {
        Vector::unit_z()
    };

    v.cross(axis).normalized()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wireframe_box_has_12_edges() {
        let camera = Camera3D::default();
        let lines = wireframe_box(&camera, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], Color32::WHITE);
        assert_eq!(lines.len(), 12);
    }

    #[test]
    fn test_coordinate_axes_has_lines() {
        let camera = Camera3D::default();
        let lines = coordinate_axes(&camera, 1.0);
        // 3 axes * (1 shaft + 3 arrowhead lines) = 12 lines
        assert!(!lines.is_empty());
    }

    #[test]
    fn test_unit_cube_vertices() {
        let vertices = unit_cube_vertices();
        assert_eq!(vertices.len(), 8);
        // Check that all vertices are within unit cube bounds
        for v in &vertices {
            assert!(v[0].abs() <= 0.5);
            assert!(v[1].abs() <= 0.5);
            assert!(v[2].abs() <= 0.5);
        }
    }

    #[test]
    fn test_perpendicular_is_orthogonal() {
        let v = Vector::new(1.0, 2.0, 3.0).normalized();
        let perp = perpendicular_vector(v);
        // Dot product should be near zero
        assert!(v.dot(perp).abs() < 1e-5);
    }
}
