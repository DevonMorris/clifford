//! 3D camera with perspective projection for visualization.
//!
//! This module provides a simple 3D camera that projects 3D points to 2D
//! coordinates suitable for rendering with egui_plot. This approach guarantees
//! WASM compatibility since it uses only standard math operations.
//!
//! # Example
//!
//! ```ignore
//! use clifford_viz::common::camera3d::{Camera3D, camera_response};
//!
//! let mut camera = Camera3D::default();
//!
//! // In render(), handle mouse interaction:
//! let response = plot_ui.response();
//! camera_response(&mut camera, &response, ui);
//!
//! // Project 3D points to 2D for drawing:
//! let point_3d = [1.0, 2.0, 3.0];
//! let point_2d = camera.project(point_3d);
//! ```

use clifford::specialized::euclidean::dim3::Vector;
use std::f32::consts::{FRAC_PI_2, PI};

/// A 3D camera with perspective projection and orbit controls.
///
/// The camera uses a spherical coordinate system for easy orbit navigation:
/// - `distance`: Distance from the target point
/// - `azimuth`: Horizontal angle (rotation around Y axis)
/// - `elevation`: Vertical angle (rotation up/down from horizontal)
///
/// The projection maps 3D points to 2D coordinates that can be drawn
/// with egui_plot.
#[derive(Debug, Clone)]
pub struct Camera3D {
    /// Point the camera is looking at (orbit center).
    pub target: [f32; 3],
    /// Distance from target.
    pub distance: f32,
    /// Horizontal angle in radians (0 = looking along +Z).
    pub azimuth: f32,
    /// Vertical angle in radians (0 = horizontal, positive = looking down).
    pub elevation: f32,
    /// Vertical field of view in radians.
    pub fov: f32,
    /// Near clipping plane distance.
    pub near: f32,
    /// Far clipping plane distance.
    pub far: f32,
}

impl Default for Camera3D {
    fn default() -> Self {
        Self {
            target: [0.0, 0.0, 0.0],
            distance: 5.0,
            azimuth: 0.4,   // Slight angle for 3D effect
            elevation: 0.3, // Slight elevation
            fov: PI / 4.0,  // 45 degrees
            near: 0.1,
            far: 100.0,
        }
    }
}

impl Camera3D {
    /// Create a camera at the default isometric-like view.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a camera looking along the X axis.
    #[must_use]
    pub fn view_x() -> Self {
        Self {
            azimuth: FRAC_PI_2,
            elevation: 0.0,
            ..Self::default()
        }
    }

    /// Create a camera looking along the Y axis (top-down view).
    #[must_use]
    pub fn view_y() -> Self {
        Self {
            azimuth: 0.0,
            elevation: FRAC_PI_2 - 0.001, // Just under 90 deg to avoid gimbal issues
            ..Self::default()
        }
    }

    /// Create a camera looking along the Z axis.
    #[must_use]
    pub fn view_z() -> Self {
        Self {
            azimuth: 0.0,
            elevation: 0.0,
            ..Self::default()
        }
    }

    /// Calculate the camera's eye position in world coordinates.
    #[must_use]
    pub fn eye_position(&self) -> [f32; 3] {
        let cos_elev = self.elevation.cos();
        let sin_elev = self.elevation.sin();
        let cos_azim = self.azimuth.cos();
        let sin_azim = self.azimuth.sin();

        [
            self.target[0] + self.distance * cos_elev * sin_azim,
            self.target[1] + self.distance * sin_elev,
            self.target[2] + self.distance * cos_elev * cos_azim,
        ]
    }

    /// Calculate the camera's up vector.
    #[must_use]
    fn up_vector(&self) -> Vector<f32> {
        // Up is always world Y, adjusted for elevation
        let cos_elev = self.elevation.cos();
        let sin_elev = self.elevation.sin();
        let cos_azim = self.azimuth.cos();
        let sin_azim = self.azimuth.sin();

        Vector::new(-sin_elev * sin_azim, cos_elev, -sin_elev * cos_azim)
    }

    /// Project a 3D point to 2D screen coordinates.
    ///
    /// Returns coordinates suitable for egui_plot, where X is horizontal
    /// and Y is vertical. The coordinate system is right-handed.
    ///
    /// Points behind the camera return coordinates at infinity.
    #[must_use]
    pub fn project(&self, point: [f32; 3]) -> [f64; 2] {
        let eye = self.eye_position();
        let up = self.up_vector();

        // Calculate view direction (from eye to target) using clifford Vector
        let forward = Vector::new(
            self.target[0] - eye[0],
            self.target[1] - eye[1],
            self.target[2] - eye[2],
        )
        .normalized();

        // Calculate right vector using GA cross product: forward x up
        let right = forward.cross(up).normalized();

        // Recalculate up to ensure orthogonality
        let up = right.cross(forward);

        // Transform point to camera space
        let rel = Vector::new(point[0] - eye[0], point[1] - eye[1], point[2] - eye[2]);

        // Camera space coordinates using GA dot product
        let cam_x = rel.dot(right);
        let cam_y = rel.dot(up);
        let cam_z = rel.dot(forward);

        // Perspective projection
        if cam_z <= self.near {
            // Point is behind camera or at near plane
            return [f64::INFINITY, f64::INFINITY];
        }

        let scale = 1.0 / (self.fov / 2.0).tan();
        let proj_x = (cam_x / cam_z) * scale;
        let proj_y = (cam_y / cam_z) * scale;

        [proj_x as f64, proj_y as f64]
    }

    /// Project a 3D point and return depth information.
    ///
    /// Returns `(screen_coords, depth)` where depth is the distance
    /// along the view direction. Useful for depth sorting.
    #[must_use]
    pub fn project_with_depth(&self, point: [f32; 3]) -> ([f64; 2], f32) {
        let eye = self.eye_position();
        let up = self.up_vector();

        let forward = Vector::new(
            self.target[0] - eye[0],
            self.target[1] - eye[1],
            self.target[2] - eye[2],
        )
        .normalized();

        let right = forward.cross(up).normalized();
        let up = right.cross(forward);

        let rel = Vector::new(point[0] - eye[0], point[1] - eye[1], point[2] - eye[2]);

        let cam_x = rel.dot(right);
        let cam_y = rel.dot(up);
        let cam_z = rel.dot(forward);

        if cam_z <= self.near {
            return ([f64::INFINITY, f64::INFINITY], f32::INFINITY);
        }

        let scale = 1.0 / (self.fov / 2.0).tan();
        let proj_x = (cam_x / cam_z) * scale;
        let proj_y = (cam_y / cam_z) * scale;

        ([proj_x as f64, proj_y as f64], cam_z)
    }

    /// Orbit the camera by the given angles (in radians).
    ///
    /// Positive `delta_azimuth` rotates counterclockwise (viewed from above).
    /// Positive `delta_elevation` tilts the camera up.
    pub fn orbit(&mut self, delta_azimuth: f32, delta_elevation: f32) {
        self.azimuth += delta_azimuth;
        self.elevation += delta_elevation;

        // Clamp elevation to avoid flipping
        self.elevation = self.elevation.clamp(-FRAC_PI_2 + 0.01, FRAC_PI_2 - 0.01);

        // Wrap azimuth to [-PI, PI]
        while self.azimuth > PI {
            self.azimuth -= 2.0 * PI;
        }
        while self.azimuth < -PI {
            self.azimuth += 2.0 * PI;
        }
    }

    /// Zoom the camera by changing the distance to target.
    ///
    /// Positive `delta` moves closer, negative moves farther.
    pub fn zoom(&mut self, delta: f32) {
        self.distance = (self.distance - delta).clamp(0.5, 50.0);
    }

    /// Pan the camera target in screen space.
    ///
    /// Moves the target point in the plane perpendicular to the view direction.
    pub fn pan(&mut self, delta_x: f32, delta_y: f32) {
        let up = self.up_vector();
        let eye = self.eye_position();

        let forward = Vector::new(
            self.target[0] - eye[0],
            self.target[1] - eye[1],
            self.target[2] - eye[2],
        )
        .normalized();

        let right = forward.cross(up).normalized();
        let up = right.cross(forward);

        // Scale pan by distance for consistent feel
        let scale = self.distance * 0.002;

        self.target[0] -= right.x() * delta_x * scale + up.x() * delta_y * scale;
        self.target[1] -= right.y() * delta_x * scale + up.y() * delta_y * scale;
        self.target[2] -= right.z() * delta_x * scale + up.z() * delta_y * scale;
    }

    /// Reset camera to default view.
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}

// =============================================================================
// UI Integration
// =============================================================================

/// Handle mouse interaction for camera controls.
///
/// Call this with the plot response to enable:
/// - Drag: Orbit camera
/// - Scroll: Zoom
/// - Shift+Drag: Pan
///
/// # Example
///
/// ```ignore
/// let response = plot_ui.response();
/// camera_response(&mut camera, &response, ui);
/// ```
pub fn camera_response(camera: &mut Camera3D, response: &egui::Response, ui: &egui::Ui) {
    // Orbit with drag
    if response.dragged() && !ui.input(|i| i.modifiers.shift) {
        let delta = response.drag_delta();
        camera.orbit(-delta.x * 0.01, -delta.y * 0.01);
    }

    // Pan with shift+drag
    if response.dragged() && ui.input(|i| i.modifiers.shift) {
        let delta = response.drag_delta();
        camera.pan(delta.x, delta.y);
    }

    // Zoom with scroll
    if response.hovered() {
        let scroll = ui.input(|i| i.raw_scroll_delta.y);
        if scroll.abs() > 0.0 {
            camera.zoom(scroll * 0.01);
        }
    }
}

/// Display camera controls in the UI.
///
/// Shows sliders for distance, azimuth, and elevation, plus view preset buttons.
pub fn camera_controls(ui: &mut egui::Ui, camera: &mut Camera3D) {
    egui::CollapsingHeader::new("Camera")
        .default_open(false)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label("Distance:");
                ui.add(
                    egui::DragValue::new(&mut camera.distance)
                        .speed(0.1)
                        .range(0.5..=50.0),
                );
            });

            ui.horizontal(|ui| {
                ui.label("Azimuth:");
                let mut deg = camera.azimuth.to_degrees();
                if ui
                    .add(
                        egui::DragValue::new(&mut deg)
                            .speed(1.0)
                            .range(-180.0..=180.0)
                            .suffix(" deg"),
                    )
                    .changed()
                {
                    camera.azimuth = deg.to_radians();
                }
            });

            ui.horizontal(|ui| {
                ui.label("Elevation:");
                let mut deg = camera.elevation.to_degrees();
                if ui
                    .add(
                        egui::DragValue::new(&mut deg)
                            .speed(1.0)
                            .range(-89.0..=89.0)
                            .suffix(" deg"),
                    )
                    .changed()
                {
                    camera.elevation = deg.to_radians();
                }
            });

            ui.separator();
            ui.label("View presets:");
            ui.horizontal(|ui| {
                if ui.button("Default").clicked() {
                    camera.reset();
                }
                if ui.button("Front").clicked() {
                    *camera = Camera3D::view_z();
                }
                if ui.button("Top").clicked() {
                    *camera = Camera3D::view_y();
                }
                if ui.button("Side").clicked() {
                    *camera = Camera3D::view_x();
                }
            });

            ui.add_space(4.0);
            ui.label(
                egui::RichText::new("Drag to orbit, scroll to zoom, shift+drag to pan")
                    .small()
                    .weak(),
            );
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_camera_eye_position() {
        let camera = Camera3D::default();
        let eye = camera.eye_position();
        // Should be somewhere in front of and above the origin
        assert!(eye[2] > 0.0, "Camera should be in front of target");
    }

    #[test]
    fn test_projection_origin() {
        let camera = Camera3D::default();
        let projected = camera.project([0.0, 0.0, 0.0]);
        // Origin should project near screen center
        assert!(projected[0].abs() < 1.0);
        assert!(projected[1].abs() < 1.0);
    }

    #[test]
    fn test_orbit_clamps_elevation() {
        let mut camera = Camera3D::default();
        camera.orbit(0.0, 10.0); // Try to rotate way up
        assert!(camera.elevation < FRAC_PI_2);
        assert!(camera.elevation > -FRAC_PI_2);
    }

    #[test]
    fn test_zoom_clamps_distance() {
        let mut camera = Camera3D::default();
        camera.zoom(100.0); // Try to zoom way in
        assert!(camera.distance >= 0.5);
        camera.zoom(-100.0); // Try to zoom way out
        assert!(camera.distance <= 50.0);
    }

    #[test]
    fn test_cross_product_handedness() {
        // Verify that X cross Y = Z (right-handed)
        let x = Vector::<f32>::unit_x();
        let y = Vector::<f32>::unit_y();
        let z = x.cross(y);
        assert!((z.x()).abs() < 1e-6);
        assert!((z.y()).abs() < 1e-6);
        assert!((z.z() - 1.0).abs() < 1e-6, "X cross Y should equal Z");
    }
}
