//! Test demo for 3D visualization infrastructure.
//!
//! This demo validates that the 3D projection and camera controls work
//! correctly, especially in WASM.

use crate::common::prelude::*;
use egui_plot::Plot;

/// Demo state for 3D visualization test.
pub struct Test3DDemo {
    /// Camera state.
    camera: Camera3D,
    /// Rotation angle for animated cube.
    angle: f32,
    /// Animation state.
    animation: Animation,
    /// Whether to show coordinate axes.
    show_axes: bool,
    /// Whether to auto-rotate the cube.
    auto_rotate: bool,
}

impl Default for Test3DDemo {
    fn default() -> Self {
        Self {
            camera: Camera3D::default(),
            angle: 0.0,
            animation: Animation::with_duration(4.0),
            show_axes: true,
            auto_rotate: false,
        }
    }
}

impl Test3DDemo {
    /// Rotate a 3D point around the Y axis.
    fn rotate_y(&self, point: [f32; 3], angle: f32) -> [f32; 3] {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        [
            point[0] * cos_a + point[2] * sin_a,
            point[1],
            -point[0] * sin_a + point[2] * cos_a,
        ]
    }
}

impl VisualizationApp for Test3DDemo {
    fn name(&self) -> &'static str {
        "3D Visualization Test"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing || self.auto_rotate {
            self.angle = self.animation.angle();
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        // Get rotated cube vertices
        let base_vertices = unit_cube_vertices();
        let rotated_vertices: [[f32; 3]; 8] =
            std::array::from_fn(|i| self.rotate_y(base_vertices[i], self.angle));

        // Build the plot
        let response = Plot::new("3d_view")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_drag(false) // We handle drag for camera
            .allow_scroll(false) // We handle scroll for zoom
            .show(ui, |plot_ui| {
                // Draw coordinate axes
                if self.show_axes {
                    for line in coordinate_axes(&self.camera, 1.5) {
                        plot_ui.line(line);
                    }
                }

                // Draw wireframe cube
                for line in wireframe_box_vertices(&self.camera, &rotated_vertices, line(&ctx)) {
                    plot_ui.line(line);
                }

                // Draw a ground plane
                for line in plane_3d(
                    &self.camera,
                    [0.0, -0.6, 0.0],
                    [0.0, 1.0, 0.0],
                    1.5,
                    with_alpha(grid(&ctx), 80),
                    4,
                ) {
                    plot_ui.line(line);
                }
            });

        // Handle camera interaction
        camera_response(&mut self.camera, &response.response, ui);
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // Camera controls
        camera_controls(ui, &mut self.camera);

        ui.separator();

        // Rotation controls
        angle_slider(ui, "Cube rotation", &mut self.angle);
        ui.checkbox(&mut self.auto_rotate, "Auto-rotate");

        if self.auto_rotate {
            animation_controls(ui, &mut self.animation);
        }

        ui.separator();

        // Display options
        ui.checkbox(&mut self.show_axes, "Show axes");

        ui.separator();

        // Info
        info_box(
            ui,
            "Drag to orbit camera\nScroll to zoom\nShift+drag to pan",
        );
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("3D visualization using custom projection via egui_plot.");
            ui.separator();
            let eye = self.camera.eye_position();
            ui.label(format!(
                "Camera: ({:.1}, {:.1}, {:.1})",
                eye[0], eye[1], eye[2]
            ));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(EducationalContent {
            title: "3D Projection",
            overview: "This demo tests the 3D visualization infrastructure. \
                       It renders a wireframe cube using perspective projection, \
                       all computed in pure Rust with no external 3D libraries.",
            math_background: "Perspective projection transforms 3D points to 2D:\n\n\
                             1. Transform to camera space (view matrix)\n\
                             2. Apply perspective division: x' = x/z, y' = y/z\n\
                             3. Scale by field of view\n\n\
                             The camera uses spherical coordinates (azimuth, elevation, distance) \
                             for intuitive orbit controls.",
            how_to_use: "- Drag the plot to orbit the camera around the cube\n\
                        - Scroll to zoom in/out\n\
                        - Hold Shift and drag to pan\n\
                        - Use the Camera section to see/edit exact values\n\
                        - Enable auto-rotate to see continuous animation",
            key_concepts: "- Custom 3D projection works in WASM (no WebGL required)\n\
                          - Wireframe rendering is efficient and clear\n\
                          - Camera orbit controls are intuitive",
            resources: &[],
        })
    }
}
