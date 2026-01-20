//! 3D Euclidean rotor visualization demonstrating gimbal-lock-free rotation.
//!
//! This demo compares GA rotors with Euler angles to show why rotors are superior
//! for 3D rotations.

use crate::common::prelude::*;
use clifford::ops::{RightComplement, Transform};
use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
use egui_plot::Plot;
use std::f32::consts::{FRAC_PI_2, TAU};

/// Rotation mode for comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RotationMode {
    /// Single rotor rotation in a plane.
    Rotor,
    /// Euler angles (yaw-pitch-roll).
    Euler,
}

/// Demo for 3D Euclidean rotors.
pub struct Euclidean3Demo {
    /// Camera state.
    camera: Camera3D,
    /// Current rotation mode.
    mode: RotationMode,

    // Rotor mode parameters
    /// Rotation angle in radians.
    rotor_angle: f32,
    /// Rotation axis (will be normalized to unit bivector).
    axis_x: f32,
    /// Rotation axis Y component.
    axis_y: f32,
    /// Rotation axis Z component.
    axis_z: f32,
    /// Animation state.
    animation: Animation,

    // Euler mode parameters
    /// Yaw (rotation around Y).
    euler_yaw: f32,
    /// Pitch (rotation around X).
    euler_pitch: f32,
    /// Roll (rotation around Z).
    euler_roll: f32,

    // Display options
    /// Show coordinate axes.
    show_axes: bool,
    /// Show the rotation plane.
    show_plane: bool,
    /// Show gimbal lock warning.
    gimbal_warning: bool,
}

impl Default for Euclidean3Demo {
    fn default() -> Self {
        Self {
            camera: Camera3D::default(),
            mode: RotationMode::Rotor,
            rotor_angle: 0.0,
            // Default to Z axis (rotation in XY plane)
            axis_x: 0.0,
            axis_y: 0.0,
            axis_z: 1.0,
            animation: Animation::with_duration(4.0),
            euler_yaw: 0.0,
            euler_pitch: 0.0,
            euler_roll: 0.0,
            show_axes: true,
            show_plane: true,
            gimbal_warning: false,
        }
    }
}

impl Euclidean3Demo {
    /// Get the rotation axis as a normalized vector.
    fn axis_vector(&self) -> Vector<f32> {
        let axis = Vector::new(self.axis_x, self.axis_y, self.axis_z);
        let norm = axis.norm();
        if norm < 1e-6 {
            // Default to Z axis if zero
            return Vector::unit_z();
        }
        axis.normalized()
    }

    /// Get the rotation plane bivector (Hodge dual of axis).
    ///
    /// In 3D Euclidean GA, the rotation plane bivector is the
    /// Hodge dual (right complement) of the rotation axis vector.
    fn axis_bivector(&self) -> Bivector<f32> {
        self.axis_vector().right_complement()
    }

    /// Get the current rotor based on mode.
    fn current_rotor(&self) -> Rotor<f32> {
        match self.mode {
            RotationMode::Rotor => Rotor::from_angle_plane(self.rotor_angle, self.axis_bivector()),
            RotationMode::Euler => {
                // Euler angles: yaw (Y) -> pitch (X) -> roll (Z)
                let r_yaw = Rotor::from_angle_plane(self.euler_yaw, Bivector::unit_ry());
                let r_pitch = Rotor::from_angle_plane(self.euler_pitch, Bivector::unit_rx());
                let r_roll = Rotor::from_angle_plane(self.euler_roll, Bivector::unit_rz());
                // Compose: roll first, then pitch, then yaw
                r_yaw * r_pitch * r_roll
            }
        }
    }

    /// Transform a 3D point using the current rotor.
    fn transform_point(&self, point: [f32; 3]) -> [f32; 3] {
        let rotor = self.current_rotor();
        let v = Vector::new(point[0], point[1], point[2]);
        let rotated = rotor.transform(&v);
        [rotated.x(), rotated.y(), rotated.z()]
    }

    /// Check if we're in gimbal lock territory.
    fn check_gimbal_lock(&mut self) {
        if self.mode == RotationMode::Euler {
            // Gimbal lock occurs when pitch is near +/- 90 degrees
            let pitch_deg = self.euler_pitch.to_degrees();
            self.gimbal_warning = (pitch_deg.abs() - 90.0).abs() < 5.0;
        } else {
            self.gimbal_warning = false;
        }
    }
}

impl VisualizationApp for Euclidean3Demo {
    fn name(&self) -> &'static str {
        "3D Euclidean Rotors"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing && self.mode == RotationMode::Rotor {
            // Animation goes 0 to 360 degrees (TAU radians)
            self.rotor_angle = self.animation.progress() * TAU;
        }
        self.check_gimbal_lock();
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        // Get transformed cube vertices
        let base_vertices = unit_cube_vertices();
        let transformed_vertices: [[f32; 3]; 8] =
            std::array::from_fn(|i| self.transform_point(base_vertices[i]));

        let response = Plot::new("euclidean3_view")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_drag(false)
            .allow_scroll(false)
            // Fix auto-resizing by setting explicit bounds
            .include_x(-3.0)
            .include_x(3.0)
            .include_y(-3.0)
            .include_y(3.0)
            .show(ui, |plot_ui| {
                // Draw world coordinate axes
                if self.show_axes {
                    for line in coordinate_axes(&self.camera, 1.8) {
                        plot_ui.line(line);
                    }
                }

                // Draw rotation plane indicator (perpendicular to axis)
                // The plane's normal is the rotation axis (dual of the bivector)
                if self.show_plane && self.mode == RotationMode::Rotor {
                    let axis = self.axis_vector();
                    for line in plane_3d(
                        &self.camera,
                        [0.0, 0.0, 0.0],
                        [axis.x(), axis.y(), axis.z()],
                        1.0,
                        with_alpha(rotor(&ctx), 60),
                        0,
                    ) {
                        plot_ui.line(line);
                    }
                }

                // Draw transformed cube
                let cube_color = if self.gimbal_warning {
                    egui::Color32::from_rgb(255, 100, 100) // Red when gimbal locked
                } else {
                    line(&ctx)
                };
                for line in wireframe_box_vertices(&self.camera, &transformed_vertices, cube_color)
                {
                    plot_ui.line(line);
                }

                // Draw local coordinate frame on the cube
                let origin = self.transform_point([0.0, 0.0, 0.0]);
                let local_x = self.transform_point([0.7, 0.0, 0.0]);
                let local_y = self.transform_point([0.0, 0.7, 0.0]);
                let local_z = self.transform_point([0.0, 0.0, 0.7]);

                for line in arrow_3d(&self.camera, origin, local_x, x_axis(&ctx)) {
                    plot_ui.line(line);
                }
                for line in arrow_3d(&self.camera, origin, local_y, y_axis(&ctx)) {
                    plot_ui.line(line);
                }
                for line in arrow_3d(&self.camera, origin, local_z, z_axis(&ctx)) {
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

        // Mode selection
        ui.heading("Rotation Mode");
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.mode, RotationMode::Rotor, "Rotor (GA)");
            ui.selectable_value(&mut self.mode, RotationMode::Euler, "Euler Angles");
        });

        ui.separator();

        match self.mode {
            RotationMode::Rotor => {
                ui.heading("Rotor Parameters");

                // Axis selection with sliders
                ui.label("Rotation Axis:");
                ui.horizontal(|ui| {
                    ui.colored_label(x_axis(ui.ctx()), "X:");
                    ui.add(egui::Slider::new(&mut self.axis_x, -1.0..=1.0).show_value(false));
                    ui.label(format!("{:.2}", self.axis_x));
                });
                ui.horizontal(|ui| {
                    ui.colored_label(y_axis(ui.ctx()), "Y:");
                    ui.add(egui::Slider::new(&mut self.axis_y, -1.0..=1.0).show_value(false));
                    ui.label(format!("{:.2}", self.axis_y));
                });
                ui.horizontal(|ui| {
                    ui.colored_label(z_axis(ui.ctx()), "Z:");
                    ui.add(egui::Slider::new(&mut self.axis_z, -1.0..=1.0).show_value(false));
                    ui.label(format!("{:.2}", self.axis_z));
                });

                // Quick axis presets
                ui.horizontal(|ui| {
                    if ui.small_button("X").clicked() {
                        self.axis_x = 1.0;
                        self.axis_y = 0.0;
                        self.axis_z = 0.0;
                    }
                    if ui.small_button("Y").clicked() {
                        self.axis_x = 0.0;
                        self.axis_y = 1.0;
                        self.axis_z = 0.0;
                    }
                    if ui.small_button("Z").clicked() {
                        self.axis_x = 0.0;
                        self.axis_y = 0.0;
                        self.axis_z = 1.0;
                    }
                    if ui.small_button("Diagonal").clicked() {
                        self.axis_x = 1.0;
                        self.axis_y = 1.0;
                        self.axis_z = 1.0;
                    }
                });

                ui.separator();

                // Angle slider 0-360 degrees
                let mut angle_deg = self.rotor_angle.to_degrees();
                ui.horizontal(|ui| {
                    ui.label("Angle:");
                    if ui
                        .add(egui::Slider::new(&mut angle_deg, 0.0..=360.0).suffix(" deg"))
                        .changed()
                    {
                        self.rotor_angle = angle_deg.to_radians();
                    }
                });

                ui.separator();
                animation_controls(ui, &mut self.animation);

                // Display rotor components
                ui.separator();
                let rotor = self.current_rotor();
                ga_value_display(
                    ui,
                    "R",
                    &[
                        ("1", rotor.s()),
                        ("e_yz", rotor.rx()),
                        ("e_xz", rotor.ry()),
                        ("e_xy", rotor.rz()),
                    ],
                );
            }
            RotationMode::Euler => {
                ui.heading("Euler Angles (YXZ)");

                // Yaw
                let mut yaw_deg = self.euler_yaw.to_degrees();
                ui.horizontal(|ui| {
                    ui.label("Yaw (Y):");
                    if ui
                        .add(egui::Slider::new(&mut yaw_deg, -180.0..=180.0).suffix(" deg"))
                        .changed()
                    {
                        self.euler_yaw = yaw_deg.to_radians();
                    }
                });

                // Pitch - with gimbal lock warning
                let mut pitch_deg = self.euler_pitch.to_degrees();
                ui.horizontal(|ui| {
                    ui.label("Pitch (X):");
                    let slider =
                        ui.add(egui::Slider::new(&mut pitch_deg, -90.0..=90.0).suffix(" deg"));
                    if slider.changed() {
                        self.euler_pitch = pitch_deg.to_radians();
                    }
                });

                // Roll
                let mut roll_deg = self.euler_roll.to_degrees();
                ui.horizontal(|ui| {
                    ui.label("Roll (Z):");
                    if ui
                        .add(egui::Slider::new(&mut roll_deg, -180.0..=180.0).suffix(" deg"))
                        .changed()
                    {
                        self.euler_roll = roll_deg.to_radians();
                    }
                });

                // Gimbal lock warning
                if self.gimbal_warning {
                    ui.separator();
                    ui.colored_label(
                        egui::Color32::from_rgb(255, 100, 100),
                        "GIMBAL LOCK! At pitch = +/-90 deg,\nyaw and roll affect the same axis.",
                    );
                    ui.label("Try changing yaw and roll - they do the same thing now!");
                }

                // Quick gimbal lock demo button
                ui.separator();
                if ui
                    .button("Set pitch to 90 deg (demo gimbal lock)")
                    .clicked()
                {
                    self.euler_pitch = FRAC_PI_2;
                }
            }
        }

        ui.separator();

        // Display options
        ui.checkbox(&mut self.show_axes, "Show world axes");
        ui.checkbox(&mut self.show_plane, "Show rotation plane");

        ui.separator();
        info_box(ui, "Drag to orbit, scroll to zoom");
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            match self.mode {
                RotationMode::Rotor => {
                    let angle_deg = self.rotor_angle.to_degrees();
                    let axis = self.axis_vector();
                    let biv = self.axis_bivector();
                    ui.label(format!(
                        "Axis: ({:.2}, {:.2}, {:.2}) | Plane: {:.2}e_yz + {:.2}e_zx + {:.2}e_xy | Angle: {:.0} deg",
                        axis.x(), axis.y(), axis.z(),
                        biv.rx(), biv.ry(), biv.rz(),
                        angle_deg
                    ));
                }
                RotationMode::Euler => {
                    if self.gimbal_warning {
                        ui.label("GIMBAL LOCK - degree of freedom lost!");
                    } else {
                        ui.label("Euler mode - yaw/pitch/roll");
                    }
                }
            }
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(EducationalContent {
            title: "3D Rotors vs Euler Angles",
            overview: "This demo compares geometric algebra rotors with traditional Euler angles \
                       for 3D rotation. Rotors are mathematically equivalent to quaternions but \
                       offer a cleaner geometric interpretation: rotation happens IN a plane, \
                       not AROUND an axis.",
            math_background: "ROTOR FORMULA\n\
                             R = cos(theta/2) + sin(theta/2) * B\n\
                             where B is the unit bivector of the rotation plane.\n\n\
                             SANDWICH PRODUCT\n\
                             v' = R v R~  (R~ is the reverse of R)\n\n\
                             EULER ANGLES PROBLEM\n\
                             At pitch = +/-90 deg, yaw and roll rotate around\n\
                             the same axis. This \"gimbal lock\" loses one\n\
                             degree of freedom.\n\n\
                             WHY ROTORS WIN\n\
                             - No gimbal lock (any orientation reachable)\n\
                             - Smooth interpolation (SLERP)\n\
                             - Easy composition (just multiply)\n\
                             - Works identically in 2D, 3D, 4D, ...",
            how_to_use: "ROTOR MODE:\n\
                        - Set any rotation axis with X/Y/Z sliders\n\
                        - Adjust the angle (0-360 deg) or animate\n\
                        - Watch the cube rotate smoothly around any axis\n\n\
                        EULER MODE:\n\
                        - Adjust yaw, pitch, and roll\n\
                        - Set pitch to 90 deg to see gimbal lock\n\
                        - Notice how yaw and roll do the same thing!\n\n\
                        CAMERA:\n\
                        - Drag to orbit\n\
                        - Scroll to zoom",
            key_concepts: "- Rotation happens IN a plane, not AROUND an axis\n\
                          - The \"axis\" is just the dual of the rotation bivector\n\
                          - Half-angle encoding (theta/2) is why rotors compose correctly\n\
                          - Gimbal lock is a fundamental flaw of Euler angles\n\
                          - Rotors and quaternions are the same thing in different notation",
            resources: &[
                (
                    "Rigid Geometric Algebra Wiki",
                    "https://rigidgeometricalgebra.org",
                ),
                (
                    "Visualizing Quaternions (video)",
                    "https://www.youtube.com/watch?v=d4EgbgTm0Bg",
                ),
            ],
        })
    }
}
