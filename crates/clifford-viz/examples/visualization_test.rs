//! Minimal test example to verify the visualization framework.
//!
//! Run with: `cargo run -p clifford-viz --example visualization_test --release`
//!
//! This example demonstrates:
//! - Basic window setup with `VisualizationApp` trait
//! - 2D grid and axis rendering
//! - Animated rotation using a rotor
//! - Control panel with sliders and animation controls

use clifford::ops::Transform;
use clifford::specialized::euclidean::dim2::{Rotor, Vector};
use clifford_viz::common::prelude::*;
use egui_plot::{Plot, Points};

/// Demo state for 2D rotor visualization.
struct RotorDemo {
    /// Current rotation angle in radians.
    angle: f32,
    /// Animation state.
    animation: Animation,
    /// Input vector (before rotation).
    input_vector: (f32, f32),
    /// Whether to show the grid.
    show_grid: bool,
}

impl Default for RotorDemo {
    fn default() -> Self {
        Self {
            angle: 0.0,
            animation: Animation::with_duration(4.0),
            input_vector: (2.0, 0.5),
            show_grid: true,
        }
    }
}

impl VisualizationApp for RotorDemo {
    fn name(&self) -> &'static str {
        "2D Rotor Demo - Clifford Visualization Framework"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            self.angle = self.animation.angle();
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        // Create the rotor from angle
        let the_rotor = Rotor::<f64>::from_angle(f64::from(self.angle));

        // Create and rotate the input vector
        let input = Vector::new(
            f64::from(self.input_vector.0),
            f64::from(self.input_vector.1),
        );
        let output = the_rotor.transform(&input);

        // Build the plot
        Plot::new("rotor_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .show(ui, |plot_ui| {
                // Draw grid
                if self.show_grid {
                    for l in grid_2d(&ctx, 5.0, 1.0) {
                        plot_ui.line(l);
                    }
                    for axis in axes_2d(&ctx, 5.0) {
                        plot_ui.line(axis);
                    }
                }

                // Draw rotation arc
                if self.angle.abs() > 0.01 {
                    let arc_radius = 0.8;
                    for arc_line in arc_with_arrow(
                        0.0,
                        0.0,
                        arc_radius,
                        0.0,
                        f64::from(self.angle),
                        rotor(&ctx),
                    ) {
                        plot_ui.line(arc_line);
                    }
                }

                // Draw input vector (blue)
                for arrow in labeled_arrow(0.0, 0.0, input.x(), input.y(), point(&ctx), "input") {
                    plot_ui.line(arrow);
                }

                // Draw output vector (green)
                for arrow in labeled_arrow(0.0, 0.0, output.x(), output.y(), plane(&ctx), "output")
                {
                    plot_ui.line(arrow);
                }

                // Draw points at vector tips
                plot_ui.points(
                    Points::new(vec![[input.x(), input.y()]])
                        .color(point(&ctx))
                        .radius(6.0)
                        .filled(true)
                        .name("Input tip"),
                );
                plot_ui.points(
                    Points::new(vec![[output.x(), output.y()]])
                        .color(plane(&ctx))
                        .radius(6.0)
                        .filled(true)
                        .name("Output tip"),
                );

                // Draw a circle showing the rotation path
                let input_len = (input.x() * input.x() + input.y() * input.y()).sqrt();
                plot_ui.line(
                    circle_2d(0.0, 0.0, input_len, with_alpha(grid(&ctx), 100), 64)
                        .name("Rotation path"),
                );
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // Angle control
        angle_slider(ui, "Rotation angle", &mut self.angle);

        ui.separator();

        // Vector input
        vector2_input(
            ui,
            "Input vector",
            &mut self.input_vector.0,
            &mut self.input_vector.1,
        );

        ui.separator();

        // Animation controls
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        ui.separator();

        // Display options
        ui.checkbox(&mut self.show_grid, "Show grid");

        ui.separator();

        // Display rotor components
        section_separator(ui, Some("Rotor Components"));
        let the_rotor = Rotor::<f64>::from_angle(f64::from(self.angle));
        ga_value_display(
            ui,
            "R",
            &[
                ("1", the_rotor.s() as f32),
                ("e\u{2081}\u{2082}", the_rotor.b() as f32),
            ],
        );

        // Display formula
        info_box(
            ui,
            &format!(
                "R = cos(\u{03b8}/2) + sin(\u{03b8}/2)e\u{2081}\u{2082}\n  = {:.3} + {:.3}e\u{2081}\u{2082}",
                the_rotor.s(),
                the_rotor.b()
            ),
        );
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("2D Euclidean rotation using rotors.");
            ui.separator();
            ui.label("Blue: input vector | Green: rotated output | Purple arc: rotation angle");
        });
    }
}

fn main() -> eframe::Result<()> {
    run_app::<RotorDemo>()
}
