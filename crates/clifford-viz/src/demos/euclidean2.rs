//! 2D Euclidean Rotor Visualization
//!
//! This demo demonstrates rotation in 2D Euclidean space using rotors from
//! Geometric Algebra. It visualizes:
//!
//! - A vector field that rotates uniformly under rotor action
//! - The bivector (rotation plane) as an oriented arc
//! - Rotor components in the `cos(θ/2) + sin(θ/2)e₁₂` form
//!
//! ## Mathematical Background
//!
//! In 2D Geometric Algebra, a rotor `R` encodes a rotation by angle `θ` as:
//!
//! ```text
//! R = cos(θ/2) + sin(θ/2)e₁₂
//! ```
//!
//! where `e₁₂` is the unit bivector (representing the oriented plane of rotation).
//!
//! Vectors transform via the sandwich product: `v' = RvR†`
//!
//! Note that the rotor uses half the angle - this is the "spinor" representation
//! that makes rotation composition natural: `R_total = R_2 R_1`.

use crate::common::prelude::*;
use clifford::ops::Transform;
use clifford::specialized::euclidean::dim2::{Rotor, Vector};
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
/// Using fixed bounds prevents disorienting viewport jumps.
const VIEWPORT_BOUNDS: f64 = 5.0;

/// Demo state for 2D Euclidean rotor visualization.
pub struct Euclidean2Demo {
    /// Current rotation angle in radians.
    angle: f32,
    /// Dilation factor (1.0 = no scaling, >1 = expand, <1 = contract).
    /// The rotor magnitude will be sqrt(dilation) so sandwich product scales by dilation.
    dilation: f32,
    /// Animation state.
    animation: Animation,
    /// Grid size for vector field (vectors go from -grid_size to +grid_size).
    grid_size: i32,
    /// Spacing between vectors in the field.
    vector_spacing: f32,
    /// Whether to show the bivector arc (half-angle visualization).
    show_bivector: bool,
    /// Whether to show the full rotation arc.
    show_rotation_arc: bool,
    /// Whether to show the coordinate grid.
    show_grid: bool,
    /// Whether to show the vector field.
    show_vector_field: bool,
    /// A single highlighted input vector for detailed visualization.
    highlight_vector: (f32, f32),
    /// Whether to show the highlighted vector and its rotation.
    show_highlight: bool,
}

impl Default for Euclidean2Demo {
    fn default() -> Self {
        Self {
            angle: 0.0,
            dilation: 1.0,
            animation: Animation::with_duration(4.0),
            grid_size: 3,
            vector_spacing: 1.0,
            show_bivector: true,
            show_rotation_arc: true,
            show_grid: true,
            show_vector_field: true,
            highlight_vector: (2.0, 0.5),
            show_highlight: true,
        }
    }
}

impl VisualizationApp for Euclidean2Demo {
    fn name(&self) -> &'static str {
        "Euclidean 2D - Rotor Animation"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            self.angle = self.animation.angle();
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        // Create the rotor from the current angle and dilation
        let the_rotor = Rotor::with_dilation(f64::from(self.angle), f64::from(self.dilation));
        let ctx = ui.ctx().clone();

        Plot::new("euclidean2_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .show(ui, |plot_ui| {
                // Draw coordinate grid (using fixed bounds for stability)
                if self.show_grid {
                    for l in grid_2d(&ctx, VIEWPORT_BOUNDS, f64::from(self.vector_spacing)) {
                        plot_ui.line(l);
                    }
                    for axis in axes_2d(&ctx, VIEWPORT_BOUNDS) {
                        plot_ui.line(axis);
                    }
                }

                // Draw bivector arc (half-angle visualization)
                if self.show_bivector && self.angle.abs() > 0.01 {
                    let bivector_radius = 0.4;
                    for arc_line in arc_with_arrow(
                        0.0,
                        0.0,
                        bivector_radius,
                        0.0,
                        f64::from(self.angle) / 2.0,
                        rotor(&ctx),
                    ) {
                        plot_ui.line(arc_line);
                    }
                }

                // Draw full rotation arc
                if self.show_rotation_arc && self.angle.abs() > 0.01 {
                    let rotation_radius = 0.7;
                    for arc_line in arc_with_arrow(
                        0.0,
                        0.0,
                        rotation_radius,
                        0.0,
                        f64::from(self.angle),
                        with_alpha(line(&ctx), 150),
                    ) {
                        plot_ui.line(arc_line);
                    }
                }

                // Draw vector field
                if self.show_vector_field {
                    let n = self.grid_size;
                    for i in -n..=n {
                        for j in -n..=n {
                            // Skip origin
                            if i == 0 && j == 0 {
                                continue;
                            }

                            // Create base vector at grid position
                            let base_x = f64::from(i) * f64::from(self.vector_spacing);
                            let base_y = f64::from(j) * f64::from(self.vector_spacing);

                            // Create a unit vector pointing outward from origin
                            let len = (base_x * base_x + base_y * base_y).sqrt();
                            let unit_x = base_x / len;
                            let unit_y = base_y / len;

                            // Apply rotor to this direction vector
                            let direction = Vector::new(unit_x, unit_y);
                            let rotated = the_rotor.transform(&direction);

                            // Scale for visibility
                            let arrow_scale = 0.4;
                            let dx = rotated.x() * arrow_scale;
                            let dy = rotated.y() * arrow_scale;

                            // Draw from grid position with softer color
                            for arrow in arrow_2d(
                                base_x - dx / 2.0,
                                base_y - dy / 2.0,
                                dx,
                                dy,
                                with_alpha(point(&ctx), 160),
                            ) {
                                plot_ui.line(arrow);
                            }
                        }
                    }
                }

                // Draw highlighted input/output vectors
                if self.show_highlight {
                    let input = Vector::new(
                        f64::from(self.highlight_vector.0),
                        f64::from(self.highlight_vector.1),
                    );
                    let output = the_rotor.transform(&input);

                    // Circle showing the rotation path (input magnitude)
                    let input_len = (input.x() * input.x() + input.y() * input.y()).sqrt();
                    plot_ui.line(
                        circle_2d(0.0, 0.0, input_len, with_alpha(grid_major(&ctx), 100), 64)
                            .name("Input magnitude"),
                    );

                    // If dilating, show the output magnitude circle too
                    if (self.dilation - 1.0).abs() > 0.01 {
                        let output_len = (output.x() * output.x() + output.y() * output.y()).sqrt();
                        plot_ui.line(
                            circle_2d(0.0, 0.0, output_len, with_alpha(plane(&ctx), 80), 64)
                                .name("Output magnitude"),
                        );
                    }

                    // Input vector (soft blue)
                    for arrow in labeled_arrow(0.0, 0.0, input.x(), input.y(), point(&ctx), "v") {
                        plot_ui.line(arrow);
                    }

                    // Output vector (soft green)
                    for arrow in labeled_arrow(0.0, 0.0, output.x(), output.y(), plane(&ctx), "Rv")
                    {
                        plot_ui.line(arrow);
                    }

                    // Point markers at tips
                    plot_ui.points(
                        Points::new(vec![[input.x(), input.y()]])
                            .color(point(&ctx))
                            .radius(6.0)
                            .filled(true)
                            .name("Input v"),
                    );
                    plot_ui.points(
                        Points::new(vec![[output.x(), output.y()]])
                            .color(plane(&ctx))
                            .radius(6.0)
                            .filled(true)
                            .name("Rotated Rv"),
                    );
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Primary Control: Rotation ===
        group_header(ui, "Rotation");
        angle_slider_range(ui, "Angle theta", &mut self.angle, -360.0, 360.0);

        // Dilation control
        ui.horizontal(|ui| {
            ui.label("Dilation");
            ui.add(
                egui::Slider::new(&mut self.dilation, 0.25..=4.0)
                    .logarithmic(true)
                    .fixed_decimals(2),
            );
            if ui.button("Reset").clicked() {
                self.dilation = 1.0;
            }
        });

        // Animation controls
        ui.add_space(spacing::XS);
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // === Rotor Components ===
        section_separator(ui, Some("Rotor Components"));
        let rotor = Rotor::with_dilation(f64::from(self.angle), f64::from(self.dilation));
        let rotor_magnitude = rotor.dilation_factor().sqrt();
        ga_value_display(
            ui,
            "R",
            &[("1", rotor.s() as f32), ("e_1_2", rotor.b() as f32)],
        );
        value_display(ui, "|R|", rotor_magnitude as f32, 4);

        // Formula box
        if (self.dilation - 1.0).abs() < 0.01 {
            info_box(
                ui,
                &format!(
                    "R = cos(theta/2) - sin(theta/2)e_1_2\n  = {:.4} - {:.4}e_1_2\n|R| = 1 (unit rotor, pure rotation)",
                    (self.angle / 2.0).cos(),
                    (self.angle / 2.0).sin()
                ),
            );
        } else {
            info_box(
                ui,
                &format!(
                    "R = k(cos(theta/2) - sin(theta/2)e_1_2)\nk = sqrt{:.2} = {:.4}\n|R| = {:.4} (scales by |R|^2 = {:.2})",
                    self.dilation,
                    self.dilation.sqrt(),
                    rotor_magnitude,
                    self.dilation
                ),
            );
        }

        // === Highlighted Vector ===
        section_separator(ui, Some("Highlighted Vector"));
        ui.checkbox(&mut self.show_highlight, "Show highlighted vector");
        if self.show_highlight {
            vector2_input_range(
                ui,
                "Vector v",
                &mut self.highlight_vector.0,
                &mut self.highlight_vector.1,
                -5.0..=5.0,
            );

            // Show input/output values
            let input = Vector::new(
                f64::from(self.highlight_vector.0),
                f64::from(self.highlight_vector.1),
            );
            let output = rotor.transform(&input);
            let input_mag = (input.x() * input.x() + input.y() * input.y()).sqrt();
            let output_mag = (output.x() * output.x() + output.y() * output.y()).sqrt();

            point2_display(ui, "v", input.x() as f32, input.y() as f32);
            point2_display(ui, "Rv", output.x() as f32, output.y() as f32);
            if (self.dilation - 1.0).abs() > 0.01 {
                ui.horizontal(|ui| {
                    ui.label(format!("|v| = {:.3}", input_mag));
                    ui.label("->");
                    ui.label(format!("|Rv| = {:.3}", output_mag));
                    ui.label(format!("(x{:.2})", output_mag / input_mag));
                });
            }
        }

        // === Display Options (compact toggle row) ===
        section_separator(ui, Some("Display Options"));
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_grid, "Grid");
            ui.checkbox(&mut self.show_vector_field, "Field");
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_bivector, "Bivector (theta/2)");
            ui.checkbox(&mut self.show_rotation_arc, "Arc (theta)");
        });

        // Vector field options
        if self.show_vector_field {
            ui.horizontal(|ui| {
                ui.label("Grid size:");
                ui.add(egui::Slider::new(&mut self.grid_size, 1..=5));
            });
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            if (self.dilation - 1.0).abs() < 0.01 {
                ui.label("2D rotation using unit rotor (pure rotation, no scaling)");
            } else {
                ui.label(format!(
                    "2D rotation + dilation (|R|={:.2}, scales by {:.2}x)",
                    self.dilation.sqrt(),
                    self.dilation
                ));
            }
            ui.separator();
            ui.colored_label(point(&ctx), "Input v");
            ui.separator();
            ui.colored_label(plane(&ctx), "Output Rv");
            ui.separator();
            ui.colored_label(rotor(&ctx), "Half-angle");
            ui.separator();
            ui.colored_label(line(&ctx), "Full angle");
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(EUCLIDEAN2_EDUCATION)
    }
}

/// Educational content for the 2D Euclidean rotor visualization.
const EUCLIDEAN2_EDUCATION: EducationalContent = EducationalContent {
    title: "2D Rotors in Geometric Algebra",

    overview: "\
This visualization demonstrates how rotations (and dilations!) work in 2D \
Geometric Algebra (GA). Unlike rotation matrices, GA uses 'rotors' - elements \
of the even subalgebra that naturally encode rotations and compose beautifully.

The key insight is that rotors use HALF the rotation angle internally. This \
enables perfect rotation composition: to combine rotations, multiply rotors.

This demo also shows what happens with non-unit rotors: they combine rotation \
with uniform scaling (dilation). Try the dilation slider to see!",

    math_background: "\
A UNIT rotor R (|R| = 1) that rotates by angle theta is:

    R = cos(theta/2) - sin(theta/2)e_1_2

where e_1_2 is the unit bivector (oriented xy-plane).

The sandwich product v' = RvRrev transforms vectors:
  - Rrev is the reverse (same scalar, negated bivector)
  - Unit rotors (|R| = 1) preserve vector length
  - Non-unit rotors scale vectors by |R|^2

SCALED ROTORS combine rotation with dilation:

    S = k*R  where k = sqrt(dilation factor)

    v' = SvSrev rotates AND scales by |S|^2 = k^2

So |R| = 2 means vectors are scaled by 4x, and |R| = 0.5 \
means vectors are scaled by 0.25x.

Composition still works: S_2S_1 gives combined rotation \
AND combined scaling!",

    how_to_use: "\
- Drag the angle slider to rotate the vector field
- Adjust the DILATION slider to see scaling effects
- Click 'Play' to animate continuous rotation
- The purple arc shows the HALF-angle (stored in rotor)
- The coral arc shows the FULL rotation angle
- When dilating, two circles show input vs output magnitude
- Watch |R| change as you adjust dilation",

    key_concepts: "\
- Rotors encode rotation using HALF the angle
- UNIT rotors (|R| = 1) perform PURE rotation
- NON-UNIT rotors combine rotation with scaling by |R|^2
- Sandwich product: v' = RvRrev
- Composition: R_total = R_2R_1 (angles add, scales multiply)
- The bivector e_1_2 represents the rotation plane
- This extends naturally to 3D and higher dimensions",

    resources: &[
        (
            "Rigid Geometric Algebra Wiki",
            "https://rigidgeometricalgebra.org/wiki/",
        ),
        (
            "Look, Ma, No Matrices!",
            "https://enkimute.github.io/LookMaNoMatrices/",
        ),
    ],
};
