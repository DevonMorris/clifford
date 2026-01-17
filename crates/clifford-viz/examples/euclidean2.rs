//! 2D Euclidean Rotor Visualization
//!
//! This demo demonstrates rotation in 2D Euclidean space using rotors from
//! Geometric Algebra. It visualizes:
//!
//! - A vector field that rotates uniformly under rotor action
//! - The bivector (rotation plane) as an oriented arc
//! - Rotor components in the `cos(θ/2) + sin(θ/2)e₁₂` form
//!
//! Run with: `cargo run -p clifford-viz --example euclidean2 --release`
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

use clifford::ops::Transform;
use clifford::specialized::euclidean::dim2::{Rotor, Vector};
use clifford_viz::common::prelude::*;
use egui_plot::{Plot, Points};

/// Demo state for 2D Euclidean rotor visualization.
struct Euclidean2Demo {
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

/// Create a scaled rotor that both rotates and dilates.
///
/// A unit rotor R with |R| = 1 preserves vector length.
/// A scaled rotor S = kR with |S| = k scales vectors by k² via sandwich product.
/// So for dilation factor d, we need |R| = sqrt(d).
fn make_dilating_rotor(angle: f64, dilation: f64) -> Rotor<f64> {
    let half = angle / 2.0;
    let scale = dilation.sqrt();
    // R = scale * (cos(θ/2) - sin(θ/2)e₁₂)
    // Note: clifford library uses negated bivector convention
    Rotor::new_unchecked(scale * half.cos(), -scale * half.sin())
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

    fn render(&self, ui: &mut egui::Ui) {
        // Create the rotor from the current angle and dilation
        let rotor = make_dilating_rotor(f64::from(self.angle), f64::from(self.dilation));

        // Calculate bounds for the plot
        let bounds = (self.grid_size as f32 * self.vector_spacing + 1.5) as f64;

        Plot::new("euclidean2_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .show(ui, |plot_ui| {
                // Draw coordinate grid
                if self.show_grid {
                    for line in grid_2d(bounds, f64::from(self.vector_spacing)) {
                        plot_ui.line(line);
                    }
                    for axis in axes_2d(bounds) {
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
                        palette::ROTOR,
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
                        with_alpha(palette::LINE, 150),
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
                            let rotated = rotor.transform(&direction);

                            // Scale for visibility
                            let arrow_scale = 0.4;
                            let dx = rotated.x() * arrow_scale;
                            let dy = rotated.y() * arrow_scale;

                            // Draw from grid position
                            for arrow in arrow_2d(
                                base_x - dx / 2.0,
                                base_y - dy / 2.0,
                                dx,
                                dy,
                                with_alpha(palette::POINT, 180),
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
                    let output = rotor.transform(&input);

                    // Circle showing the rotation path (input magnitude)
                    let input_len = (input.x() * input.x() + input.y() * input.y()).sqrt();
                    plot_ui.line(
                        circle_2d(0.0, 0.0, input_len, with_alpha(palette::GRID, 80), 64)
                            .name("Input magnitude"),
                    );

                    // If dilating, show the output magnitude circle too
                    if (self.dilation - 1.0).abs() > 0.01 {
                        let output_len = (output.x() * output.x() + output.y() * output.y()).sqrt();
                        plot_ui.line(
                            circle_2d(0.0, 0.0, output_len, with_alpha(palette::PLANE, 60), 64)
                                .name("Output magnitude"),
                        );
                    }

                    // Input vector (blue, dashed-style via thinner line)
                    for arrow in labeled_arrow(0.0, 0.0, input.x(), input.y(), palette::POINT, "v")
                    {
                        plot_ui.line(arrow);
                    }

                    // Output vector (green)
                    for arrow in
                        labeled_arrow(0.0, 0.0, output.x(), output.y(), palette::PLANE, "Rv")
                    {
                        plot_ui.line(arrow);
                    }

                    // Point markers at tips
                    plot_ui.points(
                        Points::new(vec![[input.x(), input.y()]])
                            .color(palette::POINT)
                            .radius(6.0)
                            .filled(true)
                            .name("Input v"),
                    );
                    plot_ui.points(
                        Points::new(vec![[output.x(), output.y()]])
                            .color(palette::PLANE)
                            .radius(6.0)
                            .filled(true)
                            .name("Rotated Rv"),
                    );
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // Rotation angle
        section_separator(ui, Some("Rotation"));
        angle_slider_range(ui, "Angle \u{03b8}", &mut self.angle, -360.0, 360.0);

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
        ui.add_space(4.0);
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // Rotor display
        section_separator(ui, Some("Rotor Components"));
        let rotor = make_dilating_rotor(f64::from(self.angle), f64::from(self.dilation));
        let rotor_magnitude = (rotor.s() * rotor.s() + rotor.b() * rotor.b()).sqrt();
        ga_value_display(
            ui,
            "R",
            &[
                ("1", rotor.s() as f32),
                ("e\u{2081}\u{2082}", rotor.b() as f32),
            ],
        );
        value_display(ui, "|R|", rotor_magnitude as f32, 4);

        // Formula box
        if (self.dilation - 1.0).abs() < 0.01 {
            info_box(
                ui,
                &format!(
                    "R = cos(\u{03b8}/2) - sin(\u{03b8}/2)e\u{2081}\u{2082}\n  = {:.4} - {:.4}e\u{2081}\u{2082}\n|R| = 1 (unit rotor, pure rotation)",
                    (self.angle / 2.0).cos(),
                    (self.angle / 2.0).sin()
                ),
            );
        } else {
            info_box(
                ui,
                &format!(
                    "R = k(cos(\u{03b8}/2) - sin(\u{03b8}/2)e\u{2081}\u{2082})\nk = \u{221a}{:.2} = {:.4}\n|R| = {:.4} (scales by |R|\u{00b2} = {:.2})",
                    self.dilation,
                    self.dilation.sqrt(),
                    rotor_magnitude,
                    self.dilation
                ),
            );
        }

        // Highlight vector input
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
                    ui.label("\u{2192}");
                    ui.label(format!("|Rv| = {:.3}", output_mag));
                    ui.label(format!("(\u{00d7}{:.2})", output_mag / input_mag));
                });
            }
        }

        // Display options
        section_separator(ui, Some("Display Options"));
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_vector_field, "Show vector field");
        ui.checkbox(&mut self.show_bivector, "Show bivector (\u{03b8}/2 arc)");
        ui.checkbox(&mut self.show_rotation_arc, "Show rotation arc (\u{03b8})");

        // Vector field options
        if self.show_vector_field {
            ui.horizontal(|ui| {
                ui.label("Grid size:");
                ui.add(egui::Slider::new(&mut self.grid_size, 1..=5));
            });
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if (self.dilation - 1.0).abs() < 0.01 {
                ui.label("2D rotation using unit rotor (pure rotation, no scaling)");
            } else {
                ui.label(format!(
                    "2D rotation + dilation (|R|={:.2}, scales by {:.2}\u{00d7})",
                    self.dilation.sqrt(),
                    self.dilation
                ));
            }
            ui.separator();
            ui.label("\u{1f535} Input v");
            ui.separator();
            ui.label("\u{1f7e2} Output Rv");
            ui.separator();
            ui.colored_label(palette::ROTOR, "\u{1f7e3} Half-angle");
            ui.separator();
            ui.colored_label(palette::LINE, "\u{1f534} Full angle");
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
A UNIT rotor R (|R| = 1) that rotates by angle \u{03b8} is:

    R = cos(\u{03b8}/2) - sin(\u{03b8}/2)e\u{2081}\u{2082}

where e\u{2081}\u{2082} is the unit bivector (oriented xy-plane).

The sandwich product v' = RvR\u{2020} transforms vectors:
  \u{2022} R\u{2020} is the reverse (same scalar, negated bivector)
  \u{2022} Unit rotors (|R| = 1) preserve vector length
  \u{2022} Non-unit rotors scale vectors by |R|\u{00b2}

SCALED ROTORS combine rotation with dilation:

    S = k\u{00b7}R  where k = \u{221a}(dilation factor)

    v' = SvS\u{2020} rotates AND scales by |S|\u{00b2} = k\u{00b2}

So |R| = 2 means vectors are scaled by 4\u{00d7}, and |R| = 0.5 \
means vectors are scaled by 0.25\u{00d7}.

Composition still works: S\u{2082}S\u{2081} gives combined rotation \
AND combined scaling!",

    how_to_use: "\
\u{2022} Drag the angle slider to rotate the vector field
\u{2022} Adjust the DILATION slider to see scaling effects
\u{2022} Click 'Play' to animate continuous rotation
\u{2022} The purple arc shows the HALF-angle (stored in rotor)
\u{2022} The red arc shows the FULL rotation angle
\u{2022} When dilating, two circles show input vs output magnitude
\u{2022} Watch |R| change as you adjust dilation",

    key_concepts: "\
\u{2022} Rotors encode rotation using HALF the angle
\u{2022} UNIT rotors (|R| = 1) perform PURE rotation
\u{2022} NON-UNIT rotors combine rotation with scaling by |R|\u{00b2}
\u{2022} Sandwich product: v' = RvR\u{2020}
\u{2022} Composition: R_total = R\u{2082}R\u{2081} (angles add, scales multiply)
\u{2022} The bivector e\u{2081}\u{2082} represents the rotation plane
\u{2022} This extends naturally to 3D and higher dimensions",

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

fn main() -> eframe::Result<()> {
    run_app::<Euclidean2Demo>()
}
