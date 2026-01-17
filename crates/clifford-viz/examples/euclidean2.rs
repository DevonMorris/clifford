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

    fn render(&self, ui: &mut egui::Ui) {
        // Create the rotor from the current angle
        let rotor = Rotor::<f64>::from_angle(f64::from(self.angle));

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

                    // Circle showing the rotation path
                    let input_len = (input.x() * input.x() + input.y() * input.y()).sqrt();
                    plot_ui.line(
                        circle_2d(0.0, 0.0, input_len, with_alpha(palette::GRID, 80), 64)
                            .name("Rotation path"),
                    );

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

        // Animation controls
        ui.add_space(4.0);
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // Rotor display
        section_separator(ui, Some("Rotor Components"));
        let rotor = Rotor::<f64>::from_angle(f64::from(self.angle));
        ga_value_display(
            ui,
            "R",
            &[
                ("1", rotor.s() as f32),
                ("e\u{2081}\u{2082}", rotor.b() as f32),
            ],
        );

        // Formula box
        info_box(
            ui,
            &format!(
                "R = cos(\u{03b8}/2) + sin(\u{03b8}/2)e\u{2081}\u{2082}\n  = cos({:.1}\u{00b0}) + sin({:.1}\u{00b0})e\u{2081}\u{2082}\n  = {:.4} + {:.4}e\u{2081}\u{2082}",
                self.angle.to_degrees() / 2.0,
                self.angle.to_degrees() / 2.0,
                rotor.s(),
                rotor.b()
            ),
        );

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
            point2_display(ui, "v", input.x() as f32, input.y() as f32);
            point2_display(ui, "Rv", output.x() as f32, output.y() as f32);
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
            ui.label("2D Euclidean rotation using geometric algebra rotors.");
            ui.separator();
            ui.label("\u{1f535} Input vector v");
            ui.separator();
            ui.label("\u{1f7e2} Rotated vector Rv");
            ui.separator();
            ui.colored_label(palette::ROTOR, "\u{1f7e3} Bivector (half-angle)");
            ui.separator();
            ui.colored_label(palette::LINE, "\u{1f534} Full rotation angle");
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
This visualization demonstrates how rotations work in 2D Geometric Algebra (GA). \
Unlike traditional rotation matrices or complex numbers, GA uses objects called \
'rotors' that naturally encode rotations and compose beautifully.

The key insight is that rotors use HALF the rotation angle internally. This \
might seem strange at first, but it's what makes rotation composition work \
perfectly: to combine two rotations, you simply multiply their rotors.",

    math_background: "\
A rotor R that rotates by angle \u{03b8} is defined as:

    R = cos(\u{03b8}/2) + sin(\u{03b8}/2)e\u{2081}\u{2082}

where e\u{2081}\u{2082} is the unit bivector (the oriented xy-plane).

To rotate a vector v, we use the 'sandwich product':

    v' = R v R\u{2020}

where R\u{2020} is the reverse of R (same scalar, negated bivector).

Why half-angles? Consider composing two rotations:

    R\u{2082}(R\u{2081} v R\u{2081}\u{2020})R\u{2082}\u{2020} = (R\u{2082}R\u{2081}) v (R\u{2082}R\u{2081})\u{2020}

The total rotor is simply R\u{2082}R\u{2081}, and the half-angles add up \
to give the correct total rotation!",

    how_to_use: "\
\u{2022} Drag the angle slider to rotate the vector field
\u{2022} Click 'Play' to animate continuous rotation
\u{2022} The purple arc shows the HALF-angle (what's stored in the rotor)
\u{2022} The red arc shows the FULL rotation angle
\u{2022} Adjust the highlighted vector to see specific input/output pairs
\u{2022} Toggle display options to focus on different aspects",

    key_concepts: "\
\u{2022} Rotors encode rotations using HALF the angle
\u{2022} The sandwich product v' = RvR\u{2020} applies the rotation
\u{2022} Rotations compose by multiplying rotors: R_total = R\u{2082}R\u{2081}
\u{2022} Rotors are always unit magnitude: |R| = 1
\u{2022} The bivector e\u{2081}\u{2082} represents the plane of rotation
\u{2022} This naturally extends to 3D (and higher dimensions!)",

    resources: &[
        ("Rigid Geometric Algebra Wiki", "https://rigidgeometricalgebra.org/wiki/"),
        ("Look, Ma, No Matrices!", "https://enkimute.github.io/LookMaNoMatrices/"),
    ],
};

fn main() -> eframe::Result<()> {
    run_app::<Euclidean2Demo>()
}
