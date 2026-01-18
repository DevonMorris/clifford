//! Circle from Three Points - 2D Conformal GA Visualization
//!
//! This demo demonstrates the fundamental CGA operation: creating a circle
//! from three points via the outer product (wedge).
//!
//! In CGA, three points uniquely determine a circle:
//! ```text
//! Circle = P₁ ∧ P₂ ∧ P₃
//! ```
//!
//! When the three points are collinear, the result is a line (circle through infinity).
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_circles --release`

use clifford::ops::Wedge;
use clifford::specialized::conformal::dim2::{Circle, RoundPoint};
use clifford_viz::common::prelude::*;
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
const VIEWPORT_BOUNDS: f64 = 5.0;

/// Epsilon for collinearity detection.
const COLLINEARITY_EPSILON: f64 = 1e-6;

/// Maximum radius before we treat the circle as a line.
/// This prevents viewport scaling issues as points approach collinearity.
const MAX_DRAWABLE_RADIUS: f64 = 100.0;

/// A draggable point in the visualization.
#[derive(Clone)]
struct DraggablePoint {
    /// The CGA round point.
    point: RoundPoint<f64>,
    /// Display name (e.g., "P₁", "P₂", "P₃").
    name: String,
}

impl DraggablePoint {
    /// Creates a new draggable point.
    fn new(x: f64, y: f64, name: &str) -> Self {
        Self {
            point: RoundPoint::from_euclidean(x, y),
            name: name.to_string(),
        }
    }

    /// Returns the Euclidean coordinates.
    fn euclidean(&self) -> Option<(f64, f64)> {
        self.point.to_euclidean()
    }

    /// Sets the position.
    fn set_position(&mut self, x: f64, y: f64) {
        self.point = RoundPoint::from_euclidean(x, y);
    }
}

/// Demo state for Circle from Three Points visualization.
struct Conformal2CirclesDemo {
    /// The three points that define the circle.
    points: Vec<DraggablePoint>,
    /// Index of point being dragged (if any).
    dragging_point: Option<usize>,
    /// Whether to show the grid.
    show_grid: bool,
    /// Whether to show point coordinates.
    show_coordinates: bool,
    /// Whether to show the circle's algebraic components.
    show_algebra: bool,
    /// Number of segments for circle rendering.
    circle_segments: usize,
}

impl Default for Conformal2CirclesDemo {
    fn default() -> Self {
        // Initial triangle of points that form a circle
        let points = vec![
            DraggablePoint::new(2.0, 0.0, "P\u{2081}"),
            DraggablePoint::new(0.0, 2.0, "P\u{2082}"),
            DraggablePoint::new(-2.0, 0.0, "P\u{2083}"),
        ];

        Self {
            points,
            dragging_point: None,
            show_grid: true,
            show_coordinates: true,
            show_algebra: false,
            circle_segments: 64,
        }
    }
}

impl Conformal2CirclesDemo {
    /// Computes the circle from the three points using CGA wedge product.
    fn compute_circle(&self) -> Circle<f64> {
        self.points[0]
            .point
            .wedge(&self.points[1].point)
            .wedge(&self.points[2].point)
    }

    /// Find the nearest point to mouse position within threshold.
    fn find_nearest_point(&self, mouse_x: f64, mouse_y: f64) -> Option<usize> {
        let threshold = 0.4; // Distance in plot units
        let mut nearest_idx = None;
        let mut nearest_dist = threshold;

        for (idx, point) in self.points.iter().enumerate() {
            if let Some((px, py)) = point.euclidean() {
                let dist = ((px - mouse_x).powi(2) + (py - mouse_y).powi(2)).sqrt();
                if dist < nearest_dist {
                    nearest_dist = dist;
                    nearest_idx = Some(idx);
                }
            }
        }
        nearest_idx
    }
}

impl VisualizationApp for Conformal2CirclesDemo {
    fn name(&self) -> &'static str {
        "Conformal 2D - Circle from Three Points"
    }

    fn update(&mut self, _dt: f32) {
        // No animation in this demo
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        // Compute the circle using CGA wedge product
        let circle = self.compute_circle();
        let is_line = circle.is_line(COLLINEARITY_EPSILON);

        let response = Plot::new("conformal2_circles_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .auto_bounds(egui::Vec2b::new(false, false)) // Disable auto-scaling
            .allow_zoom(false) // Fixed viewport - no zoom
            .allow_drag(false) // We handle drag ourselves
            .allow_boxed_zoom(false)
            .allow_scroll(false) // No scroll zoom
            .include_x(-VIEWPORT_BOUNDS) // Set viewport bounds
            .include_x(VIEWPORT_BOUNDS)
            .include_y(-VIEWPORT_BOUNDS)
            .include_y(VIEWPORT_BOUNDS)
            .show(ui, |plot_ui| {
                // Draw coordinate grid
                if self.show_grid {
                    for line in grid_2d(VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(line);
                    }
                    for axis in axes_2d(VIEWPORT_BOUNDS) {
                        plot_ui.line(axis);
                    }
                }

                // Determine if we should draw as line or circle
                // Treat as line if: collinear (w≈0) OR radius too large OR center too far
                let (draw_as_line, center_opt, radius_opt) = if is_line {
                    (true, None, None)
                } else if let (Some((cx, cy)), Some(r)) = (circle.center(), circle.radius()) {
                    // Check if circle is effectively a line (too large to draw sensibly)
                    let center_dist = (cx * cx + cy * cy).sqrt();
                    let too_large = r > MAX_DRAWABLE_RADIUS || center_dist > MAX_DRAWABLE_RADIUS;
                    if too_large {
                        (true, Some((cx, cy)), Some(r))
                    } else {
                        (false, Some((cx, cy)), Some(r))
                    }
                } else {
                    (true, None, None)
                };

                if draw_as_line {
                    // Draw as a line through the points
                    if let (Some((x1, y1)), Some((x2, y2))) =
                        (self.points[0].euclidean(), self.points[1].euclidean())
                    {
                        let dx = x2 - x1;
                        let dy = y2 - y1;
                        let len = (dx * dx + dy * dy).sqrt();
                        if len > 1e-10 {
                            let label = if is_line {
                                "Line (collinear)"
                            } else {
                                "Line (large radius)"
                            };
                            let plot_line = infinite_line_2d(
                                x1,
                                y1,
                                dx / len,
                                dy / len,
                                VIEWPORT_BOUNDS,
                                palette::CIRCLE,
                            )
                            .name(label);
                            plot_ui.line(plot_line);
                        }
                    }
                } else if let (Some((cx, cy)), Some(r)) = (center_opt, radius_opt) {
                    // Draw the circle
                    let circle_line =
                        circle_2d(cx, cy, r, palette::CIRCLE, self.circle_segments).name("Circle");
                    plot_ui.line(circle_line);

                    // Draw center marker (only if visible in viewport)
                    if cx.abs() < VIEWPORT_BOUNDS * 2.0 && cy.abs() < VIEWPORT_BOUNDS * 2.0 {
                        plot_ui.points(
                            Points::new(vec![[cx, cy]])
                                .color(with_alpha(palette::CIRCLE, 150))
                                .radius(4.0)
                                .filled(false)
                                .name("Center"),
                        );
                    }

                    // Draw radius line from center to first point (only if reasonably visible)
                    if let Some((px, py)) = self.points[0].euclidean() {
                        if cx.abs() < VIEWPORT_BOUNDS * 2.0 && cy.abs() < VIEWPORT_BOUNDS * 2.0 {
                            let radius_line =
                                line_segment(cx, cy, px, py, with_alpha(palette::CIRCLE, 100))
                                    .name("Radius");
                            plot_ui.line(radius_line);
                        }
                    }
                }

                // Draw the three defining points
                for point in &self.points {
                    if let Some((x, y)) = point.euclidean() {
                        plot_ui.points(
                            Points::new(vec![[x, y]])
                                .color(palette::POINT)
                                .radius(7.0)
                                .filled(true)
                                .name(&point.name),
                        );
                    }
                }
            });

        // Handle mouse interactions
        if let Some(pos) = response.response.interact_pointer_pos() {
            let plot_pos = response.transform.value_from_position(pos);
            let mouse_x = plot_pos.x;
            let mouse_y = plot_pos.y;

            // Handle drag start
            if response.response.drag_started() {
                self.dragging_point = self.find_nearest_point(mouse_x, mouse_y);
            }

            // Handle active dragging
            if response.response.dragged() {
                if let Some(idx) = self.dragging_point {
                    if idx < self.points.len() {
                        self.points[idx].set_position(mouse_x, mouse_y);
                    }
                }
            }

            // Handle drag release
            if response.response.drag_stopped() {
                self.dragging_point = None;
            }
        } else {
            self.dragging_point = None;
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // Compute circle for display
        let circle = self.compute_circle();
        let is_line = circle.is_line(COLLINEARITY_EPSILON);

        // === Point Coordinates ===
        group_header(ui, "Points (drag to move)");

        for point in &mut self.points {
            ui.horizontal(|ui| {
                ui.label(&point.name);
                if let Some((x, y)) = point.euclidean() {
                    let mut new_x = x as f32;
                    let mut new_y = y as f32;

                    ui.add(egui::DragValue::new(&mut new_x).speed(0.1).prefix("x: "));
                    ui.add(egui::DragValue::new(&mut new_y).speed(0.1).prefix("y: "));

                    if (f64::from(new_x) - x).abs() > 1e-6 || (f64::from(new_y) - y).abs() > 1e-6 {
                        point.set_position(f64::from(new_x), f64::from(new_y));
                    }
                }
            });
        }

        // === CGA Operation ===
        section_separator(ui, Some("CGA Operation"));
        ui.label("Circle = P\u{2081} \u{2227} P\u{2082} \u{2227} P\u{2083}");

        // === Result ===
        section_separator(ui, Some("Result"));

        if is_line {
            ui.colored_label(palette::CIRCLE, "Line (collinear points)");
            ui.label("w \u{2248} 0 \u{2192} circle through infinity");
            info_box(
                ui,
                "When three points are collinear,\ntheir wedge product has w = 0,\nrepresenting a line (infinite radius).",
            );
        } else if let (Some((cx, cy)), Some(r)) = (circle.center(), circle.radius()) {
            // Check if radius is too large to draw as circle
            let center_dist = (cx * cx + cy * cy).sqrt();
            let is_effectively_line = r > MAX_DRAWABLE_RADIUS || center_dist > MAX_DRAWABLE_RADIUS;

            if is_effectively_line {
                ui.colored_label(palette::CIRCLE, "Line (nearly collinear)");
                ui.label(format!("radius = {:.1} (too large to display)", r));
                info_box(
                    ui,
                    "As points approach collinearity,\nthe radius grows toward infinity.\nDrawn as a line for stability.",
                );
            } else {
                ui.colored_label(palette::CIRCLE, "Circle");
                ui.add_space(spacing::XS);
                point2_display(ui, "Center", cx as f32, cy as f32);
                value_display(ui, "Radius", r as f32, 4);
                if let Some(k) = circle.curvature() {
                    value_display(ui, "Curvature (1/r)", k as f32, 4);
                }
            }
        }

        // === Algebraic Components ===
        section_separator(ui, Some("Algebra"));
        ui.checkbox(&mut self.show_algebra, "Show circle components");

        if self.show_algebra {
            ui.add_space(spacing::XS);
            ui.label("Circle trivector components:");
            ga_value_display(
                ui,
                "C",
                &[
                    ("e\u{2081}\u{2082}\u{2083}", circle.w() as f32),
                    ("e\u{2081}\u{2082}\u{2084}", circle.cx() as f32),
                    ("e\u{2081}\u{2083}\u{2084}", circle.cy() as f32),
                    ("e\u{2082}\u{2083}\u{2084}", circle.r() as f32),
                ],
            );

            ui.add_space(spacing::XS);
            info_box(
                ui,
                "In our Cl(3,1) basis:\n\
                 \u{2022} center_x = e\u{2082}\u{2083}\u{2084} / e\u{2081}\u{2082}\u{2083}\n\
                 \u{2022} center_y = -e\u{2081}\u{2083}\u{2084} / e\u{2081}\u{2082}\u{2083}\n\
                 \u{2022} radius\u{00b2} = 2\u{00b7}e\u{2081}\u{2082}\u{2084}/w + cx\u{00b2} + cy\u{00b2}",
            );
        }

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_coordinates, "Show coordinates");

        ui.horizontal(|ui| {
            ui.label("Circle segments:");
            ui.add(egui::Slider::new(&mut self.circle_segments, 16..=128));
        });

        // === Presets ===
        section_separator(ui, Some("Presets"));

        if ui.button("Unit circle at origin").clicked() {
            self.points[0].set_position(1.0, 0.0);
            self.points[1].set_position(0.0, 1.0);
            self.points[2].set_position(-1.0, 0.0);
        }

        if ui.button("Large circle").clicked() {
            self.points[0].set_position(3.0, 0.0);
            self.points[1].set_position(0.0, 3.0);
            self.points[2].set_position(-3.0, 0.0);
        }

        if ui.button("Off-center circle").clicked() {
            self.points[0].set_position(3.0, 1.0);
            self.points[1].set_position(1.0, 3.0);
            self.points[2].set_position(-1.0, 1.0);
        }

        if ui.button("Make collinear (line)").clicked() {
            self.points[0].set_position(-2.0, 0.0);
            self.points[1].set_position(0.0, 0.0);
            self.points[2].set_position(2.0, 0.0);
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        let circle = self.compute_circle();
        let is_line = circle.is_line(COLLINEARITY_EPSILON);

        ui.horizontal(|ui| {
            if is_line {
                ui.colored_label(palette::CIRCLE, "Result: Line (collinear)");
            } else if let (Some((cx, cy)), Some(r)) = (circle.center(), circle.radius()) {
                let center_dist = (cx * cx + cy * cy).sqrt();
                if r > MAX_DRAWABLE_RADIUS || center_dist > MAX_DRAWABLE_RADIUS {
                    ui.colored_label(palette::CIRCLE, format!("Line (r={:.1})", r));
                } else {
                    ui.colored_label(
                        palette::CIRCLE,
                        format!("Circle: center=({:.2}, {:.2}) r={:.2}", cx, cy, r),
                    );
                }
            }
            ui.separator();
            ui.label(format!("w={:.4}", circle.w()));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(CONFORMAL2_CIRCLES_EDUCATION)
    }
}

/// Educational content for the Circle from Three Points visualization.
const CONFORMAL2_CIRCLES_EDUCATION: EducationalContent = EducationalContent {
    title: "Circle from Three Points in CGA",

    overview: "\
In Conformal Geometric Algebra (CGA), three points uniquely determine a circle \
through the wedge (outer) product:

    Circle = P\u{2081} \u{2227} P\u{2082} \u{2227} P\u{2083}

This is the conformal analogue of \"two points determine a line\" in projective \
geometry. The beauty of CGA is that circles and lines are treated uniformly - \
a line is simply a circle passing through the point at infinity.",

    math_background: "\
POINT EMBEDDING in CGA:
A 2D point (x, y) is embedded as a null vector:
    P = x\u{00b7}e\u{2081} + y\u{00b7}e\u{2082} + e\u{2083} + \u{00bd}(x\u{00b2}+y\u{00b2})\u{00b7}e\u{2084}

CIRCLE FROM THREE POINTS:
    C = P\u{2081} \u{2227} P\u{2082} \u{2227} P\u{2083}

The result is a trivector (grade-3 element) with components:
    C = w\u{00b7}e\u{2081}\u{2082}\u{2083} + cx\u{00b7}e\u{2081}\u{2082}\u{2084} + cy\u{00b7}e\u{2081}\u{2083}\u{2084} + r\u{00b7}e\u{2082}\u{2083}\u{2084}

EXTRACTING CENTER AND RADIUS:
In our Cl(3,1) orthonormal basis:
    center_x = r / w
    center_y = -cy / w
    radius\u{00b2} = 2\u{00b7}cx/w + center_x\u{00b2} + center_y\u{00b2}

COLLINEARITY (LINE DETECTION):
When w \u{2248} 0, the three points are collinear and the \"circle\" is \
actually a line (circle through infinity with infinite radius).",

    how_to_use: "\
\u{2022} DRAG POINTS to move them and see the circle update in real-time
\u{2022} Use PRESETS to quickly set up common configurations
\u{2022} Watch the CENTER and RADIUS update as you move points
\u{2022} Make points COLLINEAR to see the circle become a line
\u{2022} Enable SHOW ALGEBRA to see the raw trivector components",

    key_concepts: "\
\u{2022} Wedge product (\u{2227}) creates higher-grade elements
\u{2022} Three non-collinear points \u{2192} unique circle
\u{2022} Three collinear points \u{2192} line (w = 0)
\u{2022} Circle center and radius are extracted algebraically
\u{2022} CGA unifies circles and lines as the same type of object",

    resources: &[
        (
            "Conformal Geometric Algebra Wiki",
            "https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page",
        ),
        (
            "CGA Primer by Leo Dorst",
            "https://geometricalgebra.org/cga.html",
        ),
    ],
};

fn main() -> eframe::Result<()> {
    run_app::<Conformal2CirclesDemo>()
}
