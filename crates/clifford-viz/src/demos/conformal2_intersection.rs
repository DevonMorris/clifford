//! Circle-Circle Intersection - 2D Conformal GA Visualization
//!
//! This demo demonstrates the meet operation in CGA: computing the intersection
//! of two circles via the antiwedge product.
//!
//! In CGA, the intersection of two circles is computed as:
//! ```text
//! PointPair = Circle_1 ∨ Circle_2
//! ```
//!
//! The result is a point pair which can represent:
//! - Two distinct points (circles intersect at two points)
//! - A single point (circles are tangent)
//! - No real intersection (circles don't touch)

use crate::common::prelude::*;
use clifford::ops::Antiwedge;
use clifford::specialized::conformal::dim2::{Circle, PointPair};
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
const VIEWPORT_BOUNDS: f64 = 6.0;

/// Number of segments for circle rendering.
const CIRCLE_SEGMENTS: usize = 128;

/// Epsilon for numerical comparisons.
const EPSILON: f64 = 1e-6;

/// Tolerance for tangent detection.
const TANGENT_TOLERANCE: f64 = 0.1;

/// A draggable circle in the scene.
#[derive(Clone)]
struct DraggableCircle {
    /// Center x coordinate.
    cx: f64,
    /// Center y coordinate.
    cy: f64,
    /// Radius.
    radius: f64,
}

impl DraggableCircle {
    /// Creates a new draggable circle.
    fn new(cx: f64, cy: f64, radius: f64) -> Self {
        Self { cx, cy, radius }
    }

    /// Returns the CGA Circle representation.
    fn to_cga(&self) -> Circle<f64> {
        Circle::from_center_radius(self.cx, self.cy, self.radius)
    }
}

/// State of circle-circle intersection.
#[derive(Clone, Copy, PartialEq)]
enum IntersectionState {
    /// Two distinct intersection points.
    TwoPoints,
    /// Circles are tangent (one point).
    Tangent,
    /// No intersection (circles don't touch).
    NoIntersection,
    /// Circles are concentric (same center).
    Concentric,
}

impl IntersectionState {
    /// Returns a human-readable description of the intersection state.
    fn description(&self) -> &'static str {
        match self {
            Self::TwoPoints => "Two intersection points",
            Self::Tangent => "Tangent (one point)",
            Self::NoIntersection => "No intersection",
            Self::Concentric => "Concentric circles",
        }
    }
}

/// Demo state for Circle-Circle Intersection visualization.
pub struct Conformal2IntersectionDemo {
    /// First circle (blue).
    circle_a: DraggableCircle,
    /// Second circle (red).
    circle_b: DraggableCircle,

    /// Which circle is being dragged (if any).
    dragging: Option<usize>,

    /// Whether to show the grid.
    show_grid: bool,
    /// Whether to show intersection points.
    show_intersection: bool,
    /// Whether to show the algebraic details.
    show_algebra: bool,
}

impl Default for Conformal2IntersectionDemo {
    fn default() -> Self {
        Self {
            circle_a: DraggableCircle::new(-0.5, 0.0, 2.0),
            circle_b: DraggableCircle::new(2.0, 0.0, 1.5),
            dragging: None,
            show_grid: true,
            show_intersection: true,
            show_algebra: false,
        }
    }
}

impl Conformal2IntersectionDemo {
    /// Computes the intersection using CGA antiwedge and returns the PointPair.
    fn compute_intersection(&self) -> PointPair<f64> {
        let cga_a = self.circle_a.to_cga();
        let cga_b = self.circle_b.to_cga();
        cga_a.antiwedge(&cga_b)
    }

    /// Determines the intersection state from circle geometry.
    fn intersection_state(&self) -> IntersectionState {
        let dx = self.circle_b.cx - self.circle_a.cx;
        let dy = self.circle_b.cy - self.circle_a.cy;
        let d = (dx * dx + dy * dy).sqrt();
        let r1 = self.circle_a.radius;
        let r2 = self.circle_b.radius;

        // Concentric check
        if d < EPSILON {
            return IntersectionState::Concentric;
        }

        // Too far apart
        if d > r1 + r2 + TANGENT_TOLERANCE {
            return IntersectionState::NoIntersection;
        }

        // One inside the other
        if d < (r1 - r2).abs() - TANGENT_TOLERANCE {
            return IntersectionState::NoIntersection;
        }

        // Tangent (external or internal)
        if (d - (r1 + r2)).abs() < TANGENT_TOLERANCE
            || (d - (r1 - r2).abs()).abs() < TANGENT_TOLERANCE
        {
            return IntersectionState::Tangent;
        }

        IntersectionState::TwoPoints
    }

    /// Computes intersection points geometrically.
    fn compute_intersection_points(&self) -> Option<((f64, f64), (f64, f64))> {
        let c1x = self.circle_a.cx;
        let c1y = self.circle_a.cy;
        let r1 = self.circle_a.radius;
        let c2x = self.circle_b.cx;
        let c2y = self.circle_b.cy;
        let r2 = self.circle_b.radius;

        let dx = c2x - c1x;
        let dy = c2y - c1y;
        let d = (dx * dx + dy * dy).sqrt();

        if d < EPSILON {
            return None; // Concentric
        }

        if d > r1 + r2 + EPSILON || d < (r1 - r2).abs() - EPSILON {
            return None; // No intersection
        }

        // Distance from c1 to chord midpoint
        let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);

        // Height from chord midpoint to intersection
        let h_sq = r1 * r1 - a * a;

        if h_sq < -EPSILON {
            return None;
        }

        let h = if h_sq < 0.0 { 0.0 } else { h_sq.sqrt() };

        // Chord midpoint
        let mx = c1x + (a / d) * dx;
        let my = c1y + (a / d) * dy;

        // Perpendicular direction
        let perpx = -dy / d;
        let perpy = dx / d;

        let p1 = (mx + h * perpx, my + h * perpy);
        let p2 = (mx - h * perpx, my - h * perpy);

        Some((p1, p2))
    }

    /// Find which circle center is near the mouse position.
    fn find_nearest_circle(&self, mouse_x: f64, mouse_y: f64) -> Option<usize> {
        let threshold = 0.5;

        let dist_a =
            ((self.circle_a.cx - mouse_x).powi(2) + (self.circle_a.cy - mouse_y).powi(2)).sqrt();
        let dist_b =
            ((self.circle_b.cx - mouse_x).powi(2) + (self.circle_b.cy - mouse_y).powi(2)).sqrt();

        if dist_a < threshold && dist_a < dist_b {
            Some(0)
        } else if dist_b < threshold {
            Some(1)
        } else {
            None
        }
    }
}

impl VisualizationApp for Conformal2IntersectionDemo {
    fn name(&self) -> &'static str {
        "Conformal 2D - Circle Intersection"
    }

    fn update(&mut self, _dt: f32) {
        // No animation
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let state = self.intersection_state();

        let response = Plot::new("conformal2_intersection_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .auto_bounds(egui::Vec2b::new(false, false))
            .allow_zoom(false)
            .allow_drag(false)
            .allow_boxed_zoom(false)
            .allow_scroll(false)
            .include_x(-VIEWPORT_BOUNDS)
            .include_x(VIEWPORT_BOUNDS)
            .include_y(-VIEWPORT_BOUNDS)
            .include_y(VIEWPORT_BOUNDS)
            .show(ui, |plot_ui| {
                // Draw grid
                if self.show_grid {
                    for line in grid_2d(&ctx, VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(line);
                    }
                    for axis in axes_2d(&ctx, VIEWPORT_BOUNDS) {
                        plot_ui.line(axis);
                    }
                }

                // Draw Circle A (blue)
                let circle_a_line = circle_2d(
                    self.circle_a.cx,
                    self.circle_a.cy,
                    self.circle_a.radius,
                    circle(&ctx),
                    CIRCLE_SEGMENTS,
                )
                .name("Circle A");
                plot_ui.line(circle_a_line);

                // Draw Circle A center
                plot_ui.points(
                    Points::new(vec![[self.circle_a.cx, self.circle_a.cy]])
                        .color(circle(&ctx))
                        .radius(5.0)
                        .filled(true)
                        .name("Center A"),
                );

                // Draw Circle B (red/active color)
                let circle_b_line = circle_2d(
                    self.circle_b.cx,
                    self.circle_b.cy,
                    self.circle_b.radius,
                    active(&ctx),
                    CIRCLE_SEGMENTS,
                )
                .name("Circle B");
                plot_ui.line(circle_b_line);

                // Draw Circle B center
                plot_ui.points(
                    Points::new(vec![[self.circle_b.cx, self.circle_b.cy]])
                        .color(active(&ctx))
                        .radius(5.0)
                        .filled(true)
                        .name("Center B"),
                );

                // Draw intersection points
                if self.show_intersection {
                    if let Some(((p1x, p1y), (p2x, p2y))) = self.compute_intersection_points() {
                        match state {
                            IntersectionState::TwoPoints => {
                                // Two distinct points
                                plot_ui.points(
                                    Points::new(vec![[p1x, p1y], [p2x, p2y]])
                                        .color(selected(&ctx))
                                        .radius(8.0)
                                        .filled(true)
                                        .name("Intersection"),
                                );
                            }
                            IntersectionState::Tangent => {
                                // Single tangent point (both points are the same)
                                plot_ui.points(
                                    Points::new(vec![[p1x, p1y]])
                                        .color(selected(&ctx))
                                        .radius(8.0)
                                        .filled(true)
                                        .name("Tangent Point"),
                                );
                            }
                            _ => {}
                        }
                    }
                }
            });

        // Handle mouse interactions
        if let Some(pos) = response.response.interact_pointer_pos() {
            let plot_pos = response.transform.value_from_position(pos);
            let mouse_x = plot_pos.x;
            let mouse_y = plot_pos.y;

            if response.response.drag_started() {
                self.dragging = self.find_nearest_circle(mouse_x, mouse_y);
            }

            if response.response.dragged() {
                match self.dragging {
                    Some(0) => {
                        self.circle_a.cx = mouse_x;
                        self.circle_a.cy = mouse_y;
                    }
                    Some(1) => {
                        self.circle_b.cx = mouse_x;
                        self.circle_b.cy = mouse_y;
                    }
                    _ => {}
                }
            }

            if response.response.drag_stopped() {
                self.dragging = None;
            }
        } else {
            self.dragging = None;
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let state = self.intersection_state();
        let point_pair = self.compute_intersection();

        // Circle A controls
        group_header(ui, "Circle A (drag center)");
        ui.horizontal(|ui| {
            ui.label("Center:");
            let mut cx = self.circle_a.cx as f32;
            let mut cy = self.circle_a.cy as f32;
            if ui
                .add(egui::DragValue::new(&mut cx).speed(0.1).prefix("x: "))
                .changed()
            {
                self.circle_a.cx = f64::from(cx);
            }
            if ui
                .add(egui::DragValue::new(&mut cy).speed(0.1).prefix("y: "))
                .changed()
            {
                self.circle_a.cy = f64::from(cy);
            }
        });
        ui.horizontal(|ui| {
            ui.label("Radius:");
            let mut r = self.circle_a.radius as f32;
            if ui.add(egui::Slider::new(&mut r, 0.5..=4.0)).changed() {
                self.circle_a.radius = f64::from(r);
            }
        });

        // Circle B controls
        section_separator(ui, Some("Circle B (drag center)"));
        ui.horizontal(|ui| {
            ui.label("Center:");
            let mut cx = self.circle_b.cx as f32;
            let mut cy = self.circle_b.cy as f32;
            if ui
                .add(egui::DragValue::new(&mut cx).speed(0.1).prefix("x: "))
                .changed()
            {
                self.circle_b.cx = f64::from(cx);
            }
            if ui
                .add(egui::DragValue::new(&mut cy).speed(0.1).prefix("y: "))
                .changed()
            {
                self.circle_b.cy = f64::from(cy);
            }
        });
        ui.horizontal(|ui| {
            ui.label("Radius:");
            let mut r = self.circle_b.radius as f32;
            if ui.add(egui::Slider::new(&mut r, 0.5..=4.0)).changed() {
                self.circle_b.radius = f64::from(r);
            }
        });

        // Intersection result
        section_separator(ui, Some("Intersection (A v B)"));
        ui.colored_label(selected(&ctx), state.description());

        if let Some(((p1x, p1y), (p2x, p2y))) = self.compute_intersection_points() {
            match state {
                IntersectionState::TwoPoints => {
                    point2_display(ui, "P_1", p1x as f32, p1y as f32);
                    point2_display(ui, "P_2", p2x as f32, p2y as f32);
                }
                IntersectionState::Tangent => {
                    point2_display(ui, "P", p1x as f32, p1y as f32);
                }
                _ => {}
            }
        }

        // Algebraic details
        section_separator(ui, Some("Algebra"));
        ui.checkbox(&mut self.show_algebra, "Show PointPair components");

        if self.show_algebra {
            ui.add_space(spacing::XS);
            ui.label("Meet: PointPair = A v B");
            ui.add_space(spacing::XS);

            ga_value_display(
                ui,
                "PP",
                &[
                    ("m", point_pair.m() as f32),
                    ("e1ep", point_pair.e1ep() as f32),
                    ("e2ep", point_pair.e2ep() as f32),
                    ("e1em", point_pair.e1em() as f32),
                    ("e2em", point_pair.e2em() as f32),
                    ("epem", point_pair.epem() as f32),
                ],
            );

            ui.add_space(spacing::XS);
            value_display(ui, "norm^2", point_pair.norm_squared() as f32, 4);

            ui.add_space(spacing::XS);
            info_box(
                ui,
                "norm^2 > 0: two real points\n\
                 norm^2 = 0: tangent (one point)\n\
                 norm^2 < 0: no real intersection",
            );
        }

        // Display options
        section_separator(ui, Some("Display"));
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_intersection, "Show intersection points");

        // Presets
        section_separator(ui, Some("Presets"));

        if ui.button("Two intersection points").clicked() {
            self.circle_a = DraggableCircle::new(-0.5, 0.0, 2.0);
            self.circle_b = DraggableCircle::new(2.0, 0.0, 1.5);
        }

        if ui.button("Tangent (external)").clicked() {
            self.circle_a = DraggableCircle::new(-1.5, 0.0, 1.5);
            self.circle_b = DraggableCircle::new(1.5, 0.0, 1.5);
        }

        if ui.button("Tangent (internal)").clicked() {
            self.circle_a = DraggableCircle::new(0.0, 0.0, 3.0);
            self.circle_b = DraggableCircle::new(1.5, 0.0, 1.5);
        }

        if ui.button("No intersection").clicked() {
            self.circle_a = DraggableCircle::new(-2.0, 0.0, 1.0);
            self.circle_b = DraggableCircle::new(2.0, 0.0, 1.0);
        }

        if ui.button("Concentric").clicked() {
            self.circle_a = DraggableCircle::new(0.0, 0.0, 2.0);
            self.circle_b = DraggableCircle::new(0.0, 0.0, 1.0);
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let state = self.intersection_state();
        let point_pair = self.compute_intersection();

        ui.horizontal(|ui| {
            ui.colored_label(selected(&ctx), state.description());
            ui.separator();
            ui.label(format!("norm^2 = {:.4}", point_pair.norm_squared()));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(CONFORMAL2_INTERSECTION_EDUCATION)
    }
}

/// Educational content for the Circle-Circle Intersection visualization.
const CONFORMAL2_INTERSECTION_EDUCATION: EducationalContent = EducationalContent {
    title: "Circle-Circle Intersection in CGA",

    overview: "\
In Conformal Geometric Algebra, the intersection of two circles is computed using \
the meet (antiwedge) operation:

    PointPair = Circle_A v Circle_B

The result is a PointPair (grade-2 element) that encodes the intersection. \
The squared norm of this PointPair tells us the nature of the intersection.",

    math_background: "\
MEET OPERATION:
The meet (antiwedge, v) computes the common subspace of two geometric objects.
For two circles, this gives their intersection points.

POINT PAIR INTERPRETATION:
The PointPair PP = A v B has a squared norm that indicates:
- norm^2 > 0: Two distinct real intersection points
- norm^2 = 0: Single tangent point (circles touch)
- norm^2 < 0: No real intersection (imaginary points)

GEOMETRIC CASES:
- External tangent: circles touch from outside
- Internal tangent: one circle touches inside the other
- Secant: circles cross at two points
- Concentric: circles share center (no finite meet)
- Separated: circles don't touch (imaginary intersection)",

    how_to_use: "\
- DRAG circle centers to move them
- Use RADIUS SLIDERS to resize circles
- Watch the intersection points update in real-time
- Enable ALGEBRA view to see the PointPair components
- Try PRESETS for common configurations
- Observe how norm^2 changes sign at tangent configurations",

    key_concepts: "\
- Meet (v) computes intersection of geometric objects
- PointPair represents two points algebraically
- norm^2 sign determines if intersection is real
- CGA handles all cases uniformly (secant, tangent, separated)
- Tangent is the boundary between real and imaginary solutions",

    resources: &[
        (
            "Conformal Geometric Algebra Wiki",
            "https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page",
        ),
        (
            "Circle-Circle Intersection (Wikipedia)",
            "https://en.wikipedia.org/wiki/Circle%E2%80%93circle_intersection",
        ),
    ],
};
