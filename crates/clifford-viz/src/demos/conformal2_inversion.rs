//! Circle Inversion - 2D Conformal GA Visualization
//!
//! This demo demonstrates circle inversion (also known as geometric inversion),
//! a fundamental conformal transformation implemented using CGA.
//!
//! In circle inversion through a circle C with center O and radius r:
//! - A point P maps to P' such that |OP| × |OP'| = r² and P, O, P' are collinear
//! - Circles map to circles (or lines if they pass through the center)
//! - Lines map to circles (or lines if they pass through the center)
//!
//! In CGA, inversion through circle C is computed as:
//!     P' = C × P × C⁻¹
//! where × is the geometric product.

use crate::common::prelude::*;
use clifford::ops::Transform;
use clifford::specialized::conformal::dim2::{Circle, RoundPoint};
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
const VIEWPORT_BOUNDS: f64 = 6.0;

/// Epsilon for numerical comparisons.
const EPSILON: f64 = 1e-10;

/// Number of segments for circle rendering.
const CIRCLE_SEGMENTS: usize = 256;

/// A draggable point in the scene.
#[derive(Clone)]
struct ScenePoint {
    /// The CGA round point.
    point: RoundPoint<f64>,
    /// Display name.
    name: String,
}

impl ScenePoint {
    /// Creates a new scene point from Euclidean coordinates.
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

    /// Sets the position from Euclidean coordinates.
    fn set_position(&mut self, x: f64, y: f64) {
        self.point = RoundPoint::from_euclidean(x, y);
    }
}

/// A draggable circle in the scene.
#[derive(Clone)]
struct SceneCircle {
    /// The CGA circle.
    circle: Circle<f64>,
    /// Display name.
    name: String,
}

impl SceneCircle {
    /// Creates a new scene circle from center and radius.
    fn new(cx: f64, cy: f64, radius: f64, name: &str) -> Self {
        Self {
            circle: Circle::from_center_radius(cx, cy, radius),
            name: name.to_string(),
        }
    }

    /// Returns center and radius.
    fn center_radius(&self) -> Option<(f64, f64, f64)> {
        let (cx, cy) = self.circle.center()?;
        let r = self.circle.radius()?;
        Some((cx, cy, r))
    }

    /// Sets the circle from center and radius.
    fn set_center_radius(&mut self, cx: f64, cy: f64, radius: f64) {
        self.circle = Circle::from_center_radius(cx, cy, radius);
    }
}

/// Result of inverting a circle - can be a circle or a line.
#[derive(Clone)]
enum InvertedCircle {
    /// Circle with center (cx, cy) and radius.
    Circle {
        /// Center X.
        cx: f64,
        /// Center Y.
        cy: f64,
        /// Radius.
        radius: f64,
    },
    /// Line defined by nx*x + ny*y = d.
    Line {
        /// Normal X component.
        nx: f64,
        /// Normal Y component.
        ny: f64,
        /// Distance parameter.
        d: f64,
    },
}

/// Demo state for Circle Inversion visualization.
pub struct Conformal2InversionDemo {
    /// Inversion circle (CGA representation).
    inversion_circle: Circle<f64>,
    /// Inversion circle center x (for UI).
    inversion_cx: f64,
    /// Inversion circle center y (for UI).
    inversion_cy: f64,
    /// Inversion circle radius (for UI).
    inversion_radius: f64,

    /// Points in the scene.
    points: Vec<ScenePoint>,
    /// Circles in the scene.
    circles: Vec<SceneCircle>,

    /// Index of point being dragged (if any).
    dragging_point: Option<usize>,
    /// Whether dragging the inversion circle center.
    dragging_inversion_center: bool,
    /// Index of circle being dragged (if any).
    dragging_circle: Option<usize>,

    /// Whether to show the grid.
    show_grid: bool,
    /// Whether to show the inversion circle.
    show_inversion_circle: bool,
    /// Whether to show the original objects.
    show_originals: bool,
    /// Whether to show the inverted objects.
    show_inverted: bool,
    /// Whether to show transformation rules.
    show_rules: bool,
}

impl Default for Conformal2InversionDemo {
    fn default() -> Self {
        let inversion_cx = 0.0;
        let inversion_cy = 0.0;
        let inversion_radius = 2.0;

        // Initial configuration with some points and circles
        let points = vec![
            ScenePoint::new(2.5, 1.0, "P\u{2081}"),
            ScenePoint::new(3.0, -1.5, "P\u{2082}"),
            ScenePoint::new(-2.0, 2.0, "P\u{2083}"),
        ];

        let circles = vec![
            SceneCircle::new(-2.0, -1.0, 1.0, "C\u{2081}"),
            SceneCircle::new(3.5, 2.0, 0.8, "C\u{2082}"),
        ];

        Self {
            inversion_circle: Circle::from_center_radius(
                inversion_cx,
                inversion_cy,
                inversion_radius,
            ),
            inversion_cx,
            inversion_cy,
            inversion_radius,

            points,
            circles,

            dragging_point: None,
            dragging_inversion_center: false,
            dragging_circle: None,

            show_grid: true,
            show_inversion_circle: true,
            show_originals: true,
            show_inverted: true,
            show_rules: false,
        }
    }
}

impl Conformal2InversionDemo {
    /// Updates the CGA inversion circle from the UI parameters.
    fn update_inversion_circle(&mut self) {
        self.inversion_circle =
            Circle::from_center_radius(self.inversion_cx, self.inversion_cy, self.inversion_radius);
    }

    /// Transforms a point through the inversion circle using CGA.
    ///
    /// Using sandwich product (reflection): P' = C × P × C̃
    /// where × is the geometric product and C̃ is the reverse.
    fn invert_point(&self, point: &RoundPoint<f64>) -> Option<RoundPoint<f64>> {
        Some(self.inversion_circle.transform(point))
    }

    /// Transforms a circle through the inversion circle using CGA.
    ///
    /// Using sandwich product (reflection): C' = I × C × Ĩ
    /// where I is the inversion circle and × is the geometric product.
    ///
    /// Cases:
    /// - Circle NOT through center → Circle
    /// - Circle THROUGH center → Line (circle through infinity)
    fn invert_circle(&self, circle: &Circle<f64>) -> Option<InvertedCircle> {
        let result = self.inversion_circle.transform(circle);

        // Try to extract as a regular circle first
        if let (Some((cx, cy)), Some(radius)) = (result.center(), result.radius()) {
            Some(InvertedCircle::Circle { cx, cy, radius })
        } else {
            // center() or radius() returned None - this is a line
            // Extract line parameters directly from the circle components
            let nx = result.e12em();
            let ny = result.e1epem();
            let d = result.e2epem();
            Some(InvertedCircle::Line { nx, ny, d })
        }
    }

    /// Find the nearest point to mouse position within threshold.
    fn find_nearest_point(&self, mouse_x: f64, mouse_y: f64) -> Option<usize> {
        let threshold = 0.4;
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

    /// Find the nearest circle to mouse position within threshold.
    fn find_nearest_circle(&self, mouse_x: f64, mouse_y: f64) -> Option<usize> {
        let threshold = 0.4;
        let mut nearest_idx = None;
        let mut nearest_dist = threshold;

        for (idx, circle) in self.circles.iter().enumerate() {
            if let Some((cx, cy, r)) = circle.center_radius() {
                // Distance to circle edge
                let dist_to_center = ((cx - mouse_x).powi(2) + (cy - mouse_y).powi(2)).sqrt();
                let dist = (dist_to_center - r).abs();
                if dist < nearest_dist {
                    nearest_dist = dist;
                    nearest_idx = Some(idx);
                }
            }
        }
        nearest_idx
    }

    /// Check if mouse is near inversion circle center.
    fn is_near_inversion_center(&self, mouse_x: f64, mouse_y: f64) -> bool {
        let dist =
            ((self.inversion_cx - mouse_x).powi(2) + (self.inversion_cy - mouse_y).powi(2)).sqrt();
        dist < 0.4
    }
}

impl VisualizationApp for Conformal2InversionDemo {
    fn name(&self) -> &'static str {
        "Conformal 2D - Circle Inversion"
    }

    fn update(&mut self, _dt: f32) {
        // No animation in this demo
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let response = Plot::new("conformal2_inversion_plot")
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
                // Draw coordinate grid
                if self.show_grid {
                    for line in grid_2d(&ctx, VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(line);
                    }
                    for axis in axes_2d(&ctx, VIEWPORT_BOUNDS) {
                        plot_ui.line(axis);
                    }
                }

                // Draw the inversion circle
                if self.show_inversion_circle {
                    let inv_circle = circle_2d(
                        self.inversion_cx,
                        self.inversion_cy,
                        self.inversion_radius,
                        with_alpha(motor(&ctx), 200),
                        CIRCLE_SEGMENTS,
                    )
                    .name("Inversion Circle");
                    plot_ui.line(inv_circle);

                    // Draw center marker
                    plot_ui.points(
                        Points::new(vec![[self.inversion_cx, self.inversion_cy]])
                            .color(motor(&ctx))
                            .radius(6.0)
                            .filled(true)
                            .name("Inversion Center"),
                    );
                }

                // Draw original points and their inversions
                if self.show_originals {
                    for pt in &self.points {
                        if let Some((x, y)) = pt.euclidean() {
                            plot_ui.points(
                                Points::new(vec![[x, y]])
                                    .color(point(&ctx))
                                    .radius(7.0)
                                    .filled(true)
                                    .name(&pt.name),
                            );
                        }
                    }
                }

                if self.show_inverted {
                    for (idx, pt) in self.points.iter().enumerate() {
                        if let Some(inverted) = self.invert_point(&pt.point) {
                            if let Some((inv_x, inv_y)) = inverted.to_euclidean() {
                                // Only draw if within reasonable bounds
                                if inv_x.abs() < VIEWPORT_BOUNDS * 2.0
                                    && inv_y.abs() < VIEWPORT_BOUNDS * 2.0
                                {
                                    plot_ui.points(
                                        Points::new(vec![[inv_x, inv_y]])
                                            .color(active(&ctx))
                                            .radius(7.0)
                                            .filled(true)
                                            .name(format!("{}'", self.points[idx].name)),
                                    );

                                    // Draw connecting line from original to inverted
                                    if let Some((orig_x, orig_y)) = pt.euclidean() {
                                        let connect_line = line_segment(
                                            orig_x,
                                            orig_y,
                                            inv_x,
                                            inv_y,
                                            with_alpha(grid(&ctx), 80),
                                        )
                                        .name("Inversion Ray");
                                        plot_ui.line(connect_line);
                                    }
                                }
                            }
                        }
                    }
                }

                // Draw original circles and their inversions
                if self.show_originals {
                    for circ in &self.circles {
                        if let Some((cx, cy, r)) = circ.center_radius() {
                            let circle_line = circle_2d(cx, cy, r, circle(&ctx), CIRCLE_SEGMENTS)
                                .name(&circ.name);
                            plot_ui.line(circle_line);
                        }
                    }
                }

                if self.show_inverted {
                    for (idx, circ) in self.circles.iter().enumerate() {
                        if let Some(inverted) = self.invert_circle(&circ.circle) {
                            match inverted {
                                InvertedCircle::Circle { cx, cy, radius } => {
                                    // Fixed viewport - just draw the circle, let plot clip naturally
                                    let inv_circle = circle_2d(
                                        cx,
                                        cy,
                                        radius,
                                        active(&ctx),
                                        CIRCLE_SEGMENTS,
                                    )
                                    .name(format!("{}'", self.circles[idx].name));
                                    plot_ui.line(inv_circle);
                                }
                                InvertedCircle::Line { nx, ny, d } => {
                                    // Draw line nx*x + ny*y = d
                                    let len = (nx * nx + ny * ny).sqrt();
                                    if len > EPSILON {
                                        // Point on the line closest to origin
                                        let px = nx * d / (len * len);
                                        let py = ny * d / (len * len);
                                        // Direction along the line
                                        let dx = -ny / len;
                                        let dy = nx / len;

                                        let inv_line = infinite_line_2d(
                                            px,
                                            py,
                                            dx,
                                            dy,
                                            VIEWPORT_BOUNDS,
                                            active(&ctx),
                                        )
                                        .name(format!("{}' (line)", self.circles[idx].name));
                                        plot_ui.line(inv_line);
                                    }
                                }
                            }
                        }
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
                // Priority: inversion center > points > circles
                if self.is_near_inversion_center(mouse_x, mouse_y) {
                    self.dragging_inversion_center = true;
                } else if let Some(idx) = self.find_nearest_point(mouse_x, mouse_y) {
                    self.dragging_point = Some(idx);
                } else if let Some(idx) = self.find_nearest_circle(mouse_x, mouse_y) {
                    self.dragging_circle = Some(idx);
                }
            }

            // Handle active dragging
            if response.response.dragged() {
                if self.dragging_inversion_center {
                    self.inversion_cx = mouse_x;
                    self.inversion_cy = mouse_y;
                    self.update_inversion_circle();
                } else if let Some(idx) = self.dragging_point {
                    if idx < self.points.len() {
                        self.points[idx].set_position(mouse_x, mouse_y);
                    }
                } else if let Some(idx) = self.dragging_circle {
                    if idx < self.circles.len() {
                        if let Some((_, _, r)) = self.circles[idx].center_radius() {
                            self.circles[idx].set_center_radius(mouse_x, mouse_y, r);
                        }
                    }
                }
            }

            // Handle drag release
            if response.response.drag_stopped() {
                self.dragging_point = None;
                self.dragging_inversion_center = false;
                self.dragging_circle = None;
            }
        } else {
            self.dragging_point = None;
            self.dragging_inversion_center = false;
            self.dragging_circle = None;
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        // === Inversion Circle ===
        group_header(ui, "Inversion Circle (drag center)");

        ui.horizontal(|ui| {
            ui.label("Center:");
            let mut cx = self.inversion_cx as f32;
            let mut cy = self.inversion_cy as f32;
            let cx_changed = ui
                .add(egui::DragValue::new(&mut cx).speed(0.1).prefix("x: "))
                .changed();
            let cy_changed = ui
                .add(egui::DragValue::new(&mut cy).speed(0.1).prefix("y: "))
                .changed();
            if cx_changed || cy_changed {
                self.inversion_cx = f64::from(cx);
                self.inversion_cy = f64::from(cy);
                self.update_inversion_circle();
            }
        });

        ui.horizontal(|ui| {
            ui.label("Radius:");
            let mut r = self.inversion_radius as f32;
            if ui.add(egui::Slider::new(&mut r, 0.5..=4.0)).changed() {
                self.inversion_radius = f64::from(r);
                self.update_inversion_circle();
            }
        });

        // === Points ===
        section_separator(ui, Some("Points"));

        for idx in 0..self.points.len() {
            let name = self.points[idx].name.clone();
            let euclidean = self.points[idx].euclidean();

            ui.horizontal(|ui| {
                ui.label(&name);
                if let Some((x, y)) = euclidean {
                    let mut new_x = x as f32;
                    let mut new_y = y as f32;

                    let x_changed = ui
                        .add(egui::DragValue::new(&mut new_x).speed(0.1).prefix("x: "))
                        .changed();
                    let y_changed = ui
                        .add(egui::DragValue::new(&mut new_y).speed(0.1).prefix("y: "))
                        .changed();

                    if x_changed || y_changed {
                        self.points[idx].set_position(f64::from(new_x), f64::from(new_y));
                    }
                }
            });

            // Show inverted coordinates
            let inverted = self.invert_point(&self.points[idx].point);
            if let Some(inv_point) = inverted {
                if let Some((inv_x, inv_y)) = inv_point.to_euclidean() {
                    ui.horizontal(|ui| {
                        ui.label(format!("  {}' =", name));
                        ui.colored_label(active(&ctx), format!("({:.2}, {:.2})", inv_x, inv_y));
                    });
                } else {
                    ui.horizontal(|ui| {
                        ui.label(format!("  {}' =", name));
                        ui.colored_label(active(&ctx), "\u{221e} (at infinity)");
                    });
                }
            } else {
                ui.horizontal(|ui| {
                    ui.label(format!("  {}' =", name));
                    ui.colored_label(active(&ctx), "\u{221e} (at center)");
                });
            }
        }

        if ui.button("+ Add Point").clicked() && self.points.len() < 6 {
            let subscripts = [
                "\u{2081}", "\u{2082}", "\u{2083}", "\u{2084}", "\u{2085}", "\u{2086}",
            ];
            let idx = self.points.len();
            self.points.push(ScenePoint::new(
                1.0 + (idx + 1) as f64 * 0.5,
                0.5,
                &format!("P{}", subscripts[idx]),
            ));
        }

        // === Circles ===
        section_separator(ui, Some("Circles"));

        for idx in 0..self.circles.len() {
            let name = self.circles[idx].name.clone();
            if let Some((cx, cy, r)) = self.circles[idx].center_radius() {
                ui.horizontal(|ui| {
                    ui.label(&name);
                    let mut new_cx = cx as f32;
                    let mut new_cy = cy as f32;

                    let cx_changed = ui
                        .add(egui::DragValue::new(&mut new_cx).speed(0.1).prefix("x: "))
                        .changed();
                    let cy_changed = ui
                        .add(egui::DragValue::new(&mut new_cy).speed(0.1).prefix("y: "))
                        .changed();

                    if cx_changed || cy_changed {
                        self.circles[idx].set_center_radius(
                            f64::from(new_cx),
                            f64::from(new_cy),
                            r,
                        );
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("  r:");
                    let mut new_r = r as f32;
                    if ui.add(egui::Slider::new(&mut new_r, 0.2..=3.0)).changed() {
                        self.circles[idx].set_center_radius(cx, cy, f64::from(new_r));
                    }
                });

                // Check if circle passes through inversion center using CGA
                let dist_to_center =
                    ((cx - self.inversion_cx).powi(2) + (cy - self.inversion_cy).powi(2)).sqrt();
                let passes_through = (dist_to_center - r).abs() < 0.1;

                if passes_through {
                    ui.colored_label(active(&ctx), format!("  {}' = Line (through center)", name));
                } else {
                    ui.colored_label(active(&ctx), format!("  {}' = Circle", name));
                }
            }
        }

        if ui.button("+ Add Circle").clicked() && self.circles.len() < 4 {
            let idx = self.circles.len() + 1;
            let subscripts = ["\u{2081}", "\u{2082}", "\u{2083}", "\u{2084}"];
            self.circles.push(SceneCircle::new(
                2.0,
                idx as f64,
                0.7,
                &format!("C{}", subscripts[idx - 1]),
            ));
        }

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_inversion_circle, "Show inversion circle");
        ui.checkbox(&mut self.show_originals, "Show original objects");
        ui.checkbox(&mut self.show_inverted, "Show inverted objects");
        ui.checkbox(&mut self.show_rules, "Show transformation rules");

        if self.show_rules {
            ui.add_space(spacing::XS);
            info_box(
                ui,
                "Transformation Rules:\n\
                 \u{2022} Point \u{2192} Point\n\
                 \u{2022} Circle \u{2192} Circle (general)\n\
                 \u{2022} Circle thru center \u{2192} Line\n\
                 \u{2022} Line \u{2192} Circle (general)\n\
                 \u{2022} Line thru center \u{2192} Line",
            );
        }

        // === Presets ===
        section_separator(ui, Some("Presets"));

        if ui.button("Circle through center").clicked() {
            // Set up a circle that passes through the inversion center
            self.inversion_cx = 0.0;
            self.inversion_cy = 0.0;
            self.inversion_radius = 2.0;
            self.update_inversion_circle();
            self.circles.clear();
            self.circles
                .push(SceneCircle::new(2.0, 0.0, 2.0, "C\u{2081}"));
            self.points.clear();
            self.points.push(ScenePoint::new(3.0, 1.0, "P\u{2081}"));
        }

        if ui.button("Concentric circles").clicked() {
            self.inversion_cx = 0.0;
            self.inversion_cy = 0.0;
            self.inversion_radius = 2.0;
            self.update_inversion_circle();
            self.circles.clear();
            self.circles
                .push(SceneCircle::new(0.0, 0.0, 1.0, "C\u{2081}"));
            self.circles
                .push(SceneCircle::new(0.0, 0.0, 3.0, "C\u{2082}"));
            self.points.clear();
        }

        if ui.button("Reset").clicked() {
            *self = Self::default();
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            ui.colored_label(
                motor(&ctx),
                format!(
                    "Inversion: center=({:.1}, {:.1}) r={:.1}",
                    self.inversion_cx, self.inversion_cy, self.inversion_radius
                ),
            );
            ui.separator();
            ui.label(format!(
                "{} points, {} circles",
                self.points.len(),
                self.circles.len()
            ));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(CONFORMAL2_INVERSION_EDUCATION)
    }
}

/// Educational content for the Circle Inversion visualization.
const CONFORMAL2_INVERSION_EDUCATION: EducationalContent = EducationalContent {
    title: "Circle Inversion in CGA",

    overview: "\
Circle inversion (or geometric inversion) is a fundamental conformal transformation. \
Given an inversion circle with center O and radius r, every point P maps to P' such that:

    |OP| \u{00d7} |OP'| = r\u{00b2}

Points closer to O map farther away, and vice versa. The inversion circle itself \
is fixed (every point on it maps to itself).",

    math_background: "\
CGA INVERSION FORMULA:
In Conformal Geometric Algebra, inversion through a circle C is elegantly expressed as:

    P' = C \u{00d7} P \u{00d7} C\u{207b}\u{00b9}

where \u{00d7} is the geometric product and C\u{207b}\u{00b9} is the inverse of the circle.

This single formula handles ALL cases:
\u{2022} Point inversion: P' = C \u{00d7} P \u{00d7} C\u{207b}\u{00b9}
\u{2022} Circle inversion: C' = I \u{00d7} C \u{00d7} I\u{207b}\u{00b9}

PROPERTIES:
\u{2022} Circle NOT through O \u{2192} Circle
\u{2022} Circle THROUGH O \u{2192} Line (circle of infinite radius)
\u{2022} Line NOT through O \u{2192} Circle through O
\u{2022} Line THROUGH O \u{2192} Same line (self-inverse)

The CGA representation naturally handles the circle/line duality \
since lines are circles through the point at infinity.",

    how_to_use: "\
- DRAG the inversion circle center to move it
- Use the RADIUS SLIDER to change the inversion circle size
- DRAG points and circles to see their inversions update
- Try the PRESETS to see special cases like circles through the center
- Watch circles become LINES when they pass through the center",

    key_concepts: "\
- Inversion is a conformal transformation (preserves angles)
- Points inside the circle map outside, and vice versa
- Points on the inversion circle map to themselves
- Circles map to circles (or lines as special case)
- Two inversions in the same circle = identity
- CGA unifies the formula: X' = C \u{00d7} X \u{00d7} C\u{207b}\u{00b9}",

    resources: &[
        (
            "Circle Inversion (Wikipedia)",
            "https://en.wikipedia.org/wiki/Inversive_geometry",
        ),
        (
            "Conformal Geometric Algebra Wiki",
            "https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page",
        ),
    ],
};
