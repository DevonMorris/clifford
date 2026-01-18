//! Circle Inversion - 2D Conformal GA Visualization
//!
//! This demo demonstrates circle inversion (also known as geometric inversion),
//! a fundamental conformal transformation.
//!
//! In circle inversion through a circle C with center O and radius r:
//! - A point P maps to P' such that |OP| × |OP'| = r² and P, O, P' are collinear
//! - Circles map to circles (or lines if they pass through the center)
//! - Lines map to circles (or lines if they pass through the center)
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_inversion --release`

use clifford_viz::common::prelude::*;
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
const VIEWPORT_BOUNDS: f64 = 6.0;

/// Epsilon for numerical comparisons.
const EPSILON: f64 = 1e-10;

/// Maximum radius before we treat a circle as a line.
const MAX_DRAWABLE_RADIUS: f64 = 100.0;

/// Number of segments for circle rendering.
const CIRCLE_SEGMENTS: usize = 64;

/// Represents a draggable point that can be inverted.
#[derive(Clone)]
struct ScenePoint {
    /// X coordinate.
    x: f64,
    /// Y coordinate.
    y: f64,
    /// Display name.
    name: String,
}

impl ScenePoint {
    /// Creates a new scene point.
    fn new(x: f64, y: f64, name: &str) -> Self {
        Self {
            x,
            y,
            name: name.to_string(),
        }
    }
}

/// Represents a circle (original, not the inversion circle).
#[derive(Clone)]
struct SceneCircle {
    /// Center X coordinate.
    cx: f64,
    /// Center Y coordinate.
    cy: f64,
    /// Radius.
    radius: f64,
    /// Display name.
    name: String,
}

impl SceneCircle {
    /// Creates a new scene circle.
    fn new(cx: f64, cy: f64, radius: f64, name: &str) -> Self {
        Self {
            cx,
            cy,
            radius,
            name: name.to_string(),
        }
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
struct Conformal2InversionDemo {
    /// Inversion circle center x.
    inversion_cx: f64,
    /// Inversion circle center y.
    inversion_cy: f64,
    /// Inversion circle radius.
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
            inversion_cx: 0.0,
            inversion_cy: 0.0,
            inversion_radius: 2.0,

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
    /// Inverts a point through the inversion circle.
    ///
    /// For point P and inversion circle with center O and radius r:
    /// P' lies on ray OP such that |OP| × |OP'| = r²
    fn invert_point(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        // Translate to inversion circle center
        let dx = x - self.inversion_cx;
        let dy = y - self.inversion_cy;
        let dist_sq = dx * dx + dy * dy;

        // Point at center maps to infinity
        if dist_sq < EPSILON {
            return None;
        }

        // Inversion formula: P' = O + r²/|OP|² × (P - O)
        let scale = self.inversion_radius * self.inversion_radius / dist_sq;
        Some((
            self.inversion_cx + dx * scale,
            self.inversion_cy + dy * scale,
        ))
    }

    /// Inverts a circle through the inversion circle.
    ///
    /// Cases:
    /// 1. Circle through inversion center → Line
    /// 2. Circle not through center → Circle
    fn invert_circle(&self, cx: f64, cy: f64, radius: f64) -> InvertedCircle {
        // Distance from inversion center to circle center
        let dx = cx - self.inversion_cx;
        let dy = cy - self.inversion_cy;
        let d = (dx * dx + dy * dy).sqrt();

        // Check if circle passes through the inversion center
        // (distance from center to inversion center equals radius)
        let passes_through_center = (d - radius).abs() < EPSILON;

        if passes_through_center {
            // Circle through center inverts to a line
            // The line is perpendicular to the line from inversion center to circle center
            // and passes through the inverted point of the far side of the circle

            // Direction from inversion center to circle center
            if d < EPSILON {
                // Degenerate case: circle centered at inversion center
                // This should not happen if radius > 0 and passes through center
                return InvertedCircle::Line {
                    nx: 1.0,
                    ny: 0.0,
                    d: 0.0,
                };
            }

            let ux = dx / d;
            let uy = dy / d;

            // The far point of the circle from inversion center
            let far_x = cx + ux * radius;
            let far_y = cy + uy * radius;

            // Invert the far point
            if let Some((inv_x, inv_y)) = self.invert_point(far_x, far_y) {
                // Line through inv_point perpendicular to direction
                // Normal is (ux, uy), and line passes through (inv_x, inv_y)
                InvertedCircle::Line {
                    nx: ux,
                    ny: uy,
                    d: ux * inv_x + uy * inv_y,
                }
            } else {
                // Far point is at inversion center - shouldn't happen
                InvertedCircle::Line {
                    nx: ux,
                    ny: uy,
                    d: 0.0,
                }
            }
        } else {
            // Circle not through center inverts to another circle
            // Invert the two points where line from inversion center meets circle

            let r_sq = self.inversion_radius * self.inversion_radius;

            // Near and far points on the circle along the line to inversion center
            let near_dist = d - radius;
            let far_dist = d + radius;

            // Inverted distances
            let inv_near_dist = r_sq / near_dist;
            let inv_far_dist = r_sq / far_dist;

            // New radius and center distance
            let new_radius = (inv_near_dist - inv_far_dist).abs() / 2.0;
            let new_center_dist = (inv_near_dist + inv_far_dist) / 2.0;

            // Direction from inversion center to original circle center
            let (ux, uy) = if d > EPSILON {
                (dx / d, dy / d)
            } else {
                (1.0, 0.0)
            };

            // New center position
            let new_cx = self.inversion_cx + ux * new_center_dist;
            let new_cy = self.inversion_cy + uy * new_center_dist;

            InvertedCircle::Circle {
                cx: new_cx,
                cy: new_cy,
                radius: new_radius,
            }
        }
    }

    /// Find the nearest point to mouse position within threshold.
    fn find_nearest_point(&self, mouse_x: f64, mouse_y: f64) -> Option<usize> {
        let threshold = 0.4;
        let mut nearest_idx = None;
        let mut nearest_dist = threshold;

        for (idx, point) in self.points.iter().enumerate() {
            let dist = ((point.x - mouse_x).powi(2) + (point.y - mouse_y).powi(2)).sqrt();
            if dist < nearest_dist {
                nearest_dist = dist;
                nearest_idx = Some(idx);
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
            // Distance to circle edge
            let dist_to_center =
                ((circle.cx - mouse_x).powi(2) + (circle.cy - mouse_y).powi(2)).sqrt();
            let dist = (dist_to_center - circle.radius).abs();
            if dist < nearest_dist {
                nearest_dist = dist;
                nearest_idx = Some(idx);
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
                        plot_ui.points(
                            Points::new(vec![[pt.x, pt.y]])
                                .color(point(&ctx))
                                .radius(7.0)
                                .filled(true)
                                .name(&pt.name),
                        );
                    }
                }

                if self.show_inverted {
                    for (idx, pt) in self.points.iter().enumerate() {
                        if let Some((inv_x, inv_y)) = self.invert_point(pt.x, pt.y) {
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
                                let connect_line = line_segment(
                                    pt.x,
                                    pt.y,
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

                // Draw original circles and their inversions
                if self.show_originals {
                    for circ in &self.circles {
                        let circle_line =
                            circle_2d(circ.cx, circ.cy, circ.radius, circle(&ctx), CIRCLE_SEGMENTS)
                                .name(&circ.name);
                        plot_ui.line(circle_line);
                    }
                }

                if self.show_inverted {
                    for (idx, circ) in self.circles.iter().enumerate() {
                        let inverted = self.invert_circle(circ.cx, circ.cy, circ.radius);
                        match inverted {
                            InvertedCircle::Circle { cx, cy, radius } => {
                                if radius < MAX_DRAWABLE_RADIUS
                                    && cx.abs() < VIEWPORT_BOUNDS * 2.0
                                    && cy.abs() < VIEWPORT_BOUNDS * 2.0
                                {
                                    let inv_circle =
                                        circle_2d(cx, cy, radius, active(&ctx), CIRCLE_SEGMENTS)
                                            .name(format!("{}'", self.circles[idx].name));
                                    plot_ui.line(inv_circle);
                                }
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
                } else if let Some(idx) = self.dragging_point {
                    if idx < self.points.len() {
                        self.points[idx].x = mouse_x;
                        self.points[idx].y = mouse_y;
                    }
                } else if let Some(idx) = self.dragging_circle {
                    if idx < self.circles.len() {
                        self.circles[idx].cx = mouse_x;
                        self.circles[idx].cy = mouse_y;
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
            ui.add(
                egui::DragValue::new(&mut self.inversion_cx)
                    .speed(0.1)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.inversion_cy)
                    .speed(0.1)
                    .prefix("y: "),
            );
        });

        ui.horizontal(|ui| {
            ui.label("Radius:");
            ui.add(egui::Slider::new(&mut self.inversion_radius, 0.5..=4.0));
        });

        // === Points ===
        section_separator(ui, Some("Points"));

        for point in &mut self.points {
            ui.horizontal(|ui| {
                ui.label(&point.name);
                ui.add(egui::DragValue::new(&mut point.x).speed(0.1).prefix("x: "));
                ui.add(egui::DragValue::new(&mut point.y).speed(0.1).prefix("y: "));
            });

            // Show inverted coordinates
            let demo_state = Conformal2InversionDemo {
                inversion_cx: self.inversion_cx,
                inversion_cy: self.inversion_cy,
                inversion_radius: self.inversion_radius,
                ..Default::default()
            };
            if let Some((inv_x, inv_y)) = demo_state.invert_point(point.x, point.y) {
                ui.horizontal(|ui| {
                    ui.label(format!("  {}' =", point.name));
                    ui.colored_label(active(&ctx), format!("({:.2}, {:.2})", inv_x, inv_y));
                });
            } else {
                ui.horizontal(|ui| {
                    ui.label(format!("  {}' =", point.name));
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

        for circle in &mut self.circles {
            ui.horizontal(|ui| {
                ui.label(&circle.name);
                ui.add(
                    egui::DragValue::new(&mut circle.cx)
                        .speed(0.1)
                        .prefix("x: "),
                );
                ui.add(
                    egui::DragValue::new(&mut circle.cy)
                        .speed(0.1)
                        .prefix("y: "),
                );
            });
            ui.horizontal(|ui| {
                ui.label("  r:");
                ui.add(egui::Slider::new(&mut circle.radius, 0.2..=3.0));
            });

            // Check if circle passes through inversion center
            let dx = circle.cx - self.inversion_cx;
            let dy = circle.cy - self.inversion_cy;
            let dist = (dx * dx + dy * dy).sqrt();
            let passes_through = (dist - circle.radius).abs() < 0.1;

            if passes_through {
                ui.colored_label(
                    active(&ctx),
                    format!("  {}' = Line (through center)", circle.name),
                );
            } else {
                ui.colored_label(active(&ctx), format!("  {}' = Circle", circle.name));
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

    |OP| × |OP'| = r²

Points closer to O map farther away, and vice versa. The inversion circle itself \
is fixed (every point on it maps to itself).",

    math_background: "\
INVERSION FORMULA:
For a point P at distance d from center O, the inverted point P' is at distance r²/d \
along the same ray from O:

    P' = O + (r²/|OP|²) × (P - O)

CIRCLE INVERSION:
- Circle NOT through O → Circle
- Circle THROUGH O → Line (circle of infinite radius)

LINE INVERSION:
- Line NOT through O → Circle through O
- Line THROUGH O → Same line (self-inverse)

In CGA, inversion through a circle C is expressed as a reflection:
    X' = -C · X · C⁻¹
where · is the geometric product.",

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
- Composing inversions in two circles = Mobius transformation",

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

fn main() -> eframe::Result<()> {
    run_app::<Conformal2InversionDemo>()
}
