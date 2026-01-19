//! 2D Projective Geometry (PGA) Visualization
//!
//! This demo demonstrates point-line geometry in 2D Projective Geometric Algebra,
//! including:
//!
//! - Creating and manipulating points
//! - Join operation (wedge): two points → line through them
//! - Meet operation (antiwedge): two lines → intersection point
//! - Motor transformations (rotation + translation)
//!
//! ## Mathematical Background
//!
//! In 2D PGA, we represent:
//! - **Points** as grade-1 elements: `P = x·e₁ + y·e₂ + w·e₀`
//! - **Lines** as grade-2 elements: `L = d·e₁₂ + a·e₁₀ + b·e₂₀` (representing ax + by + d = 0)
//! - **Motors** transform points and lines via the antisandwich product
//!
//! The join (∧) and meet (∨) operations are fundamental:
//! - `Line = Point₁ ∧ Point₂` creates the line through two points
//! - `Point = Line₁ ∨ Line₂` finds the intersection of two lines

use crate::common::prelude::*;
use clifford::ops::{Join, Meet, Transform};
use clifford::specialized::projective::dim2::{Line, Motor, Point};
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
/// Using fixed bounds prevents disorienting viewport jumps.
const VIEWPORT_BOUNDS: f64 = 5.0;

/// A draggable point in the visualization.
#[derive(Clone)]
pub struct DraggablePoint {
    /// Unique identifier.
    id: usize,
    /// The geometric point.
    point: Point<f64>,
    /// Whether this point is currently selected.
    selected: bool,
    /// Display name (e.g., "A", "B", "C").
    name: String,
}

impl DraggablePoint {
    /// Creates a new draggable point.
    fn new(id: usize, x: f64, y: f64, name: &str) -> Self {
        Self {
            id,
            point: Point::from_cartesian(x, y),
            selected: false,
            name: name.to_string(),
        }
    }

    /// Returns the Cartesian coordinates of this point.
    fn cartesian(&self) -> Option<(f64, f64)> {
        self.point.to_cartesian()
    }
}

/// A line derived from joining two points.
#[derive(Clone)]
struct DerivedLine {
    /// Indices of the two points that define this line.
    from_points: (usize, usize),
    /// The geometric line.
    line: Line<f64>,
    /// Whether this line is currently selected.
    selected: bool,
}

/// A point derived from meeting two lines.
#[derive(Clone)]
struct DerivedPoint {
    /// Indices of the two lines that define this point.
    from_lines: (usize, usize),
    /// The geometric point.
    point: Point<f64>,
}

/// Current tool mode for the demo.
#[derive(Clone, Copy, PartialEq, Default)]
pub enum ToolMode {
    /// Select and drag points.
    #[default]
    Select,
    /// Add new points by clicking.
    AddPoint,
}

/// Demo state for 2D Projective GA visualization.
pub struct Projective2Demo {
    /// User-created points.
    points: Vec<DraggablePoint>,
    /// Lines created by joining points.
    derived_lines: Vec<DerivedLine>,
    /// Points created by meeting lines.
    derived_points: Vec<DerivedPoint>,
    /// Current tool mode.
    tool_mode: ToolMode,
    /// Next point ID for unique naming.
    next_point_id: usize,
    /// Motor for transforming all objects.
    motor: Motor<f64>,
    /// Rotation angle (radians) for motor.
    motor_rotation: f32,
    /// Translation X for motor.
    motor_tx: f32,
    /// Translation Y for motor.
    motor_ty: f32,
    /// Animation state.
    animation: Animation,
    /// Whether to show the grid.
    show_grid: bool,
    /// Whether to show coordinate labels.
    show_coordinates: bool,
    /// Whether to apply motor transformation.
    apply_motor: bool,
    /// Whether to show normal vectors on lines.
    show_normals: bool,
    /// Index of point being dragged (if any).
    dragging_point: Option<usize>,
    /// Indices of selected points for join operation.
    join_selection: Vec<usize>,
    /// Indices of selected lines for meet operation.
    meet_selection: Vec<usize>,
}

impl Default for Projective2Demo {
    fn default() -> Self {
        // Create some initial points for demonstration
        let points = vec![
            DraggablePoint::new(0, -2.0, 1.0, "A"),
            DraggablePoint::new(1, 2.0, 1.0, "B"),
            DraggablePoint::new(2, -1.0, -2.0, "C"),
            DraggablePoint::new(3, 1.0, -2.0, "D"),
        ];

        Self {
            points,
            derived_lines: Vec::new(),
            derived_points: Vec::new(),
            tool_mode: ToolMode::Select,
            next_point_id: 4,
            motor: Motor::identity(),
            motor_rotation: 0.0,
            motor_tx: 0.0,
            motor_ty: 0.0,
            animation: Animation::with_duration(4.0),
            show_grid: true,
            show_coordinates: true,
            apply_motor: false,
            show_normals: false,
            dragging_point: None,
            join_selection: Vec::new(),
            meet_selection: Vec::new(),
        }
    }
}

impl Projective2Demo {
    /// Update the motor from rotation and translation parameters.
    ///
    /// In 2D PGA, we apply rotation and translation separately since
    /// Motor × Motor = Flector. We store just the rotation motor and
    /// apply translation afterwards.
    fn update_motor(&mut self) {
        // Store the rotation motor; translation is applied separately
        self.motor = Motor::from_rotation(f64::from(self.motor_rotation));
    }

    /// Transform a point using the current motor (if enabled).
    ///
    /// Applies rotation first, then translation.
    fn transform_point(&self, p: &Point<f64>) -> Point<f64> {
        if self.apply_motor {
            // Apply rotation
            let rotated = self.motor.transform(p);
            // Apply translation
            let translation =
                Motor::from_translation(f64::from(self.motor_tx), f64::from(self.motor_ty));
            translation.transform(&rotated)
        } else {
            *p
        }
    }

    /// Transform a line using the current motor (if enabled).
    ///
    /// Applies rotation first, then translation.
    fn transform_line(&self, l: &Line<f64>) -> Line<f64> {
        if self.apply_motor {
            // Apply rotation
            let rotated = self.motor.transform(l);
            // Apply translation
            let translation =
                Motor::from_translation(f64::from(self.motor_tx), f64::from(self.motor_ty));
            translation.transform(&rotated)
        } else {
            *l
        }
    }

    /// Generate a name for a new point.
    fn next_point_name(&self) -> String {
        // Use letters A-Z, then A1, B1, etc.
        let idx = self.next_point_id;
        if idx < 26 {
            (b'A' + idx as u8) as char
        } else {
            let letter = (b'A' + (idx % 26) as u8) as char;
            let number = idx / 26;
            return format!("{}{}", letter, number);
        }
        .to_string()
    }

    /// Perform join operation on selected points.
    fn perform_join(&mut self) {
        if self.join_selection.len() == 2 {
            let p1_idx = self.join_selection[0];
            let p2_idx = self.join_selection[1];

            if p1_idx < self.points.len() && p2_idx < self.points.len() {
                let p1 = &self.points[p1_idx].point;
                let p2 = &self.points[p2_idx].point;
                let line = p1.join(p2);

                self.derived_lines.push(DerivedLine {
                    from_points: (self.points[p1_idx].id, self.points[p2_idx].id),
                    line,
                    selected: false,
                });
            }

            // Clear selection
            self.join_selection.clear();
            for point in &mut self.points {
                point.selected = false;
            }
        }
    }

    /// Perform meet operation on selected lines.
    fn perform_meet(&mut self) {
        if self.meet_selection.len() == 2 {
            let l1_idx = self.meet_selection[0];
            let l2_idx = self.meet_selection[1];

            if l1_idx < self.derived_lines.len() && l2_idx < self.derived_lines.len() {
                let l1 = &self.derived_lines[l1_idx].line;
                let l2 = &self.derived_lines[l2_idx].line;
                let point = l1.meet(l2);

                self.derived_points.push(DerivedPoint {
                    from_lines: (l1_idx, l2_idx),
                    point,
                });
            }

            // Clear selection
            self.meet_selection.clear();
            for line in &mut self.derived_lines {
                line.selected = false;
            }
        }
    }

    /// Update derived lines when points move.
    fn update_derived_lines(&mut self) {
        for derived in &mut self.derived_lines {
            // Find points by ID
            let p1 = self.points.iter().find(|p| p.id == derived.from_points.0);
            let p2 = self.points.iter().find(|p| p.id == derived.from_points.1);

            if let (Some(p1), Some(p2)) = (p1, p2) {
                derived.line = p1.point.join(&p2.point);
            }
        }
    }

    /// Update derived points when lines change.
    fn update_derived_points(&mut self) {
        use clifford::ops::Antiwedge;
        for derived in &mut self.derived_points {
            let (l1_idx, l2_idx) = derived.from_lines;
            if l1_idx < self.derived_lines.len() && l2_idx < self.derived_lines.len() {
                let l1 = &self.derived_lines[l1_idx].line;
                let l2 = &self.derived_lines[l2_idx].line;
                derived.point = l1.antiwedge(l2);
            }
        }
    }
}

impl VisualizationApp for Projective2Demo {
    fn name(&self) -> &'static str {
        "Projective 2D - Point-Line Geometry"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            self.motor_rotation = self.animation.angle();
            self.update_motor();
        }
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let response = Plot::new("projective2_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(false) // We handle drag ourselves for point manipulation
            .allow_boxed_zoom(false)
            .show(ui, |plot_ui| {
                // Draw coordinate grid (using fixed bounds for stability)
                if self.show_grid {
                    for l in grid_2d(&ctx, VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(l);
                    }
                    for axis in axes_2d(&ctx, VIEWPORT_BOUNDS) {
                        plot_ui.line(axis);
                    }
                }

                // Draw derived lines (before points so points appear on top)
                for (idx, derived) in self.derived_lines.iter().enumerate() {
                    let transformed_line = self.transform_line(&derived.line);
                    let color = if derived.selected {
                        selected(&ctx)
                    } else {
                        line(&ctx)
                    };

                    // Draw line using homogeneous coordinates
                    // Line equation: nx*x + ny*y + d = 0 (standard form ax + by + c = 0)
                    let plot_line = line_from_homogeneous(
                        transformed_line.nx(),
                        transformed_line.ny(),
                        transformed_line.d(),
                        VIEWPORT_BOUNDS,
                        color,
                    )
                    .name(format!("Line {}", idx + 1));
                    plot_ui.line(plot_line);

                    // Draw normal vector if enabled
                    if self.show_normals {
                        // Normal vector is (nx, ny) from the line equation
                        let a = transformed_line.nx();
                        let b = transformed_line.ny();
                        let c = transformed_line.d();
                        let norm_sq = a * a + b * b;
                        if norm_sq > 1e-10 {
                            let px = -a * c / norm_sq;
                            let py = -b * c / norm_sq;
                            let scale = 0.5;
                            for arrow in arrow_2d(
                                px,
                                py,
                                a * scale / norm_sq.sqrt(),
                                b * scale / norm_sq.sqrt(),
                                with_alpha(color, 150),
                            ) {
                                plot_ui.line(arrow);
                            }
                        }
                    }
                }

                // Draw derived intersection points
                for derived in &self.derived_points {
                    let transformed_pt = self.transform_point(&derived.point);
                    if let Some((x, y)) = transformed_pt.to_cartesian() {
                        // Finite intersection point
                        plot_ui.points(
                            Points::new(vec![[x, y]])
                                .color(plane(&ctx))
                                .radius(8.0)
                                .filled(true)
                                .name("Intersection"),
                        );
                    } else {
                        // Ideal point (parallel lines) - show direction indicator
                        let dx = transformed_pt.x();
                        let dy = transformed_pt.y();
                        let len = (dx * dx + dy * dy).sqrt();
                        if len > 1e-10 {
                            // Draw at edge of view to indicate direction
                            let scale = VIEWPORT_BOUNDS * 0.9 / len;
                            plot_ui.points(
                                Points::new(vec![[dx * scale, dy * scale]])
                                    .color(with_alpha(plane(&ctx), 100))
                                    .radius(6.0)
                                    .filled(false)
                                    .name("Ideal point (inf)"),
                            );
                        }
                    }
                }

                // Draw user-created points
                for point_data in &self.points {
                    let transformed_pt = self.transform_point(&point_data.point);
                    if let Some((x, y)) = transformed_pt.to_cartesian() {
                        let color = if point_data.selected {
                            selected(&ctx)
                        } else {
                            point(&ctx)
                        };
                        let radius = if point_data.selected { 8.0 } else { 6.0 };

                        plot_ui.points(
                            Points::new(vec![[x, y]])
                                .color(color)
                                .radius(radius)
                                .filled(true)
                                .name(&point_data.name),
                        );

                        // Draw coordinate label
                        if self.show_coordinates {
                            // Use a small text marker (egui_plot doesn't support text directly)
                            // We'll show coordinates in the info panel instead
                        }
                    }
                }
            });

        // Handle mouse interactions for point creation/selection/dragging
        if let Some(pos) = response.response.interact_pointer_pos() {
            let plot_pos = response.transform.value_from_position(pos);
            let mouse_x = plot_pos.x;
            let mouse_y = plot_pos.y;

            // Helper: find nearest point within threshold
            let find_nearest_point = |points: &[DraggablePoint], x: f64, y: f64| -> Option<usize> {
                let threshold = 0.4; // Distance in plot units
                let mut nearest_idx = None;
                let mut nearest_dist = threshold;

                for (idx, point_data) in points.iter().enumerate() {
                    if let Some((px, py)) = point_data.cartesian() {
                        let dist = ((px - x).powi(2) + (py - y).powi(2)).sqrt();
                        if dist < nearest_dist {
                            nearest_dist = dist;
                            nearest_idx = Some(idx);
                        }
                    }
                }
                nearest_idx
            };

            match self.tool_mode {
                ToolMode::AddPoint => {
                    // Create a new point on click
                    if response.response.clicked() {
                        let name = self.next_point_name();
                        self.points.push(DraggablePoint::new(
                            self.next_point_id,
                            mouse_x,
                            mouse_y,
                            &name,
                        ));
                        self.next_point_id += 1;
                    }
                }
                ToolMode::Select => {
                    // Handle drag start - find point to drag
                    if response.response.drag_started() {
                        self.dragging_point = find_nearest_point(&self.points, mouse_x, mouse_y);
                    }

                    // Handle active dragging - move the point
                    if response.response.dragged() {
                        if let Some(idx) = self.dragging_point {
                            if idx < self.points.len() {
                                self.points[idx].point = Point::from_cartesian(mouse_x, mouse_y);
                            }
                        }
                    }

                    // Handle drag release
                    if response.response.drag_stopped() {
                        self.dragging_point = None;
                    }

                    // Handle click (without drag) for selection toggle
                    if response.response.clicked() && self.dragging_point.is_none() {
                        if let Some(idx) = find_nearest_point(&self.points, mouse_x, mouse_y) {
                            self.points[idx].selected = !self.points[idx].selected;
                            // Update join selection
                            self.join_selection = self
                                .points
                                .iter()
                                .enumerate()
                                .filter(|(_, p)| p.selected)
                                .map(|(i, _)| i)
                                .collect();
                        } else {
                            // Click on empty space - deselect all
                            for point in &mut self.points {
                                point.selected = false;
                            }
                            self.join_selection.clear();
                        }
                    }
                }
            }
        } else {
            // Mouse left the plot area - cancel any drag
            self.dragging_point = None;
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Tool Mode Selection ===
        group_header(ui, "Tools");
        ui.horizontal(|ui| {
            if ui
                .selectable_label(self.tool_mode == ToolMode::Select, "Select")
                .clicked()
            {
                self.tool_mode = ToolMode::Select;
            }
            if ui
                .selectable_label(self.tool_mode == ToolMode::AddPoint, "Add Point")
                .clicked()
            {
                self.tool_mode = ToolMode::AddPoint;
            }
        });

        // === Point Management ===
        section_separator(ui, Some("Points"));

        // Add point button
        if ui.button("Add Point at Origin").clicked() {
            let name = self.next_point_name();
            self.points
                .push(DraggablePoint::new(self.next_point_id, 0.0, 0.0, &name));
            self.next_point_id += 1;
        }

        // List points with edit controls
        ui.add_space(spacing::XS);
        let mut points_to_remove = Vec::new();
        let mut selection_changed = false;

        for (idx, point) in self.points.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                // Selection checkbox
                if ui.checkbox(&mut point.selected, "").changed() {
                    selection_changed = true;
                }

                // Point name
                ui.label(&point.name);

                // Coordinate inputs
                if let Some((x, y)) = point.cartesian() {
                    let mut new_x = x as f32;
                    let mut new_y = y as f32;

                    ui.add(egui::DragValue::new(&mut new_x).speed(0.1).prefix("x: "));
                    ui.add(egui::DragValue::new(&mut new_y).speed(0.1).prefix("y: "));

                    if (new_x as f64 - x).abs() > 1e-6 || (new_y as f64 - y).abs() > 1e-6 {
                        point.point = Point::from_cartesian(f64::from(new_x), f64::from(new_y));
                    }
                }

                // Delete button
                if ui.button("[del]").clicked() {
                    points_to_remove.push(idx);
                }
            });
        }

        // Remove points (in reverse order to preserve indices)
        for idx in points_to_remove.into_iter().rev() {
            self.points.remove(idx);
        }

        // Update join selection from point selection
        if selection_changed {
            self.join_selection = self
                .points
                .iter()
                .enumerate()
                .filter(|(_, p)| p.selected)
                .map(|(i, _)| i)
                .collect();
        }

        // === Join Operation (wedge) ===
        section_separator(ui, Some("Join (^) - Line through Points"));
        ui.label(format!("Selected: {} points", self.join_selection.len()));

        if self.join_selection.len() == 2 {
            let p1_name = &self.points[self.join_selection[0]].name;
            let p2_name = &self.points[self.join_selection[1]].name;
            if ui
                .button(format!("Create Line {} ^ {}", p1_name, p2_name))
                .clicked()
            {
                self.perform_join();
                self.update_derived_lines();
            }
        } else {
            ui.label("Select exactly 2 points to create a line");
        }

        // Lines list
        if !self.derived_lines.is_empty() {
            ui.add_space(spacing::XS);
            ui.label("Lines:");
            let mut lines_to_remove = Vec::new();

            for (idx, line_data) in self.derived_lines.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    if ui.checkbox(&mut line_data.selected, "").changed() {
                        // Update meet selection
                    }
                    let l = &line_data.line;
                    // Line equation: nx*x + ny*y + d = 0
                    ui.label(format!(
                        "L{}: {:.2}x + {:.2}y + {:.2} = 0",
                        idx + 1,
                        l.nx(),
                        l.ny(),
                        l.d()
                    ));
                    if ui.button("[del]").clicked() {
                        lines_to_remove.push(idx);
                    }
                });
            }

            // Remove lines
            for idx in lines_to_remove.into_iter().rev() {
                self.derived_lines.remove(idx);
            }

            // Update meet selection
            self.meet_selection = self
                .derived_lines
                .iter()
                .enumerate()
                .filter(|(_, l)| l.selected)
                .map(|(i, _)| i)
                .collect();
        }

        // === Meet Operation (antiwedge) ===
        section_separator(ui, Some("Meet (v) - Intersection"));
        ui.label(format!("Selected: {} lines", self.meet_selection.len()));

        if self.meet_selection.len() == 2 {
            if ui
                .button(format!(
                    "Find Intersection L{} v L{}",
                    self.meet_selection[0] + 1,
                    self.meet_selection[1] + 1
                ))
                .clicked()
            {
                self.perform_meet();
                self.update_derived_points();
            }
        } else {
            ui.label("Select exactly 2 lines to find intersection");
        }

        // Intersection points display
        if !self.derived_points.is_empty() {
            ui.add_space(spacing::XS);
            ui.label("Intersections:");
            for (idx, derived) in self.derived_points.iter().enumerate() {
                if let Some((x, y)) = derived.point.to_cartesian() {
                    ui.label(format!("P{}: ({:.2}, {:.2})", idx + 1, x, y));
                } else {
                    ui.label(format!("P{}: ideal point (parallel lines)", idx + 1));
                }
            }
        }

        // === Motor Transform ===
        section_separator(ui, Some("Motor Transform"));
        ui.checkbox(&mut self.apply_motor, "Apply motor transformation");

        if self.apply_motor {
            angle_slider_range(
                ui,
                "Rotation theta",
                &mut self.motor_rotation,
                -360.0,
                360.0,
            );

            ui.horizontal(|ui| {
                ui.label("Translation");
                ui.add(
                    egui::DragValue::new(&mut self.motor_tx)
                        .speed(0.1)
                        .prefix("x: "),
                );
                ui.add(
                    egui::DragValue::new(&mut self.motor_ty)
                        .speed(0.1)
                        .prefix("y: "),
                );
            });

            if ui.button("Reset Motor").clicked() {
                self.motor_rotation = 0.0;
                self.motor_tx = 0.0;
                self.motor_ty = 0.0;
            }

            self.update_motor();

            // Motor display
            ga_value_display(
                ui,
                "M",
                &[
                    ("e_1", self.motor.ty() as f32),
                    ("e_2", self.motor.tx() as f32),
                    ("e_3", self.motor.r() as f32),
                    ("e_1_2_3", self.motor.ps() as f32),
                ],
            );

            ui.add_space(spacing::XS);
            animation_controls(ui, &mut self.animation);
            progress_slider(ui, &mut self.animation);
        }

        // === Display Options (compact) ===
        section_separator(ui, Some("Display Options"));
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_grid, "Grid");
            ui.checkbox(&mut self.show_normals, "Normals");
        });
        ui.checkbox(&mut self.show_coordinates, "Coordinates in list");

        // Update derived geometry when points change
        self.update_derived_lines();
        self.update_derived_points();
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            ui.colored_label(point(&ctx), format!("{} points", self.points.len()));
            ui.separator();
            ui.colored_label(line(&ctx), format!("{} lines", self.derived_lines.len()));
            ui.separator();
            ui.colored_label(
                plane(&ctx),
                format!("{} intersections", self.derived_points.len()),
            );
            ui.separator();
            if self.apply_motor {
                ui.label(format!(
                    "Motor: theta={:.1} deg t=({:.1}, {:.1})",
                    self.motor_rotation.to_degrees(),
                    self.motor_tx,
                    self.motor_ty
                ));
            } else {
                ui.label("Motor: off");
            }
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(PROJECTIVE2_EDUCATION)
    }
}

/// Educational content for the 2D Projective GA visualization.
const PROJECTIVE2_EDUCATION: EducationalContent = EducationalContent {
    title: "2D Projective Geometric Algebra (PGA)",

    overview: "\
This visualization demonstrates 2D Projective Geometric Algebra, which provides \
a unified framework for point-line geometry and rigid transformations.

Key insight: In PGA, geometric operations become algebraic products:
- The JOIN of two points (line through them) is their wedge product (^)
- The MEET of two lines (intersection point) is their antiwedge product (v)
- MOTORS encode rotation and translation as a single algebraic element

This demo lets you explore these operations interactively.",

    math_background: "\
POINTS are grade-1 elements in homogeneous coordinates:
    P = x*e_1 + y*e_2 + w*e_0

For a finite point at (x, y), we set w = 1.

LINES are grade-2 elements representing ax + by + c = 0:
    L = c*e_1_2 + a*e_1_0 + b*e_2_0

JOIN (^) - Line through two points:
    L = P_1 ^ P_2

MEET (v) - Intersection of two lines:
    P = L_1 v L_2

If lines are parallel, the result is an IDEAL POINT (w = 0), \
representing the direction at infinity.

MOTORS are elements of the odd subalgebra (grades 1 and 3):
    M = ty*e_1 + tx*e_2 + r*e_3 + ps*e_1_2_3

They transform geometry via the antisandwich product:
    P' = M^-1PM  (point transformation)
    L' = M^-1LM  (line transformation)",

    how_to_use: "\
- ADD POINTS: Select 'Add Point' tool and click on the plot
- SELECT POINTS: Select 'Select' tool and click near a point to toggle selection
- DRAG POINTS: Select 'Select' tool and drag a point to move it
- JOIN OPERATION: Select 2 points, then click 'Create Line'
- MEET OPERATION: Select 2 lines (checkboxes), click 'Find Intersection'
- MOTOR TRANSFORM: Enable 'Apply motor transformation', adjust sliders
- Lines update in real-time as you drag points
- Enable 'Normals' to visualize line orientations",

    key_concepts: "\
- Homogeneous coordinates: P = (x, y, w) with w=1 for finite points
- Ideal points (w=0) represent directions at infinity
- Wedge product (^) computes JOIN: line through two points
- Antiwedge product (v) computes MEET: intersection of two lines
- Motors compose rotation and translation into one operation
- Antisandwich product transforms geometry while preserving incidence",

    resources: &[
        (
            "Rigid Geometric Algebra Wiki - 2D PGA",
            "https://rigidgeometricalgebra.org/wiki/index.php?title=2D_projective_geometric_algebra",
        ),
        (
            "Look, Ma, No Matrices!",
            "https://enkimute.github.io/LookMaNoMatrices/",
        ),
        (
            "Steven De Keninck - PGA Tutorial",
            "https://bivector.net/PGA4CS.pdf",
        ),
    ],
};
