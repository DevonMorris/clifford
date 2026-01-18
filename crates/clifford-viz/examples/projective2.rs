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
//! Run with: `cargo run -p clifford-viz --example projective2 --release`
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

use clifford::ops::{Join, Meet, Transform};
use clifford::specialized::projective::dim2::{Line, Motor, Point};
use clifford_viz::common::prelude::*;
use egui_plot::{Plot, Points};

/// A draggable point in the visualization.
#[derive(Clone)]
struct DraggablePoint {
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
#[derive(Clone, Copy, PartialEq)]
enum ToolMode {
    /// Select and drag points.
    Select,
    /// Add new points by clicking.
    AddPoint,
}

/// Demo state for 2D Projective GA visualization.
struct Projective2Demo {
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
    /// Note: Interactive dragging requires more complex mouse handling.
    /// For now, points are edited via the control panel.
    _dragging_point: Option<usize>,
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
            _dragging_point: None,
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

    fn render(&self, ui: &mut egui::Ui) {
        let bounds = 5.0;

        let response = Plot::new("projective2_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(false) // We handle drag ourselves for point manipulation
            .allow_boxed_zoom(false)
            .show(ui, |plot_ui| {
                // Draw coordinate grid
                if self.show_grid {
                    for line in grid_2d(bounds, 1.0) {
                        plot_ui.line(line);
                    }
                    for axis in axes_2d(bounds) {
                        plot_ui.line(axis);
                    }
                }

                // Draw derived lines (before points so points appear on top)
                for (idx, derived) in self.derived_lines.iter().enumerate() {
                    let line = self.transform_line(&derived.line);
                    let color = if derived.selected {
                        palette::ROTOR
                    } else {
                        palette::LINE
                    };

                    // Draw line using homogeneous coordinates
                    let plot_line = line_from_homogeneous(
                        line.normal_x(),
                        line.normal_y(),
                        line.dist(),
                        bounds,
                        color,
                    )
                    .name(format!("Line {}", idx + 1));
                    plot_ui.line(plot_line);

                    // Draw normal vector if enabled
                    if self.show_normals {
                        // Find a point on the line for drawing normal
                        let a = line.normal_x();
                        let b = line.normal_y();
                        let c = line.dist();
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
                    let point = self.transform_point(&derived.point);
                    if let Some((x, y)) = point.to_cartesian() {
                        // Finite intersection point
                        plot_ui.points(
                            Points::new(vec![[x, y]])
                                .color(palette::PLANE)
                                .radius(8.0)
                                .filled(true)
                                .name("Intersection"),
                        );
                    } else {
                        // Ideal point (parallel lines) - show direction indicator
                        let dx = point.x();
                        let dy = point.y();
                        let len = (dx * dx + dy * dy).sqrt();
                        if len > 1e-10 {
                            // Draw at edge of view to indicate direction
                            let scale = bounds * 0.9 / len;
                            plot_ui.points(
                                Points::new(vec![[dx * scale, dy * scale]])
                                    .color(with_alpha(palette::PLANE, 100))
                                    .radius(6.0)
                                    .filled(false)
                                    .name("Ideal point (∞)"),
                            );
                        }
                    }
                }

                // Draw user-created points
                for point_data in &self.points {
                    let point = self.transform_point(&point_data.point);
                    if let Some((x, y)) = point.to_cartesian() {
                        let color = if point_data.selected {
                            palette::ROTOR
                        } else {
                            palette::POINT
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

        // Handle mouse interactions for point creation/selection
        // Note: egui_plot interactions are limited, so we use the response
        if let Some(pos) = response.response.hover_pos() {
            let _plot_pos = response.transform.value_from_position(pos);

            // Handle click
            if response.response.clicked() {
                // This would need additional state management for proper click handling
                // For now, point manipulation is done via the controls panel
            }
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // Tool mode selection
        section_separator(ui, Some("Tools"));
        ui.horizontal(|ui| {
            if ui
                .selectable_label(self.tool_mode == ToolMode::Select, "🖱 Select")
                .clicked()
            {
                self.tool_mode = ToolMode::Select;
            }
            if ui
                .selectable_label(self.tool_mode == ToolMode::AddPoint, "➕ Add Point")
                .clicked()
            {
                self.tool_mode = ToolMode::AddPoint;
            }
        });

        // Point management
        section_separator(ui, Some("Points"));

        // Add point button
        if ui.button("Add Point at Origin").clicked() {
            let name = self.next_point_name();
            self.points
                .push(DraggablePoint::new(self.next_point_id, 0.0, 0.0, &name));
            self.next_point_id += 1;
        }

        // List points with edit controls
        ui.add_space(4.0);
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
                if ui.button("🗑").clicked() {
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

        // Join operation (wedge)
        section_separator(ui, Some("Join (∧) - Line through Points"));
        ui.label(format!("Selected: {} points", self.join_selection.len()));

        if self.join_selection.len() == 2 {
            let p1_name = &self.points[self.join_selection[0]].name;
            let p2_name = &self.points[self.join_selection[1]].name;
            if ui
                .button(format!("Create Line {} ∧ {}", p1_name, p2_name))
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
            ui.add_space(4.0);
            ui.label("Lines:");
            let mut lines_to_remove = Vec::new();

            for (idx, line_data) in self.derived_lines.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    if ui.checkbox(&mut line_data.selected, "").changed() {
                        // Update meet selection
                    }
                    let l = &line_data.line;
                    ui.label(format!(
                        "L{}: {:.2}x + {:.2}y + {:.2} = 0",
                        idx + 1,
                        l.normal_x(),
                        l.normal_y(),
                        l.dist()
                    ));
                    if ui.button("🗑").clicked() {
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

        // Meet operation (antiwedge)
        section_separator(ui, Some("Meet (∨) - Intersection"));
        ui.label(format!("Selected: {} lines", self.meet_selection.len()));

        if self.meet_selection.len() == 2 {
            if ui
                .button(format!(
                    "Find Intersection L{} ∨ L{}",
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
            ui.add_space(4.0);
            ui.label("Intersections:");
            for (idx, derived) in self.derived_points.iter().enumerate() {
                if let Some((x, y)) = derived.point.to_cartesian() {
                    ui.label(format!("P{}: ({:.2}, {:.2})", idx + 1, x, y));
                } else {
                    ui.label(format!("P{}: ideal point (parallel lines)", idx + 1));
                }
            }
        }

        // Motor controls
        section_separator(ui, Some("Motor Transform"));
        ui.checkbox(&mut self.apply_motor, "Apply motor transformation");

        if self.apply_motor {
            angle_slider_range(ui, "Rotation θ", &mut self.motor_rotation, -360.0, 360.0);

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
                    ("e₁", self.motor.ty() as f32),
                    ("e₂", self.motor.tx() as f32),
                    ("e₃", self.motor.r() as f32),
                    ("e₁₂₃", self.motor.ps() as f32),
                ],
            );

            ui.add_space(4.0);
            animation_controls(ui, &mut self.animation);
            progress_slider(ui, &mut self.animation);
        }

        // Display options
        section_separator(ui, Some("Display Options"));
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_coordinates, "Show coordinates in list");
        ui.checkbox(&mut self.show_normals, "Show line normals");

        // Update derived geometry when points change
        self.update_derived_lines();
        self.update_derived_points();
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label(format!("{} points", self.points.len()));
            ui.separator();
            ui.label(format!("{} lines", self.derived_lines.len()));
            ui.separator();
            ui.label(format!("{} intersections", self.derived_points.len()));
            ui.separator();
            if self.apply_motor {
                ui.label(format!(
                    "Motor: θ={:.1}° t=({:.1}, {:.1})",
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
• The JOIN of two points (line through them) is their wedge product (∧)
• The MEET of two lines (intersection point) is their antiwedge product (∨)
• MOTORS encode rotation and translation as a single algebraic element

This demo lets you explore these operations interactively.",

    math_background: "\
POINTS are grade-1 elements in homogeneous coordinates:
    P = x·e₁ + y·e₂ + w·e₀

For a finite point at (x, y), we set w = 1.

LINES are grade-2 elements representing ax + by + c = 0:
    L = c·e₁₂ + a·e₁₀ + b·e₂₀

JOIN (∧) - Line through two points:
    L = P₁ ∧ P₂

MEET (∨) - Intersection of two lines:
    P = L₁ ∨ L₂

If lines are parallel, the result is an IDEAL POINT (w = 0), \
representing the direction at infinity.

MOTORS are elements of the odd subalgebra (grades 1 and 3):
    M = ty·e₁ + tx·e₂ + r·e₃ + ps·e₁₂₃

They transform geometry via the antisandwich product:
    P' = M⁻¹PM  (point transformation)
    L' = M⁻¹LM  (line transformation)",

    how_to_use: "\
• ADD POINTS: Click 'Add Point at Origin' and drag coordinates
• JOIN OPERATION: Select 2 points (checkboxes), click 'Create Line'
• MEET OPERATION: Select 2 lines (checkboxes), click 'Find Intersection'
• MOTOR TRANSFORM: Enable 'Apply motor transformation', adjust sliders
• Drag point coordinates to see lines update in real-time
• Enable 'Show line normals' to visualize line orientations",

    key_concepts: "\
• Homogeneous coordinates: P = (x, y, w) with w=1 for finite points
• Ideal points (w=0) represent directions at infinity
• Wedge product (∧) computes JOIN: line through two points
• Antiwedge product (∨) computes MEET: intersection of two lines
• Motors compose rotation and translation into one operation
• Antisandwich product transforms geometry while preserving incidence",

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

fn main() -> eframe::Result<()> {
    run_app::<Projective2Demo>()
}
