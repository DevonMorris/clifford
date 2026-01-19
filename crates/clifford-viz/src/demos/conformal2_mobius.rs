//! Möbius Transformations - 2D Conformal GA Visualization
//!
//! This demo demonstrates the composition of conformal transformations
//! using CGA motors: translation, rotation, dilation, and inversion.
//!
//! All Möbius transformations in 2D can be composed from these primitives.
//! CGA motors elegantly represent these transformations via the sandwich product.

use crate::common::prelude::*;
use clifford::ops::Transform;
use clifford::specialized::conformal::dim2::{Circle, Motor, RoundPoint};
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
const VIEWPORT_BOUNDS: f64 = 6.0;

/// Number of segments for shape rendering.
/// Higher values give smoother curves after conformal transformations like inversion.
const SHAPE_SEGMENTS: usize = 256;

/// Minimum scale factor for dilation.
const MIN_SCALE: f64 = 0.1;

/// Maximum scale factor for dilation.
const MAX_SCALE: f64 = 5.0;

/// Type of transformation in the chain.
#[derive(Clone, Copy, Debug, PartialEq)]
enum TransformType {
    /// Translation by a vector (dx, dy).
    Translation,
    /// Rotation around the origin by an angle.
    Rotation,
    /// Uniform scaling (dilation) from the origin.
    Dilation,
    /// Inversion through a circle.
    Inversion,
}

impl TransformType {
    /// Returns the display name of this transformation type.
    fn name(&self) -> &'static str {
        match self {
            TransformType::Translation => "Translation",
            TransformType::Rotation => "Rotation",
            TransformType::Dilation => "Dilation",
            TransformType::Inversion => "Inversion",
        }
    }

    /// Returns a short symbol for this transformation type.
    fn short_name(&self) -> &'static str {
        match self {
            TransformType::Translation => "T",
            TransformType::Rotation => "R",
            TransformType::Dilation => "D",
            TransformType::Inversion => "I",
        }
    }
}

/// A transformation in the chain.
#[derive(Clone)]
struct TransformItem {
    /// Type of transformation.
    transform_type: TransformType,
    /// Translation vector (dx, dy).
    translation: (f64, f64),
    /// Rotation angle in radians.
    rotation_angle: f64,
    /// Dilation scale factor.
    scale: f64,
    /// Inversion circle (cx, cy, r).
    inversion_circle: (f64, f64, f64),
    /// Whether this transformation is enabled.
    enabled: bool,
}

impl TransformItem {
    /// Creates a new translation transformation.
    fn new_translation(dx: f64, dy: f64) -> Self {
        Self {
            transform_type: TransformType::Translation,
            translation: (dx, dy),
            rotation_angle: 0.0,
            scale: 1.0,
            inversion_circle: (0.0, 0.0, 2.0),
            enabled: true,
        }
    }

    /// Creates a new rotation transformation.
    fn new_rotation(angle: f64) -> Self {
        Self {
            transform_type: TransformType::Rotation,
            translation: (0.0, 0.0),
            rotation_angle: angle,
            scale: 1.0,
            inversion_circle: (0.0, 0.0, 2.0),
            enabled: true,
        }
    }

    /// Creates a new dilation (scaling) transformation.
    fn new_dilation(scale: f64) -> Self {
        Self {
            transform_type: TransformType::Dilation,
            translation: (0.0, 0.0),
            rotation_angle: 0.0,
            scale,
            inversion_circle: (0.0, 0.0, 2.0),
            enabled: true,
        }
    }

    /// Creates a new inversion transformation through a circle.
    fn new_inversion(cx: f64, cy: f64, r: f64) -> Self {
        Self {
            transform_type: TransformType::Inversion,
            translation: (0.0, 0.0),
            rotation_angle: 0.0,
            scale: 1.0,
            inversion_circle: (cx, cy, r),
            enabled: true,
        }
    }

    /// Applies this transformation to a point.
    fn apply(&self, point: &RoundPoint<f64>) -> RoundPoint<f64> {
        if !self.enabled {
            return *point;
        }

        match self.transform_type {
            TransformType::Translation => {
                let motor = Motor::from_translation(self.translation.0, self.translation.1);
                motor.transform(point)
            }
            TransformType::Rotation => {
                let motor = Motor::from_rotation(self.rotation_angle);
                motor.transform(point)
            }
            TransformType::Dilation => {
                // Note: Motor::from_dilation has inverted scale semantics,
                // so we use 1/scale to get expected behavior (scale > 1 grows)
                let motor = Motor::from_dilation(1.0 / self.scale);
                motor.transform(point)
            }
            TransformType::Inversion => {
                let (cx, cy, r) = self.inversion_circle;
                let inv_circle = Circle::from_center_radius(cx, cy, r);
                inv_circle.transform(point)
            }
        }
    }

    /// Returns a description of this transformation.
    fn description(&self) -> String {
        if !self.enabled {
            return format!("[{}] (disabled)", self.transform_type.short_name());
        }

        match self.transform_type {
            TransformType::Translation => {
                format!("T({:.1}, {:.1})", self.translation.0, self.translation.1)
            }
            TransformType::Rotation => {
                format!("R({:.0} deg)", self.rotation_angle.to_degrees())
            }
            TransformType::Dilation => {
                format!("D({:.2})", self.scale)
            }
            TransformType::Inversion => {
                let (cx, cy, r) = self.inversion_circle;
                format!("I(c=({:.1},{:.1}), r={:.1})", cx, cy, r)
            }
        }
    }
}

/// Shape type for transformation.
#[derive(Clone, Copy, Debug, PartialEq)]
enum ShapeType {
    /// A circular shape.
    Circle,
    /// A square shape.
    Square,
    /// A triangular shape.
    Triangle,
}

impl ShapeType {
    /// Returns the display name of this shape type.
    fn name(&self) -> &'static str {
        match self {
            ShapeType::Circle => "Circle",
            ShapeType::Square => "Square",
            ShapeType::Triangle => "Triangle",
        }
    }
}

/// Demo state for Möbius Transformations visualization.
pub struct Conformal2MobiusDemo {
    /// Chain of transformations to apply.
    transform_chain: Vec<TransformItem>,
    /// The shape to transform.
    shape_type: ShapeType,
    /// Shape center.
    shape_center: (f64, f64),
    /// Shape radius/size.
    shape_size: f64,
    /// Animation state.
    animation: Animation,
    /// Whether to show the grid.
    show_grid: bool,
    /// Whether to show the original shape.
    show_original: bool,
    /// Whether to show intermediate shapes.
    show_intermediate: bool,
    /// Whether to show inversion circles.
    show_inversion_circles: bool,
    /// Index of transformation to add (for UI).
    add_transform_type: TransformType,
}

impl Default for Conformal2MobiusDemo {
    fn default() -> Self {
        Self {
            transform_chain: vec![
                TransformItem::new_translation(1.0, 0.5),
                TransformItem::new_rotation(std::f64::consts::FRAC_PI_4),
            ],
            shape_type: ShapeType::Circle,
            shape_center: (-2.0, 0.0),
            shape_size: 1.0,
            animation: Animation::default(),
            show_grid: true,
            show_original: true,
            show_intermediate: false,
            show_inversion_circles: true,
            add_transform_type: TransformType::Translation,
        }
    }
}

impl Conformal2MobiusDemo {
    /// Generates points on the shape boundary.
    fn shape_points(&self) -> Vec<(f64, f64)> {
        let (cx, cy) = self.shape_center;
        let r = self.shape_size;

        match self.shape_type {
            ShapeType::Circle => (0..SHAPE_SEGMENTS)
                .map(|i| {
                    let angle = 2.0 * std::f64::consts::PI * (i as f64) / (SHAPE_SEGMENTS as f64);
                    (cx + r * angle.cos(), cy + r * angle.sin())
                })
                .collect(),
            ShapeType::Square => {
                let half = r;
                let points_per_side = SHAPE_SEGMENTS / 4;
                let mut points = Vec::with_capacity(SHAPE_SEGMENTS);

                // Top side
                for i in 0..points_per_side {
                    let t = (i as f64) / (points_per_side as f64);
                    points.push((cx - half + 2.0 * half * t, cy + half));
                }
                // Right side
                for i in 0..points_per_side {
                    let t = (i as f64) / (points_per_side as f64);
                    points.push((cx + half, cy + half - 2.0 * half * t));
                }
                // Bottom side
                for i in 0..points_per_side {
                    let t = (i as f64) / (points_per_side as f64);
                    points.push((cx + half - 2.0 * half * t, cy - half));
                }
                // Left side
                for i in 0..points_per_side {
                    let t = (i as f64) / (points_per_side as f64);
                    points.push((cx - half, cy - half + 2.0 * half * t));
                }

                points
            }
            ShapeType::Triangle => {
                let points_per_side = SHAPE_SEGMENTS / 3;
                let mut points = Vec::with_capacity(SHAPE_SEGMENTS);

                // Vertices of equilateral triangle
                let v0 = (cx, cy + r);
                let v1 = (cx - r * 0.866, cy - r * 0.5);
                let v2 = (cx + r * 0.866, cy - r * 0.5);

                // Side 0 -> 1
                for i in 0..points_per_side {
                    let t = (i as f64) / (points_per_side as f64);
                    points.push((v0.0 + (v1.0 - v0.0) * t, v0.1 + (v1.1 - v0.1) * t));
                }
                // Side 1 -> 2
                for i in 0..points_per_side {
                    let t = (i as f64) / (points_per_side as f64);
                    points.push((v1.0 + (v2.0 - v1.0) * t, v1.1 + (v2.1 - v1.1) * t));
                }
                // Side 2 -> 0
                for i in 0..points_per_side {
                    let t = (i as f64) / (points_per_side as f64);
                    points.push((v2.0 + (v0.0 - v2.0) * t, v2.1 + (v0.1 - v2.1) * t));
                }

                points
            }
        }
    }

    /// Applies all enabled transformations up to a given index (for intermediate display).
    fn apply_chain_up_to(&self, points: &[(f64, f64)], up_to: usize, t: f64) -> Vec<(f64, f64)> {
        let mut result = Vec::with_capacity(points.len());

        for &(x, y) in points {
            let mut pt = RoundPoint::from_euclidean(x, y);

            for (idx, transform) in self.transform_chain.iter().enumerate() {
                if idx > up_to {
                    break;
                }
                if !transform.enabled {
                    continue;
                }

                // Apply with interpolation for the last transform
                if idx == up_to && t < 1.0 {
                    pt = self.apply_interpolated(transform, &pt, t);
                } else {
                    pt = transform.apply(&pt);
                }
            }

            if let Some((rx, ry)) = pt.to_euclidean() {
                result.push((rx, ry));
            }
        }

        result
    }

    /// Applies a transformation with interpolation (0 = identity, 1 = full).
    fn apply_interpolated(
        &self,
        transform: &TransformItem,
        point: &RoundPoint<f64>,
        t: f64,
    ) -> RoundPoint<f64> {
        if !transform.enabled || t <= 0.0 {
            return *point;
        }
        if t >= 1.0 {
            return transform.apply(point);
        }

        match transform.transform_type {
            TransformType::Translation => {
                let (dx, dy) = transform.translation;
                let motor = Motor::from_translation(dx * t, dy * t);
                motor.transform(point)
            }
            TransformType::Rotation => {
                let motor = Motor::from_rotation(transform.rotation_angle * t);
                motor.transform(point)
            }
            TransformType::Dilation => {
                // Interpolate log scale (inverted for correct semantics)
                let log_scale = (1.0 / transform.scale).ln() * t;
                let motor = Motor::from_dilation(log_scale.exp());
                motor.transform(point)
            }
            TransformType::Inversion => {
                // Inversion doesn't interpolate smoothly, apply fully after t > 0.5
                if t > 0.5 {
                    transform.apply(point)
                } else {
                    *point
                }
            }
        }
    }

    /// Applies all enabled transformations to the shape points.
    fn transform_shape(&self, points: &[(f64, f64)], t: f64) -> Vec<(f64, f64)> {
        if self.transform_chain.is_empty() {
            return points.to_vec();
        }
        self.apply_chain_up_to(points, self.transform_chain.len() - 1, t)
    }

    /// Gets the composed transformation description.
    fn composed_description(&self) -> String {
        let enabled: Vec<_> = self
            .transform_chain
            .iter()
            .filter(|t| t.enabled)
            .map(|t| t.description())
            .collect();

        if enabled.is_empty() {
            "Identity".to_string()
        } else {
            enabled.join(" o ")
        }
    }
}

impl VisualizationApp for Conformal2MobiusDemo {
    fn name(&self) -> &'static str {
        "Conformal 2D - Mobius Transformations"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let original_points = self.shape_points();

        // Compute animation progress
        let t: f64 = if self.animation.playing {
            f64::from(self.animation.progress())
        } else {
            1.0
        };

        // Transform the shape
        let transformed_points = self.transform_shape(&original_points, t);

        Plot::new("conformal2_mobius_plot")
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
                    for l in grid_2d(&ctx, VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(l);
                    }
                    for axis in axes_2d(&ctx, VIEWPORT_BOUNDS) {
                        plot_ui.line(axis);
                    }
                }

                // Draw inversion circles
                if self.show_inversion_circles {
                    for (idx, transform) in self.transform_chain.iter().enumerate() {
                        if transform.transform_type == TransformType::Inversion && transform.enabled
                        {
                            let (cx, cy, r) = transform.inversion_circle;
                            let inv_circle = circle_2d(cx, cy, r, with_alpha(motor(&ctx), 100), 64)
                                .name(format!("Inv Circle {}", idx + 1));
                            plot_ui.line(inv_circle);

                            // Draw center
                            plot_ui.points(
                                Points::new(vec![[cx, cy]])
                                    .color(with_alpha(motor(&ctx), 150))
                                    .radius(4.0)
                                    .filled(true),
                            );
                        }
                    }
                }

                // Draw original shape
                if self.show_original && !original_points.is_empty() {
                    let pts: Vec<[f64; 2]> = original_points
                        .iter()
                        .chain(std::iter::once(&original_points[0]))
                        .map(|&(x, y)| [x, y])
                        .collect();
                    let original_line = egui_plot::Line::new(pts)
                        .color(with_alpha(point(&ctx), 100))
                        .width(1.5)
                        .name("Original");
                    plot_ui.line(original_line);
                }

                // Draw intermediate shapes
                if self.show_intermediate {
                    for idx in 0..self.transform_chain.len().saturating_sub(1) {
                        let intermediate = self.apply_chain_up_to(&original_points, idx, 1.0);
                        if !intermediate.is_empty() {
                            let pts: Vec<[f64; 2]> = intermediate
                                .iter()
                                .chain(std::iter::once(&intermediate[0]))
                                .map(|&(x, y)| [x, y])
                                .collect();
                            let alpha = 50 + (idx * 30).min(100);
                            let line = egui_plot::Line::new(pts)
                                .color(with_alpha(active(&ctx), alpha as u8))
                                .width(1.0)
                                .name(format!("After {}", idx + 1));
                            plot_ui.line(line);
                        }
                    }
                }

                // Draw transformed shape
                if !transformed_points.is_empty() {
                    let pts: Vec<[f64; 2]> = transformed_points
                        .iter()
                        .chain(std::iter::once(&transformed_points[0]))
                        .map(|&(x, y)| [x, y])
                        .collect();
                    let transformed_line = egui_plot::Line::new(pts)
                        .color(active(&ctx))
                        .width(2.5)
                        .name("Transformed");
                    plot_ui.line(transformed_line);
                }

                // Draw origin marker
                plot_ui.points(
                    Points::new(vec![[0.0, 0.0]])
                        .color(with_alpha(grid(&ctx), 200))
                        .radius(3.0)
                        .filled(true)
                        .name("Origin"),
                );
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        // === Shape Configuration ===
        group_header(ui, "Shape");

        ui.horizontal(|ui| {
            ui.label("Type:");
            egui::ComboBox::from_id_salt("shape_type")
                .selected_text(self.shape_type.name())
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.shape_type, ShapeType::Circle, "Circle");
                    ui.selectable_value(&mut self.shape_type, ShapeType::Square, "Square");
                    ui.selectable_value(&mut self.shape_type, ShapeType::Triangle, "Triangle");
                });
        });

        ui.horizontal(|ui| {
            ui.label("Center:");
            let mut cx = self.shape_center.0 as f32;
            let mut cy = self.shape_center.1 as f32;
            if ui
                .add(egui::DragValue::new(&mut cx).speed(0.1).prefix("x: "))
                .changed()
                || ui
                    .add(egui::DragValue::new(&mut cy).speed(0.1).prefix("y: "))
                    .changed()
            {
                self.shape_center = (f64::from(cx), f64::from(cy));
            }
        });

        ui.horizontal(|ui| {
            ui.label("Size:");
            let mut size = self.shape_size as f32;
            if ui.add(egui::Slider::new(&mut size, 0.2..=3.0)).changed() {
                self.shape_size = f64::from(size);
            }
        });

        // === Transformation Chain ===
        section_separator(ui, Some("Transformation Chain"));

        ui.label(format!("Composed: {}", self.composed_description()));
        ui.add_space(spacing::XS);

        // List transformations
        let mut to_remove: Option<usize> = None;
        let mut to_move_up: Option<usize> = None;
        let mut to_move_down: Option<usize> = None;
        let chain_len = self.transform_chain.len();

        for (idx, transform) in self.transform_chain.iter_mut().enumerate() {
            ui.push_id(idx, |ui| {
                ui.group(|ui| {
                    ui.horizontal(|ui| {
                        ui.checkbox(&mut transform.enabled, "");
                        ui.label(format!("{}. {}", idx + 1, transform.transform_type.name()));

                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui.small_button("X").clicked() {
                                to_remove = Some(idx);
                            }
                            if idx < chain_len - 1 && ui.small_button("v").clicked() {
                                to_move_down = Some(idx);
                            }
                            if idx > 0 && ui.small_button("^").clicked() {
                                to_move_up = Some(idx);
                            }
                        });
                    });

                    // Transform-specific controls
                    match transform.transform_type {
                        TransformType::Translation => {
                            ui.horizontal(|ui| {
                                let mut dx = transform.translation.0 as f32;
                                let mut dy = transform.translation.1 as f32;
                                ui.label("  ");
                                if ui
                                    .add(egui::DragValue::new(&mut dx).speed(0.1).prefix("dx: "))
                                    .changed()
                                {
                                    transform.translation.0 = f64::from(dx);
                                }
                                if ui
                                    .add(egui::DragValue::new(&mut dy).speed(0.1).prefix("dy: "))
                                    .changed()
                                {
                                    transform.translation.1 = f64::from(dy);
                                }
                            });
                        }
                        TransformType::Rotation => {
                            ui.horizontal(|ui| {
                                ui.label("  ");
                                let mut degrees = transform.rotation_angle.to_degrees() as f32;
                                if ui
                                    .add(
                                        egui::Slider::new(&mut degrees, -180.0..=180.0)
                                            .suffix(" deg"),
                                    )
                                    .changed()
                                {
                                    transform.rotation_angle = (degrees as f64).to_radians();
                                }
                            });
                        }
                        TransformType::Dilation => {
                            ui.horizontal(|ui| {
                                ui.label("  Scale:");
                                let mut scale = transform.scale as f32;
                                if ui
                                    .add(
                                        egui::Slider::new(
                                            &mut scale,
                                            MIN_SCALE as f32..=MAX_SCALE as f32,
                                        )
                                        .logarithmic(true),
                                    )
                                    .changed()
                                {
                                    transform.scale = f64::from(scale);
                                }
                            });
                        }
                        TransformType::Inversion => {
                            let (mut cx, mut cy, mut r) = transform.inversion_circle;
                            let mut cx_f = cx as f32;
                            let mut cy_f = cy as f32;
                            let mut r_f = r as f32;

                            ui.horizontal(|ui| {
                                ui.label("  Center:");
                                if ui
                                    .add(egui::DragValue::new(&mut cx_f).speed(0.1).prefix("x: "))
                                    .changed()
                                {
                                    cx = f64::from(cx_f);
                                }
                                if ui
                                    .add(egui::DragValue::new(&mut cy_f).speed(0.1).prefix("y: "))
                                    .changed()
                                {
                                    cy = f64::from(cy_f);
                                }
                            });
                            ui.horizontal(|ui| {
                                ui.label("  Radius:");
                                if ui.add(egui::Slider::new(&mut r_f, 0.5..=4.0)).changed() {
                                    r = f64::from(r_f);
                                }
                            });

                            transform.inversion_circle = (cx, cy, r);
                        }
                    }
                });
            });
        }

        // Apply modifications
        if let Some(idx) = to_remove {
            self.transform_chain.remove(idx);
        }
        if let Some(idx) = to_move_up {
            self.transform_chain.swap(idx, idx - 1);
        }
        if let Some(idx) = to_move_down {
            self.transform_chain.swap(idx, idx + 1);
        }

        // Add new transformation
        ui.add_space(spacing::XS);
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("add_transform")
                .selected_text(self.add_transform_type.name())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.add_transform_type,
                        TransformType::Translation,
                        "Translation",
                    );
                    ui.selectable_value(
                        &mut self.add_transform_type,
                        TransformType::Rotation,
                        "Rotation",
                    );
                    ui.selectable_value(
                        &mut self.add_transform_type,
                        TransformType::Dilation,
                        "Dilation",
                    );
                    ui.selectable_value(
                        &mut self.add_transform_type,
                        TransformType::Inversion,
                        "Inversion",
                    );
                });

            if ui.button("+ Add").clicked() {
                let new_transform = match self.add_transform_type {
                    TransformType::Translation => TransformItem::new_translation(1.0, 0.0),
                    TransformType::Rotation => {
                        TransformItem::new_rotation(std::f64::consts::FRAC_PI_4)
                    }
                    TransformType::Dilation => TransformItem::new_dilation(1.5),
                    TransformType::Inversion => TransformItem::new_inversion(0.0, 0.0, 2.0),
                };
                self.transform_chain.push(new_transform);
            }
        });

        // === Animation ===
        section_separator(ui, Some("Animation"));
        animation_controls(ui, &mut self.animation);

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_original, "Show original shape");
        ui.checkbox(&mut self.show_intermediate, "Show intermediate shapes");
        ui.checkbox(&mut self.show_inversion_circles, "Show inversion circles");

        // === Presets ===
        section_separator(ui, Some("Presets"));

        if ui.button("Translation + Rotation").clicked() {
            self.transform_chain = vec![
                TransformItem::new_translation(2.0, 1.0),
                TransformItem::new_rotation(std::f64::consts::FRAC_PI_3),
            ];
            self.shape_center = (-2.0, 0.0);
        }

        if ui.button("Rotation + Dilation").clicked() {
            self.transform_chain = vec![
                TransformItem::new_rotation(std::f64::consts::FRAC_PI_4),
                TransformItem::new_dilation(1.5),
            ];
            self.shape_center = (1.5, 0.0);
        }

        if ui.button("Circle Inversion").clicked() {
            self.transform_chain = vec![TransformItem::new_inversion(0.0, 0.0, 2.0)];
            self.shape_center = (3.0, 0.0);
            self.shape_size = 1.0;
        }

        if ui.button("Double Inversion").clicked() {
            self.transform_chain = vec![
                TransformItem::new_inversion(0.0, 0.0, 2.0),
                TransformItem::new_inversion(0.0, 0.0, 2.0),
            ];
            self.shape_center = (3.0, 0.0);
            self.shape_size = 1.0;
        }

        if ui.button("Full Möbius").clicked() {
            self.transform_chain = vec![
                TransformItem::new_translation(1.0, 0.0),
                TransformItem::new_rotation(std::f64::consts::FRAC_PI_6),
                TransformItem::new_dilation(1.3),
                TransformItem::new_inversion(0.0, 0.0, 3.0),
            ];
            self.shape_center = (-3.0, -1.0);
        }

        if ui.button("Clear All").clicked() {
            self.transform_chain.clear();
        }

        // === Math Info ===
        section_separator(ui, Some("Mathematics"));
        info_box(
            ui,
            "Mobius transformations:\n\
             - Translation: T = 1 - (1/2)*d*e_inf\n\
             - Rotation: R = cos(a/2) - sin(a/2)*e12\n\
             - Dilation: D = cosh(ln k/2) + sinh(ln k/2)*e0_inf\n\
             - Inversion: P' = C*P*C^-1",
        );

        if !self.transform_chain.is_empty() {
            ui.add_space(spacing::XS);
            ui.colored_label(
                active(&ctx),
                "Double inversion = identity\n(two inversions through same circle cancel)",
            );
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            ui.label(format!("{}: ", self.shape_type.name()));
            ui.colored_label(point(&ctx), "original");
            ui.label(" -> ");
            ui.colored_label(active(&ctx), "transformed");
            ui.separator();
            ui.label(format!(
                "{} transforms",
                self.transform_chain.iter().filter(|t| t.enabled).count()
            ));
            if self.animation.playing {
                ui.separator();
                ui.label(format!("t = {:.2}", self.animation.progress()));
            }
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(CONFORMAL2_MOBIUS_EDUCATION)
    }
}

/// Educational content for the Mobius Transformations visualization.
const CONFORMAL2_MOBIUS_EDUCATION: EducationalContent = EducationalContent {
    title: "Mobius Transformations in CGA",

    overview: "\
Mobius transformations are conformal (angle-preserving) mappings of the plane \
that generalize translations, rotations, scalings, and inversions. Every Mobius \
transformation can be written as a composition of these four primitive operations.

In Conformal Geometric Algebra, these transformations are elegantly represented as \
\"motors\" (even-grade multivectors) that act via the sandwich product: X' = M*X*M^-1",

    math_background: "\
TRANSFORMATION MOTORS in CGA:

TRANSLATION by vector d = (dx, dy):
    T = 1 - (1/2)*d*e_inf
    = 1 - (1/2)*dx*e1*e_inf - (1/2)*dy*e2*e_inf

ROTATION by angle a around origin:
    R = cos(a/2) - sin(a/2)*e1*e2

DILATION (scaling) by factor k around origin:
    D = cosh(ln k / 2) + sinh(ln k / 2)*e0*e_inf

INVERSION through circle C:
    P' = C*P*C^-1

COMPOSITION:
Motors compose by multiplication: M_total = M_n * ... * M_2 * M_1
(applied right-to-left)

KEY PROPERTY:
Two inversions through the same circle = identity (they cancel out).",

    how_to_use: "\
- Add transformations using the dropdown and + Add button
- REORDER transformations with ^/v arrows
- TOGGLE transformations on/off with checkboxes
- Adjust parameters for each transformation
- Use PRESETS for common combinations
- Play ANIMATION to see transformation smoothly applied
- Enable INTERMEDIATE SHAPES to see step-by-step",

    key_concepts: "\
- Mobius = composition of T, R, D, and I
- All Mobius transforms preserve angles (conformal)
- Circles map to circles (including lines as infinite-radius circles)
- Two inversions through same circle = identity
- Order of composition matters (not commutative)
- CGA unifies all transforms as sandwich products",

    resources: &[
        (
            "Mobius Transformations (Wikipedia)",
            "https://en.wikipedia.org/wiki/M%C3%B6bius_transformation",
        ),
        (
            "Conformal Geometric Algebra Wiki",
            "https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page",
        ),
    ],
};
