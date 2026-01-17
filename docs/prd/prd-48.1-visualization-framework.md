# PRD-48.1: Visualization Framework Setup

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-49 (workspace restructure)
**Goal**: Establish common infrastructure for all visualization demos

## Overview

Create the shared utilities, rendering primitives, and egui boilerplate that all visualization demos will use.

## Crate Setup

Create a new crate `clifford-viz` in the workspace:

```toml
# crates/clifford-viz/Cargo.toml
[package]
name = "clifford-viz"
version = "0.1.0"
edition = "2024"
description = "Interactive visualizations for clifford geometric algebra"
license = "MIT OR Apache-2.0"

[dependencies]
clifford = { path = "../clifford" }
eframe = "0.30"
egui = "0.30"
egui_plot = "0.30"

[dev-dependencies]
image = "0.25"  # For visual testing
```

Add to workspace in root `Cargo.toml`:
```toml
[workspace]
members = [
    ".",
    "crates/clifford-codegen",
    "crates/clifford-viz",  # Add this
]
```

## Directory Structure

```
crates/clifford-viz/
  src/
    lib.rs            # Public API, re-exports
    common/
      mod.rs          # Re-exports
      app.rs          # Base app trait and runner
      colors.rs       # Color palette
      grid.rs         # 2D/3D grid rendering
      shapes.rs       # Primitive shape drawing
      widgets.rs      # Reusable UI components
      animation.rs    # Animation timing utilities
    euclidean.rs      # Reusable Euclidean widgets
    projective.rs     # Reusable PGA widgets
    conformal.rs      # Reusable CGA widgets
    quaternion.rs     # Reusable quaternion widgets
    spacetime.rs      # Reusable Minkowski widgets
  examples/           # Interactive demos
  tests/visual/       # Visual regression tests
```

## Components

### 1. Base Application Trait

```rust
// examples/visualization/common/app.rs

/// Base trait for visualization demos
pub trait VisualizationApp {
    /// Display name shown in window title
    fn name(&self) -> &'static str;

    /// Update logic (called every frame)
    fn update(&mut self, dt: f32);

    /// Render the main visualization area
    fn render(&self, ui: &mut egui::Ui);

    /// Render the control panel (sliders, buttons, etc.)
    fn controls(&mut self, ui: &mut egui::Ui);

    /// Optional: info panel with explanatory text
    fn info(&self, ui: &mut egui::Ui) {
        // Default: no info panel
    }
}

/// Run a visualization app
pub fn run_app<T: VisualizationApp + Default + 'static>() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1024.0, 768.0]),
        ..Default::default()
    };
    eframe::run_native(
        T::default().name(),
        options,
        Box::new(|_cc| Ok(Box::new(AppWrapper::<T>::default()))),
    )
}
```

### 2. Color Palette

```rust
// examples/visualization/common/colors.rs

use egui::Color32;

/// Consistent color palette across all demos
pub mod palette {
    use super::*;

    // Primary geometric objects
    pub const POINT: Color32 = Color32::from_rgb(66, 133, 244);    // Blue
    pub const LINE: Color32 = Color32::from_rgb(234, 67, 53);      // Red
    pub const PLANE: Color32 = Color32::from_rgb(52, 168, 83);     // Green
    pub const CIRCLE: Color32 = Color32::from_rgb(251, 188, 5);    // Yellow
    pub const SPHERE: Color32 = Color32::from_rgb(154, 160, 166);  // Gray

    // Transformations
    pub const ROTOR: Color32 = Color32::from_rgb(171, 71, 188);    // Purple
    pub const MOTOR: Color32 = Color32::from_rgb(0, 172, 193);     // Cyan

    // Axes
    pub const X_AXIS: Color32 = Color32::from_rgb(244, 67, 54);    // Red
    pub const Y_AXIS: Color32 = Color32::from_rgb(76, 175, 80);    // Green
    pub const Z_AXIS: Color32 = Color32::from_rgb(33, 150, 243);   // Blue
    pub const T_AXIS: Color32 = Color32::from_rgb(255, 193, 7);    // Amber (time)

    // UI states
    pub const SELECTED: Color32 = Color32::from_rgb(255, 235, 59); // Highlight
    pub const HOVERED: Color32 = Color32::from_rgb(255, 255, 255); // White
    pub const GRID: Color32 = Color32::from_rgb(66, 66, 66);       // Dark gray
    pub const BACKGROUND: Color32 = Color32::from_rgb(30, 30, 30); // Near black
}

/// Generate a color from hue (0-1) for domain coloring
pub fn hue_to_color(hue: f32, saturation: f32, lightness: f32) -> Color32 {
    // HSL to RGB conversion
    // ...
}
```

### 3. Grid Rendering

```rust
// examples/visualization/common/grid.rs

use egui_plot::{Line, PlotPoints};

/// Draw a 2D Cartesian grid
pub fn grid_2d(bounds: f64, step: f64) -> Vec<Line> {
    let mut lines = Vec::new();
    let n = (bounds / step).ceil() as i32;

    for i in -n..=n {
        let v = i as f64 * step;
        // Vertical line
        lines.push(Line::new(PlotPoints::new(vec![[v, -bounds], [v, bounds]]))
            .color(palette::GRID)
            .width(if i == 0 { 2.0 } else { 1.0 }));
        // Horizontal line
        lines.push(Line::new(PlotPoints::new(vec![[-bounds, v], [bounds, v]]))
            .color(palette::GRID)
            .width(if i == 0 { 2.0 } else { 1.0 }));
    }
    lines
}

/// Draw coordinate axes with labels
pub fn axes_2d(bounds: f64) -> Vec<Line> {
    vec![
        Line::new(PlotPoints::new(vec![[-bounds, 0.0], [bounds, 0.0]]))
            .color(palette::X_AXIS).width(2.0).name("x"),
        Line::new(PlotPoints::new(vec![[0.0, -bounds], [0.0, bounds]]))
            .color(palette::Y_AXIS).width(2.0).name("y"),
    ]
}
```

### 4. Shape Primitives

```rust
// examples/visualization/common/shapes.rs

use egui_plot::{Points, Line, Polygon, PlotPoints};

/// Draw a point marker
pub fn point_marker(x: f64, y: f64, color: Color32) -> Points {
    Points::new(vec![[x, y]])
        .color(color)
        .radius(6.0)
        .filled(true)
}

/// Draw an arrow from origin to (dx, dy)
pub fn arrow_2d(ox: f64, oy: f64, dx: f64, dy: f64, color: Color32) -> Vec<Line> {
    let len = (dx * dx + dy * dy).sqrt();
    let head_size = len * 0.1;
    let angle = dy.atan2(dx);

    // Main shaft
    let shaft = Line::new(PlotPoints::new(vec![[ox, oy], [ox + dx, oy + dy]]))
        .color(color)
        .width(2.0);

    // Arrowhead
    let head_angle = 0.4; // radians
    let h1 = Line::new(PlotPoints::new(vec![
        [ox + dx, oy + dy],
        [ox + dx - head_size * (angle + head_angle).cos(),
         oy + dy - head_size * (angle + head_angle).sin()],
    ])).color(color).width(2.0);
    let h2 = Line::new(PlotPoints::new(vec![
        [ox + dx, oy + dy],
        [ox + dx - head_size * (angle - head_angle).cos(),
         oy + dy - head_size * (angle - head_angle).sin()],
    ])).color(color).width(2.0);

    vec![shaft, h1, h2]
}

/// Draw a circle (as line segments)
pub fn circle_2d(cx: f64, cy: f64, radius: f64, color: Color32, segments: usize) -> Line {
    let points: Vec<[f64; 2]> = (0..=segments)
        .map(|i| {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / segments as f64;
            [cx + radius * angle.cos(), cy + radius * angle.sin()]
        })
        .collect();
    Line::new(PlotPoints::new(points)).color(color).width(2.0)
}

/// Draw a line segment
pub fn line_segment(x1: f64, y1: f64, x2: f64, y2: f64, color: Color32) -> Line {
    Line::new(PlotPoints::new(vec![[x1, y1], [x2, y2]]))
        .color(color)
        .width(2.0)
}

/// Draw an infinite line (clipped to bounds) given point and direction
pub fn infinite_line_2d(px: f64, py: f64, dx: f64, dy: f64, bounds: f64, color: Color32) -> Line {
    // Extend in both directions to bounds
    let t_max = bounds * 2.0 / (dx.abs().max(dy.abs()).max(0.001));
    Line::new(PlotPoints::new(vec![
        [px - t_max * dx, py - t_max * dy],
        [px + t_max * dx, py + t_max * dy],
    ])).color(color).width(2.0)
}

/// Draw an arc (portion of circle) for visualizing rotations
pub fn arc_2d(cx: f64, cy: f64, radius: f64, start_angle: f64, end_angle: f64, color: Color32) -> Line {
    let segments = ((end_angle - start_angle).abs() * 20.0) as usize + 1;
    let points: Vec<[f64; 2]> = (0..=segments)
        .map(|i| {
            let t = i as f64 / segments as f64;
            let angle = start_angle + t * (end_angle - start_angle);
            [cx + radius * angle.cos(), cy + radius * angle.sin()]
        })
        .collect();
    Line::new(PlotPoints::new(points)).color(color).width(2.0)
}
```

### 5. Reusable Widgets

```rust
// examples/visualization/common/widgets.rs

/// Angle slider with degree display
pub fn angle_slider(ui: &mut egui::Ui, label: &str, radians: &mut f32) {
    let mut degrees = radians.to_degrees();
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::Slider::new(&mut degrees, -180.0..=180.0).suffix("°"));
    });
    *radians = degrees.to_radians();
}

/// Vector2 input
pub fn vector2_input(ui: &mut egui::Ui, label: &str, x: &mut f32, y: &mut f32) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(x).prefix("x: ").speed(0.1));
        ui.add(egui::DragValue::new(y).prefix("y: ").speed(0.1));
    });
}

/// Vector3 input
pub fn vector3_input(ui: &mut egui::Ui, label: &str, x: &mut f32, y: &mut f32, z: &mut f32) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(x).prefix("x: ").speed(0.1));
        ui.add(egui::DragValue::new(y).prefix("y: ").speed(0.1));
        ui.add(egui::DragValue::new(z).prefix("z: ").speed(0.1));
    });
}

/// Display a GA value with proper notation
pub fn ga_value_display(ui: &mut egui::Ui, label: &str, components: &[(&str, f32)]) {
    ui.horizontal(|ui| {
        ui.label(format!("{} =", label));
        ui.monospace(
            components
                .iter()
                .map(|(basis, val)| format!("{:.3}{}", val, basis))
                .collect::<Vec<_>>()
                .join(" + ")
        );
    });
}

/// Info tooltip with markdown-like formatting
pub fn info_box(ui: &mut egui::Ui, text: &str) {
    egui::Frame::none()
        .fill(Color32::from_rgba_unmultiplied(0, 0, 0, 200))
        .inner_margin(8.0)
        .rounding(4.0)
        .show(ui, |ui| {
            ui.label(egui::RichText::new(text).color(Color32::WHITE).size(12.0));
        });
}
```

### 6. Animation Utilities

```rust
// examples/visualization/common/animation.rs

/// Animation state for continuous playback
pub struct Animation {
    pub playing: bool,
    pub time: f32,
    pub speed: f32,
    pub loop_duration: f32,
}

impl Default for Animation {
    fn default() -> Self {
        Self {
            playing: false,
            time: 0.0,
            speed: 1.0,
            loop_duration: 4.0, // 4 second loop
        }
    }
}

impl Animation {
    pub fn update(&mut self, dt: f32) {
        if self.playing {
            self.time += dt * self.speed;
            if self.time > self.loop_duration {
                self.time -= self.loop_duration;
            }
        }
    }

    /// Normalized progress (0.0 to 1.0)
    pub fn progress(&self) -> f32 {
        self.time / self.loop_duration
    }

    /// Progress as angle (0 to 2π)
    pub fn angle(&self) -> f32 {
        self.progress() * std::f32::consts::TAU
    }

    pub fn toggle(&mut self) {
        self.playing = !self.playing;
    }

    pub fn reset(&mut self) {
        self.time = 0.0;
    }
}

/// Animation controls widget
pub fn animation_controls(ui: &mut egui::Ui, anim: &mut Animation) {
    ui.horizontal(|ui| {
        if ui.button(if anim.playing { "⏸ Pause" } else { "▶ Play" }).clicked() {
            anim.toggle();
        }
        if ui.button("⏮ Reset").clicked() {
            anim.reset();
        }
        ui.add(egui::Slider::new(&mut anim.speed, 0.1..=3.0).text("Speed"));
    });
}
```

## 3D Rendering Infrastructure

Many demos (euclidean3, projective3, conformal3, elliptic2, minkowski3) require 3D visualization. We use `three-d` crate integrated with egui.

### Dependencies

```toml
[dependencies]
three-d = "0.18"
three-d-asset = "0.18"
```

### 3D Viewport Component

```rust
// crates/clifford-viz/src/common/viewport3d.rs

use three_d::*;

/// 3D viewport that integrates with egui
pub struct Viewport3D {
    context: Context,
    camera: Camera,
    orbit_control: OrbitControl,
}

impl Viewport3D {
    pub fn new(width: u32, height: u32) -> Self {
        let context = Context::new().unwrap();
        let camera = Camera::new_perspective(
            Viewport::new_at_origo(width, height),
            vec3(3.0, 3.0, 3.0),  // position
            vec3(0.0, 0.0, 0.0),  // target
            vec3(0.0, 1.0, 0.0),  // up
            degrees(45.0),        // fov
            0.1,                  // near
            100.0,                // far
        );
        let orbit_control = OrbitControl::new(camera.target(), 1.0, 10.0);

        Self { context, camera, orbit_control }
    }

    pub fn handle_events(&mut self, events: &[egui::Event]) {
        // Convert egui events to three-d events
        // Update orbit control
    }

    pub fn render(&mut self, objects: &[&dyn Object]) -> egui::TextureId {
        // Render to texture, return for egui display
        todo!()
    }
}
```

### 3D Shape Primitives

```rust
// crates/clifford-viz/src/common/shapes3d.rs

use three_d::*;

/// Create a coordinate frame (RGB axes)
pub fn coordinate_frame(context: &Context, size: f32) -> Gm<Mesh, ColorMaterial> {
    // X axis (red), Y axis (green), Z axis (blue)
    todo!()
}

/// Create a wireframe box
pub fn wireframe_box(context: &Context, size: Vec3) -> Gm<Mesh, ColorMaterial> {
    todo!()
}

/// Create a line segment in 3D
pub fn line_3d(context: &Context, start: Vec3, end: Vec3, color: Color) -> Gm<Mesh, ColorMaterial> {
    todo!()
}

/// Create a point marker (small sphere)
pub fn point_3d(context: &Context, position: Vec3, radius: f32, color: Color) -> Gm<Mesh, ColorMaterial> {
    todo!()
}

/// Create a plane (bounded quad)
pub fn plane_3d(context: &Context, center: Vec3, normal: Vec3, size: f32, color: Color) -> Gm<Mesh, ColorMaterial> {
    todo!()
}

/// Create a sphere (wireframe or solid)
pub fn sphere_3d(context: &Context, center: Vec3, radius: f32, color: Color, wireframe: bool) -> Gm<Mesh, ColorMaterial> {
    todo!()
}

/// Create a circle in 3D (ring)
pub fn circle_3d(context: &Context, center: Vec3, normal: Vec3, radius: f32, color: Color) -> Gm<Mesh, ColorMaterial> {
    todo!()
}

/// Create an arrow in 3D
pub fn arrow_3d(context: &Context, origin: Vec3, direction: Vec3, color: Color) -> Gm<Mesh, ColorMaterial> {
    todo!()
}
```

### Camera Controls

```rust
// crates/clifford-viz/src/common/camera.rs

/// Camera state for 3D demos
pub struct Camera3D {
    pub position: [f32; 3],
    pub target: [f32; 3],
    pub up: [f32; 3],
    pub fov: f32,

    // Orbit control state
    orbit_radius: f32,
    orbit_theta: f32,  // Horizontal angle
    orbit_phi: f32,    // Vertical angle
}

impl Camera3D {
    pub fn orbit(&mut self, delta_theta: f32, delta_phi: f32) {
        self.orbit_theta += delta_theta;
        self.orbit_phi = (self.orbit_phi + delta_phi).clamp(-89.0_f32.to_radians(), 89.0_f32.to_radians());
        self.update_position();
    }

    pub fn zoom(&mut self, delta: f32) {
        self.orbit_radius = (self.orbit_radius - delta).clamp(1.0, 20.0);
        self.update_position();
    }

    pub fn pan(&mut self, delta_x: f32, delta_y: f32) {
        // Move target in camera-relative directions
        todo!()
    }

    fn update_position(&mut self) {
        self.position = [
            self.target[0] + self.orbit_radius * self.orbit_phi.cos() * self.orbit_theta.sin(),
            self.target[1] + self.orbit_radius * self.orbit_phi.sin(),
            self.target[2] + self.orbit_radius * self.orbit_phi.cos() * self.orbit_theta.cos(),
        ];
    }
}

/// Camera controls widget
pub fn camera_controls(ui: &mut egui::Ui, camera: &mut Camera3D) {
    ui.collapsing("Camera", |ui| {
        ui.horizontal(|ui| {
            ui.label("Distance:");
            ui.add(egui::DragValue::new(&mut camera.orbit_radius).speed(0.1).range(1.0..=20.0));
        });
        if ui.button("Reset View").clicked() {
            *camera = Camera3D::default();
        }
    });
}
```

### Integration with egui

```rust
// Render 3D content to a texture, display in egui panel

pub fn show_3d_viewport(ui: &mut egui::Ui, viewport: &mut Viewport3D, objects: &[&dyn Object]) {
    let available_size = ui.available_size();
    let (rect, response) = ui.allocate_exact_size(available_size, egui::Sense::drag());

    // Handle mouse input for camera control
    if response.dragged() {
        let delta = response.drag_delta();
        viewport.orbit_control.handle_events(/* convert delta */);
    }

    // Render 3D scene to texture
    let texture_id = viewport.render(objects);

    // Display texture in egui
    ui.painter().image(
        texture_id,
        rect,
        egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
        egui::Color32::WHITE,
    );
}
```

### Directory Structure Update

```
crates/clifford-viz/
  src/
    lib.rs
    common/
      mod.rs
      app.rs          # Base app trait
      colors.rs       # Color palette
      grid.rs         # 2D grid rendering
      shapes.rs       # 2D shape primitives
      shapes3d.rs     # 3D shape primitives (NEW)
      viewport3d.rs   # 3D viewport component (NEW)
      camera.rs       # 3D camera controls (NEW)
      widgets.rs      # UI components
      animation.rs    # Animation utilities
```

## Implementation Tasks

1. [ ] Create `crates/clifford-viz/src/common/mod.rs` with re-exports
2. [ ] Implement `app.rs` with base trait and runner
3. [ ] Implement `colors.rs` color palette
4. [ ] Implement `grid.rs` for 2D grids and axes
5. [ ] Implement `shapes.rs` 2D primitives (point, arrow, circle, line, arc)
6. [ ] Implement `shapes3d.rs` 3D primitives (box, sphere, plane, line, arrow)
7. [ ] Implement `viewport3d.rs` 3D viewport with three-d
8. [ ] Implement `camera.rs` 3D camera with orbit controls
9. [ ] Implement `widgets.rs` reusable UI components
10. [ ] Implement `animation.rs` timing utilities
11. [ ] Create minimal 2D test app
12. [ ] Create minimal 3D test app

## Verification

```bash
# Create a minimal test example
cargo run --example visualization_test --release
```

Test example should display:
- Window with title
- 2D plot area with grid
- Control panel with sliders
- Animation controls that work

## Dependencies on Other PRDs

None - this is the foundation.

## Dependent PRDs

All other PRD-48.x sub-PRDs depend on this.
