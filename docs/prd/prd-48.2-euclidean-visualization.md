# PRD-48.2: Euclidean Geometry Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for euclidean2 and euclidean3 algebras

## Overview

Demonstrate rotors, vectors, and bivectors in 2D and 3D Euclidean space. Focus on showing the elegance of GA rotations vs traditional approaches.

## Demo 1: euclidean2 - 2D Rotor Animation

**File**: `examples/visualization/euclidean2.rs`

### Features

1. **Vector Field Display**
   - Grid of arrows showing vector directions
   - All vectors rotate together under rotor action

2. **Rotor Controls**
   - Angle slider (0° to 360°)
   - Animation toggle for continuous rotation
   - Speed control

3. **Bivector Visualization**
   - Show bivector as oriented arc at origin
   - Arc angle corresponds to rotor angle/2

4. **Component Display**
   - Show rotor as `R = cos(θ/2) + sin(θ/2)e₁₂`
   - Real-time update as angle changes

### Implementation

```rust
use clifford::specialized::euclidean::dim2::{Vector, Bivector, Rotor};
use crate::visualization::common::*;

pub struct Euclidean2Demo {
    angle: f32,
    animation: Animation,
    show_bivector: bool,
    vector_grid_size: usize,
}

impl VisualizationApp for Euclidean2Demo {
    fn name(&self) -> &'static str { "Euclidean 2D - Rotor Animation" }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            self.angle = self.animation.angle();
        }
    }

    fn render(&self, ui: &mut egui::Ui) {
        let rotor = Rotor::from_angle(self.angle);

        egui_plot::Plot::new("euclidean2")
            .data_aspect(1.0)
            .show(ui, |plot_ui| {
                // Draw grid
                for line in grid_2d(5.0, 1.0) {
                    plot_ui.line(line);
                }

                // Draw rotated vector field
                for i in -2..=2 {
                    for j in -2..=2 {
                        if i == 0 && j == 0 { continue; }
                        let base = Vector::new(i as f32, j as f32);
                        let rotated = rotor.sandwich(&base);
                        for arrow in arrow_2d(0.0, 0.0, rotated.x() as f64, rotated.y() as f64, palette::POINT) {
                            plot_ui.line(arrow);
                        }
                    }
                }

                // Draw bivector arc
                if self.show_bivector {
                    plot_ui.line(arc_2d(0.0, 0.0, 0.5, 0.0, self.angle as f64 / 2.0, palette::ROTOR));
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        angle_slider(ui, "Rotation angle", &mut self.angle);
        ui.checkbox(&mut self.show_bivector, "Show bivector (half-angle)");
        ui.separator();
        animation_controls(ui, &mut self.animation);
    }

    fn info(&self, ui: &mut egui::Ui) {
        let rotor = Rotor::from_angle(self.angle);
        info_box(ui, &format!(
            "Rotor: R = {:.3} + {:.3}e₁₂\n\
             Angle: {:.1}°\n\n\
             The rotor encodes rotation as R = cos(θ/2) + sin(θ/2)e₁₂.\n\
             Vectors transform via v' = RvR†",
            rotor.s(), rotor.b(),
            self.angle.to_degrees()
        ));
    }
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Euclidean 2D - Rotor Animation                        [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Controls              │
│     ↗  ↗  ↗  ↗  ↗                  │                       │
│                                     │ Rotation: [====●====] │
│     ↗  ↗  ↗  ↗  ↗                  │           135°        │
│           ╭─╮                       │                       │
│     ↗  ↗ ╱   ↗  ↗                  │ [✓] Show bivector     │
│          ╰──                        │                       │
│     ↗  ↗  ↗  ↗  ↗                  │ ─────────────────     │
│                                     │ [▶ Play] [⏮] Speed:1x │
│     ↗  ↗  ↗  ↗  ↗                  │                       │
│                                     ├───────────────────────┤
│                                     │ Info                  │
│                                     │                       │
│                                     │ Rotor: 0.707 + 0.707e₁₂│
│                                     │ Angle: 90.0°          │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 2: euclidean3 - 3D Rotor vs Euler Angles

**File**: `examples/visualization/euclidean3.rs`

### Features

1. **3D Object Display**
   - Wireframe cube or coordinate frame
   - Smooth rotation via rotor

2. **Dual Mode Comparison**
   - Toggle between Euler angles and Rotor
   - Euler mode shows gimbal lock
   - Rotor mode is always smooth

3. **Bivector Plane Selection**
   - Buttons for XY, XZ, YZ planes
   - Custom bivector input
   - Visualize rotation plane as disk

4. **Rotor Composition**
   - Apply multiple rotations
   - Show that order matters: R₂R₁ ≠ R₁R₂
   - Compare with matrix multiplication

### Implementation

```rust
use clifford::specialized::euclidean::dim3::{Vector, Bivector, Rotor};

pub struct Euclidean3Demo {
    mode: RotationMode,
    // Rotor mode
    rotor_angle: f32,
    rotor_plane: BivectorPlane,
    // Euler mode
    euler_x: f32,
    euler_y: f32,
    euler_z: f32,
    // Visualization
    show_rotation_plane: bool,
    show_gimbal: bool,
    animation: Animation,
    camera: Camera3D,
}

#[derive(PartialEq)]
enum RotationMode {
    Rotor,
    Euler,
}

#[derive(PartialEq, Clone, Copy)]
enum BivectorPlane {
    XY,
    XZ,
    YZ,
    Custom(f32, f32, f32),
}

impl BivectorPlane {
    fn to_bivector(&self) -> Bivector<f32> {
        match self {
            Self::XY => Bivector::unit_xy(),
            Self::XZ => Bivector::unit_xz(),
            Self::YZ => Bivector::unit_yz(),
            Self::Custom(xy, xz, yz) => Bivector::new(*xy, *xz, *yz).normalize(),
        }
    }
}
```

### Gimbal Lock Demonstration

```
Euler Angles Mode:
┌─────────────────────────────────────┐
│                                     │
│    Pitch (Y): 90°  ← GIMBAL LOCK!   │
│                                     │
│    At pitch=90°, roll and yaw       │
│    rotate around the SAME axis.    │
│    One degree of freedom is lost.   │
│                                     │
│         ┌───┐                       │
│         │ × │ ← stuck!              │
│         └───┘                       │
│                                     │
└─────────────────────────────────────┘

Rotor Mode:
┌─────────────────────────────────────┐
│                                     │
│    Same orientation, but rotor      │
│    works perfectly!                 │
│                                     │
│    R = cos(45°) + sin(45°)e₁₃      │
│                                     │
│         ┌───┐                       │
│         │   │ ← smooth rotation     │
│         └───┘                       │
│                                     │
└─────────────────────────────────────┘
```

### Rotor Composition Demo

Show that `R_total = R_2 * R_1` (right-to-left application):

```
┌─────────────────────────────────────┐
│ Rotation Composition                │
│                                     │
│ R₁: 90° in XY plane                │
│ R₂: 90° in XZ plane                │
│                                     │
│ R₂R₁ result:   [3D view]           │
│ R₁R₂ result:   [3D view]           │
│                                     │
│ Note: R₂R₁ ≠ R₁R₂ (non-commutative)│
└─────────────────────────────────────┘
```

---

## Implementation Tasks

### euclidean2.rs
1. [x] Basic app structure with `VisualizationApp` trait
2. [x] Vector field rendering (grid of arrows)
3. [x] Rotor application to all vectors
4. [x] Angle slider with degree display
5. [x] Bivector arc visualization (half-angle)
6. [x] Animation support (play/pause, speed, progress)
7. [x] Component display panel (rotor components with formulas)
8. [x] **Educational content** - "Learn About This" popup with:
   - Mathematical background (rotor formula, sandwich product)
   - How to use the visualization
   - Key concepts (half-angles, composition, unit magnitude)
   - Links to external resources

### euclidean3.rs
1. [ ] 3D rendering setup (camera, projection)
2. [ ] Wireframe cube rendering
3. [ ] Rotor mode with bivector plane selection
4. [ ] Euler angles mode with gimbal visualization
5. [ ] Side-by-side comparison view
6. [ ] Rotation composition demo
7. [ ] Animation support
8. [ ] **Educational content** - gimbal lock explanation, rotor advantages

## Verification

```bash
cargo run -p clifford-viz --example euclidean2 --release
cargo run -p clifford-viz --example euclidean3 --release
```

### euclidean2 Checklist
- [x] Vectors rotate smoothly with angle slider
- [x] Animation loops continuously
- [x] Bivector arc shows half-angle correctly
- [x] Component display updates in real-time
- [x] Educational popup explains the math

### euclidean3 Checklist
- [ ] 3D cube renders correctly
- [ ] Rotor rotation is smooth for all angles
- [ ] Euler mode shows gimbal lock at pitch=90°
- [ ] Rotor composition demo shows non-commutativity
- [ ] Camera controls work (drag to orbit)
- [ ] Educational popup explains gimbal lock vs rotors

## API Extraction

**IMPORTANT**: While developing visualizations, identify useful methods that should be promoted to the main algebra API.

After completing each visualization, review for methods that would benefit users:

### Candidates from euclidean2
- [ ] `Rotor::with_scale(angle, scale)` - create dilating rotor
- [ ] `Rotor::scale(&self)` / `Rotor::dilation_factor(&self)` - extract scale from rotor

### Candidates from euclidean3
- [ ] TBD after implementation

### Process
1. Identify helper functions written in visualization code
2. Evaluate if they have general utility beyond visualization
3. If useful, add to `extensions.rs` for the corresponding algebra
4. Follow CLAUDE.md guidelines (no shadowing traits, semantic naming)
5. Add tests and documentation
