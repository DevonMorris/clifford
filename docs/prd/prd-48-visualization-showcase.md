# PRD-48: Visualization Showcase

**Status**: Draft
**Depends on**: [PRD-49](prd-49-workspace-restructure.md) (workspace restructure)
**Goal**: Create interactive visualizations demonstrating each algebra's unique geometric capabilities

## Sub-PRDs

| PRD | Title | Algebras | Phase |
|-----|-------|----------|-------|
| [PRD-48.1](prd-48.1-visualization-framework.md) | Visualization Framework | (infrastructure) | 1 |
| [PRD-48.10](prd-48.10-visual-testing.md) | Visual Regression Testing | (infrastructure) | 1 |
| [PRD-48.2](prd-48.2-euclidean-visualization.md) | Euclidean Visualization | euclidean2, euclidean3 | 2 |
| [PRD-48.3](prd-48.3-projective2-visualization.md) | Projective 2D Visualization | projective2 | 3 |
| [PRD-48.4](prd-48.4-projective3-visualization.md) | Projective 3D Visualization | projective3 | 3 |
| [PRD-48.5](prd-48.5-conformal3-visualization.md) | Conformal 3D Visualization | conformal3 | 4 |
| [PRD-48.6](prd-48.6-quaternion-dualquat-visualization.md) | Quaternion & Dual Quaternion | quaternion, dualquat | 5 |
| [PRD-48.7](prd-48.7-complex-dual-visualization.md) | Complex & Dual Numbers | complex, dual | 5 |
| [PRD-48.8](prd-48.8-non-euclidean-visualization.md) | Non-Euclidean Geometry | hyperbolic2, elliptic2 | 6 |
| [PRD-48.9](prd-48.9-minkowski-visualization.md) | Minkowski Spacetime | minkowski2, minkowski3 | 6 |

**Implementation Order**: PRD-48.1 and PRD-48.10 are foundational and should be implemented first. Visual testing infrastructure (48.10) enables regression testing for all subsequent demos.

## Motivation

The clifford library now supports 14 distinct algebras with powerful code generation. However, the mathematical beauty and practical utility of these algebras is hard to appreciate without visual demonstrations. This PRD proposes a comprehensive visualization suite that:

1. **Showcases library capabilities** - Demonstrates what users can build
2. **Educational value** - Helps users understand GA concepts visually
3. **Validates correctness** - Visual output catches errors that unit tests miss
4. **Marketing/adoption** - Eye-catching demos attract users to the library

## Design Principles

### Framework Choice: egui + eframe

Use [egui](https://github.com/emilk/egui) for cross-platform, immediate-mode GUI:

```toml
[dev-dependencies]
eframe = "0.30"
egui = "0.30"
egui_plot = "0.30"  # 2D plotting
three-d = "0.18"    # 3D rendering (optional, for complex scenes)
```

**Rationale**:
- Pure Rust, minimal dependencies
- Works on desktop and WASM (web demos!)
- Immediate mode pairs well with real-time animation
- Strong ecosystem (egui_plot for 2D, three-d for 3D)

### Crate Structure: `clifford-viz`

A separate crate for visualization, keeping the main `clifford` crate dependency-free:

```
crates/
  clifford-viz/
    Cargo.toml          # Depends on clifford, egui, eframe
    src/
      lib.rs            # Public API for reusable widgets
      common/           # Shared rendering utilities
        mod.rs
        app.rs          # Base app trait and runner
        camera.rs       # 3D camera controls
        colors.rs       # Consistent color scheme
        grid.rs         # Background grids
        shapes.rs       # Primitive shape drawing
        widgets.rs      # Reusable UI components
        animation.rs    # Animation utilities
      euclidean.rs      # Euclidean viz widgets (RotorWidget, etc.)
      projective.rs     # PGA viz widgets (MotorWidget, etc.)
      conformal.rs      # CGA viz widgets
      quaternion.rs     # Quaternion viz widgets
      spacetime.rs      # Minkowski viz widgets
    examples/           # Interactive demos
      euclidean2.rs
      euclidean3.rs
      projective2.rs
      projective3.rs
      conformal3.rs
      quaternion.rs
      dualquat.rs
      complex.rs
      dual.rs
      hyperbolic2.rs
      elliptic2.rs
      minkowski2.rs
      minkowski3.rs
    tests/
      visual/           # Visual regression tests
        mod.rs
        golden/         # Baseline images
```

Run with: `cargo run -p clifford-viz --example projective3 --release`

### Cargo.toml

```toml
[package]
name = "clifford-viz"
version = "0.1.0"
edition = "2024"

[dependencies]
clifford = { path = "../clifford" }
eframe = "0.30"
egui = "0.30"
egui_plot = "0.30"

[dev-dependencies]
image = "0.25"
```

---

## Visualization Catalog

### 1. Euclidean 2D (`euclidean2`)

**Demo: Vector Field Rotor Animation**

Visualize how rotors smoothly rotate a field of vectors:

```
┌─────────────────────────────────────┐
│  ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗    │
│  ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗    │
│  ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗    │  [Angle: 45°]
│  ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗    │  [▶ Play/Pause]
│  ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗ ↗    │
└─────────────────────────────────────┘
```

**Features**:
- Slider for rotation angle (0° to 360°)
- Animation mode: auto-rotate at constant angular velocity
- Show bivector as oriented arc at origin
- Display rotor components: `R = cos(θ/2) + sin(θ/2)e₁₂`

---

### 2. Euclidean 3D (`euclidean3`)

**Demo: Gimbal-Free 3D Rotation**

Interactive 3D rotation demonstrating rotor superiority over Euler angles:

```
┌─────────────────────────────────────┐
│                                     │
│         ┌───────────┐               │  Rotation Mode:
│        ╱           ╱│               │  ○ Euler Angles
│       ╱───────────╱ │               │  ● Rotor (GA)
│      │           │  │               │
│      │    3D     │  │               │  Bivector Plane:
│      │   Cube    │ ╱                │  [XY] [XZ] [YZ]
│      │           │╱                 │
│      └───────────┘                  │  Angle: [====●===]
│                                     │
└─────────────────────────────────────┘
```

**Features**:
- Drag to rotate object
- Compare Euler angles (shows gimbal lock) vs rotor (smooth)
- Visualize rotation plane (bivector) as a disk
- Compose multiple rotations: `R_total = R_2 R_1`
- Show rotor interpolation (SLERP-like)

---

### 3. Projective 2D (`projective2`)

**Demo: Interactive Line-Point Geometry**

Demonstrate join (wedge) and meet (antiwedge) operations:

```
┌─────────────────────────────────────┐
│                                     │
│     A ●─────────────────────● B     │  Points: A, B, C
│              ╲    ╱                 │  Lines: AB, CD
│               ╲  ╱                  │
│                ╳ ← intersection     │  A ∧ B = Line AB
│               ╱  ╲                  │  AB ∨ CD = Point P
│              ╱    ╲                 │
│     C ●─────────────────────● D     │  [Drag points]
│                                     │
└─────────────────────────────────────┘
```

**Features**:
- Drag points to see lines update in real-time
- Click to create new points
- Join any two points → line appears
- Meet any two lines → intersection point appears
- Motor animation: rigid motion of entire figure

**Advanced: Animated 2D Robot**

Simple 2D robot arm with motor-based kinematics:

```
┌─────────────────────────────────────┐
│                                     │
│     ○────────●────────○             │  Joint 1: [===●===]
│     │         ╲                     │  Joint 2: [===●===]
│   base         ╲                    │
│                 ● end effector      │  Forward Kinematics
│                                     │  using Motor chain
└─────────────────────────────────────┘
```

---

### 4. Projective 3D (`projective3`)

**Demo: Plücker Line Coordinates**

Visualize the 6 components of a 3D line:

```
┌─────────────────────────────────────┐
│                                     │
│            ↗ direction (d)          │  Direction: [dx, dy, dz]
│           ╱                         │  Moment:    [mx, my, mz]
│     ─────●─────────                 │
│         ╱  ╲                        │  d × p₀ = m
│        ╱    ╲ moment (m)            │  (moment encodes position)
│       ○ origin                      │
│                                     │  [Drag line endpoints]
└─────────────────────────────────────┘
```

**Features**:
- Interactive line placement (two-point definition)
- Show direction vector and moment vector
- Meet operations: line ∨ line = point, line ∨ plane = point
- Join operations: point ∧ point = line, point ∧ line = plane

**Advanced: 3D Rigid Body Animation**

```
┌─────────────────────────────────────┐
│                                     │
│    ┌─────┐                          │  Motor: M = e^(d·L/2)
│    │     │  ──→ screw motion        │
│    │     │      along axis          │  L = line (screw axis)
│    └─────┘                          │  d = displacement
│                                     │
│    ════════════════ axis            │  [Animate] [Step]
│                                     │
└─────────────────────────────────────┘
```

---

### 5. Conformal 3D (`conformal3`)

**Demo: Sphere Inversion**

The crown jewel - demonstrate conformal transformations:

```
┌─────────────────────────────────────┐
│                                     │
│         ╭─────╮                     │  Inversion Sphere:
│        ╱       ╲     → →            │  Center: [x, y, z]
│       │    ●    │   transforms      │  Radius: [====●===]
│        ╲       ╱     ← ←            │
│         ╰─────╯                     │  Objects:
│              ↓                      │  ○ → ○  (sphere → sphere)
│                                     │  ● → ●  (point → point)
│    ────────────────────────         │  ─ → ○  (line → circle)
│    line through origin              │  ○ → ─  (circle → line)
│    becomes circle!                  │
└─────────────────────────────────────┘
```

**Features**:
- Interactive inversion sphere (drag center, resize)
- Place spheres, circles, lines, planes
- Watch objects transform under inversion
- Möbius transformations on circles
- Ray-circle/sphere intersection demos

**Advanced: Circle Packing**

```
┌─────────────────────────────────────┐
│                                     │
│     ○○○○○○○○○○○○○○○○○○○○○○○○○○○    │  Apollonian Gasket
│    ○   ○   ○   ○   ○   ○   ○   ○   │
│   ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○    │  Each circle tangent
│    ○   ○   ○   ○   ○   ○   ○   ○   │  to three others
│     ○○○○○○○○○○○○○○○○○○○○○○○○○○○    │
│                                     │  [Iterations: 5]
└─────────────────────────────────────┘
```

Use CGA to compute tangent circles efficiently via meet operations.

---

### 6. Quaternion (`quaternion`)

**Demo: Rotation Interpolation**

Compare quaternion SLERP with naive linear interpolation:

```
┌─────────────────────────────────────┐
│                                     │
│   Start ●─────────────────● End     │  Interpolation t: [●═══]
│          ╲               ╱          │
│           ╲─────────────╱           │  ○ Linear (broken)
│            ╲           ╱            │  ● SLERP (smooth)
│             ╲─────────╱             │
│              ╲       ╱              │  Current rotation:
│               ╲─────╱ arc path      │  w=0.9 i=0.2 j=0.1 k=0.3
│                                     │
└─────────────────────────────────────┘
```

**Features**:
- Set start and end rotations
- Slider for interpolation parameter t
- Compare SLERP (spherical) vs LERP (straight line - wrong!)
- Show path on unit 4-sphere (projected to 3D)
- Visualize quaternion double-cover (360° vs 720° rotation)

---

### 7. Dual Quaternion (`dualquat`)

**Demo: Screw Motion Animation**

Demonstrate the unified rotation+translation of dual quaternions:

```
┌─────────────────────────────────────┐
│                                     │
│    ═══════════════════════════      │  Screw Axis
│         ╱                           │
│   ┌────┼────┐                       │  Pitch: [====●===]
│   │    ●    │  ──spiral──→          │  (translation per rotation)
│   └────┼────┘                       │
│        │                            │  Total rotation: 180°
│        ↓                            │  Total translation: 2.0
│    ═══════════════════════════      │
│                                     │  [▶ Animate]
└─────────────────────────────────────┘
```

**Features**:
- Animate object along screw trajectory
- Adjust pitch (translation per rotation)
- Compare with separate rotation + translation
- Show dual quaternion interpolation (ScLERP)
- Robot arm inverse kinematics demo

---

### 8. Complex Numbers (`complex`)

**Demo: Domain Coloring**

Visualize complex functions using domain coloring:

```
┌─────────────────────────────────────┐
│                                     │
│  ████████████████████████████████   │  f(z) = z²
│  ████████████████████████████████   │
│  ████████████████████████████████   │  Functions:
│  ████████████████████████████████   │  ○ z²   ○ z³
│  ████████████████████████████████   │  ○ 1/z  ○ e^z
│  ████████████████████████████████   │  ○ sin(z)
│  ████████████████████████████████   │
│                                     │  Color = arg(f(z))
│                                     │  Brightness = |f(z)|
└─────────────────────────────────────┘
```

**Features**:
- Complex function selection dropdown
- Hue encodes argument (angle)
- Brightness encodes magnitude
- Zoom and pan
- Optional: Mandelbrot/Julia set explorer

---

### 9. Dual Numbers (`dual`)

**Demo: Automatic Differentiation Visualizer**

Show how dual numbers compute derivatives:

```
┌─────────────────────────────────────┐
│                                     │
│         ╭                           │  f(x) = sin(x²)
│        ╱                            │
│    ───●───────────────────          │  At x = 1.5:
│      ╱│╲                            │    f(x)  = 0.778
│     ╱ │ ╲ tangent line              │    f'(x) = 1.942
│    ╱  │                             │
│       │                             │  Dual: sin((1.5 + ε)²)
│       x                             │       = 0.778 + 1.942ε
│                                     │
│  x = [=========●========]           │  [Edit function]
└─────────────────────────────────────┘
```

**Features**:
- Plot function f(x)
- Slider to move evaluation point
- Show tangent line (computed via dual numbers)
- Display dual number calculation breakdown
- Compare with numerical differentiation

---

### 10. Hyperbolic Plane (`hyperbolic2`)

**Demo: Poincaré Disk Model**

Visualize hyperbolic geometry in the Poincaré disk:

```
┌─────────────────────────────────────┐
│                                     │
│        ╭───────────────╮            │  Hyperbolic:
│       ╱   ╱───────╲     ╲           │  • Lines are circular arcs
│      │   │    ●    │     │          │  • Parallel lines exist!
│      │   │ center  │     │          │  • Sum of triangle angles
│       ╲   ╲───────╱     ╱           │    < 180°
│        ╰───────────────╯            │
│              boundary               │  [Add Point] [Add Line]
│         (points at infinity)        │  [Translate] [Rotate]
└─────────────────────────────────────┘
```

**Features**:
- Draw hyperbolic geodesics (circular arcs)
- Demonstrate parallel lines (Lobachevsky geometry)
- Hyperbolic rotations (rotor animation)
- Hyperbolic translations
- Tessellation patterns (optional)

---

### 11. Elliptic/Spherical (`elliptic2`)

**Demo: Spherical Geometry**

Geometry on the sphere surface:

```
┌─────────────────────────────────────┐
│                                     │
│           ╱ ╲                       │  Spherical Geometry:
│          ╱   ╲ great circle         │  • Lines are great circles
│         ╱  ●  ╲                     │  • All lines intersect!
│        ╱───────╲                    │  • Sum of triangle angles
│       │         │                   │    > 180°
│       │    ●────│────●              │
│        ╲       ╱                    │  [Add Point] [Add Line]
│         ╲─────╱                     │  [Rotate Globe]
│                                     │
└─────────────────────────────────────┘
```

**Features**:
- 3D sphere visualization
- Great circles as "lines"
- Spherical triangles with angle sum > 180°
- Parallel transport demonstration
- Rotor-based rotations on sphere

---

### 12. Minkowski Plane (`minkowski2`)

**Demo: 1+1 Spacetime Diagram**

Interactive spacetime visualization:

```
┌─────────────────────────────────────┐
│  t                                  │
│  ↑    ╲ light    light ╱            │  Boost velocity:
│  │     ╲  cone    cone╱             │  v = [====●===] 0.5c
│  │      ╲     │     ╱               │
│  │       ╲    │    ╱                │  Events:
│  │        ╲   │   ╱                 │  ● A: (0, 0)
│  │         ╲  ●  ╱                  │  ● B: (2, 1)
│  │──────────╲│╱──────────→ x        │
│             ╱│╲                     │  Proper time A→B:
│            ╱ │ ╲                    │  τ = 1.73
│                                     │
└─────────────────────────────────────┘
```

**Features**:
- Draw worldlines
- Apply Lorentz boosts (hyperbolic rotations)
- Visualize light cones
- Calculate proper time/length
- Time dilation demonstration

---

### 13. Minkowski 3+1 (`minkowski3`)

**Demo: 3D Light Cones**

Full spacetime visualization:

```
┌─────────────────────────────────────┐
│                                     │
│              ╲     ╱                │  3D space + 1D time
│               ╲   ╱ future          │  (time shown as vertical)
│                ╲ ╱  light           │
│        ─────────●─────────          │  Lorentz Transform:
│                ╱ ╲  cone            │  Boost: [x] [y] [z]
│               ╱   ╲                 │  Rotation: [xy][xz][yz]
│              ╱     ╲ past           │
│                                     │  [▶ Animate Boost]
└─────────────────────────────────────┘
```

**Features**:
- 3D visualization with time as color/height
- Lorentz boosts in any direction
- Spatial rotations
- Length contraction visualization
- Simultaneity demonstration

---

## Implementation Phases

### Phase 1: Framework Setup (Foundation)
1. Create `examples/visualization/common/` utilities
2. Set up egui/eframe boilerplate
3. Implement 2D grid, axis drawing
4. Implement basic 3D camera
5. Define consistent color palette

### Phase 2: Euclidean Demos (Core GA)
1. `euclidean2.rs` - 2D rotor animation
2. `euclidean3.rs` - 3D rotor, gimbal lock demo
3. Both with vector field visualization

### Phase 3: Projective Demos (PGA Showcase)
1. `projective2.rs` - Point-line geometry, motor animation
2. `projective3.rs` - Plücker coordinates, rigid body motion
3. Interactive drag-and-drop for points/lines

### Phase 4: Conformal Demo (Advanced)
1. `conformal3.rs` - Sphere inversion, Möbius transforms
2. Circle/sphere visualization
3. Optional: Apollonian gasket

### Phase 5: Specialty Algebras
1. `quaternion.rs` - SLERP comparison
2. `dualquat.rs` - Screw motion
3. `complex.rs` - Domain coloring
4. `dual.rs` - Autodiff visualization

### Phase 6: Non-Euclidean Geometries
1. `hyperbolic2.rs` - Poincaré disk
2. `elliptic2.rs` - Spherical geometry
3. `minkowski2.rs` - 1+1 spacetime
4. `minkowski3.rs` - 3+1 spacetime (advanced)

---

## Cross-Cutting Features

### Animation Framework

```rust
/// Trait for animatable parameters
pub trait Animatable {
    fn update(&mut self, dt: f32);
    fn is_playing(&self) -> bool;
    fn toggle_play(&mut self);
    fn reset(&mut self);
}

/// Common animation controls widget
pub fn animation_controls(ui: &mut egui::Ui, anim: &mut impl Animatable) {
    ui.horizontal(|ui| {
        if ui.button(if anim.is_playing() { "⏸" } else { "▶" }).clicked() {
            anim.toggle_play();
        }
        if ui.button("⏮").clicked() {
            anim.reset();
        }
    });
}
```

### Value Display

Show GA values with proper notation:

```rust
/// Format a rotor for display
pub fn format_rotor_2d<T: Float + Display>(r: &Rotor<T>) -> String {
    format!("{:.3} + {:.3}e₁₂", r.s(), r.b())
}

/// Format with subscript unicode
pub fn format_bivector_3d<T: Float + Display>(b: &Bivector<T>) -> String {
    format!("{:.3}e₁₂ + {:.3}e₁₃ + {:.3}e₂₃", b.rz(), b.ry(), b.rx())
}
```

### Interactive Object Selection

```rust
pub struct Selectable<T> {
    pub value: T,
    pub selected: bool,
    pub hovered: bool,
    pub dragging: bool,
}
```

---

## WASM Support

Enable web deployment for maximum accessibility:

```toml
[target.'cfg(target_arch = "wasm32")'.dependencies]
eframe = { version = "0.30", default-features = false, features = ["wgpu"] }
```

Build with: `trunk build examples/visualization/projective3.rs`

Host on GitHub Pages for easy sharing.

---

## Files to Create

### New Files
- `examples/visualization/common/mod.rs`
- `examples/visualization/common/camera.rs`
- `examples/visualization/common/colors.rs`
- `examples/visualization/common/grid.rs`
- `examples/visualization/common/widgets.rs`
- `examples/visualization/euclidean2.rs`
- `examples/visualization/euclidean3.rs`
- `examples/visualization/projective2.rs`
- `examples/visualization/projective3.rs`
- `examples/visualization/conformal3.rs`
- `examples/visualization/quaternion.rs`
- `examples/visualization/dualquat.rs`
- `examples/visualization/complex.rs`
- `examples/visualization/dual.rs`
- `examples/visualization/hyperbolic2.rs`
- `examples/visualization/elliptic2.rs`
- `examples/visualization/minkowski2.rs`
- `examples/visualization/minkowski3.rs`

### Modified Files
- `Cargo.toml` - Add dev-dependencies for egui, eframe, egui_plot, three-d

---

## Success Criteria

- [ ] All 14 algebras have at least one visualization
- [ ] Demos run at 60fps on mid-range hardware
- [ ] WASM builds work for web deployment
- [ ] Each demo has explanatory text for educational value
- [ ] Interactive controls are intuitive
- [ ] Visual output matches mathematical expectations

---

## Future Enhancements

1. **Guided tutorials** - Step-by-step walkthroughs of GA concepts
2. **Export to Rerun** - Log visualization data for replay
3. **VR support** - Immersive 3D GA exploration
4. **Benchmark visualizations** - Performance graphs
5. **Comparison mode** - Side-by-side matrix vs GA implementations
