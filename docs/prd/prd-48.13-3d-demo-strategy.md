# PRD-48.13: 3D Demo Strategy and Lessons Learned

**Status**: Implemented (Phase 1-2 complete)
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Strategy for implementing 3D demos based on lessons from 2D development

## Lessons Learned from 2D Demos

### What Worked Well

#### 1. Unified Framework Architecture
The `VisualizationApp` trait pattern proved excellent:
- Consistent structure across all demos (render, controls, info, educational_content)
- Easy to add new demos without boilerplate
- Educational content popup encourages mathematical explanations
- Responsive design (mobile/tablet/desktop) handled centrally

#### 2. egui_plot for 2D Visualization
- Simple API for points, lines, polygons
- Built-in pan/zoom/data-aspect ratio
- Good performance for real-time updates
- Mouse interaction via PlotResponse

#### 3. Demo Complexity Progression
Building simpler demos first revealed patterns:
- Start with basic rendering, then add interactivity
- Educational content written alongside implementation (not after)
- Test interactive features manually during development

#### 4. Reusable Shape Primitives
`common/shapes.rs` helpers were valuable:
- `arrow_2d`, `circle_2d`, `arc_2d`, `line_segment`
- Color palette in `colors.rs` ensures consistency
- Grid rendering in `grid.rs`

### Challenges Encountered

#### 1. Interactive Object Selection
Each demo reinvented selection/dragging:
- `projective2.rs` has DraggablePoint, Selection enum
- `conformal2_circles.rs` has its own point selection
- **Lesson**: Need a shared interaction abstraction

#### 2. Coordinate Mapping
Plot coordinates vs screen coordinates required careful handling:
- `world_to_screen` / `screen_to_world` transforms
- Plot bounds management for consistent views
- **Lesson**: Standardize coordinate helpers

#### 3. Animation Timing
Each demo managed animation differently:
- Some use `Animation` struct from common
- Some have custom timing
- **Lesson**: Animation utilities need enhancement for 3D

#### 4. Large Demo Files
Some demos grew to 20-35KB:
- `conformal2_inversion.rs`: 28KB
- `conformal2_mobius.rs`: 34KB
- `projective2.rs`: 31KB
- **Lesson**: Consider splitting demos or extracting shared code

### Missing Infrastructure for 3D - RESOLVED

We now have TWO approaches for 3D rendering:

#### Approach A: Custom Projection via egui_plot (WASM-compatible)
- `camera3d.rs` - perspective projection, orbit controls
- `shapes3d.rs` - wireframe primitives
- Works in WASM, simpler, sufficient for wireframes
- Used by: `euclidean3.rs` demo

#### Approach B: three-d Library (Real 3D Rendering)
- PBR materials with proper lighting
- three-d manages its own window + OpenGL context
- egui integration via `GUI` struct
- Used by: `euclidean3_three_d.rs` demo

**Decision**: Use both approaches:
- 2D demos + simple 3D wireframes: eframe + egui_plot
- Full 3D demos with lighting: three-d + egui

---

## 3D Demo Strategy

### Dual Architecture

**2D Demos**: eframe + egui + egui_plot (current)
- All conformal2, projective2, euclidean2 demos
- Works in WASM via trunk
- Simple and reliable

**3D Demos**: three-d + egui
- Real 3D rendering with PBR lighting
- three-d's `Window` + `GUI` for egui overlay
- WASM support via wasm-pack (needs validation)

### Phase 1: 3D Infrastructure - COMPLETE

#### Custom Projection (egui_plot approach)

Given that WASM compatibility is critical, we use **Option B: Custom egui-based 3D**.

**Why not three-d?**
- three-d uses wgpu which has complex WASM requirements
- Adds significant dependency weight
- Risk of browser compatibility issues

**Why custom projection works:**
- egui_plot already proven in WASM (all 12 2D demos work)
- Projecting 3D to 2D is straightforward math
- Wireframe rendering is sufficient for our educational demos
- No new dependencies required

#### Implementation Tasks

1. **Create camera3d.rs** (`common/camera3d.rs`)
   - `Camera3D` struct with position, target, up, fov
   - `project(point: [f32; 3]) -> [f64; 2]` - perspective projection
   - Orbit controls: `orbit(delta_theta, delta_phi)`
   - Zoom: `zoom(delta)`
   - Pan: `pan(dx, dy)`
   - `camera_controls_widget(ui, camera)` for settings panel

2. **Create shapes3d.rs** (`common/shapes3d.rs`)
   - `wireframe_box(camera, center, size, color) -> Vec<Line>`
   - `coordinate_axes(camera, length) -> Vec<Line>` (RGB)
   - `point_3d(camera, pos, radius, color) -> (circle approximation)`
   - `line_3d(camera, start, end, color) -> Line`
   - `arrow_3d(camera, origin, direction, color) -> Vec<Line>`
   - `circle_3d(camera, center, normal, radius, color) -> Line` (projected ellipse)
   - `sphere_wireframe(camera, center, radius, color) -> Vec<Line>`

3. **Create test_3d.rs example**
   - Render wireframe cube with coordinate axes
   - Mouse drag to orbit camera
   - Scroll to zoom
   - Validate works in WASM via `trunk serve`

4. **Integrate with existing app framework**
   - Modify `render()` to handle 3D plot setup
   - Data aspect ratio considerations for 3D projection

### Phase 1.5: WASM Validation (CRITICAL)

Before building any real demos, create a minimal 3D test that:
1. Renders a simple wireframe cube
2. Has basic mouse orbit controls
3. Compiles to WASM and runs in browser

**If three-d + WASM fails**: Fall back to custom 3D projection using egui_plot:
- Project 3D points to 2D using perspective/orthographic projection
- Draw edges as 2D line segments
- This is simpler but sufficient for our needs (wireframes, not photorealistic)

### Phase 2: euclidean3 Demo (PRD-48.2 completion)

The euclidean3 demo is the ideal first 3D demo because:
- Conceptually simple (rotors rotating a cube)
- Tests all basic 3D infrastructure
- Establishes patterns for later demos
- **WASM compatibility validated early**

#### Features
1. Wireframe cube rendering
2. Rotor-based rotation with angle slider
3. Bivector plane selection (XY, XZ, YZ, custom)
4. Gimbal lock comparison mode:
   - Split view: Euler angles (left) vs Rotor (right)
   - At pitch=90deg, Euler loses a degree of freedom
   - Rotor always works smoothly
5. Rotation composition demo (R2 * R1 != R1 * R2)
6. Educational content explaining:
   - Why rotors avoid gimbal lock
   - Half-angle encoding
   - Sandwich product

#### Implementation Approach

Given WASM requirement, prefer **custom 3D projection** approach:

```rust
// Simple 3D projection - works everywhere
pub struct Camera3D {
    pub eye: [f32; 3],      // Camera position
    pub target: [f32; 3],   // Look-at point
    pub up: [f32; 3],       // Up vector
    pub fov: f32,           // Field of view
}

impl Camera3D {
    /// Project 3D point to 2D screen coordinates
    pub fn project(&self, point: [f32; 3], aspect: f32) -> [f64; 2] {
        // View matrix (look-at)
        // Projection matrix (perspective)
        // Return 2D coordinates for egui_plot
    }
}

// Render cube as 12 line segments
fn render_wireframe_cube(plot_ui: &mut PlotUi, camera: &Camera3D, vertices: &[[f32; 3]; 8]) {
    let edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),  // Front face
        (4, 5), (5, 6), (6, 7), (7, 4),  // Back face
        (0, 4), (1, 5), (2, 6), (3, 7),  // Connecting edges
    ];
    for (a, b) in edges {
        let p1 = camera.project(vertices[a], aspect);
        let p2 = camera.project(vertices[b], aspect);
        plot_ui.line(Line::new(PlotPoints::new(vec![p1, p2])));
    }
}
```

This approach:
- Uses only egui_plot (already works in WASM)
- No additional 3D library dependencies
- Simple to understand and maintain
- Sufficient for wireframe visualizations

### Phase 3: projective3 Demos (PRD-48.4)

After euclidean3 validates the 3D framework, build PGA demos:

#### Demo Order
1. **projective3_geometry.rs** - Point-line-plane operations
   - Join and meet operations
   - Interactive point/line/plane placement
   - Establishes interaction patterns for 3D

2. **projective3_lines.rs** - Plucker coordinates
   - Line visualization with direction/moment
   - Line-line and line-plane meets

3. **projective3_motor.rs** - Rigid body motion
   - Motor visualization (screw axis)
   - Screw motion animation
   - Motor interpolation

4. **projective3_robot.rs** - Robot arm
   - Multi-joint robot
   - Forward kinematics via motors
   - Most complex, built last

### Phase 4: conformal3 Demos (PRD-48.5)

CGA demos require sphere/circle rendering:

#### Demo Order
1. **conformal3_circles.rs** - Circle from 3 points
   - Circle rendering in 3D (ring)
   - Tests circle extraction from CGA

2. **conformal3_inversion.rs** - Sphere inversion
   - Sphere rendering (wireframe)
   - Inversion transformations

3. **conformal3_mobius.rs** - Mobius transformations
   - Transformation composition UI
   - Animation

4. **conformal3_apollonian.rs** - Apollonian gasket
   - Performance-critical (many spheres)
   - Optional/advanced

### Phase 5: Specialty 3D Demos

1. **quaternion.rs** (PRD-48.6)
   - SLERP visualization
   - Uses euclidean3 infrastructure

2. **dualquat.rs** (PRD-48.6)
   - Screw motion (similar to projective3_motor)

3. **elliptic2.rs** (PRD-48.8)
   - Spherical geometry (needs sphere rendering)

4. **minkowski3.rs** (PRD-48.9)
   - 3D+time visualization
   - Light cones in 3D

---

## Updated Success Criteria

### Infrastructure
- [x] 3D viewport integrated with egui (camera3d.rs)
- [x] Camera controls work (orbit, zoom, pan) (camera_response, camera_controls)
- [x] All 3D shape primitives render correctly (shapes3d.rs)
- [x] 3D demos work on native desktop (tested)
- [ ] 3D demos work in WASM (needs validation with `trunk serve`)

### Demos
- [x] euclidean3 demonstrates rotor rotation
- [ ] projective3 demonstrates PGA operations
- [ ] conformal3 demonstrates sphere/circle operations
- [ ] quaternion demonstrates SLERP
- [x] euclidean3 has educational content

### Performance
- [x] 3D demos run at 60fps with reasonable object counts
- [x] Camera controls are responsive

---

## Risk Assessment

### ~~High Risk: three-d + WASM Compatibility~~ RESOLVED
~~three-d uses wgpu which has WASM support but it's more complex than 2D.~~

**Update**: three-d uses WebGL (not wgpu), and the library compiles for WASM.
We now have TWO approaches:
1. Custom projection via egui_plot - proven WASM compatible
2. three-d library - compiles for WASM, needs runtime validation

Both approaches available. Use custom projection for simple wireframes,
three-d for demos that need real lighting/materials.

### Medium Risk: 3D Interaction Complexity
Selecting/dragging objects in 3D is harder than 2D.
**Mitigation**: Start with simple interactions (no 3D dragging), add ray-casting later.

### Medium Risk: Demo File Size
3D demos will likely be even larger than 2D.
**Mitigation**: Extract shared code to common modules aggressively.

### Low Risk: Learning Curve
three-d API is well-documented.
**Mitigation**: Start with their examples, adapt to our needs.

---

## Implementation Timeline

| Phase | Deliverable | Depends On |
|-------|-------------|------------|
| 1a | three-d dependency + viewport3d | - |
| 1b | camera.rs + shapes3d.rs | 1a |
| 1c | test_3d.rs validation | 1a, 1b |
| 2 | euclidean3.rs | Phase 1 |
| 3a | projective3_geometry.rs | Phase 2 |
| 3b | projective3_lines.rs | 3a |
| 3c | projective3_motor.rs | 3b |
| 3d | projective3_robot.rs | 3c |
| 4a | conformal3_circles.rs | Phase 2 |
| 4b | conformal3_inversion.rs | 4a |
| 4c | conformal3_mobius.rs | 4b |
| 5 | quaternion, dualquat, etc. | Phase 2 |

---

## Candidates for API Extraction (from 2D development)

Methods discovered during 2D visualization that should be promoted:

### Already Useful (Add Now)
- `Circle::from_center_radius()` - Used in conformal2
- `Circle::center()` / `Circle::radius()` - Used in conformal2
- `PointPair::to_points()` - Used in conformal2_intersection

### 3D Candidates (Extract During 3D Development)
- Motor decomposition (rotation angle, translation vector)
- Motor interpolation (slerp-like)
- Plucker line distance calculations
- Sphere extraction methods for CGA

---

## Open Questions

1. ~~**3D Library Choice**: Should we prototype with three-d before committing?~~ RESOLVED: Use custom projection
2. ~~**WASM Priority**: Should WASM work for 3D be deferred entirely?~~ RESOLVED: Custom projection works in WASM
3. **apollonian demos**: Worth the complexity or defer indefinitely?
4. **Non-Euclidean 2D**: hyperbolic2 and elliptic2 need Poincare disk model - prioritize?

---

## Lessons Learned from 3D Implementation (Phase 1-2)

### What Worked Well

#### 1. Custom 3D Projection via egui_plot
The decision to use custom projection instead of three-d was correct:
- No new dependencies, no WASM compatibility issues
- Perspective projection math is straightforward
- Wireframe rendering is sufficient for educational demos
- Same plot API as 2D demos - consistent developer experience

#### 2. Spherical Coordinates for Camera
Using `(distance, azimuth, elevation)` instead of raw eye position:
- Natural orbit controls (drag = change azimuth/elevation)
- Easy to clamp elevation to avoid gimbal lock
- Easy to implement zoom (just change distance)
- Reset to default view is trivial

#### 3. Using clifford Library Types
Refactoring camera3d.rs and shapes3d.rs to use `clifford::specialized::euclidean::dim3::Vector`:
- `forward.cross(up)` instead of manual cross product
- `vector.dot(other)` instead of manual dot product
- `vector.normalized()` instead of manual normalization
- Code is cleaner and validates our library API

#### 4. Hodge Dual for Axis/Bivector Conversion
In euclidean3.rs, using `right_complement()` to convert between rotation axis and rotation plane:
```rust
fn axis_bivector(&self) -> Bivector<f32> {
    self.axis_vector().right_complement()
}
```
This demonstrates the axis/plane duality that is fundamental to GA rotations.

### Issues Encountered and Fixes

#### 1. Right-Handed Coordinate System
**Problem**: Initial implementation appeared left-handed (X axis on wrong side).

**Root Cause**: Considered swapping `cross(forward, up)` to `cross(up, forward)`, but this was WRONG.

**Solution**: The correct formula for a right-handed system is:
```rust
let right = forward.cross(up).normalized();  // Correct!
```
This gives +X when looking down -Z with Y up. The issue was in camera interpretation, not the math.

#### 2. Auto-Resizing Plot Bounds
**Problem**: Plot auto-resized when cube rotated to edges, causing jarring visual jumps.

**Solution**: Set explicit fixed bounds on the plot:
```rust
Plot::new("euclidean3_view")
    .include_x(-3.0)
    .include_x(3.0)
    .include_y(-3.0)
    .include_y(3.0)
```

#### 3. Arbitrary Rotation Axis UX
**Problem**: Dropdown for axis selection (X, Y, Z, Diagonal) was limiting.

**Solution**: Three sliders (X, Y, Z) for arbitrary axis + preset buttons:
- Allows any rotation axis
- Normalizes internally
- Preset buttons for common cases (X, Y, Z, Diagonal)

#### 4. Angle Range
**Problem**: Initial angle range -PI to PI (-180 to 180) was confusing for rotation.

**Solution**: Use 0 to 360 degrees for full rotation visualization.

### Design Patterns Established

#### 1. Camera3D + camera_response Pattern
```rust
// In render():
let response = Plot::new("my_plot").show(ui, |plot_ui| {
    // Draw things using camera.project()
});
camera_response(&mut self.camera, &response.response, ui);

// In controls():
camera_controls(ui, &mut self.camera);
```

#### 2. 3D Shape Primitives Pattern
All shape functions take `camera: &Camera3D` as first parameter and return `Vec<Line>`:
```rust
fn wireframe_box_vertices(camera: &Camera3D, vertices: &[[f32; 3]; 8], color: Color32) -> Vec<Line>
fn arrow_3d(camera: &Camera3D, origin: [f32; 3], tip: [f32; 3], color: Color32) -> Vec<Line>
fn coordinate_axes(camera: &Camera3D, length: f32) -> Vec<Line>
```

#### 3. Color Convention
RGB -> XYZ consistently throughout all visualizations:
- Red = X axis
- Green = Y axis
- Blue = Z axis

### Recommendations for Future 3D Demos

1. **Always use fixed plot bounds** to prevent auto-resizing artifacts

2. **Use clifford types** for all vector math - validates our API and keeps code clean

3. **Demonstrate GA operations** like Hodge dual where appropriate

4. **Test coordinate handedness early** - rotate a known shape and verify X/Y/Z directions

5. **Keep camera controls consistent** - orbit, zoom, pan should work the same across all 3D demos

---

## Lessons Learned from three-d Integration

### three-d Library Setup

1. **Feature flag**: Add `three-d = ["dep:three-d", "three-d?/egui-gui"]` for egui integration
2. **WASM getrandom**: Add `getrandom = { version = "0.3", features = ["wasm_js"] }` for WASM target
3. **Separate window**: three-d manages its own window, not compatible with eframe

### Flat Shading for Cubes

**Problem**: `CpuMesh::cube()` has shared vertices at corners, so `compute_normals()`
averages them into useless values.

**Solution**: Create cube mesh with 24 vertices (4 per face) with explicit per-face normals:
```rust
let faces = [
    (4, 5, 6, 7, Vector::new(0.0, 0.0, 1.0)),   // Front +Z
    // ... 6 faces total
];
```

### Transformed Normals

**Critical**: When transforming vertices with a rotor, **normals must also be transformed**:
```rust
let rotated_normal = rotor.transform(&base_normal);
```
Without this, shading appears flat/wrong because normals point in original directions.

### Showcasing Clifford Types

Keep all geometry as clifford types until render time:
- Cube vertices: `[Vector<f32>; 8]`
- Axis direction: `Vector<f32>`
- Rotation plane: `Bivector<f32>` (Hodge dual of axis)
- All transforms: `rotor.transform(&vector)`

Only convert to three-d's `Vec3`/`Mat4` in rendering helper functions.

### UI for GA Education

Display GA values prominently:
- Axis vector: `0.00e1 + 0.00e2 + 1.00e3`
- Plane bivector: `0.00e23 + 0.00e31 + 1.00e12` (Hodge dual)
- Rotor: `0.707 + 0.00e23 + 0.00e31 + 0.707e12`
- Rotor norm: `|R| = 1.0000`

### Lighting Setup

Good defaults for educational demos:
```rust
let ambient = AmbientLight::new(&context, 0.3, Srgba::WHITE);
let key_light = DirectionalLight::new(&context, 2.5, Srgba::WHITE, vec3(-1.0, -1.0, -1.0));
let fill_light = DirectionalLight::new(&context, 1.0, Srgba::new_opaque(200, 200, 255), vec3(1.0, 0.5, 1.0));
```

Material with some specularity:
```rust
PhysicalMaterial::new_opaque(context, &CpuMaterial {
    albedo: color,
    roughness: 0.4,
    metallic: 0.1,
    ..Default::default()
})
