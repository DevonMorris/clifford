# PRD-48.3: Projective 2D (PGA) Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for 2D projective geometric algebra

## Overview

Demonstrate the power of 2D PGA for point-line geometry and rigid transformations. Show join (wedge) and meet (antiwedge) operations visually.

## Demo: Interactive Point-Line Geometry

**File**: `examples/visualization/projective2.rs`

### Features

1. **Point Creation and Manipulation**
   - Click to create points
   - Drag to move points
   - Delete points with right-click

2. **Join Operation (Wedge)**
   - Select two points → shows line through them
   - `Line = Point₁ ∧ Point₂`
   - Line updates in real-time as points move

3. **Meet Operation (Antiwedge)**
   - Select two lines → shows intersection point
   - `Point = Line₁ ∨ Line₂`
   - Handles parallel lines (point at infinity)

4. **Motor Animation**
   - Apply rigid transformation to all objects
   - Decompose motor into rotation + translation
   - Animate motor interpolation

5. **Line Representation Display**
   - Show line in form `ax + by + c = 0`
   - Show distance from origin
   - Show normal vector

### Implementation

```rust
use clifford::specialized::projective::dim2::{Point, Line, Motor};
use clifford::ops::{Wedge, Antiwedge};

pub struct Projective2Demo {
    points: Vec<DraggablePoint>,
    derived_lines: Vec<DerivedLine>,
    derived_points: Vec<DerivedPoint>,
    selection: Selection,
    motor: Motor<f32>,
    motor_animation: Animation,
    show_coordinates: bool,
}

struct DraggablePoint {
    id: usize,
    point: Point<f32>,
    selected: bool,
}

struct DerivedLine {
    from_points: (usize, usize),
    line: Line<f32>,
}

struct DerivedPoint {
    from_lines: (usize, usize),
    point: Point<f32>,
}

enum Selection {
    None,
    Point(usize),
    TwoPoints(usize, usize),  // Ready for join
    Line(usize),
    TwoLines(usize, usize),   // Ready for meet
}

impl Projective2Demo {
    fn compute_join(&mut self, p1_idx: usize, p2_idx: usize) {
        let p1 = &self.points[p1_idx].point;
        let p2 = &self.points[p2_idx].point;
        let line = p1.wedge(p2);
        self.derived_lines.push(DerivedLine {
            from_points: (p1_idx, p2_idx),
            line,
        });
    }

    fn compute_meet(&mut self, l1_idx: usize, l2_idx: usize) {
        let l1 = &self.derived_lines[l1_idx].line;
        let l2 = &self.derived_lines[l2_idx].line;
        let point = l1.antiwedge(l2);
        self.derived_points.push(DerivedPoint {
            from_lines: (l1_idx, l2_idx),
            point,
        });
    }

    fn apply_motor(&self, point: &Point<f32>) -> Point<f32> {
        self.motor.transform(point)
    }
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Projective 2D - Point-Line Geometry                   [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Tools                 │
│                                     │ ○ Select  ● Add Point │
│     A ●───────────────────● B       │                       │
│              ╲    ╱                 │ Operations            │
│               ╲  ╱                  │ [Join A∧B] [Meet L∨M] │
│                ╳ P                  │                       │
│               ╱  ╲                  │ Selected: Point A     │
│              ╱    ╲                 │ x: 1.50  y: 2.00      │
│     C ●───────────────────● D       │ w: 1.00 (normalized)  │
│                                     │                       │
│                                     │ ─────────────────     │
│                                     │ Motor Controls        │
│                                     │ Rotation: [===●===]   │
│                                     │ Trans X:  [===●===]   │
│                                     │ Trans Y:  [===●===]   │
│                                     │ [▶ Animate Motor]     │
│                                     │                       │
└─────────────────────────────────────┴───────────────────────┘
```

### Interactive Workflow

```
Step 1: Create points A and B
┌─────────────────┐
│                 │
│  A ●      ● B   │
│                 │
└─────────────────┘

Step 2: Select both → Join creates line
┌─────────────────┐
│                 │
│  A ●──────● B   │  Line AB = A ∧ B
│                 │
└─────────────────┘

Step 3: Create points C and D, join them
┌─────────────────┐
│  A ●──────● B   │
│       ╲ ╱       │  Line CD = C ∧ D
│  C ●───╳───● D  │
│       intersection
└─────────────────┘

Step 4: Meet lines → intersection point P
┌─────────────────┐
│  A ●──────● B   │
│       ╲ ╱       │  P = AB ∨ CD
│  C ●───●───● D  │
│        P        │
└─────────────────┘
```

### Motor Animation

Show how motors smoothly transform all geometry:

```
Frame 0:          Frame 1:          Frame 2:
●───●             ●───●               ●───●
 ╲ ╱     Motor     ╲ ╱      Motor      ╲ ╱
  ●     ────→       ●      ────→        ●
                  (rotated)         (rotated+translated)
```

Motor decomposition display:
```
Motor M = e^(d·L/2)

Rotation component: 45.0°
Translation: (2.0, 1.0)
Screw axis: through (0.5, 0.3)

M = [s: 0.924, tx: 0.383, ty: 0.924, r: 0.383]
```

---

## Advanced Demo: 2D Robot Arm

**File**: `examples/visualization/projective2_robot.rs`

### Features

1. **Two-Link Robot Arm**
   - Base joint (rotation only)
   - Elbow joint (rotation only)
   - End effector position display

2. **Forward Kinematics via Motors**
   - `M_total = M_base · M_link1 · M_elbow · M_link2`
   - Each joint is a motor

3. **Interactive Control**
   - Sliders for joint angles
   - Drag end effector (inverse kinematics hint)

### Implementation

```rust
pub struct Robot2D {
    joint1_angle: f32,  // Base rotation
    joint2_angle: f32,  // Elbow rotation
    link1_length: f32,
    link2_length: f32,
}

impl Robot2D {
    fn forward_kinematics(&self) -> Point<f32> {
        // Motor for base rotation
        let m1 = Motor::from_rotation(self.joint1_angle);

        // Motor for translation along link1
        let t1 = Motor::from_translation(self.link1_length, 0.0);

        // Motor for elbow rotation
        let m2 = Motor::from_rotation(self.joint2_angle);

        // Motor for translation along link2
        let t2 = Motor::from_translation(self.link2_length, 0.0);

        // Compose: M_total = M1 · T1 · M2 · T2
        let m_total = m1.geometric_product(&t1)
                       .geometric_product(&m2)
                       .geometric_product(&t2);

        // Transform origin to get end effector
        m_total.transform(&Point::origin())
    }
}
```

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  2D Robot Arm - Motor Kinematics                       [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Joint Controls        │
│                                     │                       │
│             ○ end effector          │ Joint 1: [====●====]  │
│            ╱                        │          45.0°        │
│           ╱ link2                   │                       │
│          ●───────────               │ Joint 2: [====●====]  │
│         ╱   elbow                   │          -30.0°       │
│        ╱                            │                       │
│       ╱ link1                       │ Link lengths          │
│      ●                              │ Link 1: [====●====]   │
│      │ base                         │ Link 2: [====●====]   │
│   ═══╪═══                           │                       │
│                                     │ End effector:         │
│                                     │ x: 3.21  y: 2.45      │
│                                     │                       │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Implementation Tasks

### projective2.rs (Main Demo)
1. [ ] Point creation on click
2. [ ] Point dragging
3. [ ] Point selection (single and multi)
4. [ ] Join operation (wedge) visualization
5. [ ] Meet operation (antiwedge) visualization
6. [ ] Line rendering (infinite, clipped to view)
7. [ ] Motor controls (rotation + translation sliders)
8. [ ] Motor application to all objects
9. [ ] Motor animation
10. [ ] Coordinate display panel

### projective2_robot.rs (Advanced)
1. [ ] Robot arm structure
2. [ ] Forward kinematics via motor composition
3. [ ] Joint angle sliders
4. [ ] Link visualization
5. [ ] End effector position display
6. [ ] Workspace boundary visualization (optional)

## Verification

```bash
cargo run --example projective2 --release
cargo run --example projective2_robot --release
```

### Checklist
- [ ] Click creates points at correct coordinates
- [ ] Dragging points updates derived lines in real-time
- [ ] Join creates correct line through two points
- [ ] Meet finds correct intersection (or point at infinity for parallel)
- [ ] Motor smoothly transforms all geometry
- [ ] Robot arm end effector matches manual calculation
- [ ] Motor composition order is correct (right-to-left)
