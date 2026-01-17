# PRD-48.4: Projective 3D (PGA) Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for 3D projective geometric algebra

## Overview

Demonstrate 3D PGA with Plücker line coordinates, plane geometry, and motor-based rigid body transformations. This is the flagship demo for robotics and graphics applications.

## Demo 1: Plücker Line Coordinates

**File**: `examples/visualization/projective3_lines.rs`

### Features

1. **Line Definition**
   - Create lines via two points
   - Display Plücker coordinates: direction (d) and moment (m)
   - Show geometric meaning: `m = d × p₀` where p₀ is point on line

2. **Interactive Manipulation**
   - Drag line endpoints
   - See Plücker coordinates update in real-time
   - Visualize direction vector and moment vector

3. **Line-Line Meet**
   - Two lines in general position → no intersection
   - Coplanar lines → intersection point
   - Show distance between skew lines

4. **Line-Plane Meet**
   - Line intersects plane at a point
   - Parallel line → point at infinity

### Plücker Coordinate Display

```
Line L through points A(1,0,0) and B(0,1,1):

Direction: d = B - A = (-1, 1, 1)  [normalized]
Moment:    m = A × d = (1, 1, 1)

Plücker coords: L = [dx, dy, dz, mx, my, mz]
                  = [-0.58, 0.58, 0.58, 0.58, 0.58, 0.58]

Constraint: d · m = 0  ✓ (lines satisfy this)
```

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  Plücker Line Coordinates                              [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Line Definition       │
│         B●                          │                       │
│          ╲                          │ Point A: [x][y][z]    │
│           ╲  ← line L               │ Point B: [x][y][z]    │
│            ╲                        │                       │
│             ●A                      │ Plücker Coordinates   │
│              ↘ direction d          │ Direction:            │
│               ○→ moment m           │   d = (-0.58, 0.58, 0.58)
│                                     │ Moment:               │
│                                     │   m = (0.58, 0.58, 0.58)
│      ════════════════ plane        │                       │
│              ●P intersection        │ d·m = 0.000 ✓         │
│                                     │                       │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 2: Rigid Body Motion

**File**: `examples/visualization/projective3_motor.rs`

### Features

1. **Motor Visualization**
   - Show screw axis (the line around which rotation occurs)
   - Show pitch (translation per rotation)
   - Animate object along screw trajectory

2. **Motor Decomposition**
   - Extract rotation angle
   - Extract translation vector
   - Show rotation axis

3. **Motor Interpolation**
   - Smooth interpolation between two poses
   - Compare with naive matrix interpolation

4. **Object Transformation**
   - Transform coordinate frame (three orthogonal arrows)
   - Transform wireframe box
   - Trail showing motion path

### Screw Motion Visualization

```
                    screw axis
                       ║
                       ║
    ┌─────┐           ║           ┌─────┐
    │     │           ║           │     │
    │  1  │ ─────────►║──────────►│  2  │  final
    │     │           ║           │     │
    └─────┘           ║           └─────┘
     start            ║
                      ║
                      ↓ pitch = translation/rotation

Motor: M = exp(d·L/2) where L is screw axis line, d is displacement
```

### Motor Controls

```
┌───────────────────────────────────────┐
│ Motor Parameters                      │
│                                       │
│ Screw Axis:                          │
│   Point: [0.0] [0.0] [0.0]           │
│   Direction: [0.0] [0.0] [1.0]       │
│                                       │
│ Displacement: [========●=====] 2.5   │
│                                       │
│ ─── Decomposition ───                │
│ Rotation: 143.2°                     │
│ Translation: (0.0, 0.0, 1.2)         │
│                                       │
│ [▶ Animate] [Reset] Speed: [●==]     │
└───────────────────────────────────────┘
```

---

## Demo 3: Point-Line-Plane Operations

**File**: `examples/visualization/projective3_geometry.rs`

### Features

1. **Join Operations**
   - Point ∧ Point = Line
   - Point ∧ Line = Plane
   - Line ∧ Plane = (degenerate)

2. **Meet Operations**
   - Plane ∨ Plane = Line
   - Line ∨ Plane = Point
   - Line ∨ Line = Point (if coplanar)

3. **Interactive Construction**
   - Click to place points
   - Join selected points to create lines
   - Join point and line to create planes
   - Meet planes to find intersection lines

### Operation Reference

```
JOIN (∧) - creates higher-grade element
┌─────────────────────────────────────┐
│ Point ∧ Point = Line                │
│   ●────────●                        │
│                                     │
│ Point ∧ Line = Plane                │
│      ●                              │
│     ╱│╲                             │
│ ───●─┼─●───                         │
│     ╲│╱                             │
└─────────────────────────────────────┘

MEET (∨) - finds intersection
┌─────────────────────────────────────┐
│ Plane ∨ Plane = Line                │
│    ╲ ╱                              │
│     ╳────── line                    │
│    ╱ ╲                              │
│                                     │
│ Line ∨ Plane = Point                │
│      │                              │
│ ─────●───── plane                   │
│      │                              │
└─────────────────────────────────────┘
```

---

## Demo 4: 3D Robot Arm

**File**: `examples/visualization/projective3_robot.rs`

### Features

1. **Multi-Joint Robot**
   - 3+ joints with configurable axes
   - Each joint is a motor (rotation or translation)
   - Forward kinematics via motor chain

2. **Interactive Control**
   - Sliders for each joint angle
   - Real-time end effector update
   - Joint limit visualization

3. **Workspace Visualization**
   - Sample end effector positions
   - Show reachable workspace boundary

### Robot Structure

```
            ┌───┐
            │ J3│ ← wrist
            └─┬─┘
              │
            ┌─┴─┐
            │ J2│ ← elbow
            └─┬─┘
              │
              │
            ┌─┴─┐
            │ J1│ ← shoulder
            └─┬─┘
              │
           ═══╪═══ base

Joint motors:
  M₁ = rotation about Z axis
  M₂ = rotation about Y axis
  M₃ = rotation about Y axis

End effector: M_total = M₁ · T₁ · M₂ · T₂ · M₃ · T₃
```

---

## Implementation Tasks

### projective3_lines.rs
1. [ ] 3D viewport with camera controls
2. [ ] Line rendering (as cylinder or thick line)
3. [ ] Point rendering (as sphere)
4. [ ] Two-point line creation
5. [ ] Plücker coordinate display
6. [ ] Direction/moment vector visualization
7. [ ] Line-line distance calculation
8. [ ] Line-plane intersection

### projective3_motor.rs
1. [ ] Wireframe box rendering
2. [ ] Coordinate frame rendering (RGB arrows)
3. [ ] Screw axis line rendering
4. [ ] Motor parameter controls
5. [ ] Screw motion animation
6. [ ] Motion trail rendering
7. [ ] Motor decomposition display
8. [ ] Motor interpolation demo

### projective3_geometry.rs
1. [ ] Interactive point placement (click on plane)
2. [ ] Point selection and multi-select
3. [ ] Join operation (point∧point, point∧line)
4. [ ] Meet operation (plane∨plane, line∨plane)
5. [ ] Plane rendering (as bounded quad)
6. [ ] Operation history/undo

### projective3_robot.rs
1. [ ] Joint definition structure
2. [ ] Forward kinematics via motor chain
3. [ ] Joint angle sliders
4. [ ] Link and joint rendering
5. [ ] End effector position display
6. [ ] Workspace sampling (optional)

## Verification

```bash
cargo run --example projective3_lines --release
cargo run --example projective3_motor --release
cargo run --example projective3_geometry --release
cargo run --example projective3_robot --release
```

### Checklist
- [ ] Lines display correct Plücker coordinates
- [ ] Plücker constraint d·m=0 satisfied for all lines
- [ ] Motor animation follows screw trajectory
- [ ] Join/meet operations produce correct results
- [ ] Robot forward kinematics matches expected positions
- [ ] Camera controls are intuitive (orbit, pan, zoom)
