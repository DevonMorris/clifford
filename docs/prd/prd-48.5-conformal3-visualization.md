# PRD-48.5: Conformal 3D (CGA) Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for 3D conformal geometric algebra

## Overview

Conformal GA is the most powerful algebra in the clifford library, capable of representing spheres, circles, and planes as first-class objects, with transformations that preserve angles. This demo showcases the "magic" of CGA.

## Demo 1: Sphere Inversion

**File**: `examples/visualization/conformal3_inversion.rs`

### The Key Insight

In CGA, inversion through a sphere is a simple reflection:
- Spheres map to spheres
- Planes map to spheres (unless through inversion center → plane)
- Lines map to circles (unless through center → line)
- Points map to points

### Features

1. **Inversion Sphere Control**
   - Drag to move center
   - Slider for radius
   - Toggle visibility

2. **Object Placement**
   - Add spheres, planes, lines, circles
   - Each object shows its inverted image

3. **Real-Time Transformation**
   - Move inversion sphere → all images update
   - Animate inversion sphere movement

4. **Special Cases**
   - Object through center → maps to infinity
   - Plane through center → stays a plane
   - Line through center → stays a line

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  Sphere Inversion                                      [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Inversion Sphere      │
│                                     │ Center: [x][y][z]     │
│         ╭─────╮                     │ Radius: [====●====]   │
│        ╱   I   ╲    inversion       │                       │
│       │    ●    │   sphere          │ Objects               │
│        ╲       ╱                    │ [+Sphere] [+Plane]    │
│         ╰─────╯                     │ [+Line] [+Circle]     │
│              ↓                      │                       │
│         ╭───────╮                   │ Transformations:      │
│        ╱         ╲  inverted        │ ● Sphere → Sphere     │
│       │     ○     │ sphere          │ ○ Plane → Sphere      │
│        ╲         ╱  (larger!)       │ ○ Line → Circle       │
│         ╰───────╯                   │                       │
│                                     │ [▶ Animate]           │
└─────────────────────────────────────┴───────────────────────┘
```

### Transformation Rules Display

```
┌────────────────────────────────────────┐
│ Inversion Transformations              │
│                                        │
│ Object          → Image               │
│ ─────────────────────────────────────  │
│ Sphere (outside) → Sphere (inside)    │
│ Sphere (inside)  → Sphere (outside)   │
│ Sphere (through) → Sphere (through)   │
│ Plane (not thru) → Sphere (thru ctr)  │
│ Plane (thru ctr) → Plane (same)       │
│ Line (not thru)  → Circle (thru ctr)  │
│ Line (thru ctr)  → Line (same)        │
│ Circle           → Circle             │
│ Point            → Point              │
└────────────────────────────────────────┘
```

---

## Demo 2: Circle Operations

**File**: `examples/visualization/conformal3_circles.rs`

### Features

1. **Circle from Three Points**
   - Place three points → unique circle
   - `Circle = P₁ ∧ P₂ ∧ P₃`

2. **Circle-Circle Intersection**
   - Two circles in general position → two points
   - Coplanar circles → two points
   - Tangent circles → one point

3. **Circle-Sphere Intersection**
   - Circle intersecting sphere → two points
   - Circle on sphere → circle (degenerate)

4. **Circle Visualization**
   - Render as 3D ring
   - Show center and normal
   - Display radius

### Implementation

```rust
use clifford::specialized::conformal::dim3::{RoundPoint, Circle, Sphere};

pub struct CircleDemo {
    points: Vec<RoundPoint<f32>>,
    circles: Vec<DerivedCircle>,
}

struct DerivedCircle {
    from_points: (usize, usize, usize),
    circle: Circle<f32>,
}

impl CircleDemo {
    fn create_circle(&mut self, p1: usize, p2: usize, p3: usize) {
        let circle = self.points[p1]
            .wedge(&self.points[p2])
            .wedge(&self.points[p3]);
        self.circles.push(DerivedCircle {
            from_points: (p1, p2, p3),
            circle,
        });
    }
}
```

---

## Demo 3: Möbius Transformations

**File**: `examples/visualization/conformal3_mobius.rs`

### Features

1. **Transformation Gallery**
   - Translation (move everything)
   - Rotation (rotate everything)
   - Dilation (scale from point)
   - Inversion (reflect through sphere)
   - General Möbius (composition)

2. **Interactive Composition**
   - Build up transformation from primitives
   - See transformation matrix equivalent
   - Apply to scene objects

3. **Animation**
   - Smooth transformation interpolation
   - Show transformation path

### Transformation Types

```
┌────────────────────────────────────────┐
│ Möbius Transformations                 │
│                                        │
│ Translation T(v):                      │
│   Moves all objects by vector v        │
│   Preserves: all shapes and sizes      │
│                                        │
│ Rotation R(θ, axis):                   │
│   Rotates around axis through origin   │
│   Preserves: all shapes and sizes      │
│                                        │
│ Dilation D(k, center):                 │
│   Scales by factor k from center       │
│   Preserves: angles, not sizes         │
│                                        │
│ Inversion I(sphere):                   │
│   Reflects through sphere              │
│   Preserves: angles (conformal!)       │
│                                        │
│ General Möbius = composition of above  │
└────────────────────────────────────────┘
```

---

## Demo 4: Apollonian Gasket (Advanced)

**File**: `examples/visualization/conformal3_apollonian.rs`

### Features

1. **Initial Configuration**
   - Four mutually tangent circles
   - Computed via CGA tangency conditions

2. **Recursive Construction**
   - For any three mutually tangent circles, find the two circles tangent to all three
   - Use CGA meet operations
   - Recursively fill gaps

3. **Interactive Depth**
   - Slider for recursion depth
   - Color by generation
   - Zoom to see detail

### Apollonian Gasket

```
         ╭─────────────────────────╮
        ╱  ○   ○   ○   ○   ○   ○   ╲
       ╱ ○   ○   ○   ○   ○   ○   ○  ╲
      │ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ │
      │○   ╭───╮       ╭───╮       ○│
      │ ○ ╱     ╲     ╱     ╲  ○  ○ │
      │○ │   ○   │   │   ○   │ ○  ○│
      │ ○ ╲     ╱     ╲     ╱  ○  ○ │
      │○   ╰───╯       ╰───╯       ○│
      │ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ │
       ╲ ○   ○   ○   ○   ○   ○   ○  ╱
        ╲  ○   ○   ○   ○   ○   ○   ╱
         ╰─────────────────────────╯

Each circle tangent to its neighbors!
Computed via CGA tangency conditions.
```

---

## Implementation Tasks

### conformal3_inversion.rs
1. [ ] Sphere rendering (wireframe or solid)
2. [ ] Inversion sphere controls
3. [ ] Object placement UI
4. [ ] Inversion computation for each object type
5. [ ] Inverted object rendering (different color)
6. [ ] Special case handling (through center)
7. [ ] Animation support

### conformal3_circles.rs
1. [ ] Point placement in 3D
2. [ ] Circle from three points
3. [ ] Circle rendering (3D ring)
4. [ ] Circle-circle intersection
5. [ ] Circle-sphere intersection
6. [ ] Circle properties display (center, radius, normal)

### conformal3_mobius.rs
1. [ ] Translation motor
2. [ ] Rotation motor
3. [ ] Dilation motor
4. [ ] Inversion transformation
5. [ ] Composition UI
6. [ ] Apply to scene objects
7. [ ] Transformation animation

### conformal3_apollonian.rs
1. [ ] Initial four-circle configuration
2. [ ] Tangent circle computation
3. [ ] Recursive circle placement
4. [ ] Depth control slider
5. [ ] Generation-based coloring
6. [ ] Zoom controls

## Verification

```bash
cargo run --example conformal3_inversion --release
cargo run --example conformal3_circles --release
cargo run --example conformal3_mobius --release
cargo run --example conformal3_apollonian --release
```

### Checklist
- [ ] Sphere inversion produces correct transformed objects
- [ ] Circle from three points passes through all three
- [ ] Circle-circle intersection finds correct points
- [ ] Möbius transformations compose correctly
- [ ] Apollonian gasket circles are all mutually tangent
- [ ] Transformations preserve angles (conformal property)
