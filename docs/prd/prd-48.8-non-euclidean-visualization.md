# PRD-48.8: Non-Euclidean Geometry Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for hyperbolic and elliptic/spherical geometries

## Overview

Visualize the two non-Euclidean geometries: hyperbolic (negative curvature, Lobachevsky plane) and elliptic (positive curvature, spherical geometry). Show how parallel postulate differs in each.

---

## Demo 1: Poincaré Disk (Hyperbolic)

**File**: `examples/visualization/hyperbolic2_poincare.rs`

### Features

1. **Geodesic Lines**
   - Lines are circular arcs perpendicular to boundary
   - Lines through center are straight
   - Click to create geodesics

2. **Parallel Lines Demo**
   - Given line L and point P not on L
   - Infinitely many lines through P not intersecting L!
   - Show limiting parallel lines

3. **Hyperbolic Triangles**
   - Sum of angles < 180°
   - Display angle sum
   - Compare with Euclidean triangle

4. **Transformations**
   - Hyperbolic rotations (about a point)
   - Hyperbolic translations (along a geodesic)
   - Animate transformations

### The Parallel Postulate

```
Euclidean (flat):              Hyperbolic (curved):
┌────────────────────┐         ╭──────────────────╮
│                    │        ╱                    ╲
│      ● P           │       │        ● P          │
│      │             │       │       ╱│╲           │
│      │ one         │       │      ╱ │ ╲          │
│ ─────┼───── L      │       │  ╭──╱──┼──╲──╮      │
│      │ parallel    │       │  │ ╱   │   ╲ │      │
│                    │        ╲╱     │     ╲      ╱
└────────────────────┘         ╰─────┴─────╯
                              infinitely many
                              parallels!
```

### Implementation

```rust
use clifford::specialized::hyperbolic::dim2::{Point, Line, Rotor};

pub struct PoincareDiskDemo {
    points: Vec<HyperbolicPoint>,
    geodesics: Vec<HyperbolicGeodesic>,
    triangles: Vec<HyperbolicTriangle>,
    mode: InteractionMode,
}

struct HyperbolicPoint {
    // Poincaré disk coordinates (x, y) with x² + y² < 1
    x: f64,
    y: f64,
}

impl HyperbolicPoint {
    /// Convert to hyperboloid model point
    fn to_hyperboloid(&self) -> Point<f64> {
        let r2 = self.x * self.x + self.y * self.y;
        let t = (1.0 + r2) / (1.0 - r2);
        let x = 2.0 * self.x / (1.0 - r2);
        let y = 2.0 * self.y / (1.0 - r2);
        Point::new(x, y, t)
    }
}

/// Geodesic in Poincaré disk is a circular arc perpendicular to boundary
fn draw_geodesic(p1: &HyperbolicPoint, p2: &HyperbolicPoint) -> impl Iterator<Item = [f64; 2]> {
    // Calculate circle through p1, p2 perpendicular to unit circle
    // Special case: if p1, p2, origin collinear → straight line
    todo!()
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Poincaré Disk - Hyperbolic Geometry                   [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Tools                 │
│         ╭───────────────╮           │ ○ Add Point           │
│        ╱     ╭───╮       ╲          │ ● Add Geodesic        │
│       │     ╱  A  ╲       │         │ ○ Add Triangle        │
│       │    ╱───────╲      │         │                       │
│       │   │    ╲    │     │         │ Selected Triangle:    │
│       │   │ B   ╲   │ C   │         │  Angle A: 35°         │
│       │   │      ╲  │     │         │  Angle B: 42°         │
│       │    ╲──────╲╱      │         │  Angle C: 38°         │
│        ╲    geodesics    ╱          │  Sum: 115° < 180°!    │
│         ╰───────────────╯           │                       │
│              boundary               │ Parallel Demo:        │
│           (points at ∞)             │ [Show Parallels]      │
│                                     │                       │
│                                     │ [▶ Animate Transform] │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 2: Hyperbolic Tessellation

**File**: `examples/visualization/hyperbolic2_tessellation.rs`

### Features

1. **Regular Tessellations**
   - {p, q} tessellation: p-gons, q meeting at each vertex
   - Hyperbolic when (p-2)(q-2) > 4
   - Examples: {5,4}, {7,3}, {4,5}

2. **Interactive Selection**
   - Choose p and q
   - See tessellation fill disk

3. **Transformation**
   - Drag to apply hyperbolic translation
   - See tessellation move

### Hyperbolic Tessellations

```
{5, 4} tessellation:          {7, 3} tessellation:
5-gons, 4 at each vertex      7-gons, 3 at each vertex

    ╱─────╲                       ╱───────╲
   ╱ ╱───╲ ╲                     ╱  ╱───╲  ╲
  │ │ ╱─╲ │ │                   │  │ ╱─╲ │  │
  │ │ │ │ │ │                   │  │ │ │ │  │
  │ │ ╲─╱ │ │                   │  │ ╲─╱ │  │
   ╲ ╲───╱ ╱                     ╲  ╲───╱  ╱
    ╲─────╱                       ╲───────╱

In Euclidean: only {3,6}, {4,4}, {6,3} work
In Hyperbolic: infinitely many tessellations!
```

---

## Demo 3: Elliptic/Spherical Geometry

**File**: `examples/visualization/elliptic2_sphere.rs`

### Features

1. **Spherical View**
   - 3D sphere with camera controls
   - Points on sphere surface
   - Great circles as "lines"

2. **Geodesic Lines**
   - All lines are great circles
   - Any two lines intersect (no parallels!)
   - Create geodesics through two points

3. **Spherical Triangles**
   - Sum of angles > 180°
   - Display angle sum and spherical excess
   - Area = spherical excess × R²

4. **Spherical Transformations**
   - Rotations (rotor-based)
   - Animate rotations

### The Parallel Postulate (Elliptic)

```
Euclidean (flat):              Elliptic (spherical):
┌────────────────────┐              ╭───────╮
│                    │             ╱    ●    ╲
│      ● P           │            │   ╱ P╲    │
│      │             │            │  ╱    ╲   │
│      │ one         │            │ ╱      ╲  │
│ ─────┼───── L      │           ─●────────●─ L
│      │ parallel    │            │          │
│                    │             ╲        ╱
└────────────────────┘              ╰──────╯
                               All lines intersect!
                               No parallels exist!
```

### Implementation

```rust
use clifford::specialized::elliptic::dim2::{Point, Line, Rotor};

pub struct SphericalDemo {
    points: Vec<SphericalPoint>,
    great_circles: Vec<GreatCircle>,
    triangles: Vec<SphericalTriangle>,
    camera: Camera3D,
}

struct SphericalPoint {
    // Homogeneous coordinates on unit sphere
    x: f64,
    y: f64,
    z: f64,
}

impl SphericalPoint {
    fn to_point(&self) -> Point<f64> {
        Point::new(self.x, self.y, self.z).normalize()
    }

    fn latitude(&self) -> f64 {
        self.z.asin()
    }

    fn longitude(&self) -> f64 {
        self.y.atan2(self.x)
    }
}

struct SphericalTriangle {
    vertices: [SphericalPoint; 3],
}

impl SphericalTriangle {
    fn angle_sum(&self) -> f64 {
        // Calculate angles using GA
        // Will be > π
        todo!()
    }

    fn spherical_excess(&self) -> f64 {
        self.angle_sum() - std::f64::consts::PI
    }

    fn area(&self, radius: f64) -> f64 {
        self.spherical_excess() * radius * radius
    }
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Spherical Geometry                                    [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Tools                 │
│           ╱───────╲                 │ ○ Add Point           │
│          ╱    A    ╲                │ ● Add Great Circle    │
│         │  ╱─────╲  │               │ ○ Add Triangle        │
│         │ ╱   ●   ╲ │               │                       │
│         │╱    │    ╲│               │ Selected Triangle:    │
│     ────●─────┼─────●────           │  Angle A: 72°         │
│         │╲    │    ╱│               │  Angle B: 68°         │
│         │ ╲   │   ╱ │               │  Angle C: 75°         │
│         │  ╲──┼──╱  │               │  Sum: 215° > 180°!    │
│          ╲   B│C   ╱                │                       │
│           ╲───────╱                 │ Spherical Excess:     │
│              3D sphere              │  35° = 0.61 rad       │
│        [drag to rotate]             │  Area: 0.61 R²        │
│                                     │                       │
│                                     │ [▶ Animate Rotation]  │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 4: Geometry Comparison

**File**: `examples/visualization/geometry_comparison.rs`

### Features

1. **Side-by-Side Views**
   - Euclidean | Hyperbolic | Elliptic
   - Same operations in each geometry

2. **Parallel Lines**
   - Euclidean: exactly one parallel
   - Hyperbolic: infinitely many
   - Elliptic: none

3. **Triangle Angle Sums**
   - Euclidean: = 180°
   - Hyperbolic: < 180°
   - Elliptic: > 180°

4. **Interactive Comparison**
   - Place points in one view
   - See corresponding construction in others

### Comparison Table

```
┌──────────────────┬─────────────┬─────────────┬─────────────┐
│ Property         │ Euclidean   │ Hyperbolic  │ Elliptic    │
├──────────────────┼─────────────┼─────────────┼─────────────┤
│ Curvature        │ 0           │ -1          │ +1          │
│ Parallel lines   │ exactly 1   │ infinitely  │ none        │
│                  │             │ many        │             │
│ Triangle angles  │ = 180°      │ < 180°      │ > 180°      │
│ Circles          │ constant    │ exp growth  │ sine growth │
│                  │ circumf.    │             │             │
│ Model            │ plane       │ Poincaré    │ sphere      │
│                  │             │ disk        │             │
└──────────────────┴─────────────┴─────────────┴─────────────┘
```

---

## Implementation Tasks

### hyperbolic2_poincare.rs
1. [ ] Poincaré disk boundary rendering
2. [ ] Point placement (constrained to disk)
3. [ ] Geodesic arc calculation
4. [ ] Geodesic rendering
5. [ ] Parallel line demonstration
6. [ ] Triangle angle calculation
7. [ ] Hyperbolic transformation animation

### hyperbolic2_tessellation.rs
1. [ ] {p, q} tessellation generation
2. [ ] p and q parameter selection
3. [ ] Tessellation rendering
4. [ ] Hyperbolic translation (drag)
5. [ ] Validation of hyperbolic condition

### elliptic2_sphere.rs
1. [ ] 3D sphere rendering
2. [ ] Camera orbit controls
3. [ ] Point placement on sphere
4. [ ] Great circle calculation and rendering
5. [ ] Spherical triangle angles
6. [ ] Spherical excess calculation
7. [ ] Rotation animation

### geometry_comparison.rs
1. [ ] Three-panel layout
2. [ ] Synchronized point placement
3. [ ] Parallel line demonstration
4. [ ] Triangle angle comparison
5. [ ] Property table display

## Verification

```bash
cargo run --example hyperbolic2_poincare --release
cargo run --example hyperbolic2_tessellation --release
cargo run --example elliptic2_sphere --release
cargo run --example geometry_comparison --release
```

### Checklist
- [ ] Poincaré disk geodesics are circular arcs ⊥ to boundary
- [ ] Lines through center are straight diameters
- [ ] Hyperbolic triangle angle sum < 180°
- [ ] {p,q} tessellation fills disk correctly
- [ ] Spherical great circles render correctly
- [ ] Spherical triangle angle sum > 180°
- [ ] Spherical excess × R² = area
- [ ] Comparison view shows all three geometries
