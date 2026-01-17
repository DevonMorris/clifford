# PRD-48.12: Conformal 2D (CGA2D) Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1, PRD-50
**Goal**: Interactive demos for 2D conformal geometric algebra

## Overview

2D Conformal Geometric Algebra (CGA2D) offers the elegance of conformal transformations in a more accessible 2D setting. This visualization showcases circles, lines, and points as first-class objects with interactive conformal transformations.

CGA2D is pedagogically valuable as a stepping stone before tackling the more complex 3D CGA, while still demonstrating the core conformal concepts: circle inversion, Mobius transformations, and unified treatment of circles and lines.

## Demo 1: Circle from Points

**File**: `crates/clifford-viz/examples/conformal2_circles.rs`

### The Key Insight

In CGA, three points uniquely determine a circle via the outer product:
```
Circle = P₁ ∧ P₂ ∧ P₃
```

This is the conformal analogue of "two points determine a line" - but now we can handle circles just as naturally.

### Features

1. **Interactive Point Placement**
   - Click to place points
   - Drag points to move them
   - Delete points with right-click

2. **Automatic Circle Generation**
   - Any three points → unique circle
   - Visual feedback showing which points define which circle
   - Color coding for multiple circles

3. **Circle Properties Display**
   - Center coordinates
   - Radius
   - Curvature (1/r)
   - Algebraic representation

4. **Collinear Detection**
   - When three points become collinear, circle → line
   - Visual transition as points approach collinearity
   - Display "Circle at infinity" when degenerate

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  Circle from Three Points                              [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Points                │
│                                     │ ● P₁: (2.0, 1.5)      │
│         ╭──────────────╮            │ ● P₂: (4.0, 3.0)      │
│        ╱                ╲           │ ● P₃: (3.0, 0.5)      │
│       │      P₂ ●        │          │                       │
│       │                  │          │ Circle Properties     │
│       │     ○ center     │          │ Center: (3.0, 1.5)    │
│       │                  │          │ Radius: 1.58          │
│        ╲   P₁ ●    ● P₃ ╱           │ Curvature: 0.63       │
│         ╰──────────────╯            │                       │
│                                     │ [Clear All]           │
│                                     │ [+ Add Point]         │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 2: Circle Inversion

**File**: `crates/clifford-viz/examples/conformal2_inversion.rs`

### The Key Insight

Inversion through a circle is a fundamental conformal transformation. In CGA, it's elegantly expressed as a reflection:
```
P' = C · P · C⁻¹
```

Inversion maps:
- Points → Points (reflected through circle)
- Circles → Circles (in general)
- Lines → Circles (through inversion center)
- Circles through center → Lines

### Features

1. **Inversion Circle Control**
   - Drag to move center
   - Slider/scroll for radius
   - Toggle visibility

2. **Object Library**
   - Add circles (click center + drag for radius)
   - Add lines (click two points)
   - Add points
   - Each object shows its inverted image

3. **Real-Time Transformation**
   - Move inversion circle → all images update
   - Animate inversion circle movement
   - Morph animation between original and inverted

4. **Special Cases Highlighted**
   - Line through center → stays a line
   - Circle through center → becomes a line
   - Point at center → maps to infinity (shown visually)

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  Circle Inversion                                      [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Inversion Circle      │
│    original        inverted         │ Center: [x]──●──[x]   │
│    circle          circle           │         [-2]   [2]    │
│    (blue)          (orange)         │         [y]──●──[y]   │
│         ╭──╮                        │         [-2]   [2]    │
│        ╱    ╲      ╭────────╮       │ Radius: [====●====]   │
│       │      │    ╱          ╲      │         0.5    2.0    │
│       │   ╭──┼──╮│            │     │                       │
│       │  │   ●  ││            │     │ Add Objects           │
│       │   ╰──┼──╯│            │     │ [+ Circle] [+ Line]   │
│        ╲    ╱     ╲          ╱      │ [+ Point]             │
│         ╰──╯       ╰────────╯       │                       │
│      inversion                      │ Transformation Rules: │
│      circle                         │ ○→○  Circle→Circle    │
│      (dashed)                       │ ─→○  Line→Circle      │
│                                     │ ○→─  Circle(thru)→Line│
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 3: Mobius Transformations

**File**: `crates/clifford-viz/examples/conformal2_mobius.rs`

### The Key Insight

All Mobius transformations in 2D can be composed from:
- Translations
- Rotations
- Dilations (uniform scaling)
- Inversions

CGA motors represent these transformations elegantly.

### Features

1. **Transformation Primitives**
   - Translation: Click and drag vector
   - Rotation: Angle slider around origin
   - Dilation: Scale factor slider
   - Inversion: Toggle and set circle

2. **Composition Chain**
   - Build up transformation from primitives
   - Reorder transformations (drag and drop)
   - See intermediate results

3. **Before/After Visualization**
   - Original shape (semi-transparent)
   - Transformed shape (solid)
   - Animation between states

4. **Fixed Points**
   - Highlight fixed points of transformation
   - Explain why they're fixed

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  Mobius Transformations                                [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Transformation Chain  │
│                                     │ ┌─────────────────┐   │
│      before       after             │ │ 1. Translate    │   │
│      (ghost)      (solid)           │ │    dx=1, dy=0.5 │   │
│                                     │ ├─────────────────┤   │
│        ╭──╮    ╭─────╮              │ │ 2. Rotate       │   │
│       │    │  │       │             │ │    θ = 45°      │   │
│       │    │  │       │             │ ├─────────────────┤   │
│        ╰──╯    ╰─────╯              │ │ 3. Dilate       │   │
│                                     │ │    k = 1.5      │   │
│                                     │ └─────────────────┘   │
│                                     │ [+ Translation]       │
│                                     │ [+ Rotation]          │
│                                     │ [+ Dilation]          │
│                                     │ [+ Inversion]         │
│                                     │                       │
│                                     │ [▶ Animate] [Reset]   │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 4: Circle-Circle Intersection

**File**: `crates/clifford-viz/examples/conformal2_intersection.rs`

### The Key Insight

In CGA, the intersection of two circles is computed via the meet (antiwedge):
```
PointPair = Circle₁ ∨ Circle₂
```

The result is a point pair (two points), a single point (tangent), or empty (no intersection).

### Features

1. **Two Circle Placement**
   - Drag circles to position
   - Resize with scroll/slider
   - Color coding

2. **Intersection Display**
   - Two points: show both
   - Tangent: single highlighted point
   - No intersection: indicate this state

3. **Intersection Algebra**
   - Show PointPair components
   - Explain how to extract points

4. **Special Cases**
   - Concentric circles (no intersection)
   - Tangent circles (one point)
   - Same circle (degenerate)

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  Circle-Circle Intersection                            [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Circle A (Blue)       │
│                                     │ Center: (1.0, 1.0)    │
│                                     │ Radius: 2.0           │
│         ╭───────╮                   │                       │
│        ╱    A    ╲   ╭────╮         │ Circle B (Red)        │
│       │          ●━━━━━● B │        │ Center: (3.5, 1.0)    │
│       │           ╲ ╱│     │        │ Radius: 1.5           │
│        ╲          ●━━━━━━━╱         │                       │
│         ╰───────╯   ╰────╯          │ Intersection          │
│                                     │ Points: 2             │
│    ● = intersection points          │ P₁: (2.8, 2.1)        │
│                                     │ P₂: (2.8, -0.1)       │
│                                     │                       │
│                                     │ Meet: A ∨ B           │
│                                     │ [Show algebra]        │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 5: Apollonian Gasket (2D)

**File**: `crates/clifford-viz/examples/conformal2_apollonian.rs`

### The Key Insight

The Apollonian gasket demonstrates CGA's power: given three mutually tangent circles, we can find the two circles tangent to all three using CGA operations.

In 2D, this is more visually accessible than the 3D version.

### Features

1. **Initial Configuration**
   - Three mutually tangent circles
   - Optional: outer bounding circle

2. **Recursive Generation**
   - For each "gap" (three mutually tangent circles), compute the inscribed circle
   - Recursively fill smaller gaps
   - Control recursion depth

3. **Interactive Controls**
   - Modify initial configuration
   - Adjust recursion depth (slider)
   - Color by generation/curvature

4. **Mathematical Display**
   - Show Descartes Circle Theorem: (k₁+k₂+k₃+k₄)² = 2(k₁²+k₂²+k₃²+k₄²)
   - Display curvatures (1/radius)

### Visualization

```
┌─────────────────────────────────────────────────────────────┐
│  Apollonian Gasket                                     [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Configuration         │
│         ╭─────────────────────╮     │ ○ Three circles       │
│        ╱ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○  ╲    │ ● With outer bound    │
│       │ ○   ╭───╮   ╭───╮   ○ │    │                       │
│      │ ○  ╱       ╲╱       ╲  ○│    │ Recursion Depth       │
│      │○  │    ○    │    ○   │ ○│    │ [●═══════════]  6     │
│      │ ○ │   ╱╲   │   ╱╲   │ ○│    │                       │
│      │  ○ ╲ │  │ ╱ ╲ │  │ ╱ ○ │    │ Color By              │
│       │ ○  ╲╰──╯╱   ╲╰──╯╱  ○│     │ ○ Generation          │
│        ╲ ○ ○ ○ ○ ○ ○ ○ ○ ○  ╱      │ ● Curvature           │
│         ╰─────────────────────╯     │                       │
│                                     │ Circles: 127          │
│                                     │ [Reset] [Randomize]   │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Implementation Tasks

### conformal2_circles.rs
1. [ ] Point placement and dragging
2. [ ] Circle-from-three-points computation
3. [ ] Circle rendering with fill/stroke options
4. [ ] Center and radius extraction
5. [ ] Collinearity detection and line display
6. [ ] Multi-circle support with color coding

### conformal2_inversion.rs
1. [ ] Inversion circle controls
2. [ ] Object library (circles, lines, points)
3. [ ] Inversion computation for each object type
4. [ ] Inverted object rendering
5. [ ] Special case handling (through center)
6. [ ] Animation support

### conformal2_mobius.rs
1. [ ] Translation motor
2. [ ] Rotation motor
3. [ ] Dilation motor
4. [ ] Inversion transformation
5. [ ] Transformation chain UI
6. [ ] Composition and reordering
7. [ ] Animation between states

### conformal2_intersection.rs
1. [ ] Two-circle placement
2. [ ] Meet operation for intersection
3. [ ] Point pair extraction
4. [ ] Tangent detection
5. [ ] No-intersection state handling
6. [ ] Algebraic display

### conformal2_apollonian.rs
1. [ ] Initial three-circle configuration
2. [ ] Tangent circle computation
3. [ ] Recursive generation algorithm
4. [ ] Depth control
5. [ ] Color schemes (generation, curvature)
6. [ ] Performance optimization for deep recursion

## Candidates for API Extraction

Methods discovered during visualization development that may be useful in the main algebra:

### RoundPoint Extensions
- `from_euclidean(x, y)` - Embed 2D point into CGA
- `to_euclidean() -> Option<(T, T)>` - Extract 2D coordinates
- `is_null()` - Check if on null cone (actual point)

### Circle Extensions
- `from_center_radius(cx, cy, r)` - Create from center and radius
- `from_three_points(p1, p2, p3)` - Create from three points
- `center() -> Option<(T, T)>` - Extract center
- `radius() -> Option<T>` - Extract radius
- `curvature() -> T` - Return 1/radius
- `is_line() -> bool` - Check if degenerate (infinite radius)

### Line Extensions
- `from_points(p1, p2)` - Create line through two points
- `from_normal_distance(nx, ny, d)` - Create from normal form

### PointPair Extensions
- `to_points() -> Option<((T, T), (T, T))>` - Extract both points
- `is_real() -> bool` - Check if intersection exists
- `is_tangent() -> bool` - Check if single point (tangent)

### Motor Extensions
- `translation(dx, dy)` - Translation motor
- `rotation(angle)` - Rotation motor about origin
- `dilation(factor)` - Uniform scaling motor
- `inversion(circle)` - Inversion through circle

## Verification

```bash
cargo run -p clifford-viz --example conformal2_circles --release
cargo run -p clifford-viz --example conformal2_inversion --release
cargo run -p clifford-viz --example conformal2_mobius --release
cargo run -p clifford-viz --example conformal2_intersection --release
cargo run -p clifford-viz --example conformal2_apollonian --release
```

### Checklist
- [ ] Circle from three points passes through all three
- [ ] Inversion produces correct transformed objects
- [ ] Lines through inversion center remain lines
- [ ] Mobius transformations compose correctly
- [ ] Circle-circle intersection finds correct points
- [ ] Tangent circles produce single intersection point
- [ ] Apollonian gasket circles are all mutually tangent
- [ ] All transformations preserve angles (conformal property)
- [ ] Demo runs at 60fps on mid-range hardware

## Educational Value

Each demo includes explanatory panels:
1. **Mathematical foundation** - What operation is being performed
2. **CGA representation** - How objects are encoded
3. **Algebraic computation** - What products are used
4. **Geometric interpretation** - What it means visually

This makes CGA2D an ideal "gateway" to understanding conformal geometric algebra before moving to the more complex 3D case.
