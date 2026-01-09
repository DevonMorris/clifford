# PRD-5: Projective GA (PGA)

**Status**: Pending
**Goal**: PGA for points, lines, planes, and rigid transforms

## Background

Projective Geometric Algebra adds a degenerate (null) basis vector to represent
points at infinity, enabling elegant representation of:
- Points, lines, planes as blades
- Rigid body transforms (rotation + translation) as motors
- Meet (intersection) and join (span) operations

## Deliverables

### 1. PGA Signatures (`src/signature/pga.rs`)

```rust
/// 2D Projective GA: Cl(2,0,1)
/// Basis: e1, e2, e0 where e0² = 0
#[derive(Clone, Copy, Debug, Default)]
pub struct PGA2D;

impl Signature for PGA2D {
    const P: usize = 2;
    const Q: usize = 0;
    const R: usize = 1;
    const DIM: usize = 3;
    const NUM_BLADES: usize = 8;

    fn metric(i: usize) -> i8 {
        if i < 2 { 1 } else { 0 }
    }
}

/// 3D Projective GA: Cl(3,0,1)
/// Basis: e1, e2, e3, e0 where e0² = 0
#[derive(Clone, Copy, Debug, Default)]
pub struct PGA3D;

impl Signature for PGA3D {
    const P: usize = 3;
    const Q: usize = 0;
    const R: usize = 1;
    const DIM: usize = 4;
    const NUM_BLADES: usize = 16;

    fn metric(i: usize) -> i8 {
        if i < 3 { 1 } else { 0 }
    }
}
```

### 2. PGA3D Types (`src/specialized/peuclidean::dim3/`)

```rust
/// A point in 3D PGA (grade 3 trivector)
/// P = x*e032 + y*e013 + z*e021 + w*e123
#[derive(Clone, Copy, Debug)]
pub struct Point<T: Float> {
    pub x: T,
    pub y: T,
    pub z: T,
    pub w: T,  // homogeneous coordinate
}

/// A plane in 3D PGA (grade 1 vector)
/// π = d*e0 + a*e1 + b*e2 + c*e3
/// Represents plane ax + by + cz + d = 0
#[derive(Clone, Copy, Debug)]
pub struct Plane<T: Float> {
    pub d: T,  // e0 coefficient
    pub a: T,  // e1 coefficient (normal x)
    pub b: T,  // e2 coefficient (normal y)
    pub c: T,  // e3 coefficient (normal z)
}

/// A line in 3D PGA (grade 2 bivector)
/// 6 components: e01, e02, e03, e23, e31, e12
#[derive(Clone, Copy, Debug)]
pub struct Line<T: Float> {
    // Direction (moment)
    pub dx: T,  // e01
    pub dy: T,  // e02
    pub dz: T,  // e03
    // Position (direction)
    pub mx: T,  // e23
    pub my: T,  // e31
    pub mz: T,  // e12
}

/// Motor: rigid body transformation
/// Even subalgebra element: scalar + bivector + quadvector
#[derive(Clone, Copy, Debug)]
pub struct Motor<T: Float> {
    // Scalar
    pub s: T,
    // Bivector (6 components)
    pub b01: T, pub b02: T, pub b03: T,
    pub b23: T, pub b31: T, pub b12: T,
    // Pseudoscalar
    pub q: T,
}
```

### 3. Motor Operations (`src/transforms/motor.rs`)

```rust
impl<T: Float> Motor<T> {
    /// Identity motor (no transformation)
    pub fn identity() -> Self;

    /// Pure rotation around origin
    pub fn from_rotor(angle: T, axis: [T; 3]) -> Self;

    /// Pure translation
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self;

    /// Rotation around arbitrary line
    pub fn from_line_angle(line: &Line<T>, angle: T) -> Self;

    /// Apply to point: M * P * M̃
    pub fn transform_point(&self, p: &Point<T>) -> Point<T>;

    /// Apply to line: M * L * M̃
    pub fn transform_line(&self, l: &Line<T>) -> Line<T>;

    /// Apply to plane: M * π * M̃
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T>;

    /// Compose: M2 * M1
    pub fn compose(&self, other: &Self) -> Self;

    /// Inverse transformation
    pub fn inverse(&self) -> Self;

    /// Interpolate motors (screw interpolation)
    pub fn lerp(&self, other: &Self, t: T) -> Self;
}
```

### 4. Meet and Join

```rust
impl<T: Float> Point<T> {
    /// Join two points -> line through them
    pub fn join(&self, other: &Point<T>) -> Line<T>;
}

impl<T: Float> Line<T> {
    /// Join line and point -> plane through them
    pub fn join_point(&self, p: &Point<T>) -> Plane<T>;

    /// Meet two lines -> point (if they intersect)
    pub fn meet(&self, other: &Line<T>) -> Option<Point<T>>;
}

impl<T: Float> Plane<T> {
    /// Meet two planes -> line of intersection
    pub fn meet(&self, other: &Plane<T>) -> Line<T>;

    /// Meet plane and line -> point
    pub fn meet_line(&self, l: &Line<T>) -> Point<T>;
}
```

## Files to Create

- `src/signature/pga.rs`
- `src/specialized/peuclidean::dim3/mod.rs`
- `src/specialized/peuclidean::dim3/types.rs`
- `src/specialized/peuclidean::dim3/ops.rs`
- `src/transforms/motor.rs`

## Testing (proptest)

```rust
proptest! {
    #[test]
    fn motor_preserves_distance(
        m in arb_motor::<f64>(),
        p1 in arb_point::<f64>(),
        p2 in arb_point::<f64>(),
    ) {
        let d_before = p1.distance(&p2);
        let d_after = m.transform_point(&p1).distance(&m.transform_point(&p2));
        prop_assert!((d_before - d_after).abs() < 1e-10);
    }

    #[test]
    fn motor_composition_associative(
        m1 in arb_motor::<f64>(),
        m2 in arb_motor::<f64>(),
        m3 in arb_motor::<f64>(),
    ) {
        let lhs = m1.compose(&m2).compose(&m3);
        let rhs = m1.compose(&m2.compose(&m3));
        prop_assert!(lhs.approx_eq(&rhs, 1e-10));
    }

    #[test]
    fn join_meet_duality(
        p1 in arb_point::<f64>(),
        p2 in arb_point::<f64>(),
        plane in arb_plane::<f64>(),
    ) {
        let line = p1.join(&p2);
        let intersection = plane.meet_line(&line);
        // intersection lies on the plane
        prop_assert!(plane.contains(&intersection, 1e-10));
    }
}
```

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - all proptest properties pass
- [ ] `cargo clippy` - no warnings
- [ ] Comprehensive geometric documentation
