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

## nalgebra Interoperability

When the `nalgebra-0_33` or `nalgebra-0_34` feature is enabled, PGA types should provide
conversions to/from nalgebra equivalents:

### Point Conversions

```rust
// PGA Point <-> nalgebra Point3 (homogeneous <-> Cartesian)
impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T> {
    fn from(p: na::Point3<T>) -> Self {
        Point { x: p.x, y: p.y, z: p.z, w: T::ONE }
    }
}

impl<T: Float + na::Scalar> TryFrom<Point<T>> for na::Point3<T> {
    type Error = NalgebraConversionError;

    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.w.abs() < T::EPSILON {
            return Err(NalgebraConversionError::PointAtInfinity);
        }
        Ok(na::Point3::new(p.x / p.w, p.y / p.w, p.z / p.w))
    }
}
```

### Plane Conversions

```rust
// Plane normal + distance representation
impl<T: Float + na::Scalar> From<Plane<T>> for (na::Unit<na::Vector3<T>>, T) {
    fn from(p: Plane<T>) -> Self {
        let normal = na::Unit::new_normalize(na::Vector3::new(p.a, p.b, p.c));
        let n = (p.a * p.a + p.b * p.b + p.c * p.c).sqrt();
        (normal, p.d / n)
    }
}
```

### Motor <-> Isometry3

```rust
// Motor represents rigid transform, equivalent to nalgebra's Isometry3
impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry3<T> {
    fn from(m: Motor<T>) -> Self {
        let rotation: na::UnitQuaternion<T> = m.rotation_part().into();
        let translation = m.translation_part();
        na::Isometry3::from_parts(
            na::Translation3::new(translation[0], translation[1], translation[2]),
            rotation,
        )
    }
}

impl<T: Float + na::RealField> From<na::Isometry3<T>> for Motor<T> {
    fn from(iso: na::Isometry3<T>) -> Self {
        let rotor = Rotor::from(iso.rotation);
        let t = iso.translation.vector;
        Motor::from_rotation_translation(rotor, [t.x, t.y, t.z])
    }
}
```

### Error Type Extension

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NalgebraConversionError {
    PointAtInfinity,  // PGA point has w ≈ 0
    InvalidLine,      // Line is ideal (at infinity)
}
```

## Benchmarks

PGA operations should be benchmarked against both:
1. Generic `Multivector` implementation
2. nalgebra equivalents (when feature enabled)

### Benchmark File (`benches/pga.rs`)

```rust
//! Benchmarks for PGA operations.
//!
//! Run with: `cargo bench --bench pga`

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use clifford::specialized::pga3d::{Point, Line, Plane, Motor};

// Motor operations
fn bench_motor_transform_point(c: &mut Criterion) {
    let motor = Motor::from_translation(1.0, 2.0, 3.0);
    let point = Point::new(1.0, 0.0, 0.0);

    c.bench_function("pga3d/motor_transform_point", |bencher| {
        bencher.iter(|| black_box(motor).transform_point(&black_box(point)))
    });
}

fn bench_motor_compose(c: &mut Criterion) {
    let m1 = Motor::from_translation(1.0, 0.0, 0.0);
    let m2 = Motor::from_rotor(0.5, [0.0, 0.0, 1.0]);

    c.bench_function("pga3d/motor_compose", |bencher| {
        bencher.iter(|| black_box(m1).compose(&black_box(m2)))
    });
}

// Meet and join operations
fn bench_point_join(c: &mut Criterion) {
    let p1 = Point::new(0.0, 0.0, 0.0);
    let p2 = Point::new(1.0, 1.0, 1.0);

    c.bench_function("pga3d/point_join_line", |bencher| {
        bencher.iter(|| black_box(p1).join(&black_box(p2)))
    });
}

fn bench_plane_meet(c: &mut Criterion) {
    let plane1 = Plane::new(0.0, 1.0, 0.0, 0.0); // y = 0
    let plane2 = Plane::new(0.0, 0.0, 1.0, 0.0); // z = 0

    c.bench_function("pga3d/plane_meet_line", |bencher| {
        bencher.iter(|| black_box(plane1).meet(&black_box(plane2)))
    });
}

fn bench_plane_meet_line(c: &mut Criterion) {
    let plane = Plane::new(0.0, 0.0, 1.0, -1.0); // z = 1
    let line = Line::through_points(
        &Point::new(0.0, 0.0, 0.0),
        &Point::new(0.0, 0.0, 2.0),
    );

    c.bench_function("pga3d/plane_meet_line_point", |bencher| {
        bencher.iter(|| black_box(plane).meet_line(&black_box(line)))
    });
}

criterion_group!(
    benches,
    bench_motor_transform_point,
    bench_motor_compose,
    bench_point_join,
    bench_plane_meet,
    bench_plane_meet_line,
);
criterion_main!(benches);
```

### nalgebra Comparison (`benches/pga_nalgebra.rs`)

When `nalgebra-0_33` or `nalgebra-0_34` is enabled, compare PGA motors with `na::Isometry3`:

```rust
#![cfg(any(feature = "nalgebra-0_33", feature = "nalgebra-0_34"))]

// Motor vs Isometry3 for rigid transforms
fn bench_motor_transform_clifford(c: &mut Criterion) {
    let motor = Motor::from_translation(1.0, 2.0, 3.0)
        .compose(&Motor::from_rotor(0.5, [0.0, 0.0, 1.0]));
    let point = Point::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rigid_transform/pga_motor", |bencher| {
        bencher.iter(|| black_box(motor).transform_point(&black_box(point)))
    });
}

fn bench_isometry_transform_nalgebra(c: &mut Criterion) {
    let iso = na::Isometry3::new(
        na::Vector3::new(1.0, 2.0, 3.0),
        na::Vector3::z() * 0.5,
    );
    let point = na::Point3::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rigid_transform/nalgebra_isometry", |bencher| {
        bencher.iter(|| black_box(iso) * black_box(point))
    });
}

// Motor composition vs Isometry composition
fn bench_motor_chain(c: &mut Criterion) {
    let m1 = Motor::from_translation(1.0, 0.0, 0.0);
    let m2 = Motor::from_rotor(0.3, [1.0, 0.0, 0.0]);
    let m3 = Motor::from_translation(0.0, 1.0, 0.0);

    c.bench_function("comparison/compose_3/pga_motor", |bencher| {
        bencher.iter(|| {
            black_box(m1).compose(&black_box(m2)).compose(&black_box(m3))
        })
    });
}

fn bench_isometry_chain(c: &mut Criterion) {
    let i1 = na::Isometry3::translation(1.0, 0.0, 0.0);
    let i2 = na::Isometry3::rotation(na::Vector3::x() * 0.3);
    let i3 = na::Isometry3::translation(0.0, 1.0, 0.0);

    c.bench_function("comparison/compose_3/nalgebra_isometry", |bencher| {
        bencher.iter(|| {
            black_box(i1) * black_box(i2) * black_box(i3)
        })
    });
}
```

### Expected Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| Motor transform point | < 10 ns | Sandwich product optimization |
| Motor compose | < 15 ns | Even subalgebra product |
| Point join (line) | < 5 ns | Regressive product |
| Plane meet (line) | < 5 ns | Outer product |
| Motor vs Isometry3 | Comparable | Similar underlying math |

### Benchmark README Section

Add to `benches/README.md`:

```markdown
## PGA Benchmarks

Operations on specialized PGA types.

### Motor Operations
![pga3d_motor_transform_point](reports/pga3d_motor_transform_point_pdf.svg)
![pga3d_motor_compose](reports/pga3d_motor_compose_pdf.svg)

### Meet/Join Operations
![pga3d_point_join_line](reports/pga3d_point_join_line_pdf.svg)
![pga3d_plane_meet_line](reports/pga3d_plane_meet_line_pdf.svg)

### Comparison with nalgebra (with nalgebra-0_34 feature)
| Operation | PGA Motor | nalgebra Isometry3 |
|-----------|-----------|-------------------|
| Transform point | ~X ns | ~Y ns |
| Compose 3 | ~X ns | ~Y ns |
```

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - all proptest properties pass
- [ ] `cargo clippy` - no warnings
- [ ] Comprehensive geometric documentation
- [ ] nalgebra conversions tested with feature flags
- [ ] Benchmarks run and documented
