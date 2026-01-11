# PRD-6.7: Operations and Integration

**Status**: Pending
**Parent**: PRD-6 (Conformal Geometric Algebra)
**Depends On**: PRD-6.6 (Dilator and Inversor Versors)
**Goal**: Implement meet/join operations, carriers, nalgebra integration, and benchmarks

## Reference

- https://conformalgeometricalgebra.org/wiki/index.php?title=Expansion
- https://conformalgeometricalgebra.org/wiki/index.php?title=Projections

## Background

### Meet and Join

In CGA, **meet** (∧) and **join** (∨) operations compute geometric intersections and spans:

| Operation | Result |
|-----------|--------|
| Sphere ∧ Sphere | Circle (intersection) |
| Sphere ∧ Plane | Circle |
| Plane ∧ Plane | Line |
| Plane ∧ Line | Point |
| Point ∨ Point | Line (through both) |
| Point ∨ Line | Plane (containing both) |
| Point ∨ Point ∨ Point | Circle (through all three) |

### Carriers and Attitudes

- **Carrier**: The flat element that carries a round element (e.g., the plane containing a circle)
- **Attitude**: The direction/orientation component of an element

### Projections

- **Projection onto flat**: Project a round element onto a plane/line
- **Orthogonal complement**: The perpendicular component

## Deliverables

### 1. Meet Operations (`dim3/meet.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere, Plane, Circle, Line, Dipole};

/// Meet (intersection) of two spheres: a circle.
pub fn sphere_meet_sphere<T: Float>(s1: &Sphere<T>, s2: &Sphere<T>) -> Circle<T> {
    // Circle = S₁ ∧ S₂
    todo!("Derive from SymPy")
}

/// Meet of sphere and plane: a circle.
pub fn sphere_meet_plane<T: Float>(s: &Sphere<T>, p: &Plane<T>) -> Circle<T> {
    todo!("Derive from SymPy")
}

/// Meet of two planes: a line.
pub fn plane_meet_plane<T: Float>(p1: &Plane<T>, p2: &Plane<T>) -> Line<T> {
    todo!("Derive from SymPy")
}

/// Meet of plane and line: a point.
pub fn plane_meet_line<T: Float>(p: &Plane<T>, l: &Line<T>) -> Point<T> {
    todo!("Derive from SymPy")
}

/// Meet of three planes: a point.
pub fn planes_meet<T: Float>(p1: &Plane<T>, p2: &Plane<T>, p3: &Plane<T>) -> Point<T> {
    plane_meet_line(&p1, &plane_meet_plane(&p2, &p3))
}

/// Meet of sphere and line: a dipole (point pair).
pub fn sphere_meet_line<T: Float>(s: &Sphere<T>, l: &Line<T>) -> Dipole<T> {
    todo!("Derive from SymPy")
}

/// Add meet method to Sphere
impl<T: Float> Sphere<T> {
    /// Meet with another sphere: their intersection circle.
    pub fn meet_sphere(&self, other: &Sphere<T>) -> Circle<T> {
        sphere_meet_sphere(self, other)
    }

    /// Meet with a plane: intersection circle.
    pub fn meet_plane(&self, plane: &Plane<T>) -> Circle<T> {
        sphere_meet_plane(self, plane)
    }

    /// Meet with a line: intersection points (dipole).
    pub fn meet_line(&self, line: &Line<T>) -> Dipole<T> {
        sphere_meet_line(self, line)
    }
}
```

### 2. Join Operations (`dim3/join.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Circle, Line, Plane, Sphere};

/// Join of two points: line through them.
pub fn point_join_point<T: Float>(p1: &Point<T>, p2: &Point<T>) -> Line<T> {
    // L = P₁ ∧ P₂ ∧ e∞
    todo!("Derive from SymPy")
}

/// Join of three points: circle through them.
pub fn points_join_circle<T: Float>(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> Circle<T> {
    // C = P₁ ∧ P₂ ∧ P₃
    Circle::from_three_points(p1, p2, p3)
}

/// Join of four points: sphere through them.
pub fn points_join_sphere<T: Float>(
    p1: &Point<T>, p2: &Point<T>, p3: &Point<T>, p4: &Point<T>
) -> Sphere<T> {
    // S = P₁ ∧ P₂ ∧ P₃ ∧ P₄
    Sphere::from_four_points(p1, p2, p3, p4)
}

/// Join of point and line: plane containing both.
pub fn point_join_line<T: Float>(p: &Point<T>, l: &Line<T>) -> Plane<T> {
    todo!("Derive from SymPy")
}

/// Join of two lines: plane containing both (if coplanar).
pub fn line_join_line<T: Float>(l1: &Line<T>, l2: &Line<T>) -> Option<Plane<T>> {
    // Only works if lines are coplanar
    todo!("Derive from SymPy")
}

/// Add join method to Point
impl<T: Float> Point<T> {
    /// Join with another point: line through both.
    pub fn join_line(&self, other: &Point<T>) -> Line<T> {
        point_join_point(self, other)
    }

    /// Join with two other points: circle through all three.
    pub fn join_circle(&self, p2: &Point<T>, p3: &Point<T>) -> Circle<T> {
        points_join_circle(self, p2, p3)
    }

    /// Join with three other points: sphere through all four.
    pub fn join_sphere(&self, p2: &Point<T>, p3: &Point<T>, p4: &Point<T>) -> Sphere<T> {
        points_join_sphere(self, p2, p3, p4)
    }

    /// Join with a line: plane containing both.
    pub fn join_plane(&self, line: &Line<T>) -> Plane<T> {
        point_join_line(self, line)
    }
}
```

### 3. Carrier and Attitude (`dim3/carrier.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Circle, Line, Plane, Sphere};

/// Extracts the carrier (flat element) of a round element.
pub trait Carrier<T: Float> {
    type CarrierType;

    /// Returns the flat element that carries this round element.
    fn carrier(&self) -> Self::CarrierType;
}

impl<T: Float> Carrier<T> for Circle<T> {
    type CarrierType = Plane<T>;

    /// Returns the plane containing this circle.
    fn carrier(&self) -> Plane<T> {
        self.carrier_plane()
    }
}

impl<T: Float> Carrier<T> for Sphere<T> {
    type CarrierType = Point<T>;

    /// Returns the center point (carrier of a sphere).
    fn carrier(&self) -> Point<T> {
        let (cx, cy, cz) = self.center();
        Point::new(cx, cy, cz)
    }
}

/// Extracts the attitude (directional component) of an element.
pub trait Attitude<T: Float> {
    type AttitudeType;

    /// Returns the attitude of this element.
    fn attitude(&self) -> Self::AttitudeType;
}

impl<T: Float> Attitude<T> for Plane<T> {
    type AttitudeType = (T, T, T);  // Normal vector

    fn attitude(&self) -> (T, T, T) {
        self.normal()
    }
}

impl<T: Float> Attitude<T> for Line<T> {
    type AttitudeType = (T, T, T);  // Direction vector

    fn attitude(&self) -> (T, T, T) {
        self.direction()
    }
}
```

### 4. Projections (`dim3/projection.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Line, Plane, Circle, Sphere};

/// Projects a point onto a plane.
pub fn project_point_onto_plane<T: Float>(p: &Point<T>, plane: &Plane<T>) -> Point<T> {
    plane.project(p)
}

/// Projects a point onto a line.
pub fn project_point_onto_line<T: Float>(p: &Point<T>, line: &Line<T>) -> Point<T> {
    line.closest_point(p)
}

/// Projects a point onto a sphere (closest point on surface).
pub fn project_point_onto_sphere<T: Float>(p: &Point<T>, sphere: &Sphere<T>) -> Point<T> {
    let (cx, cy, cz) = sphere.center();
    let r = sphere.radius();

    // Direction from center to point
    let dx = p.x() - cx;
    let dy = p.y() - cy;
    let dz = p.z() - cz;
    let dist = (dx*dx + dy*dy + dz*dz).sqrt();

    if dist < T::epsilon() {
        // Point at center - return arbitrary point on surface
        Point::new(cx + r, cy, cz)
    } else {
        // Scale to sphere surface
        let scale = r / dist;
        Point::new(cx + dx * scale, cy + dy * scale, cz + dz * scale)
    }
}

/// Rejection: the component of a point perpendicular to a plane.
pub fn reject_point_from_plane<T: Float>(p: &Point<T>, plane: &Plane<T>) -> Point<T> {
    let projected = plane.project(p);
    Point::new(
        p.x() - projected.x(),
        p.y() - projected.y(),
        p.z() - projected.z(),
    )
}
```

### 5. nalgebra Integration (`dim3/nalgebra.rs`)

Comprehensive nalgebra conversions:

```rust
#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra_impl {
    use super::*;

    cfg_if::cfg_if! { /* version selection */ }

    // ================================================================
    // Point conversions
    // ================================================================

    impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T> {
        fn from(p: na::Point3<T>) -> Self {
            Point::new(p.x, p.y, p.z)
        }
    }

    impl<T: Float + na::Scalar> From<Point<T>> for na::Point3<T> {
        fn from(p: Point<T>) -> Self {
            na::Point3::new(p.x(), p.y(), p.z())
        }
    }

    // ================================================================
    // Sphere conversions
    // ================================================================

    /// Convert sphere to nalgebra representation (center, radius).
    impl<T: Float + na::Scalar> From<Sphere<T>> for (na::Point3<T>, T) {
        fn from(s: Sphere<T>) -> Self {
            let (cx, cy, cz) = s.center();
            (na::Point3::new(cx, cy, cz), s.radius())
        }
    }

    impl<T: Float + na::Scalar> From<(na::Point3<T>, T)> for Sphere<T> {
        fn from((center, radius): (na::Point3<T>, T)) -> Self {
            Sphere::from_center_radius(center.x, center.y, center.z, radius)
        }
    }

    // ================================================================
    // Plane conversions
    // ================================================================

    /// Plane to (normal, distance) representation.
    impl<T: Float + na::Scalar> From<Plane<T>> for (na::Unit<na::Vector3<T>>, T) {
        fn from(p: Plane<T>) -> Self {
            let (a, b, c) = p.normal();
            let d = p.distance_from_origin();
            (na::Unit::new_normalize(na::Vector3::new(a, b, c)), d)
        }
    }

    // ================================================================
    // Translator conversions
    // ================================================================

    impl<T: Float + na::RealField> From<na::Translation3<T>> for Translator<T> {
        fn from(t: na::Translation3<T>) -> Self {
            Translator::new(t.x, t.y, t.z)
        }
    }

    impl<T: Float + na::RealField> From<Translator<T>> for na::Translation3<T> {
        fn from(t: Translator<T>) -> Self {
            let (dx, dy, dz) = t.displacement();
            na::Translation3::new(dx, dy, dz)
        }
    }

    // ================================================================
    // Rotor conversions
    // ================================================================

    impl<T: Float + na::RealField> From<na::UnitQuaternion<T>> for Rotor<T> {
        fn from(q: na::UnitQuaternion<T>) -> Self {
            // nalgebra quaternion: (x, y, z, w) where w is scalar
            let (axis, angle) = q.axis_angle()
                .map(|(a, ang)| ((a.x, a.y, a.z), ang))
                .unwrap_or(((T::one(), T::zero(), T::zero()), T::zero()));
            Rotor::from_axis_angle(axis, angle)
        }
    }

    impl<T: Float + na::RealField> From<Rotor<T>> for na::UnitQuaternion<T> {
        fn from(r: Rotor<T>) -> Self {
            if let Some(axis) = r.axis() {
                let na_axis = na::Unit::new_normalize(
                    na::Vector3::new(axis.0, axis.1, axis.2)
                );
                na::UnitQuaternion::from_axis_angle(&na_axis, r.angle())
            } else {
                na::UnitQuaternion::identity()
            }
        }
    }

    // ================================================================
    // Motor conversions
    // ================================================================

    impl<T: Float + na::RealField> From<na::Isometry3<T>> for Motor<T> {
        fn from(iso: na::Isometry3<T>) -> Self {
            let t = Translator::from(iso.translation);
            let r = Rotor::from(iso.rotation);
            Motor::from_translator_rotor(&t, &r)
        }
    }

    impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry3<T> {
        fn from(m: Motor<T>) -> Self {
            let (t, r) = m.decompose();
            let na_t: na::Translation3<T> = t.into();
            let na_r: na::UnitQuaternion<T> = r.into();
            na::Isometry3::from_parts(na_t, na_r)
        }
    }

    // ================================================================
    // Dilator conversions
    // ================================================================

    impl<T: Float + na::RealField> From<Dilator<T>> for na::Similarity3<T> {
        fn from(d: Dilator<T>) -> Self {
            na::Similarity3::from_scaling(d.scale_factor())
        }
    }
}
```

### 6. Benchmarks (`benches/cga.rs`)

```rust
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use clifford::specialized::conformal::dim3::*;

fn bench_point_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("cga_point");

    group.bench_function("embed", |b| {
        b.iter(|| Point::new(black_box(1.0), black_box(2.0), black_box(3.0)))
    });

    let p1 = Point::new(1.0, 2.0, 3.0);
    let p2 = Point::new(4.0, 5.0, 6.0);

    group.bench_function("distance", |b| {
        b.iter(|| p1.distance(black_box(&p2)))
    });

    group.finish();
}

fn bench_sphere_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("cga_sphere");

    let s1 = Sphere::from_center_radius(0.0, 0.0, 0.0, 1.0);
    let s2 = Sphere::from_center_radius(1.0, 0.0, 0.0, 1.0);
    let p = Point::new(1.0, 0.0, 0.0);

    group.bench_function("contains", |b| {
        b.iter(|| s1.contains(black_box(&p), 1e-10))
    });

    group.bench_function("meet_sphere", |b| {
        b.iter(|| s1.meet_sphere(black_box(&s2)))
    });

    group.finish();
}

fn bench_translator(c: &mut Criterion) {
    let mut group = c.benchmark_group("cga_translator");

    let t = Translator::new(1.0, 2.0, 3.0);
    let p = Point::new(0.0, 0.0, 0.0);

    group.bench_function("transform_point", |b| {
        b.iter(|| t.transform_point(black_box(&p)))
    });

    let t2 = Translator::new(4.0, 5.0, 6.0);
    group.bench_function("compose", |b| {
        b.iter(|| t.compose(black_box(&t2)))
    });

    group.finish();
}

fn bench_motor(c: &mut Criterion) {
    let mut group = c.benchmark_group("cga_motor");

    let r = Rotor::from_axis_angle((0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_4);
    let t = Translator::new(1.0, 2.0, 3.0);
    let m = Motor::from_translator_rotor(&t, &r);
    let p = Point::new(1.0, 0.0, 0.0);

    group.bench_function("transform_point", |b| {
        b.iter(|| m.transform_point(black_box(&p)))
    });

    let m2 = Motor::from_translator_rotor(&Translator::new(0.0, 1.0, 0.0), &r);
    group.bench_function("compose", |b| {
        b.iter(|| m.compose(black_box(&m2)))
    });

    group.finish();
}

fn bench_cga_vs_nalgebra(c: &mut Criterion) {
    let mut group = c.benchmark_group("cga_vs_nalgebra");

    // Translation comparison
    let cga_t = Translator::new(1.0, 2.0, 3.0);
    let cga_p = Point::new(0.0, 0.0, 0.0);
    let na_t = na::Translation3::new(1.0, 2.0, 3.0);
    let na_p = na::Point3::new(0.0, 0.0, 0.0);

    group.bench_function("translate_cga", |b| {
        b.iter(|| cga_t.transform_point(black_box(&cga_p)))
    });

    group.bench_function("translate_nalgebra", |b| {
        b.iter(|| na_t * black_box(na_p))
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_point_operations,
    bench_sphere_operations,
    bench_translator,
    bench_motor,
    bench_cga_vs_nalgebra,
);
criterion_main!(benches);
```

## Property Tests

```rust
proptest! {
    // ================================================================
    // Meet/Join duality
    // ================================================================

    #[test]
    fn sphere_meet_contains_common_points(
        c1x in -5.0f64..5.0, c1y in -5.0f64..5.0, c1z in -5.0f64..5.0,
        c2x in -5.0f64..5.0, c2y in -5.0f64..5.0, c2z in -5.0f64..5.0,
        r1 in 1.0f64..5.0,
        r2 in 1.0f64..5.0,
    ) {
        let s1 = Sphere::from_center_radius(c1x, c1y, c1z, r1);
        let s2 = Sphere::from_center_radius(c2x, c2y, c2z, r2);
        let circle = s1.meet_sphere(&s2);

        // Any point on the circle should be on both spheres
        // (test by checking center is equidistant from both sphere centers)
        let center = circle.center();
        let d1 = center.distance(&Point::new(c1x, c1y, c1z));
        let d2 = center.distance(&Point::new(c2x, c2y, c2z));

        // Circle center should be equidistant from both sphere centers
        // (only if spheres intersect)
    }

    #[test]
    fn point_join_is_line_through_points(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
    ) {
        let line = p1.join_line(&p2);
        prop_assert!(line.contains(&p1, ABS_DIFF_EQ_EPS));
        prop_assert!(line.contains(&p2, ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // nalgebra consistency
    // ================================================================

    #[test]
    fn nalgebra_point_roundtrip(
        x in -100.0f64..100.0,
        y in -100.0f64..100.0,
        z in -100.0f64..100.0,
    ) {
        let cga = Point::new(x, y, z);
        let na: na::Point3<f64> = cga.into();
        let back: Point<f64> = na.into();

        prop_assert!(abs_diff_eq!(cga.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(cga.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(cga.z(), back.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn nalgebra_motor_roundtrip(
        axis_x in -1.0f64..1.0, axis_y in -1.0f64..1.0, axis_z in -1.0f64..1.0,
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        dx in -10.0f64..10.0, dy in -10.0f64..10.0, dz in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let len = (axis_x*axis_x + axis_y*axis_y + axis_z*axis_z).sqrt();
        if len < ABS_DIFF_EQ_EPS {
            return Ok(());
        }
        let axis = (axis_x/len, axis_y/len, axis_z/len);

        let cga_motor = Motor::from_translator_rotor(
            &Translator::new(dx, dy, dz),
            &Rotor::from_axis_angle(axis, angle)
        );

        // Convert to nalgebra and back
        let na_iso: na::Isometry3<f64> = cga_motor.into();
        let back: Motor<f64> = na_iso.into();

        // Transform a point with both and compare
        let p = Point::new(px, py, pz);
        let q1 = cga_motor.transform_point(&p);
        let q2 = back.transform_point(&p);

        prop_assert!(abs_diff_eq!(q1.x(), q2.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q1.y(), q2.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q1.z(), q2.z(), epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Files to Create/Modify

### New Files
- `src/specialized/conformal/dim3/meet.rs`
- `src/specialized/conformal/dim3/join.rs`
- `src/specialized/conformal/dim3/carrier.rs`
- `src/specialized/conformal/dim3/projection.rs`
- `benches/cga.rs`
- `benches/cga_nalgebra_comparison.rs`

### Modified Files
- `src/specialized/conformal/dim3/mod.rs` - Export new modules
- `src/specialized/conformal/dim3/nalgebra.rs` - Complete conversions
- `Cargo.toml` - Add CGA benchmark
- `derivations/src/clifford_derivations/cga.py` - Final derivations

## Verification Checklist

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] `cargo bench` - CGA benchmarks run
- [ ] Meet operations produce correct intersections
- [ ] Join operations produce correct spans
- [ ] nalgebra conversions roundtrip correctly
- [ ] Performance comparable to nalgebra for equivalent operations

## Dependencies

- PRD-6.6 (Dilator and Inversor Versors) - must be complete
- All previous PRD-6.x phases

## Completion Criteria

When all sub-PRDs (6.1-6.7) are complete:

1. Update `docs/prd/prd-6-cga.md` status to **Complete**
2. Update `docs/prd/README.md` to show PRD-6 as complete
3. Add CGA section to `README.md`
4. Update `CHANGELOG.md` with CGA features
5. Consider creating examples:
   - `examples/cga_visualization.rs` (with rerun)
   - `examples/cga_tutorial.rs`

## Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| Point embed | < 10 ns | x² + y² + z² computation |
| Point extract | < 5 ns | Division by weight |
| Distance | < 15 ns | Inner product |
| Translator transform | < 20 ns | Simple addition |
| Rotor transform | < 30 ns | Quaternion-like |
| Motor transform | < 50 ns | Combined operation |
| Sphere meet | < 30 ns | Outer product |
| CGA vs nalgebra | Within 2x | CGA has embed/extract overhead |
