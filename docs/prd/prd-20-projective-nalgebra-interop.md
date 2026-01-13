# PRD-20: Projective GA nalgebra and Rerun Interoperability

## Overview

This PRD defines the nalgebra and Rerun interoperability for 2D and 3D Projective Geometric Algebra (PGA) types. The goal is seamless bidirectional conversion between clifford's PGA types and nalgebra's equivalent types, with comprehensive testing to verify operational equivalence.

## Problem Statement

Currently, the projective modules (`projective::dim2` and `projective::dim3`) have placeholder comments for nalgebra integration but no implementation. Users working with PGA need to:

1. Convert between clifford and nalgebra types for interop with existing code
2. Trust that equivalent operations produce equivalent results
3. Visualize PGA elements using Rerun (dim3 has partial support)

Without nalgebra interop, users cannot easily integrate clifford's PGA implementation into existing robotics, graphics, or physics codebases that use nalgebra.

## Goals

1. **Complete nalgebra conversions** for all PGA types with nalgebra equivalents
2. **Comprehensive property-based testing** verifying operational equivalence
3. **Complete Rerun support** for dim2 (dim3 already has partial support)
4. **Document mathematical correspondences** clearly in rustdoc

## Non-Goals

- Supporting nalgebra versions prior to 0.32
- Creating new visualization primitives in Rerun
- Optimizing for performance over correctness (correctness first)

## Type Mappings

### 3D PGA (projective::dim3)

| clifford | nalgebra | Conversion | Notes |
|----------|----------|------------|-------|
| `Point<T>` | `na::Point3<T>` | Bidirectional | Homogeneous ↔ Cartesian |
| `Motor<T>` | `na::Isometry3<T>` | Bidirectional | Rigid transformation |
| `Motor<T>` | `na::UnitDualQuaternion<T>` | Bidirectional | Direct representation |
| `Plane<T>` | `(na::Unit<na::Vector3<T>>, T)` | Extract only | Normal + distance |
| `Line<T>` | Plücker coords | Extract only | Direction + moment |
| `Flector<T>` | N/A | No direct equiv | Improper isometry |

### 2D PGA (projective::dim2)

| clifford | nalgebra | Conversion | Notes |
|----------|----------|------------|-------|
| `Point<T>` | `na::Point2<T>` | Bidirectional | Homogeneous ↔ Cartesian |
| `Motor<T>` | `na::Isometry2<T>` | Bidirectional | Rigid transformation |
| `Motor<T>` | `na::UnitComplex<T>` + translation | Bidirectional | Rotation + translation |
| `Line<T>` | `(na::Unit<na::Vector2<T>>, T)` | Extract only | Normal + distance |

### Rerun Mappings

**dim3** (existing, to be enhanced):
| clifford | rerun | Notes |
|----------|-------|-------|
| `Point<f32>` | `Position3D` | Already implemented |
| `Motor<f32>` | `Transform3D` | Already implemented |
| `Plane<f32>` | `Vec3D` | Already implemented (normal) |
| `Line<f32>` | `LineStrips3D` | New: two points on line |

**dim2** (new):
| clifford | rerun | Notes |
|----------|-------|-------|
| `Point<f32>` | `Position2D` | Cartesian extraction |
| `Motor<f32>` | `Transform3D` | As 2D embedded in 3D |
| `Line<f32>` | `LineStrips2D` | Two points on line |

## Design

### 1. Point Conversions

```rust
// dim3
impl<T: Float + na::Scalar> From<Point<T>> for na::Point3<T> {
    fn from(p: Point<T>) -> Self {
        // Extract Cartesian coordinates (divide by w)
        na::Point3::new(p.x(), p.y(), p.z())
    }
}

impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T> {
    fn from(p: na::Point3<T>) -> Self {
        // Create homogeneous point with w=1
        Point::from_cartesian(p.x, p.y, p.z)
    }
}

// dim2 - similar pattern
```

### 2. Motor <-> Isometry Conversions

The key insight is that PGA motors and nalgebra isometries both represent rigid transformations (rotation + translation), but with different internal representations.

**Motor representation** (3D):
- Scalar `s` and bivector `(e23, e31, e12)` encode rotation
- Null bivector `(e01, e02, e03)` encodes translation

**Isometry representation** (3D):
- `UnitQuaternion` for rotation
- `Translation3` for translation

```rust
impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry3<T> {
    fn from(motor: Motor<T>) -> Self {
        let m = motor.unitized();

        // Extract rotation quaternion
        // Motor: (s, e23, e31, e12) -> Quaternion: (w, x, y, z)
        let rotation = na::UnitQuaternion::from_quaternion(
            na::Quaternion::new(m.s(), m.e23(), m.e31(), m.e12())
        );

        // Extract translation
        // t = 2 * (s*d - b×d) where d is the null bivector part
        let tx = T::TWO * (m.s() * m.e01() + m.e12() * m.e02() - m.e31() * m.e03());
        let ty = T::TWO * (m.s() * m.e02() + m.e23() * m.e03() - m.e12() * m.e01());
        let tz = T::TWO * (m.s() * m.e03() + m.e31() * m.e01() - m.e23() * m.e02());
        let translation = na::Translation3::new(tx, ty, tz);

        na::Isometry3::from_parts(translation, rotation)
    }
}

impl<T: Float + na::RealField> From<na::Isometry3<T>> for Motor<T> {
    fn from(iso: na::Isometry3<T>) -> Self {
        let q = iso.rotation.quaternion();
        let t = iso.translation.vector;

        // Rotation part: quaternion -> motor bivector
        let s = q.w;
        let e23 = q.i;
        let e31 = q.j;
        let e12 = q.k;

        // Translation part: t -> null bivector (with 1/2 factor)
        // d = (t/2) encoded in the null bivector, but modified by rotation
        let half = T::from_f64(0.5);
        let e01 = half * (s * t.x - e12 * t.y + e31 * t.z);
        let e02 = half * (s * t.y + e12 * t.x - e23 * t.z);
        let e03 = half * (s * t.z - e31 * t.x + e23 * t.y);

        // e0123 = 0 for valid motor
        Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, T::zero())
    }
}
```

### 3. Comprehensive Testing Strategy

The testing strategy focuses on **operational equivalence**: if we perform the same operation using clifford and nalgebra, we should get the same result.

#### 3.1 Round-trip Tests

Verify that converting to nalgebra and back preserves the value:

```rust
proptest! {
    #[test]
    fn point_roundtrip(p in any::<Point<f64>>()) {
        // Skip ideal points (w ≈ 0)
        prop_assume!(p.e0().abs() > 1e-6);

        let na_p: na::Point3<f64> = p.into();
        let back: Point<f64> = na_p.into();

        // Compare Cartesian coordinates
        prop_assert!(relative_eq!(p.x(), back.x(), epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(p.y(), back.y(), epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(p.z(), back.z(), epsilon = EPS, max_relative = EPS));
    }

    #[test]
    fn motor_roundtrip(m in any::<BulkMotor<f64>>()) {
        let iso: na::Isometry3<f64> = (*m).into();
        let back: Motor<f64> = iso.into();

        // Motors have double cover, test by applying to a point
        let test_p = Point::from_cartesian(1.0, 2.0, 3.0);
        let result1 = m.transform_point(&test_p);
        let result2 = back.transform_point(&test_p);

        prop_assert!(relative_eq!(result1.x(), result2.x(), epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result1.y(), result2.y(), epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result1.z(), result2.z(), epsilon = EPS, max_relative = EPS));
    }
}
```

#### 3.2 Operational Equivalence Tests

**Critical**: Verify that the same operation produces the same result in both representations.

```rust
proptest! {
    /// Verify motor.transform_point() == isometry.transform_point()
    #[test]
    fn transform_point_equivalence(
        m in any::<BulkMotor<f64>>(),
        p in finite_point_strategy(),  // Only finite points
    ) {
        // Transform with clifford motor
        let result_ga = m.transform_point(&p);

        // Transform with nalgebra isometry
        let iso: na::Isometry3<f64> = (*m).into();
        let na_p: na::Point3<f64> = p.into();
        let na_result = iso.transform_point(&na_p);

        // Compare results
        prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
    }

    /// Verify motor composition == isometry composition
    #[test]
    fn composition_equivalence(
        m1 in any::<BulkMotor<f64>>(),
        m2 in any::<BulkMotor<f64>>(),
        p in finite_point_strategy(),
    ) {
        // Compose with clifford
        let composed_ga = m1.compose(&*m2);
        let result_ga = composed_ga.transform_point(&p);

        // Compose with nalgebra
        let iso1: na::Isometry3<f64> = (*m1).into();
        let iso2: na::Isometry3<f64> = (*m2).into();
        let composed_na = iso1 * iso2;
        let na_p: na::Point3<f64> = p.into();
        let na_result = composed_na.transform_point(&na_p);

        // Compare results
        prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
    }

    /// Verify motor inverse == isometry inverse
    #[test]
    fn inverse_equivalence(
        m in any::<BulkMotor<f64>>(),
        p in finite_point_strategy(),
    ) {
        // Inverse with clifford
        let inv_ga = m.inverse();
        let result_ga = inv_ga.transform_point(&p);

        // Inverse with nalgebra
        let iso: na::Isometry3<f64> = (*m).into();
        let inv_na = iso.inverse();
        let na_p: na::Point3<f64> = p.into();
        let na_result = inv_na.transform_point(&na_p);

        // Compare results
        prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
    }
}
```

#### 3.3 Factory Method Equivalence Tests

Verify that factory methods produce equivalent transforms:

```rust
proptest! {
    /// Verify from_translation produces equivalent results
    #[test]
    fn translation_factory_equivalence(
        tx in -100.0f64..100.0,
        ty in -100.0f64..100.0,
        tz in -100.0f64..100.0,
        p in finite_point_strategy(),
    ) {
        let motor = Motor::from_translation(tx, ty, tz);
        let iso = na::Isometry3::from_parts(
            na::Translation3::new(tx, ty, tz),
            na::UnitQuaternion::identity(),
        );

        let result_ga = motor.transform_point(&p);
        let na_p: na::Point3<f64> = p.into();
        let na_result = iso.transform_point(&na_p);

        prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
    }

    /// Verify from_rotation_z produces equivalent results
    #[test]
    fn rotation_z_factory_equivalence(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        p in finite_point_strategy(),
    ) {
        let motor = Motor::from_rotation_z(angle);
        let iso = na::Isometry3::from_parts(
            na::Translation3::identity(),
            na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), angle),
        );

        let result_ga = motor.transform_point(&p);
        let na_p: na::Point3<f64> = p.into();
        let na_result = iso.transform_point(&na_p);

        prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
    }
}
```

#### 3.4 2D-Specific Tests

Similar tests for 2D, plus verification that 2D operations match 3D operations when embedded:

```rust
proptest! {
    /// 2D motor transformation equivalence
    #[test]
    fn transform_point_equivalence_2d(
        m in any::<dim2::BulkMotor<f64>>(),
        p in finite_point_2d_strategy(),
    ) {
        let result_ga = m.transform_point(&p);

        let iso: na::Isometry2<f64> = (*m).into();
        let na_p: na::Point2<f64> = p.into();
        let na_result = iso.transform_point(&na_p);

        prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
        prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
    }
}
```

### 4. Error Handling

For conversions that can fail (e.g., ideal points to Cartesian):

```rust
/// Error when converting from clifford types to nalgebra types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NalgebraConversionError {
    /// Point is ideal (at infinity), cannot convert to finite coordinates.
    IdealPoint,
    /// Motor is not normalized.
    NotNormalized,
}

impl<T: Float + na::Scalar> TryFrom<Point<T>> for na::Point3<T> {
    type Error = NalgebraConversionError;

    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0().abs() < T::epsilon() {
            return Err(NalgebraConversionError::IdealPoint);
        }
        Ok(na::Point3::new(p.x(), p.y(), p.z()))
    }
}
```

Note: We also provide infallible `From` implementations that use the accessors `x()`, `y()`, `z()` which already handle the division by `w` and may produce NaN for ideal points. Users who need error handling can use `TryFrom`.

## Implementation Plan

### Phase 1: dim3 nalgebra (Priority)

1. Add `nalgebra.rs` module to `projective::dim3`
2. Implement Point <-> Point3 conversions
3. Implement Motor <-> Isometry3 conversions
4. Implement Motor <-> UnitDualQuaternion conversions
5. Add comprehensive property-based tests
6. Wire up module in `mod.rs`

### Phase 2: dim2 nalgebra

1. Add `nalgebra.rs` module to `projective::dim2`
2. Implement Point <-> Point2 conversions
3. Implement Motor <-> Isometry2 conversions
4. Add comprehensive property-based tests
5. Wire up module in `mod.rs`

### Phase 3: Rerun dim2

1. Add `rerun.rs` module to `projective::dim2`
2. Implement Point -> Position2D
3. Implement Motor -> Transform3D (2D embedded in 3D)
4. Implement Line -> LineStrips2D
5. Wire up module in `mod.rs`

### Phase 4: Enhanced Testing

1. Add cross-library equivalence tests to integration tests
2. Add edge case tests (identity, near-singular, etc.)
3. Benchmark conversion overhead

## Files to Create/Modify

### New Files
- `src/specialized/projective/dim3/nalgebra.rs`
- `src/specialized/projective/dim2/nalgebra.rs`
- `src/specialized/projective/dim2/rerun.rs`

### Modified Files
- `src/specialized/projective/dim3/mod.rs` - Enable nalgebra module
- `src/specialized/projective/dim2/mod.rs` - Enable nalgebra and rerun modules

## Success Criteria

1. All nalgebra conversions compile and pass tests
2. Property-based tests demonstrate operational equivalence:
   - `motor.transform_point(&p)` ≈ `isometry.transform_point(&p)`
   - `m1.compose(&m2)` ≈ `iso1 * iso2`
   - `motor.inverse()` ≈ `isometry.inverse()`
3. Factory methods produce equivalent transforms
4. Rerun visualization works for dim2
5. Documentation clearly explains mathematical correspondences
6. No performance regression in existing benchmarks

## Testing Checklist

- [ ] Point3 round-trip preserves coordinates
- [ ] Point2 round-trip preserves coordinates
- [ ] Motor3 round-trip preserves transformation behavior
- [ ] Motor2 round-trip preserves transformation behavior
- [ ] `transform_point` equivalence (3D)
- [ ] `transform_point` equivalence (2D)
- [ ] `compose` equivalence (3D)
- [ ] `compose` equivalence (2D)
- [ ] `inverse` equivalence (3D)
- [ ] `inverse` equivalence (2D)
- [ ] `from_translation` equivalence
- [ ] `from_rotation_x/y/z` equivalence
- [ ] `from_axis_angle` equivalence
- [ ] Ideal point handling (error or NaN)
- [ ] Identity motor/isometry equivalence
- [ ] Large translation values
- [ ] Small rotation angles (near identity)

## References

- [nalgebra Isometry documentation](https://docs.rs/nalgebra/latest/nalgebra/geometry/type.Isometry3.html)
- [Rigid Geometric Algebra Wiki - Motors](https://rigidgeometricalgebra.org/wiki/index.php?title=Motor)
- [Geometric Algebra for Computer Science - Chapter 11](http://www.geometricalgebra.net/)
- Existing euclidean nalgebra interop: `src/specialized/euclidean/dim3/nalgebra.rs`
