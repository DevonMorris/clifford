# PRD-15.3: Projective 3D Codegen Migration

**Status**: Proposed
**Parent**: PRD-15 (Codegen Migration)
**Depends On**: PRD-15.1, PRD-15.2 (pattern validation)
**Goal**: Migrate `specialized/projective/dim3` from hand-rolled to generated code

## Overview

Projective 3D is the most complex algebra with 5 types and geometric constraints (Study condition, Plücker condition). This migration tests constraint handling in the codegen framework.

### Current Types

| Type | Grades | Fields | Constraints | Key Methods |
|------|--------|--------|-------------|-------------|
| `Point` | [1] | e1, e2, e3, e0 | None | `origin()`, `ideal()`, `is_ideal()`, `to_cartesian()`, `join()`, `distance()`, `midpoint()` |
| `Line` | [2] | e01, e02, e03, e23, e31, e12 | Plücker | `from_point_and_direction()`, `direction()`, `moment()`, `satisfies_plucker_condition()` |
| `Plane` | [3] | e032, e013, e021, e123 | None | `from_normal_and_distance()`, `normal()`, `unitized()` |
| `Motor` | [0,2,4] | s, e23, e31, e12, e01, e02, e03, e0123 | Study | `from_translation()`, `from_axis_angle()`, `transform_*()`, `compose()`, `satisfies_study_condition()` |
| `Flector` | [1,3] | e1, e2, e3, e0, e032, e013, e021, e123 | None | `from_plane()`, `transform_*()`, `is_pure_reflection()` |

### Key Challenges

1. **Study condition**: Motors must satisfy `s·e₀₁₂₃ - e₂₃·e₀₁ - e₃₁·e₀₂ - e₁₂·e₀₃ = 0`
2. **Plücker condition**: Lines must satisfy `d·m = e₀₁·e₂₃ + e₀₂·e₃₁ + e₀₃·e₁₂ = 0`
3. **Complex products**: 8-component motor composition
4. **Sandwich products**: Motor transformations of points, lines, planes

## Lessons Learned from PRD-15.1

The following lessons were learned during the Euclidean 3D migration:

1. **Flat field constructors**: Generated types use flat fields, not nested types. For example:
   - Old: `Rotor::new(s, Bivector::new(xy, xz, yz))`
   - New: `Rotor::new(s, xy, xz, yz)`

2. **Remove wrapper types**: Don't create wrapper types like `UnitVector`, `UnitRotor`, `NonZeroVector`. Instead:
   - Use `any::<Type<f64>>()` with `.normalized()` in tests
   - The generated Arbitrary implementations are sufficient

3. **Generated conversions handle From traits**: The generated `conversions.rs` includes:
   - `From<Type> for Multivector<T, Signature>`
   - `From<Multivector<T, Signature>> for Type` (via `from_multivector_unchecked`)
   - Don't create a separate hand-written `conversions.rs`

4. **Doc tests need explicit type annotations**: Generic types need turbofish syntax:
   - ✓ `Vector::<f64>::unit_x()`
   - ✗ `Vector::unit_x()` (fails with "type annotations needed")

5. **Update dependent modules**: Modules that import types from the migrated module need their API calls updated to match the new constructor signatures. For example, projective/dim3/conversions.rs needed updates when euclidean/dim3 Rotor changed.

6. **extensions.rs pattern**: Domain-specific methods (transform_point, compose, from_axis_angle, etc.) go in `extensions.rs`, importing from `generated/products` and `generated/types`.

7. **No arbitrary.rs needed**: The generated code includes Arbitrary implementations for all types. Delete any hand-written arbitrary.rs.

8. **Imports may conflict**: If both hand-written and generated code define the same From impls or Arbitrary impls, you'll get conflicts. Remove the hand-written versions entirely.

9. **Unit constraint as type-level invariant**: For types that must satisfy a unit constraint (like Rotor, Motor), use `solve_for` in the TOML constraint to make the constraint a type-level invariant:
   ```toml
   constraints = [
       { name = "unit", expression = "s*s + xy*xy + xz*xz + yz*yz = 1", solve_for = "s" }
   ]
   ```
   This generates:
   - `new(xy, xz, yz) -> Option<Self>` - computes `s` from constraint, returns `None` if invalid
   - `new_unchecked(s, xy, xz, yz) -> Self` - bypasses check for when constraint is guaranteed

10. **Use new_unchecked for internal operations**: When operations mathematically preserve a constraint (like motor composition, conversions from Euclidean types), use `new_unchecked` since the constraint is guaranteed by construction.

11. **nalgebra interop tests in nalgebra.rs**: Keep nalgebra conversion tests in the nalgebra.rs file using generated Arbitrary impls, not in a separate arbitrary.rs with wrapper types.

12. **Motor-to-Rotor conversions**: When converting a unit Motor's rotation part to a Euclidean Rotor, normalization is unnecessary since the rotation part of a unit motor is already unit

13. **Use `normalize()` from generated code**: Don't add a `normalized()` method to extensions.rs - the generated code already provides `normalize()` and `try_normalize()` methods. Use these consistently

## Phase 1: Create TOML Specification

### Deliverable: `algebras/projective3.toml`

```toml
# 3D Projective Geometric Algebra (Point-based)
# Signature: Cl(3,0,1)
# Basis: e1, e2, e3, e0 where e0² = 0

[algebra]
name = "projective3"
module_path = "projective::dim3"
description = "3D Projective Geometric Algebra (Point-based PGA)"

[signature]
positive = ["e1", "e2", "e3"]
zero = ["e0"]

[blades]
# Grade 1
e1 = "e1"
e2 = "e2"
e3 = "e3"
e0 = "e0"
# Grade 2 (Euclidean bivectors)
e23 = "e23"
e31 = "e31"
e12 = "e12"
# Grade 2 (Translation bivectors)
e01 = "e01"
e02 = "e02"
e03 = "e03"
# Grade 3
e032 = "e032"
e013 = "e013"
e021 = "e021"
e123 = "e123"
# Grade 4
e0123 = "e0123"

[types.Point]
grades = [1]
description = "3D point in homogeneous coordinates"
fields = ["e1", "e2", "e3", "e0"]

[types.Line]
grades = [2]
description = "3D Plücker line"
fields = ["e01", "e02", "e03", "e23", "e31", "e12"]

# Plücker condition: direction · moment = 0
[types.Line.constraints.plucker]
expression = "e23*e01 + e31*e02 + e12*e03 = 0"
# Note: This is a geometric constraint, not used for construction

[types.Plane]
grades = [3]
description = "3D plane"
fields = ["e032", "e013", "e021", "e123"]

[types.Motor]
grades = [0, 2, 4]
description = "3D rigid transformation (rotation + translation)"
fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03", "e0123"]

# Study condition: s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0
[types.Motor.constraints.study]
expression = "s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0"
solve_for = "e0123"
enforce = "project"

[types.Flector]
grades = [1, 3]
description = "3D reflection/glide transformation"
fields = ["e1", "e2", "e3", "e0", "e032", "e013", "e021", "e123"]

[products.geometric]
Motor_Motor = "Motor"

[products.regressive]
Point_Point = "Line"
Point_Line = "Plane"

[products.outer]
Plane_Plane = "Line"
Plane_Line = "Point"
Line_Line = "Point"

[products.sandwich]
Motor_Point = "Point"
Motor_Line = "Line"
Motor_Plane = "Plane"
Flector_Point = "Point"
Flector_Line = "Line"
Flector_Plane = "Plane"

[options]
generate_serde = true
generate_arbitrary = true
generate_tests = true
```

## Phase 2: Generate Code

### Command

```bash
cargo run --package clifford-codegen -- generate algebras/projective3.toml \
    -o src/specialized/projective/dim3/generated/ --force
```

### Expected Output

```
src/specialized/projective/dim3/generated/
├── mod.rs
├── types.rs      # Point, Line, Plane, Motor, Flector
├── products.rs   # Product functions including sandwich
└── traits.rs     # Trait implementations
```

### Generated Motor Type (with Study constraint)

```rust
// types.rs (generated)

/// 3D rigid transformation (rotation + translation)
///
/// Motors satisfy the Study condition:
/// `s·e₀₁₂₃ - e₂₃·e₀₁ - e₃₁·e₀₂ - e₁₂·e₀₃ = 0`
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Motor<T: Float> {
    s: T,
    e23: T,
    e31: T,
    e12: T,
    e01: T,
    e02: T,
    e03: T,
    e0123: T,
}

impl<T: Float> Motor<T> {
    /// Create motor from components without validation.
    ///
    /// # Safety
    ///
    /// Caller must ensure Study condition is satisfied.
    #[inline]
    pub fn new_unchecked(
        s: T, e23: T, e31: T, e12: T,
        e01: T, e02: T, e03: T, e0123: T,
    ) -> Self {
        Self { s, e23, e31, e12, e01, e02, e03, e0123 }
    }

    /// Create motor with e0123 computed from Study condition.
    ///
    /// The e0123 component is determined by:
    /// `e0123 = (e23*e01 + e31*e02 + e12*e03) / s`
    ///
    /// Returns `None` if `s ≈ 0` (pure translation case requires different handling).
    pub fn try_from_components(
        s: T, e23: T, e31: T, e12: T,
        e01: T, e02: T, e03: T,
    ) -> Option<Self> {
        if s.abs() < T::epsilon() {
            return None;
        }
        let e0123 = (e23 * e01 + e31 * e02 + e12 * e03) / s;
        Some(Self::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123))
    }

    /// Project motor onto Study manifold.
    ///
    /// Adjusts e0123 to satisfy Study condition.
    pub fn from_components_projected(
        s: T, e23: T, e31: T, e12: T,
        e01: T, e02: T, e03: T, _e0123: T,
    ) -> Self {
        if s.abs() < T::epsilon() {
            // Pure translation: e0123 = 0
            Self::new_unchecked(s, e23, e31, e12, e01, e02, e03, T::zero())
        } else {
            let e0123 = (e23 * e01 + e31 * e02 + e12 * e03) / s;
            Self::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
        }
    }

    // Accessors
    pub fn s(&self) -> T { self.s }
    pub fn e23(&self) -> T { self.e23 }
    pub fn e31(&self) -> T { self.e31 }
    pub fn e12(&self) -> T { self.e12 }
    pub fn e01(&self) -> T { self.e01 }
    pub fn e02(&self) -> T { self.e02 }
    pub fn e03(&self) -> T { self.e03 }
    pub fn e0123(&self) -> T { self.e0123 }

    // Norm operations
    pub fn weight_norm_squared(&self) -> T {
        self.s * self.s + self.e23 * self.e23 + self.e31 * self.e31 + self.e12 * self.e12
    }

    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            self.s / wn,
            self.e23 / wn,
            self.e31 / wn,
            self.e12 / wn,
            self.e01 / wn,
            self.e02 / wn,
            self.e03 / wn,
            self.e0123 / wn,
        )
    }
}
```

### Generated Sandwich Products

```rust
// products.rs (generated)

/// Sandwich product: Motor * Point * reverse(Motor)
#[inline]
pub fn sandwich_motor_point<T: Float>(m: &Motor<T>, p: &Point<T>) -> Point<T> {
    // Generated from symbolic computation of M P M̃
    let s = m.s();
    let e23 = m.e23();
    let e31 = m.e31();
    let e12 = m.e12();
    let e01 = m.e01();
    let e02 = m.e02();
    let e03 = m.e03();

    let px = p.e1();
    let py = p.e2();
    let pz = p.e3();
    let pw = p.e0();

    // Rotation part (3x3 orthogonal transformation)
    let r00 = s*s + e23*e23 - e31*e31 - e12*e12;
    let r01 = T::TWO * (e23*e31 - s*e12);
    let r02 = T::TWO * (e23*e12 + s*e31);
    let r10 = T::TWO * (e23*e31 + s*e12);
    let r11 = s*s - e23*e23 + e31*e31 - e12*e12;
    let r12 = T::TWO * (e31*e12 - s*e23);
    let r20 = T::TWO * (e23*e12 - s*e31);
    let r21 = T::TWO * (e31*e12 + s*e23);
    let r22 = s*s - e23*e23 - e31*e31 + e12*e12;

    // Translation part
    let tx = T::TWO * (s*e01 + e23*e0123 + e31*e03 - e12*e02);
    let ty = T::TWO * (s*e02 + e31*e0123 + e12*e01 - e23*e03);
    let tz = T::TWO * (s*e03 + e12*e0123 + e23*e02 - e31*e01);

    Point::new(
        r00*px + r01*py + r02*pz + tx*pw,
        r10*px + r11*py + r12*pz + ty*pw,
        r20*px + r21*py + r22*pz + tz*pw,
        pw,
    )
}

/// Sandwich product: Motor * Line * reverse(Motor)
#[inline]
pub fn sandwich_motor_line<T: Float>(m: &Motor<T>, l: &Line<T>) -> Line<T> {
    // Generated from symbolic computation
    // ... (complex 6x6 transformation)
}

/// Sandwich product: Motor * Plane * reverse(Motor)
#[inline]
pub fn sandwich_motor_plane<T: Float>(m: &Motor<T>, p: &Plane<T>) -> Plane<T> {
    // Generated from symbolic computation
    // ... (4x4 transformation)
}

/// Geometric product of two motors
#[inline]
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // Generated 8x8 product
    // Result automatically satisfies Study condition
    // ...
}

/// Regressive product: Point ∨ Point = Line
#[inline]
pub fn regressive_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    // Join of two points gives the line through them
    Line::new_unchecked(
        a.e1()*b.e0() - a.e0()*b.e1(),
        a.e2()*b.e0() - a.e0()*b.e2(),
        a.e3()*b.e0() - a.e0()*b.e3(),
        a.e2()*b.e3() - a.e3()*b.e2(),
        a.e3()*b.e1() - a.e1()*b.e3(),
        a.e1()*b.e2() - a.e2()*b.e1(),
    )
}
```

## Phase 3: Create Extensions Module

### File: `src/specialized/projective/dim3/extensions.rs`

```rust
//! Domain-specific extensions for 3D PGA types.

use super::generated::products;
use super::generated::types::{Flector, Line, Motor, Plane, Point};
use crate::scalar::Float;
use crate::specialized::euclidean::dim3::Vector as EuclideanVector;

// ============================================================================
// Point extensions
// ============================================================================

impl<T: Float> Point<T> {
    /// Point at the origin.
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Ideal point (point at infinity) in direction (dx, dy, dz).
    #[inline]
    pub fn ideal(dx: T, dy: T, dz: T) -> Self {
        Self::new(dx, dy, dz, T::zero())
    }

    /// Create from homogeneous coordinates.
    #[inline]
    pub fn from_homogeneous(e1: T, e2: T, e3: T, e0: T) -> Self {
        Self::new(e1, e2, e3, e0)
    }

    /// Returns true if this is an ideal point.
    #[inline]
    pub fn is_ideal(&self) -> bool {
        self.e0().abs() < T::epsilon()
    }

    /// Returns true if this is a finite point.
    #[inline]
    pub fn is_finite(&self) -> bool {
        !self.is_ideal()
    }

    /// Cartesian x-coordinate (requires finite point).
    #[inline]
    pub fn x(&self) -> T {
        self.e1() / self.e0()
    }

    /// Cartesian y-coordinate (requires finite point).
    #[inline]
    pub fn y(&self) -> T {
        self.e2() / self.e0()
    }

    /// Cartesian z-coordinate (requires finite point).
    #[inline]
    pub fn z(&self) -> T {
        self.e3() / self.e0()
    }

    /// Homogeneous weight.
    #[inline]
    pub fn w(&self) -> T {
        self.e0()
    }

    /// Convert to Cartesian coordinates.
    pub fn to_cartesian(&self) -> Option<(T, T, T)> {
        if self.is_ideal() {
            None
        } else {
            Some((self.x(), self.y(), self.z()))
        }
    }

    /// Join with another point to create a line.
    #[inline]
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::regressive_point_point(self, other)
    }

    /// Euclidean distance to another finite point.
    pub fn distance(&self, other: &Point<T>) -> T {
        self.distance_squared(other).sqrt()
    }

    /// Squared Euclidean distance.
    pub fn distance_squared(&self, other: &Point<T>) -> T {
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        let dz = self.z() - other.z();
        dx * dx + dy * dy + dz * dz
    }

    /// Midpoint between two finite points.
    pub fn midpoint(&self, other: &Point<T>) -> Point<T> {
        let two = T::TWO;
        Point::new(
            (self.e1() * other.e0() + other.e1() * self.e0()) / two,
            (self.e2() * other.e0() + other.e2() * self.e0()) / two,
            (self.e3() * other.e0() + other.e3() * self.e0()) / two,
            self.e0() * other.e0(),
        )
    }

    /// Weight norm (absolute value of e0).
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.e0().abs()
    }

    /// Bulk norm (Euclidean norm of direction part).
    #[inline]
    pub fn bulk_norm(&self) -> T {
        (self.e1() * self.e1() + self.e2() * self.e2() + self.e3() * self.e3()).sqrt()
    }

    /// Normalize to unit weight.
    pub fn normalize(&self) -> Option<Self> {
        let w = self.e0();
        if w.abs() < T::epsilon() {
            None
        } else {
            Some(Self::new(
                self.e1() / w,
                self.e2() / w,
                self.e3() / w,
                T::one(),
            ))
        }
    }
}

// ============================================================================
// Line extensions
// ============================================================================

impl<T: Float> Line<T> {
    /// Create line through two points.
    #[inline]
    pub fn through_points(p1: &Point<T>, p2: &Point<T>) -> Self {
        p1.join(p2)
    }

    /// Create line from point and direction.
    pub fn from_point_and_direction(point: &Point<T>, direction: &EuclideanVector<T>) -> Self {
        let ideal = Point::ideal(direction.x(), direction.y(), direction.z());
        point.join(&ideal)
    }

    /// Create from Plücker coordinates without validation.
    #[inline]
    pub fn from_plucker_unchecked(
        e01: T, e02: T, e03: T,
        e23: T, e31: T, e12: T,
    ) -> Self {
        Self::new_unchecked(e01, e02, e03, e23, e31, e12)
    }

    /// Create from Plücker coordinates with validation.
    pub fn from_plucker(
        e01: T, e02: T, e03: T,
        e23: T, e31: T, e12: T,
        tolerance: T,
    ) -> Result<Self, super::super::errors::LineConstraintError> {
        let line = Self::new_unchecked(e01, e02, e03, e23, e31, e12);
        if line.plucker_residual().abs() > tolerance {
            return Err(super::super::errors::LineConstraintError::PluckerConditionViolation);
        }
        Ok(line)
    }

    /// Direction vector (e23, e31, e12).
    #[inline]
    pub fn direction(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.e23(), self.e31(), self.e12())
    }

    /// Moment vector (e01, e02, e03).
    #[inline]
    pub fn moment(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.e01(), self.e02(), self.e03())
    }

    /// X-axis (line through origin along x).
    #[inline]
    pub fn x_axis() -> Self {
        Self::new_unchecked(T::zero(), T::zero(), T::zero(), T::one(), T::zero(), T::zero())
    }

    /// Y-axis.
    #[inline]
    pub fn y_axis() -> Self {
        Self::new_unchecked(T::zero(), T::zero(), T::zero(), T::zero(), T::one(), T::zero())
    }

    /// Z-axis.
    #[inline]
    pub fn z_axis() -> Self {
        Self::new_unchecked(T::zero(), T::zero(), T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Weight norm (direction magnitude).
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.direction().norm()
    }

    /// Unitize to unit direction.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            self.e01() / wn,
            self.e02() / wn,
            self.e03() / wn,
            self.e23() / wn,
            self.e31() / wn,
            self.e12() / wn,
        )
    }

    /// Plücker condition residual: d · m.
    #[inline]
    pub fn plucker_residual(&self) -> T {
        self.e23() * self.e01() + self.e31() * self.e02() + self.e12() * self.e03()
    }

    /// Check if line satisfies Plücker condition.
    #[inline]
    pub fn satisfies_plucker_condition(&self, tolerance: T) -> bool {
        self.plucker_residual().abs() < tolerance
    }

    /// Check if line passes through origin.
    pub fn through_origin(&self) -> bool {
        self.moment().norm() < T::epsilon()
    }

    /// Check if two lines are parallel.
    pub fn is_parallel(&self, other: &Line<T>) -> bool {
        self.direction().cross(&other.direction()).norm() < T::epsilon()
    }

    /// Check if two lines intersect.
    pub fn intersects(&self, other: &Line<T>) -> bool {
        // Lines intersect if their join is zero (they're coplanar and not parallel)
        let d1 = self.direction();
        let d2 = other.direction();
        let m1 = self.moment();
        let m2 = other.moment();

        // Reciprocal product: d1·m2 + d2·m1 = 0 for intersecting lines
        let reciprocal = d1.dot(&m2) + d2.dot(&m1);
        reciprocal.abs() < T::epsilon() && !self.is_parallel(other)
    }

    /// Reverse.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new_unchecked(
            -self.e01(), -self.e02(), -self.e03(),
            -self.e23(), -self.e31(), -self.e12(),
        )
    }
}

// ============================================================================
// Plane extensions
// ============================================================================

impl<T: Float> Plane<T> {
    /// Create plane from normal and distance.
    ///
    /// The plane equation is `n·x + d = 0`.
    pub fn from_normal_and_distance(nx: T, ny: T, nz: T, d: T) -> Self {
        // In point-based PGA:
        // e032 corresponds to x coefficient
        // e013 corresponds to y coefficient
        // e021 corresponds to z coefficient
        // e123 corresponds to constant term
        Self::new(nx, ny, nz, d)
    }

    /// XY plane (z = 0).
    #[inline]
    pub fn xy() -> Self {
        Self::from_normal_and_distance(T::zero(), T::zero(), T::one(), T::zero())
    }

    /// XZ plane (y = 0).
    #[inline]
    pub fn xz() -> Self {
        Self::from_normal_and_distance(T::zero(), T::one(), T::zero(), T::zero())
    }

    /// YZ plane (x = 0).
    #[inline]
    pub fn yz() -> Self {
        Self::from_normal_and_distance(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Normal vector.
    #[inline]
    pub fn normal(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.e032(), self.e013(), self.e021())
    }

    /// Distance from origin (signed).
    #[inline]
    pub fn distance(&self) -> T {
        self.e123()
    }

    /// Weight norm (normal magnitude).
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.normal().norm()
    }

    /// Unitize to unit normal.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new(
            self.e032() / wn,
            self.e013() / wn,
            self.e021() / wn,
            self.e123() / wn,
        )
    }

    /// Reverse.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(-self.e032(), -self.e013(), -self.e021(), -self.e123())
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Identity motor (no transformation).
    #[inline]
    pub fn identity() -> Self {
        Self::new_unchecked(
            T::one(),
            T::zero(), T::zero(), T::zero(),
            T::zero(), T::zero(), T::zero(),
            T::zero(),
        )
    }

    /// Pure translation motor.
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self {
        let half = T::one() / T::TWO;
        Self::new_unchecked(
            T::one(),
            T::zero(), T::zero(), T::zero(),
            dx * half, dy * half, dz * half,
            T::zero(),
        )
    }

    /// Pure rotation around x-axis through origin.
    pub fn from_rotation_x(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            half.cos(),
            half.sin(), T::zero(), T::zero(),
            T::zero(), T::zero(), T::zero(),
            T::zero(),
        )
    }

    /// Pure rotation around y-axis through origin.
    pub fn from_rotation_y(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            half.cos(),
            T::zero(), half.sin(), T::zero(),
            T::zero(), T::zero(), T::zero(),
            T::zero(),
        )
    }

    /// Pure rotation around z-axis through origin.
    pub fn from_rotation_z(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            half.cos(),
            T::zero(), T::zero(), half.sin(),
            T::zero(), T::zero(), T::zero(),
            T::zero(),
        )
    }

    /// Rotation around arbitrary axis through origin.
    pub fn from_axis_angle(axis: &EuclideanVector<T>, angle: T) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let axis_norm = axis.normalized();
        Self::new_unchecked(
            cos_half,
            sin_half * axis_norm.x(),
            sin_half * axis_norm.y(),
            sin_half * axis_norm.z(),
            T::zero(), T::zero(), T::zero(),
            T::zero(),
        )
    }

    /// Screw motion along a line.
    pub fn from_line(line: &Line<T>, angle: T, distance: T) -> Self {
        let line_unit = line.unitized();
        let half_angle = angle / T::TWO;
        let half_dist = distance / T::TWO;

        let (sin_a, cos_a) = (half_angle.sin(), half_angle.cos());

        let d = line_unit.direction();
        let m = line_unit.moment();

        Self::new_unchecked(
            cos_a,
            sin_a * d.x(),
            sin_a * d.y(),
            sin_a * d.z(),
            sin_a * m.x() + half_dist * cos_a * d.x(),
            sin_a * m.y() + half_dist * cos_a * d.y(),
            sin_a * m.z() + half_dist * cos_a * d.z(),
            -half_dist * sin_a,
        )
    }

    /// Transform a point.
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::sandwich_motor_point(self, p)
    }

    /// Transform a line.
    #[inline]
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        products::sandwich_motor_line(self, l)
    }

    /// Transform a plane.
    #[inline]
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T> {
        products::sandwich_motor_plane(self, p)
    }

    /// Compose motors: self then other.
    #[inline]
    pub fn compose(&self, other: &Motor<T>) -> Motor<T> {
        products::geometric_motor_motor(other, self)
    }

    /// Reverse.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new_unchecked(
            self.s(),
            -self.e23(), -self.e31(), -self.e12(),
            -self.e01(), -self.e02(), -self.e03(),
            self.e0123(),
        )
    }

    /// Inverse motor.
    pub fn inverse(&self) -> Self {
        let wn_sq = self.weight_norm_squared();
        let rev = self.reverse();
        Self::new_unchecked(
            rev.s() / wn_sq,
            rev.e23() / wn_sq,
            rev.e31() / wn_sq,
            rev.e12() / wn_sq,
            rev.e01() / wn_sq,
            rev.e02() / wn_sq,
            rev.e03() / wn_sq,
            rev.e0123() / wn_sq,
        )
    }

    /// Check if motor is unitized.
    #[inline]
    pub fn is_unitized(&self, tolerance: T) -> bool {
        (self.weight_norm_squared() - T::one()).abs() < tolerance
    }

    /// Study condition residual.
    #[inline]
    pub fn study_residual(&self) -> T {
        self.s() * self.e0123()
            - self.e23() * self.e01()
            - self.e31() * self.e02()
            - self.e12() * self.e03()
    }

    /// Check if motor satisfies Study condition.
    #[inline]
    pub fn satisfies_study_condition(&self, tolerance: T) -> bool {
        self.study_residual().abs() < tolerance
    }

    /// Extract rotation angle.
    pub fn rotation_angle(&self) -> T {
        let cos_half = self.s().min(T::one()).max(-T::one());
        cos_half.acos() * T::TWO
    }

    /// Extract translation vector.
    pub fn translation(&self) -> EuclideanVector<T> {
        let two = T::TWO;
        EuclideanVector::new(
            two * (self.s() * self.e01() + self.e31() * self.e03() - self.e12() * self.e02() + self.e23() * self.e0123()),
            two * (self.s() * self.e02() + self.e12() * self.e01() - self.e23() * self.e03() + self.e31() * self.e0123()),
            two * (self.s() * self.e03() + self.e23() * self.e02() - self.e31() * self.e01() + self.e12() * self.e0123()),
        )
    }
}

// ============================================================================
// Flector extensions
// ============================================================================

impl<T: Float> Flector<T> {
    /// Create flector from reflection plane.
    pub fn from_plane(plane: &Plane<T>) -> Self {
        let p = plane.unitized();
        Self::new(
            T::zero(), T::zero(), T::zero(), T::zero(),
            p.e032(), p.e013(), p.e021(), p.e123(),
        )
    }

    /// Create reflection through plane at origin.
    pub fn from_plane_through_origin(nx: T, ny: T, nz: T) -> Self {
        let norm = (nx*nx + ny*ny + nz*nz).sqrt();
        Self::new(
            T::zero(), T::zero(), T::zero(), T::zero(),
            nx / norm, ny / norm, nz / norm, T::zero(),
        )
    }

    /// Reflect through XY plane.
    #[inline]
    pub fn reflect_xy() -> Self {
        Self::from_plane(&Plane::xy())
    }

    /// Reflect through XZ plane.
    #[inline]
    pub fn reflect_xz() -> Self {
        Self::from_plane(&Plane::xz())
    }

    /// Reflect through YZ plane.
    #[inline]
    pub fn reflect_yz() -> Self {
        Self::from_plane(&Plane::yz())
    }

    /// Point part (grade 1).
    #[inline]
    pub fn point_part(&self) -> Point<T> {
        Point::new(self.e1(), self.e2(), self.e3(), self.e0())
    }

    /// Plane part (grade 3).
    #[inline]
    pub fn plane_part(&self) -> Plane<T> {
        Plane::new(self.e032(), self.e013(), self.e021(), self.e123())
    }

    /// Check if this is a pure reflection (no point part).
    #[inline]
    pub fn is_pure_reflection(&self) -> bool {
        let pt = self.point_part();
        pt.bulk_norm() < T::epsilon() && pt.weight_norm() < T::epsilon()
    }

    /// Transform a point.
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::sandwich_flector_point(self, p)
    }

    /// Transform a line.
    #[inline]
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        products::sandwich_flector_line(self, l)
    }

    /// Transform a plane.
    #[inline]
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T> {
        products::sandwich_flector_plane(self, p)
    }

    /// Reverse.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(
            self.e1(), self.e2(), self.e3(), self.e0(),
            -self.e032(), -self.e013(), -self.e021(), -self.e123(),
        )
    }
}
```

## Phase 4: Tests

### File: `src/specialized/projective/dim3/tests.rs`

```rust
//! Tests for 3D PGA codegen migration.

use super::*;
use crate::multivector::Multivector;
use crate::signature::Projective3;
use crate::test_utils::ABS_DIFF_EQ_EPS;
use approx::abs_diff_eq;
use proptest::prelude::*;

proptest! {
    // ========================================================================
    // Motor constraints
    // ========================================================================

    /// All constructed motors satisfy Study condition.
    #[test]
    fn motor_satisfies_study(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0,
        ty in -10.0f64..10.0,
        tz in -10.0f64..10.0,
    ) {
        let rotation = Motor::from_rotation_z(angle);
        let translation = Motor::from_translation(tx, ty, tz);
        let combined = translation.compose(&rotation);

        prop_assert!(rotation.satisfies_study_condition(ABS_DIFF_EQ_EPS));
        prop_assert!(translation.satisfies_study_condition(ABS_DIFF_EQ_EPS));
        prop_assert!(combined.satisfies_study_condition(ABS_DIFF_EQ_EPS));
    }

    /// Motor composition preserves Study condition.
    #[test]
    fn motor_compose_preserves_study(
        a1 in -std::f64::consts::PI..std::f64::consts::PI,
        a2 in -std::f64::consts::PI..std::f64::consts::PI,
        t1x in -10.0f64..10.0, t1y in -10.0f64..10.0, t1z in -10.0f64..10.0,
        t2x in -10.0f64..10.0, t2y in -10.0f64..10.0, t2z in -10.0f64..10.0,
    ) {
        let m1 = Motor::from_translation(t1x, t1y, t1z).compose(&Motor::from_rotation_z(a1));
        let m2 = Motor::from_translation(t2x, t2y, t2z).compose(&Motor::from_rotation_x(a2));
        let composed = m1.compose(&m2);

        prop_assert!(composed.satisfies_study_condition(ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Line constraints
    // ========================================================================

    /// Lines from point join satisfy Plücker condition.
    #[test]
    fn line_from_join_satisfies_plucker(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
    ) {
        let line = p1.join(&p2);
        prop_assert!(line.satisfies_plucker_condition(ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Transformation consistency
    // ========================================================================

    /// Motor sandwich matches generic Multivector sandwich.
    #[test]
    fn motor_transform_point_matches_multivector(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0, ty in -10.0f64..10.0, tz in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let motor = Motor::from_translation(tx, ty, tz)
            .compose(&Motor::from_rotation_z(angle));
        let point = Point::new(px, py, pz);

        let spec = motor.transform_point(&point);

        let mv_motor = Multivector::<f64, Projective3>::from(motor);
        let mv_point = Multivector::<f64, Projective3>::from(point);
        let gen = mv_motor.sandwich(&mv_point);
        let gen_point = Point::from_multivector_unchecked(&gen);

        if let (Some(spec_cart), Some(gen_cart)) = (spec.to_cartesian(), gen_point.to_cartesian()) {
            prop_assert!(abs_diff_eq!(spec_cart.0, gen_cart.0, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(spec_cart.1, gen_cart.1, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(spec_cart.2, gen_cart.2, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    /// Motor preserves distance between points.
    #[test]
    fn motor_preserves_distance(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0, ty in -10.0f64..10.0, tz in -10.0f64..10.0,
        p1x in -10.0f64..10.0, p1y in -10.0f64..10.0, p1z in -10.0f64..10.0,
        p2x in -10.0f64..10.0, p2y in -10.0f64..10.0, p2z in -10.0f64..10.0,
    ) {
        let motor = Motor::from_translation(tx, ty, tz)
            .compose(&Motor::from_rotation_z(angle))
            .unitized();
        let p1 = Point::new(p1x, p1y, p1z);
        let p2 = Point::new(p2x, p2y, p2z);

        let d_before = p1.distance(&p2);
        let t1 = motor.transform_point(&p1);
        let t2 = motor.transform_point(&p2);
        let d_after = t1.distance(&t2);

        prop_assert!(abs_diff_eq!(d_before, d_after, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Verification Checklist

- [ ] TOML specification with Study and Plücker constraints
- [ ] Code generation succeeds
- [ ] Generated products verified symbolically
- [ ] Extensions module complete with all methods
- [ ] Operators updated
- [ ] Constraint verification tests pass
- [ ] Transformation consistency tests pass
- [ ] All existing tests pass
- [ ] nalgebra conversions work
- [ ] `cargo test --all-features` passes

## Files Changed

| File | Action |
|------|--------|
| `algebras/projective3.toml` | Create |
| `src/specialized/projective/dim3/mod.rs` | Update |
| `src/specialized/projective/dim3/generated/` | Create |
| `src/specialized/projective/dim3/extensions.rs` | Create |
| `src/specialized/projective/dim3/nalgebra.rs` | Update (new constructor API) |
| `src/specialized/projective/dim3/rerun.rs` | Update (new Arbitrary API) |
| `src/specialized/projective/dim3/types.rs` | Delete (replaced by generated) |
| `src/specialized/projective/dim3/ops.rs` | Delete (replaced by generated traits) |
| `src/specialized/projective/dim3/arbitrary.rs` | Delete (generated handles this) |
| `src/specialized/projective/dim3/conversions.rs` | Delete (generated handles this) |
