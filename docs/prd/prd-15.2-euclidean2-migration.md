# PRD-15.2: Euclidean 2D Codegen Migration

**Status**: Proposed
**Parent**: PRD-15 (Codegen Migration)
**Depends On**: PRD-15.1 (pattern validation)
**Goal**: Migrate `specialized/euclidean/dim2` from hand-rolled to generated code

## Overview

Euclidean 2D is simpler than 3D with 3 types (Vector, Bivector, Rotor). This migration validates the pattern established in PRD-15.1.

### Current Types

| Type | Grades | Fields | Key Methods |
|------|--------|--------|-------------|
| `Vector` | [1] | x, y | `dot()`, `wedge()`, `geometric()`, `perp()` |
| `Bivector` | [2] | xy | `reverse()` |
| `Rotor` | [0,2] | s, xy | `from_angle()`, `from_vectors()`, `rotate()`, `compose()`, `inverse()`, `lerp()`, `slerp()`, `angle()` |

### Key Differences from 3D

- No `Trivector` (bivector is the pseudoscalar)
- No `cross()` product (2D specific)
- No `dual()` on bivector (trivial in 2D)
- Simpler rotor formulas

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

5. **Update dependent modules**: Modules that import types from the migrated module need their API calls updated to match the new constructor signatures.

6. **extensions.rs pattern**: Domain-specific methods (dot, wedge, cross, from_angle, rotate, etc.) go in `extensions.rs`, importing from `generated/products` and `generated/types`.

7. **No arbitrary.rs needed**: The generated code includes Arbitrary implementations for all types. Delete any hand-written arbitrary.rs.

8. **Unit constraint as type-level invariant**: For types that must satisfy a unit constraint (like Rotor), use `solve_for` in the TOML constraint to make the constraint a type-level invariant:
   ```toml
   constraints = [
       { name = "unit", expression = "s*s + xy*xy + xz*xz + yz*yz = 1", solve_for = "s" }
   ]
   ```
   This generates:
   - `new(xy, xz, yz) -> Option<Self>` - computes `s` from constraint, returns `None` if invalid
   - `new_unchecked(s, xy, xz, yz) -> Self` - bypasses check for when constraint is guaranteed

9. **Use new_unchecked for internal operations**: When operations mathematically preserve the unit constraint (like rotor composition, slerp, from_angle_plane), use `new_unchecked` since normalization is unnecessary.

10. **nalgebra interop tests in nalgebra.rs**: Keep nalgebra conversion tests in the nalgebra.rs file using generated Arbitrary impls, not in a separate arbitrary.rs with wrapper types

11. **Use `normalize()` from generated code**: Don't add a `normalized()` method to extensions.rs - the generated code already provides `normalize()` and `try_normalize()` methods. Use these consistently

## Phase 1: Update TOML Specification

### Deliverable: `algebras/euclidean2.toml`

```toml
# 2D Euclidean Geometric Algebra
# Signature: (2, 0, 0)
# Basis: e1, e2
# Blades: 1, e1, e2, e12

[algebra]
name = "euclidean2"
module_path = "euclidean::dim2"
description = "2D Euclidean Geometric Algebra"

[signature]
positive = ["e1", "e2"]

[blades]
e1 = "x"
e2 = "y"
e12 = "xy"

[types.Vector]
grades = [1]
description = "2D vector (grade-1 element)"
fields = ["x", "y"]

[types.Bivector]
grades = [2]
description = "2D bivector (pseudoscalar, grade-2 element)"
fields = ["xy"]

[types.Rotor]
grades = [0, 2]
description = "2D rotor (rotation element)"
fields = ["s", "xy"]

[types.Even]
grades = [0, 2]
description = "Even subalgebra element"
fields = ["s", "xy"]
alias_of = "Rotor"

[products.geometric]
Vector_Vector = "Rotor"
Rotor_Rotor = "Rotor"

[products.outer]
Vector_Vector = "Bivector"

[products.inner]
Vector_Vector = "T"

[products.scalar]
Vector_Vector = "T"

[products.sandwich]
Rotor_Vector = "Vector"

[options]
generate_serde = true
generate_arbitrary = true
generate_tests = true
```

## Phase 2: Generate Code

### Command

```bash
cargo run --package clifford-codegen -- generate algebras/euclidean2.toml \
    -o src/specialized/euclidean/dim2/generated/ --force
```

### Expected Output

```
src/specialized/euclidean/dim2/generated/
├── mod.rs
├── types.rs      # Vector, Bivector, Rotor, Even
├── products.rs   # Product functions
└── traits.rs     # Trait implementations
```

### Generated Types

```rust
// types.rs (generated)

/// 2D vector (grade-1 element)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Vector<T: Float> {
    x: T,
    y: T,
}

impl<T: Float> Vector<T> {
    pub fn new(x: T, y: T) -> Self { Self { x, y } }
    pub fn x(&self) -> T { self.x }
    pub fn y(&self) -> T { self.y }
    pub fn zero() -> Self { Self::new(T::zero(), T::zero()) }
    pub fn unit_x() -> Self { Self::new(T::one(), T::zero()) }
    pub fn unit_y() -> Self { Self::new(T::zero(), T::one()) }
    pub fn norm_squared(&self) -> T { self.x * self.x + self.y * self.y }
    pub fn norm(&self) -> T { self.norm_squared().sqrt() }
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.x / n, self.y / n)
    }
    pub fn scale(self, s: T) -> Self { Self::new(self.x * s, self.y * s) }
}

/// 2D bivector (pseudoscalar)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Bivector<T: Float> {
    xy: T,
}

impl<T: Float> Bivector<T> {
    pub fn new(xy: T) -> Self { Self { xy } }
    pub fn xy(&self) -> T { self.xy }
    pub fn zero() -> Self { Self::new(T::zero()) }
    pub fn unit() -> Self { Self::new(T::one()) }
    pub fn norm_squared(&self) -> T { self.xy * self.xy }
    pub fn norm(&self) -> T { self.xy.abs() }
    pub fn normalized(&self) -> Self {
        Self::new(self.xy / self.xy.abs())
    }
    pub fn scale(self, s: T) -> Self { Self::new(self.xy * s) }
}

/// 2D rotor (rotation element)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Rotor<T: Float> {
    s: T,
    xy: T,
}

impl<T: Float> Rotor<T> {
    pub fn new(s: T, xy: T) -> Self { Self { s, xy } }
    pub fn s(&self) -> T { self.s }
    pub fn xy(&self) -> T { self.xy }
    pub fn norm_squared(&self) -> T { self.s * self.s + self.xy * self.xy }
    pub fn norm(&self) -> T { self.norm_squared().sqrt() }
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.s / n, self.xy / n)
    }
    pub fn scale(self, s: T) -> Self { Self::new(self.s * s, self.xy * s) }
}
```

### Generated Products

```rust
// products.rs (generated)

/// Geometric product of two 2D vectors.
#[inline]
pub fn geometric_vector_vector<T: Float>(a: &Vector<T>, b: &Vector<T>) -> Rotor<T> {
    Rotor::new(
        a.x() * b.x() + a.y() * b.y(),  // scalar: dot product
        a.x() * b.y() - a.y() * b.x(),  // xy: wedge product
    )
}

/// Outer product of two 2D vectors.
#[inline]
pub fn outer_vector_vector<T: Float>(a: &Vector<T>, b: &Vector<T>) -> Bivector<T> {
    Bivector::new(a.x() * b.y() - a.y() * b.x())
}

/// Scalar product of two 2D vectors.
#[inline]
pub fn scalar_vector_vector<T: Float>(a: &Vector<T>, b: &Vector<T>) -> T {
    a.x() * b.x() + a.y() * b.y()
}

/// Geometric product of two 2D rotors.
#[inline]
pub fn geometric_rotor_rotor<T: Float>(a: &Rotor<T>, b: &Rotor<T>) -> Rotor<T> {
    Rotor::new(
        a.s() * b.s() - a.xy() * b.xy(),
        a.s() * b.xy() + a.xy() * b.s(),
    )
}

/// Sandwich product: rotor * vector * reverse(rotor)
#[inline]
pub fn sandwich_rotor_vector<T: Float>(r: &Rotor<T>, v: &Vector<T>) -> Vector<T> {
    let s = r.s();
    let xy = r.xy();
    let x = v.x();
    let y = v.y();

    // R v R̃ for 2D rotor
    let s2_minus_xy2 = s * s - xy * xy;
    let two_s_xy = T::TWO * s * xy;

    Vector::new(
        x * s2_minus_xy2 + y * two_s_xy,
        y * s2_minus_xy2 - x * two_s_xy,
    )
}
```

## Phase 3: Create Extensions Module

### File: `src/specialized/euclidean/dim2/extensions.rs`

```rust
//! Domain-specific extensions for 2D Euclidean GA types.

use super::generated::products;
use super::generated::types::{Bivector, Rotor, Vector};
use crate::scalar::Float;

// ============================================================================
// Vector extensions
// ============================================================================

impl<T: Float> Vector<T> {
    /// Dot product (inner product).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let a = Vector::new(3.0, 4.0);
    /// let b = Vector::new(1.0, 0.0);
    /// assert_eq!(a.dot(&b), 3.0);
    /// ```
    #[inline]
    pub fn dot(&self, other: &Self) -> T {
        products::scalar_vector_vector(self, other)
    }

    /// Wedge (outer) product.
    ///
    /// Returns the bivector (pseudoscalar in 2D) representing the
    /// signed area of the parallelogram spanned by the vectors.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let a = Vector::unit_x();
    /// let b = Vector::unit_y();
    /// assert_eq!(a.wedge(&b).xy(), 1.0);
    /// ```
    #[inline]
    pub fn wedge(&self, other: &Self) -> Bivector<T> {
        products::outer_vector_vector(self, other)
    }

    /// Geometric product.
    ///
    /// Returns a rotor: `a * b = a·b + a∧b`.
    #[inline]
    pub fn geometric(&self, other: &Self) -> Rotor<T> {
        products::geometric_vector_vector(self, other)
    }

    /// Perpendicular vector (90° counterclockwise rotation).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let v = Vector::new(1.0, 0.0);
    /// let perp = v.perp();
    /// assert_eq!(perp.x(), 0.0);
    /// assert_eq!(perp.y(), 1.0);
    /// ```
    #[inline]
    pub fn perp(&self) -> Self {
        Self::new(-self.y(), self.x())
    }
}

// ============================================================================
// Bivector extensions
// ============================================================================

impl<T: Float> Bivector<T> {
    /// Reverse (reversion).
    ///
    /// For grade-2, reverse negates: `B̃ = -B`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(-self.xy())
    }

    /// Alternative accessor for the single component.
    #[inline]
    pub fn value(&self) -> T {
        self.xy()
    }
}

// ============================================================================
// Rotor extensions
// ============================================================================

impl<T: Float> Rotor<T> {
    /// Identity rotor (no rotation).
    #[inline]
    pub fn identity() -> Self {
        Self::new(T::one(), T::zero())
    }

    /// Create from rotation angle.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::{Rotor, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let r = Rotor::from_angle(FRAC_PI_2);
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// ```
    pub fn from_angle(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new(half.cos(), half.sin())
    }

    /// Create rotation that takes vector `a` to vector `b`.
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        let a_norm = a.normalized();
        let b_norm = b.normalized();
        let ba = b_norm.geometric(&a_norm);

        let unnorm = Self::new(T::one() + ba.s(), ba.xy());

        if unnorm.norm_squared() < T::epsilon() {
            // 180° rotation
            return Self::new(T::zero(), T::one());
        }

        unnorm.normalized()
    }

    /// Reverse (conjugate).
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(self.s(), -self.xy())
    }

    /// Inverse rotor.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        let rev = self.reverse();
        Self::new(rev.s() / norm_sq, rev.xy() / norm_sq)
    }

    /// Apply rotation to vector.
    #[inline]
    pub fn rotate(&self, v: Vector<T>) -> Vector<T> {
        products::sandwich_rotor_vector(self, &v)
    }

    /// Compose rotations.
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        products::geometric_rotor_rotor(other, self)
    }

    /// Linear interpolation.
    #[inline]
    pub fn lerp(&self, other: &Self, t: T) -> Self {
        let one_minus_t = T::one() - t;
        Self::new(
            self.s() * one_minus_t + other.s() * t,
            self.xy() * one_minus_t + other.xy() * t,
        )
    }

    /// Spherical linear interpolation.
    pub fn slerp(&self, other: &Self, t: T) -> Self {
        let rel = other.compose(&self.inverse());
        let angle = rel.angle();

        if angle.abs() < T::epsilon() {
            return *self;
        }

        let interp_rel = Self::from_angle(angle * t);
        interp_rel.compose(self)
    }

    /// Extract rotation angle in radians.
    #[inline]
    pub fn angle(&self) -> T {
        self.xy().atan2(self.s()) * T::TWO
    }
}
```

## Phase 4: Update Module Structure

### File: `src/specialized/euclidean/dim2/mod.rs`

```rust
//! 2D Euclidean Geometric Algebra.
//!
//! This module provides optimized types for 2D Euclidean GA (Cl(2,0,0)):
//!
//! - [`Vector`] - Grade-1 elements (directions, positions)
//! - [`Bivector`] - Grade-2 element (pseudoscalar, oriented area)
//! - [`Rotor`] - Even subalgebra elements (rotations)
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim2::{Vector, Rotor};
//! use std::f64::consts::FRAC_PI_2;
//!
//! let v = Vector::new(1.0, 0.0);
//! let r = Rotor::from_angle(FRAC_PI_2);
//! let rotated = r.rotate(v);
//! ```

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions
mod extensions;

// Note: conversions.rs and arbitrary.rs are NOT needed - generated code handles these

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra;

#[cfg(feature = "rerun-0_28")]
mod rerun;

// Re-export generated types
pub use generated::types::{Bivector, Rotor, Scalar, Vector};

// Re-export Even alias from extensions
pub use extensions::Even;

#[cfg(test)]
mod tests;
```

## Phase 5: Update Operators

### File: `src/specialized/euclidean/dim2/ops.rs`

```rust
//! Operator overloads for 2D Euclidean GA types.

use super::generated::products;
use super::generated::types::{Bivector, Rotor, Vector};
use crate::scalar::Float;
use std::ops::{Add, BitXor, Mul, Neg, Sub};

// ============================================================================
// Vector operators
// ============================================================================

impl<T: Float> Neg for Vector<T> {
    type Output = Self;
    fn neg(self) -> Self { Self::new(-self.x(), -self.y()) }
}

impl<T: Float> Add for Vector<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x() + rhs.x(), self.y() + rhs.y())
    }
}

impl<T: Float> Sub for Vector<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x() - rhs.x(), self.y() - rhs.y())
    }
}

impl<T: Float> Mul<T> for Vector<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self { self.scale(rhs) }
}

impl<T: Float> BitXor<Vector<T>> for Vector<T> {
    type Output = Bivector<T>;
    fn bitxor(self, rhs: Vector<T>) -> Bivector<T> {
        products::outer_vector_vector(&self, &rhs)
    }
}

impl Mul<Vector<f64>> for f64 {
    type Output = Vector<f64>;
    fn mul(self, rhs: Vector<f64>) -> Vector<f64> { rhs.scale(self) }
}

impl Mul<Vector<f32>> for f32 {
    type Output = Vector<f32>;
    fn mul(self, rhs: Vector<f32>) -> Vector<f32> { rhs.scale(self) }
}

// ============================================================================
// Bivector operators
// ============================================================================

impl<T: Float> Neg for Bivector<T> {
    type Output = Self;
    fn neg(self) -> Self { Self::new(-self.xy()) }
}

impl<T: Float> Add for Bivector<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self { Self::new(self.xy() + rhs.xy()) }
}

impl<T: Float> Sub for Bivector<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self { Self::new(self.xy() - rhs.xy()) }
}

impl<T: Float> Mul<T> for Bivector<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self { self.scale(rhs) }
}

// ============================================================================
// Rotor operators
// ============================================================================

impl<T: Float> Neg for Rotor<T> {
    type Output = Self;
    fn neg(self) -> Self { Self::new(-self.s(), -self.xy()) }
}

impl<T: Float> Add for Rotor<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.s() + rhs.s(), self.xy() + rhs.xy())
    }
}

impl<T: Float> Sub for Rotor<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.s() - rhs.s(), self.xy() - rhs.xy())
    }
}

impl<T: Float> Mul<T> for Rotor<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self { self.scale(rhs) }
}

impl<T: Float> Mul<Rotor<T>> for Rotor<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        products::geometric_rotor_rotor(&self, &rhs)
    }
}
```

## Phase 6: Tests

### File: `src/specialized/euclidean/dim2/tests.rs`

```rust
//! Tests for 2D Euclidean GA codegen migration.

use super::*;
use crate::multivector::Multivector;
use crate::signature::Euclidean2;
use crate::test_utils::ABS_DIFF_EQ_EPS;
use approx::abs_diff_eq;
use proptest::prelude::*;

proptest! {
    /// Generated geometric product matches generic Multivector.
    #[test]
    fn vector_geometric_matches_multivector(
        a in any::<Vector<f64>>(),
        b in any::<Vector<f64>>(),
    ) {
        let spec = a.geometric(&b);

        let mv_a = Multivector::<f64, Euclidean2>::from(a);
        let mv_b = Multivector::<f64, Euclidean2>::from(b);
        let gen = &mv_a * &mv_b;
        let gen_rotor = Rotor::from_multivector_unchecked(&gen);

        prop_assert!(abs_diff_eq!(spec, gen_rotor, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Generated sandwich product matches explicit R * v * R̃.
    #[test]
    fn rotor_sandwich_matches_explicit(
        r in any::<Rotor<f64>>(),
        v in any::<Vector<f64>>(),
    ) {
        let r = r.normalized();  // Use normalized() instead of wrapper type
        let sandwich = r.rotate(v);

        let mv_r = Multivector::<f64, Euclidean2>::from(r);
        let mv_v = Multivector::<f64, Euclidean2>::from(v);
        let explicit = &(&mv_r * &mv_v) * &mv_r.reverse();
        let explicit_vec = Vector::from_multivector_unchecked(&explicit);

        prop_assert!(abs_diff_eq!(sandwich, explicit_vec, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Rotation preserves vector norm.
    #[test]
    fn rotation_preserves_norm(
        r in any::<Rotor<f64>>(),
        v in any::<Vector<f64>>(),
    ) {
        let r = r.normalized();  // Use normalized() instead of wrapper type
        let rotated = r.rotate(v);
        prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));
    }

    /// from_angle round-trip.
    #[test]
    fn from_angle_roundtrip(angle in -std::f64::consts::PI..std::f64::consts::PI) {
        let r = Rotor::from_angle(angle);
        let extracted = r.angle();
        prop_assert!(abs_diff_eq!(angle, extracted, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// perp is perpendicular and same length.
    #[test]
    fn perp_properties(v in any::<Vector<f64>>()) {
        let perp = v.perp();
        prop_assert!(abs_diff_eq!(v.dot(&perp), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(v.norm(), perp.norm(), epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Verification Checklist

- [ ] TOML specification updated
- [ ] Code generation succeeds
- [ ] Extensions module complete
- [ ] Operators updated
- [ ] All existing tests pass
- [ ] Consistency tests pass
- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo test --all-features` passes

## Files Changed

| File | Action |
|------|--------|
| `algebras/euclidean2.toml` | Update |
| `src/specialized/euclidean/dim2/mod.rs` | Update |
| `src/specialized/euclidean/dim2/generated/` | Create |
| `src/specialized/euclidean/dim2/extensions.rs` | Create |
| `src/specialized/euclidean/dim2/nalgebra.rs` | Update (new constructor API) |
| `src/specialized/euclidean/dim2/rerun.rs` | Update (new Arbitrary API) |
| `src/specialized/euclidean/dim2/types.rs` | Delete (replaced by generated) |
| `src/specialized/euclidean/dim2/ops.rs` | Delete (replaced by generated traits) |
| `src/specialized/euclidean/dim2/arbitrary.rs` | Delete (generated handles this) |
| `src/specialized/euclidean/dim2/conversions.rs` | Delete (generated handles this) |
