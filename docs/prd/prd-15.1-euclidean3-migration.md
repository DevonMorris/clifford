# PRD-15.1: Euclidean 3D Codegen Migration

**Status**: Proposed
**Parent**: PRD-15 (Codegen Migration)
**Goal**: Migrate `specialized/euclidean/dim3` from hand-rolled to generated code

## Overview

Euclidean 3D is the most complex of the euclidean algebras with 5 types (Vector, Bivector, Trivector, Rotor, Even). We migrate this first to validate the pattern and catch issues early.

### Current Types

| Type | Grades | Fields | Key Methods |
|------|--------|--------|-------------|
| `Vector` | [1] | x, y, z | `dot()`, `wedge()`, `cross()`, `geometric()`, `perp()` |
| `Bivector` | [2] | xy, xz, yz | `reverse()`, `dual()` |
| `Trivector` | [3] | xyz | `reverse()` |
| `Rotor` | [0,2] | s, xy, xz, yz | `from_angle_plane()`, `from_vectors()`, `rotate()`, `compose()`, `inverse()`, `slerp()` |
| `Even` | [0,2] | s, xy, xz, yz | `to_rotor()` (alias of Rotor) |

## Phase 1: Update TOML Specification

### Deliverable: `algebras/euclidean3.toml`

Verify and update the existing specification:

```toml
# 3D Euclidean Geometric Algebra
# Signature: (3, 0, 0)
# Basis: e1, e2, e3
# Blades: 1, e1, e2, e3, e12, e13, e23, e123

[algebra]
name = "euclidean3"
module_path = "euclidean::dim3"
description = "3D Euclidean Geometric Algebra"

[signature]
positive = ["e1", "e2", "e3"]

[blades]
# Grade 1
e1 = "x"
e2 = "y"
e3 = "z"
# Grade 2
e12 = "xy"
e13 = "xz"
e23 = "yz"
# Grade 3
e123 = "xyz"

[types.Vector]
grades = [1]
description = "3D vector (grade-1 element)"
fields = ["x", "y", "z"]

[types.Bivector]
grades = [2]
description = "3D bivector (oriented plane, grade-2 element)"
fields = ["xy", "xz", "yz"]

[types.Trivector]
grades = [3]
description = "3D trivector (pseudoscalar, grade-3 element)"
fields = ["xyz"]

[types.Rotor]
grades = [0, 2]
description = "3D rotor (rotation element in even subalgebra)"
fields = ["s", "xy", "xz", "yz"]

[types.Even]
grades = [0, 2]
description = "Even subalgebra element"
fields = ["s", "xy", "xz", "yz"]
alias_of = "Rotor"

[products.geometric]
Vector_Vector = "Rotor"
Rotor_Rotor = "Rotor"
Bivector_Vector = "Odd"
Vector_Bivector = "Odd"

[products.outer]
Vector_Vector = "Bivector"
Vector_Bivector = "Trivector"
Bivector_Vector = "Trivector"

[products.inner]
Vector_Vector = "T"
Bivector_Bivector = "T"
Vector_Bivector = "Vector"
Bivector_Vector = "Vector"

[products.scalar]
Vector_Vector = "T"
Rotor_Rotor = "T"
Bivector_Bivector = "T"

[products.sandwich]
Rotor_Vector = "Vector"
Rotor_Bivector = "Bivector"

[options]
generate_serde = true
generate_arbitrary = true
generate_tests = true
```

### Verification

- [ ] Run `cargo run --package clifford-codegen -- blades algebras/euclidean3.toml` to verify blade mappings
- [ ] Run `cargo run --package clifford-codegen -- products algebras/euclidean3.toml --table` to verify products

## Phase 2: Generate Code

### Command

```bash
cargo run --package clifford-codegen -- generate algebras/euclidean3.toml \
    -o src/specialized/euclidean/dim3/generated/ --force
```

### Expected Output Files

```
src/specialized/euclidean/dim3/generated/
├── mod.rs          # Module exports
├── types.rs        # Vector, Bivector, Trivector, Rotor, Even structs
├── products.rs     # Product functions
└── traits.rs       # Clone, Copy, Debug, PartialEq, approx traits
```

### Generated Types Structure

```rust
// types.rs (generated)

/// 3D vector (grade-1 element)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Vector<T: Float> {
    x: T,
    y: T,
    z: T,
}

impl<T: Float> Vector<T> {
    /// Creates a new vector.
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    /// Returns the x component.
    #[inline]
    pub fn x(&self) -> T { self.x }

    /// Returns the y component.
    #[inline]
    pub fn y(&self) -> T { self.y }

    /// Returns the z component.
    #[inline]
    pub fn z(&self) -> T { self.z }

    /// Zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Unit vector along x-axis.
    #[inline]
    pub fn unit_x() -> Self {
        Self::new(T::one(), T::zero(), T::zero())
    }

    /// Unit vector along y-axis.
    #[inline]
    pub fn unit_y() -> Self {
        Self::new(T::zero(), T::one(), T::zero())
    }

    /// Unit vector along z-axis.
    #[inline]
    pub fn unit_z() -> Self {
        Self::new(T::zero(), T::zero(), T::one())
    }

    /// Squared norm.
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Norm (magnitude).
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Returns normalized vector.
    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.x / n, self.y / n, self.z / n)
    }

    /// Scale by scalar.
    #[inline]
    pub fn scale(self, s: T) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}
```

### Generated Products

```rust
// products.rs (generated)

use super::types::{Vector, Bivector, Trivector, Rotor};
use crate::scalar::Float;

/// Geometric product of two vectors.
#[inline]
pub fn geometric_vector_vector<T: Float>(a: &Vector<T>, b: &Vector<T>) -> Rotor<T> {
    // s = a.x*b.x + a.y*b.y + a.z*b.z
    // xy = a.x*b.y - a.y*b.x
    // xz = a.x*b.z - a.z*b.x
    // yz = a.y*b.z - a.z*b.y
    Rotor::new(
        a.x() * b.x() + a.y() * b.y() + a.z() * b.z(),
        a.x() * b.y() - a.y() * b.x(),
        a.x() * b.z() - a.z() * b.x(),
        a.y() * b.z() - a.z() * b.y(),
    )
}

/// Outer product of two vectors.
#[inline]
pub fn outer_vector_vector<T: Float>(a: &Vector<T>, b: &Vector<T>) -> Bivector<T> {
    Bivector::new(
        a.x() * b.y() - a.y() * b.x(),
        a.x() * b.z() - a.z() * b.x(),
        a.y() * b.z() - a.z() * b.y(),
    )
}

/// Scalar (inner) product of two vectors.
#[inline]
pub fn scalar_vector_vector<T: Float>(a: &Vector<T>, b: &Vector<T>) -> T {
    a.x() * b.x() + a.y() * b.y() + a.z() * b.z()
}

/// Geometric product of two rotors.
#[inline]
pub fn geometric_rotor_rotor<T: Float>(a: &Rotor<T>, b: &Rotor<T>) -> Rotor<T> {
    // Generated from symbolic computation
    Rotor::new(
        a.s() * b.s() - a.xy() * b.xy() - a.xz() * b.xz() - a.yz() * b.yz(),
        a.s() * b.xy() + a.xy() * b.s() - a.xz() * b.yz() + a.yz() * b.xz(),
        a.s() * b.xz() + a.xy() * b.yz() + a.xz() * b.s() - a.yz() * b.xy(),
        a.s() * b.yz() - a.xy() * b.xz() + a.xz() * b.xy() + a.yz() * b.s(),
    )
}

/// Sandwich product: rotor * vector * reverse(rotor)
#[inline]
pub fn sandwich_rotor_vector<T: Float>(r: &Rotor<T>, v: &Vector<T>) -> Vector<T> {
    // Generated from symbolic computation
    let s = r.s();
    let xy = r.xy();
    let xz = r.xz();
    let yz = r.yz();
    let x = v.x();
    let y = v.y();
    let z = v.z();

    // Full sandwich expansion R * v * R̃
    Vector::new(
        x * (s*s + xy*xy - xz*xz - yz*yz)
            + y * T::TWO * (xy*s + xz*yz)
            + z * T::TWO * (xz*s - xy*yz),
        x * T::TWO * (xy*s - xz*yz)
            + y * (s*s - xy*xy + xz*xz - yz*yz)
            + z * T::TWO * (yz*s + xy*xz),
        x * T::TWO * (xz*s + xy*yz)
            + y * T::TWO * (yz*s - xy*xz)
            + z * (s*s - xy*xy - xz*xz + yz*yz),
    )
}
```

### Verification

- [ ] Generated files compile without errors
- [ ] Generated types have correct field count
- [ ] Generated products have correct signatures

## Phase 3: Create Extensions Module

### File: `src/specialized/euclidean/dim3/extensions.rs`

```rust
//! Domain-specific extensions for 3D Euclidean GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types.

use super::generated::products;
use super::generated::types::{Bivector, Rotor, Trivector, Vector};
use crate::scalar::Float;

// ============================================================================
// Vector extensions
// ============================================================================

impl<T: Float> Vector<T> {
    /// Dot product (inner product).
    ///
    /// Returns the scalar part of the geometric product `a * b`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let a = Vector::new(1.0, 2.0, 3.0);
    /// let b = Vector::new(4.0, 5.0, 6.0);
    /// assert_eq!(a.dot(&b), 32.0); // 1*4 + 2*5 + 3*6
    /// ```
    #[inline]
    pub fn dot(&self, other: &Self) -> T {
        products::scalar_vector_vector(self, other)
    }

    /// Wedge (outer) product.
    ///
    /// Returns the bivector representing the oriented plane spanned by
    /// the two vectors.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Vector, Bivector};
    ///
    /// let a = Vector::unit_x();
    /// let b = Vector::unit_y();
    /// let ab = a.wedge(&b);
    /// assert_eq!(ab, Bivector::unit_xy());
    /// ```
    #[inline]
    pub fn wedge(&self, other: &Self) -> Bivector<T> {
        products::outer_vector_vector(self, other)
    }

    /// Cross product.
    ///
    /// Computed as the Hodge dual of the wedge product.
    /// Returns a vector perpendicular to both inputs.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let x = Vector::unit_x();
    /// let y = Vector::unit_y();
    /// let z = x.cross(&y);
    /// assert_eq!(z, Vector::unit_z());
    /// ```
    #[inline]
    pub fn cross(&self, other: &Self) -> Self {
        self.wedge(other).dual()
    }

    /// Geometric product with another vector.
    ///
    /// Returns a rotor (scalar + bivector) representing `a * b = a·b + a∧b`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let a = Vector::new(1.0, 0.0, 0.0);
    /// let b = Vector::new(0.0, 1.0, 0.0);
    /// let ab = a.geometric(&b);
    /// assert_eq!(ab.s(), 0.0);  // perpendicular, no dot product
    /// ```
    #[inline]
    pub fn geometric(&self, other: &Self) -> Rotor<T> {
        products::geometric_vector_vector(self, other)
    }

    /// Returns a perpendicular vector in the xy-plane.
    ///
    /// Rotates by 90 degrees counterclockwise when viewed from +z.
    /// Only uses x and y components.
    #[inline]
    pub fn perp(&self) -> Self {
        Self::new(-self.y(), self.x(), T::zero())
    }
}

// ============================================================================
// Bivector extensions
// ============================================================================

impl<T: Float> Bivector<T> {
    /// Reverse (reversion).
    ///
    /// For a bivector, reverse negates the sign: `B̃ = -B`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(-self.xy(), -self.xz(), -self.yz())
    }

    /// Hodge dual.
    ///
    /// Maps bivector to vector: `*(e₁₂) = e₃`, `*(e₁₃) = -e₂`, `*(e₂₃) = e₁`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Vector, Bivector};
    ///
    /// let b = Bivector::unit_xy();
    /// assert_eq!(b.dual(), Vector::unit_z());
    /// ```
    #[inline]
    pub fn dual(&self) -> Vector<T> {
        // In 3D: *(e12) = e3, *(e13) = -e2, *(e23) = e1
        Vector::new(self.yz(), -self.xz(), self.xy())
    }
}

// ============================================================================
// Trivector extensions
// ============================================================================

impl<T: Float> Trivector<T> {
    /// Reverse (reversion).
    ///
    /// For a trivector (grade 3), reverse negates: `T̃ = -T`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(-self.xyz())
    }
}

// ============================================================================
// Rotor extensions
// ============================================================================

impl<T: Float> Rotor<T> {
    /// Identity rotor (no rotation).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Vector};
    ///
    /// let r = Rotor::identity();
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// assert_eq!(r.rotate(v), v);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self::new(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Create from rotation angle and plane (unit bivector).
    ///
    /// The rotation is by `angle` radians in the plane defined by the bivector,
    /// following the right-hand rule.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Bivector, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // 90° rotation in xy-plane
    /// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// // x rotates to y
    /// ```
    pub fn from_angle_plane(angle: T, plane: Bivector<T>) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let plane_norm = plane.normalized();
        Self::new(
            cos_half,
            -sin_half * plane_norm.xy(),
            -sin_half * plane_norm.xz(),
            -sin_half * plane_norm.yz(),
        )
    }

    /// Create rotation that takes vector `a` to vector `b`.
    ///
    /// Both vectors should be non-zero. The rotation is the shortest
    /// arc from `a` to `b`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Vector};
    /// use approx::abs_diff_eq;
    ///
    /// let a = Vector::unit_x();
    /// let b = Vector::unit_y();
    /// let r = Rotor::from_vectors(a, b);
    /// let rotated = r.rotate(a);
    /// assert!(abs_diff_eq!(rotated.x(), b.x(), epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.y(), b.y(), epsilon = 1e-10));
    /// ```
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        // R = (1 + ba) / |1 + ba|
        // This gives the rotor that rotates a to b
        let a_norm = a.normalized();
        let b_norm = b.normalized();
        let ba = b_norm.geometric(&a_norm);

        // Add 1 to scalar part
        let unnorm = Self::new(
            T::one() + ba.s(),
            ba.xy(),
            ba.xz(),
            ba.yz(),
        );

        // Handle anti-parallel case (180° rotation)
        if unnorm.norm_squared() < T::epsilon() {
            // Find perpendicular axis
            let perp = if a_norm.x().abs() < T::from_f64(0.9) {
                a_norm.cross(&Vector::unit_x())
            } else {
                a_norm.cross(&Vector::unit_y())
            };
            let plane = a_norm.wedge(&perp.normalized());
            return Self::from_angle_plane(T::PI, plane);
        }

        unnorm.normalized()
    }

    /// Reverse (conjugate) of rotor.
    ///
    /// For a unit rotor, `R̃ = R⁻¹`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(self.s(), -self.xy(), -self.xz(), -self.yz())
    }

    /// Inverse rotor.
    ///
    /// For any rotor: `R⁻¹ = R̃ / |R|²`.
    /// For unit rotors: `R⁻¹ = R̃`.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        let rev = self.reverse();
        Self::new(
            rev.s() / norm_sq,
            rev.xy() / norm_sq,
            rev.xz() / norm_sq,
            rev.yz() / norm_sq,
        )
    }

    /// Apply rotation to vector.
    ///
    /// Computes the sandwich product `R v R̃`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Bivector, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn rotate(&self, v: Vector<T>) -> Vector<T> {
        products::sandwich_rotor_vector(self, &v)
    }

    /// Compose two rotations.
    ///
    /// Returns a rotor that applies `self` first, then `other`.
    /// Mathematically: `other * self`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Bivector, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let r1 = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let r2 = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let r3 = r1.compose(&r2); // 180° rotation
    /// ```
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        products::geometric_rotor_rotor(other, self)
    }

    /// Linear interpolation between rotors.
    ///
    /// Note: This does not produce constant angular velocity.
    /// Use `slerp` for geodesic interpolation.
    #[inline]
    pub fn lerp(&self, other: &Self, t: T) -> Self {
        let one_minus_t = T::one() - t;
        Self::new(
            self.s() * one_minus_t + other.s() * t,
            self.xy() * one_minus_t + other.xy() * t,
            self.xz() * one_minus_t + other.xz() * t,
            self.yz() * one_minus_t + other.yz() * t,
        )
    }

    /// Spherical linear interpolation.
    ///
    /// Interpolates along the geodesic (great arc) between two rotors.
    /// Produces constant angular velocity.
    ///
    /// # Arguments
    ///
    /// * `other` - Target rotor
    /// * `t` - Interpolation parameter in [0, 1]
    pub fn slerp(&self, other: &Self, t: T) -> Self {
        // Compute relative rotation: other * self.reverse()
        let rel = other.compose(&self.inverse());

        // Extract angle from relative rotation
        let cos_half = rel.s().min(T::one()).max(-T::one());
        let half_angle = cos_half.acos();
        let angle = half_angle * T::TWO;

        if angle.abs() < T::epsilon() {
            return *self;
        }

        // Interpolate angle
        let t_angle = angle * t;

        // Reconstruct interpolated rotation
        let sin_half = (T::one() - rel.s() * rel.s()).sqrt();
        if sin_half.abs() < T::epsilon() {
            return self.lerp(other, t).normalized();
        }

        let axis_scale = T::one() / sin_half;
        let plane = Bivector::new(
            rel.xy() * axis_scale,
            rel.xz() * axis_scale,
            rel.yz() * axis_scale,
        );

        // Apply interpolated rotation to self
        let interp_rel = Self::from_angle_plane(t_angle, plane);
        interp_rel.compose(self)
    }

    /// Extract rotation angle in radians.
    ///
    /// Returns the angle in [0, 2π].
    #[inline]
    pub fn angle(&self) -> T {
        let cos_half = self.s().min(T::one()).max(-T::one());
        cos_half.acos() * T::TWO
    }
}

// ============================================================================
// Even extensions
// ============================================================================

impl<T: Float> super::generated::types::Even<T> {
    /// Convert to rotor (Even is alias of Rotor).
    #[inline]
    pub fn to_rotor(&self) -> Rotor<T> {
        Rotor::new(self.s(), self.xy(), self.xz(), self.yz())
    }
}
```

### Verification

- [ ] All methods from current `types.rs` are covered
- [ ] Documentation matches current style
- [ ] Examples compile as doctests

## Phase 4: Update Module Structure

### File: `src/specialized/euclidean/dim3/mod.rs`

```rust
//! 3D Euclidean Geometric Algebra.
//!
//! This module provides optimized types for 3D Euclidean GA (Cl(3,0,0)):
//!
//! - [`Vector`] - Grade-1 elements (directions, positions)
//! - [`Bivector`] - Grade-2 elements (oriented planes, rotation generators)
//! - [`Trivector`] - Grade-3 elements (pseudoscalar, oriented volumes)
//! - [`Rotor`] - Even subalgebra elements (rotations)
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim3::{Vector, Rotor, Bivector};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // Create vectors
//! let v = Vector::new(1.0, 0.0, 0.0);
//!
//! // Create a 90° rotation in the xy-plane
//! let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
//!
//! // Apply rotation
//! let rotated = r.rotate(v);
//! ```

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions
mod extensions;

// Operator overloads
mod ops;

// Multivector conversions
mod conversions;

// Re-export generated types
pub use generated::types::{Bivector, Even, Rotor, Trivector, Vector};

// Re-export conversion utilities
pub use conversions::CONVERSION_TOLERANCE;

// Proptest support
#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;

// nalgebra integration
#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra;

// Rerun visualization
#[cfg(feature = "rerun")]
mod rerun;

#[cfg(test)]
mod tests;
```

### File: `src/specialized/euclidean/dim3/generated/mod.rs`

```rust
//! Generated code for 3D Euclidean GA.
//!
//! **Do not edit manually.** This code is generated by `clifford-codegen`.
//! To regenerate:
//!
//! ```bash
//! cargo run --package clifford-codegen -- generate algebras/euclidean3.toml \
//!     -o src/specialized/euclidean/dim3/generated/ --force
//! ```

pub mod products;
pub mod traits;
pub mod types;
```

## Phase 5: Update Operator Overloads

### File: `src/specialized/euclidean/dim3/ops.rs`

Update to use generated products:

```rust
//! Operator overloads for 3D Euclidean GA types.

use super::generated::products;
use super::generated::types::{Bivector, Rotor, Trivector, Vector};
use crate::scalar::Float;
use std::ops::{Add, BitXor, Mul, Neg, Sub};

// ============================================================================
// Vector operators
// ============================================================================

impl<T: Float> Neg for Vector<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x(), -self.y(), -self.z())
    }
}

impl<T: Float> Add for Vector<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x() + rhs.x(), self.y() + rhs.y(), self.z() + rhs.z())
    }
}

impl<T: Float> Sub for Vector<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x() - rhs.x(), self.y() - rhs.y(), self.z() - rhs.z())
    }
}

impl<T: Float> Mul<T> for Vector<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        self.scale(rhs)
    }
}

/// Wedge product via `^` operator.
impl<T: Float> BitXor<Vector<T>> for Vector<T> {
    type Output = Bivector<T>;

    #[inline]
    fn bitxor(self, rhs: Vector<T>) -> Bivector<T> {
        products::outer_vector_vector(&self, &rhs)
    }
}

// Scalar * Vector (for f64)
impl Mul<Vector<f64>> for f64 {
    type Output = Vector<f64>;

    #[inline]
    fn mul(self, rhs: Vector<f64>) -> Vector<f64> {
        rhs.scale(self)
    }
}

// Scalar * Vector (for f32)
impl Mul<Vector<f32>> for f32 {
    type Output = Vector<f32>;

    #[inline]
    fn mul(self, rhs: Vector<f32>) -> Vector<f32> {
        rhs.scale(self)
    }
}

// ============================================================================
// Bivector operators
// ============================================================================

impl<T: Float> Neg for Bivector<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.xy(), -self.xz(), -self.yz())
    }
}

impl<T: Float> Add for Bivector<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(
            self.xy() + rhs.xy(),
            self.xz() + rhs.xz(),
            self.yz() + rhs.yz(),
        )
    }
}

impl<T: Float> Sub for Bivector<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(
            self.xy() - rhs.xy(),
            self.xz() - rhs.xz(),
            self.yz() - rhs.yz(),
        )
    }
}

impl<T: Float> Mul<T> for Bivector<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        self.scale(rhs)
    }
}

// ============================================================================
// Trivector operators
// ============================================================================

impl<T: Float> Neg for Trivector<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.xyz())
    }
}

impl<T: Float> Add for Trivector<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.xyz() + rhs.xyz())
    }
}

impl<T: Float> Sub for Trivector<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.xyz() - rhs.xyz())
    }
}

impl<T: Float> Mul<T> for Trivector<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        self.scale(rhs)
    }
}

// ============================================================================
// Rotor operators
// ============================================================================

impl<T: Float> Neg for Rotor<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.s(), -self.xy(), -self.xz(), -self.yz())
    }
}

impl<T: Float> Add for Rotor<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(
            self.s() + rhs.s(),
            self.xy() + rhs.xy(),
            self.xz() + rhs.xz(),
            self.yz() + rhs.yz(),
        )
    }
}

impl<T: Float> Sub for Rotor<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(
            self.s() - rhs.s(),
            self.xy() - rhs.xy(),
            self.xz() - rhs.xz(),
            self.yz() - rhs.yz(),
        )
    }
}

impl<T: Float> Mul<T> for Rotor<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        self.scale(rhs)
    }
}

/// Rotor composition via `*` operator.
impl<T: Float> Mul<Rotor<T>> for Rotor<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        products::geometric_rotor_rotor(&self, &rhs)
    }
}
```

## Phase 6: Update Conversions

### File: `src/specialized/euclidean/dim3/conversions.rs`

Update to work with generated types (keep existing validation logic).

## Phase 7: Tests

### New Test File: `src/specialized/euclidean/dim3/tests.rs`

```rust
//! Tests for 3D Euclidean GA codegen migration.

use super::*;
use crate::multivector::Multivector;
use crate::signature::Euclidean3;
use crate::test_utils::ABS_DIFF_EQ_EPS;
use approx::abs_diff_eq;
use proptest::prelude::*;

#[cfg(any(test, feature = "proptest-support"))]
use super::arbitrary::{NonZeroVector, UnitRotor};

proptest! {
    // ========================================================================
    // Consistency with Multivector
    // ========================================================================

    /// Generated geometric product matches generic Multivector.
    #[test]
    fn vector_geometric_matches_multivector(
        a in any::<Vector<f64>>(),
        b in any::<Vector<f64>>(),
    ) {
        let spec = a.geometric(&b);

        let mv_a = Multivector::<f64, Euclidean3>::from(a);
        let mv_b = Multivector::<f64, Euclidean3>::from(b);
        let gen = &mv_a * &mv_b;
        let gen_rotor = Rotor::from_multivector_unchecked(&gen);

        prop_assert!(abs_diff_eq!(spec, gen_rotor, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Generated outer product matches generic Multivector.
    #[test]
    fn vector_outer_matches_multivector(
        a in any::<Vector<f64>>(),
        b in any::<Vector<f64>>(),
    ) {
        let spec = a.wedge(&b);

        let mv_a = Multivector::<f64, Euclidean3>::from(a);
        let mv_b = Multivector::<f64, Euclidean3>::from(b);
        let gen = mv_a.outer(&mv_b);
        let gen_bivec = Bivector::from_multivector_unchecked(&gen);

        prop_assert!(abs_diff_eq!(spec, gen_bivec, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Generated scalar product matches generic inner product.
    #[test]
    fn vector_dot_matches_multivector(
        a in any::<Vector<f64>>(),
        b in any::<Vector<f64>>(),
    ) {
        let spec = a.dot(&b);

        let mv_a = Multivector::<f64, Euclidean3>::from(a);
        let mv_b = Multivector::<f64, Euclidean3>::from(b);
        let gen = mv_a.inner(&mv_b).scalar_part();

        prop_assert!(abs_diff_eq!(spec, gen, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Generated sandwich product matches explicit R * v * R̃.
    #[test]
    fn rotor_sandwich_matches_explicit(
        r in any::<UnitRotor<f64>>(),
        v in any::<Vector<f64>>(),
    ) {
        let sandwich = r.rotate(v);

        let mv_r = Multivector::<f64, Euclidean3>::from(**r);
        let mv_v = Multivector::<f64, Euclidean3>::from(v);
        let explicit = &(&mv_r * &mv_v) * &mv_r.reverse();
        let explicit_vec = Vector::from_multivector_unchecked(&explicit);

        prop_assert!(abs_diff_eq!(sandwich, explicit_vec, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Generated rotor composition matches generic product.
    #[test]
    fn rotor_compose_matches_multivector(
        r1 in any::<Rotor<f64>>(),
        r2 in any::<Rotor<f64>>(),
    ) {
        let spec = r1.compose(&r2);

        let mv_r1 = Multivector::<f64, Euclidean3>::from(r1);
        let mv_r2 = Multivector::<f64, Euclidean3>::from(r2);
        let gen = &mv_r2 * &mv_r1;
        let gen_rotor = Rotor::from_multivector_unchecked(&gen);

        prop_assert!(abs_diff_eq!(spec, gen_rotor, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Algebraic properties
    // ========================================================================

    /// Rotation preserves vector norm.
    #[test]
    fn rotation_preserves_norm(
        r in any::<UnitRotor<f64>>(),
        v in any::<Vector<f64>>(),
    ) {
        let rotated = r.rotate(v);
        prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Rotor composition is associative.
    #[test]
    fn rotor_compose_associative(
        r1 in any::<Rotor<f64>>(),
        r2 in any::<Rotor<f64>>(),
        r3 in any::<Rotor<f64>>(),
    ) {
        let lhs = r1.compose(&r2).compose(&r3);
        let rhs = r1.compose(&r2.compose(&r3));
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Inverse rotor undoes rotation.
    #[test]
    fn rotor_inverse(
        r in any::<UnitRotor<f64>>(),
        v in any::<Vector<f64>>(),
    ) {
        let rotated = r.rotate(v);
        let back = r.inverse().rotate(rotated);
        prop_assert!(abs_diff_eq!(v, back, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// from_vectors produces rotor that maps a to b.
    #[test]
    fn from_vectors_correctness(
        a in any::<NonZeroVector<f64>>(),
        b in any::<NonZeroVector<f64>>(),
    ) {
        let r = Rotor::from_vectors(**a, **b);
        let rotated = r.rotate(a.normalized());
        let expected = b.normalized();
        prop_assert!(abs_diff_eq!(rotated, expected, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Conversion round-trips
    // ========================================================================

    /// Vector round-trip through Multivector.
    #[test]
    fn vector_roundtrip(v in any::<Vector<f64>>()) {
        let mv = Multivector::<f64, Euclidean3>::from(v);
        let back = Vector::try_from(mv).unwrap();
        prop_assert!(abs_diff_eq!(v, back, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Bivector round-trip through Multivector.
    #[test]
    fn bivector_roundtrip(b in any::<Bivector<f64>>()) {
        let mv = Multivector::<f64, Euclidean3>::from(b);
        let back = Bivector::try_from(mv).unwrap();
        prop_assert!(abs_diff_eq!(b, back, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Rotor round-trip through Multivector.
    #[test]
    fn rotor_roundtrip(r in any::<Rotor<f64>>()) {
        let mv = Multivector::<f64, Euclidean3>::from(r);
        let back = Rotor::try_from(mv).unwrap();
        prop_assert!(abs_diff_eq!(r, back, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Verification Checklist

### Phase 1: TOML Specification
- [ ] TOML parses without errors
- [ ] Blade mappings verified
- [ ] Product outputs verified

### Phase 2: Code Generation
- [ ] `cargo run --package clifford-codegen -- generate ...` succeeds
- [ ] Generated files compile
- [ ] No warnings in generated code

### Phase 3: Extensions
- [ ] All current methods implemented
- [ ] Documentation complete
- [ ] Doctests pass

### Phase 4: Module Structure
- [ ] Re-exports work correctly
- [ ] Feature flags respected

### Phase 5: Operators
- [ ] All operators implemented
- [ ] Operators use generated products

### Phase 6: Conversions
- [ ] All conversions work
- [ ] Validation logic preserved

### Phase 7: Tests
- [ ] All existing tests pass
- [ ] New consistency tests pass
- [ ] Benchmarks show no regression

### Final Verification
- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] Public API unchanged (semver compatible)

## Files Changed

| File | Action |
|------|--------|
| `algebras/euclidean3.toml` | Update |
| `src/specialized/euclidean/dim3/mod.rs` | Update |
| `src/specialized/euclidean/dim3/generated/` | Create (generated) |
| `src/specialized/euclidean/dim3/extensions.rs` | Create |
| `src/specialized/euclidean/dim3/ops.rs` | Update |
| `src/specialized/euclidean/dim3/conversions.rs` | Update |
| `src/specialized/euclidean/dim3/types.rs` | Delete (moved to generated) |
| `src/specialized/euclidean/dim3/tests.rs` | Create |
