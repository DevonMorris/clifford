# PRD-15.4: Projective 2D Codegen Migration

**Status**: Proposed
**Parent**: PRD-15 (Codegen Migration)
**Depends On**: PRD-15.1, PRD-15.2, PRD-15.3 (pattern validation)
**Goal**: Migrate `specialized/projective/dim2` from hand-rolled to generated code

## Overview

Projective 2D is simpler than 3D with 3 types (Point, Line, Motor). This migration validates the PGA pattern established in PRD-15.3.

### Current Types

| Type | Grades | Fields | Constraints | Key Methods |
|------|--------|--------|-------------|-------------|
| `Point` | [1] | e1, e2, e0 | None | `origin()`, `ideal()`, `is_ideal()`, `to_cartesian()`, `join()`, `distance()` |
| `Line` | [2] | e12, e20, e01 | None | `from_implicit()`, `x_axis()`, `y_axis()`, `normal()`, `meet()`, `distance_to_point()`, `project()`, `reflect()` |
| `Motor` | [0,2] | s, e12, e20, e01 | Study (trivial) | `from_rotation()`, `from_translation()`, `from_rotation_around()`, `transform_point()`, `transform_line()`, `compose()` |

### Key Differences from 3D

- No `Plane` type (line is the boundary element)
- No `Flector` type (simpler reflections)
- Motor has only 4 components (vs 8 in 3D)
- Study condition is trivial: `s·e₀₁₂ = 0` (e012 is always 0)
- No Plücker condition (2D lines are always valid)

## Phase 1: Create TOML Specification

### Deliverable: `algebras/projective2.toml`

```toml
# 2D Projective Geometric Algebra (Point-based)
# Signature: Cl(2,0,1)
# Basis: e1, e2, e0 where e0² = 0

[algebra]
name = "projective2"
module_path = "projective::dim2"
description = "2D Projective Geometric Algebra (Point-based PGA)"

[signature]
positive = ["e1", "e2"]
zero = ["e0"]

[blades]
# Grade 1
e1 = "e1"
e2 = "e2"
e0 = "e0"
# Grade 2
e12 = "e12"
e20 = "e20"
e01 = "e01"
# Grade 3 (pseudoscalar)
e012 = "e012"

[types.Point]
grades = [1]
description = "2D point in homogeneous coordinates"
fields = ["e1", "e2", "e0"]

[types.Line]
grades = [2]
description = "2D line (ax + by + c = 0)"
fields = ["e12", "e20", "e01"]

[types.Motor]
grades = [0, 2]
description = "2D rigid transformation (rotation + translation)"
fields = ["s", "e12", "e20", "e01"]

# In 2D, motors always satisfy the Study condition trivially
# (no e012 component in the motor)

[products.geometric]
Motor_Motor = "Motor"

[products.regressive]
Point_Point = "Line"

[products.outer]
Line_Line = "Point"

[products.sandwich]
Motor_Point = "Point"
Motor_Line = "Line"

[options]
generate_serde = true
generate_arbitrary = true
generate_tests = true
```

## Phase 2: Generate Code

### Command

```bash
cargo run --package clifford-codegen -- generate algebras/projective2.toml \
    -o src/specialized/projective/dim2/generated/ --force
```

### Expected Output

```
src/specialized/projective/dim2/generated/
├── mod.rs
├── types.rs      # Point, Line, Motor
├── products.rs   # Product functions
└── traits.rs     # Trait implementations
```

### Generated Types

```rust
// types.rs (generated)

/// 2D point in homogeneous coordinates.
///
/// A point is represented as `P = x·e₁ + y·e₂ + w·e₀` where
/// `(x/w, y/w)` are the Cartesian coordinates for `w ≠ 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Point<T: Float> {
    e1: T,
    e2: T,
    e0: T,
}

impl<T: Float> Point<T> {
    pub fn new(e1: T, e2: T, e0: T) -> Self { Self { e1, e2, e0 } }
    pub fn e1(&self) -> T { self.e1 }
    pub fn e2(&self) -> T { self.e2 }
    pub fn e0(&self) -> T { self.e0 }

    pub fn zero() -> Self { Self::new(T::zero(), T::zero(), T::zero()) }

    pub fn bulk_norm_squared(&self) -> T {
        self.e1 * self.e1 + self.e2 * self.e2
    }

    pub fn bulk_norm(&self) -> T {
        self.bulk_norm_squared().sqrt()
    }

    pub fn weight_norm(&self) -> T {
        self.e0.abs()
    }
}

/// 2D line in implicit form.
///
/// A line is represented as `L = d·e₁₂ + a·e₂₀ + b·e₀₁`
/// corresponding to the equation `ax + by + d = 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Line<T: Float> {
    e12: T,  // distance from origin (scaled)
    e20: T,  // x-component of normal
    e01: T,  // y-component of normal
}

impl<T: Float> Line<T> {
    pub fn new(e12: T, e20: T, e01: T) -> Self { Self { e12, e20, e01 } }
    pub fn e12(&self) -> T { self.e12 }
    pub fn e20(&self) -> T { self.e20 }
    pub fn e01(&self) -> T { self.e01 }

    pub fn zero() -> Self { Self::new(T::zero(), T::zero(), T::zero()) }

    pub fn weight_norm_squared(&self) -> T {
        self.e20 * self.e20 + self.e01 * self.e01
    }

    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    pub fn bulk_norm(&self) -> T {
        self.e12.abs()
    }
}

/// 2D motor (rotation + translation).
///
/// A motor is `M = s + d·e₁₂ + tx·e₂₀ + ty·e₀₁` where:
/// - `s = cos(θ/2)` for rotation angle θ
/// - `d = sin(θ/2)` for rotation
/// - `(tx, ty)` encode translation
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Motor<T: Float> {
    s: T,
    e12: T,
    e20: T,
    e01: T,
}

impl<T: Float> Motor<T> {
    pub fn new(s: T, e12: T, e20: T, e01: T) -> Self {
        Self { s, e12, e20, e01 }
    }

    /// Create motor without validation (alias for `new` in 2D).
    pub fn new_unchecked(s: T, e12: T, e20: T, e01: T) -> Self {
        Self::new(s, e12, e20, e01)
    }

    pub fn s(&self) -> T { self.s }
    pub fn e12(&self) -> T { self.e12 }
    pub fn e20(&self) -> T { self.e20 }
    pub fn e01(&self) -> T { self.e01 }

    pub fn weight_norm_squared(&self) -> T {
        self.s * self.s + self.e12 * self.e12
    }

    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    pub fn normalized(&self) -> Self {
        let wn = self.weight_norm();
        Self::new(
            self.s / wn,
            self.e12 / wn,
            self.e20 / wn,
            self.e01 / wn,
        )
    }
}
```

### Generated Products

```rust
// products.rs (generated)

/// Geometric product of two 2D motors.
#[inline]
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    Motor::new(
        a.s() * b.s() - a.e12() * b.e12(),
        a.s() * b.e12() + a.e12() * b.s(),
        a.s() * b.e20() + a.e12() * b.e01() + a.e20() * b.s() - a.e01() * b.e12(),
        a.s() * b.e01() - a.e12() * b.e20() + a.e20() * b.e12() + a.e01() * b.s(),
    )
}

/// Regressive product: Point ∨ Point = Line
#[inline]
pub fn regressive_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    Line::new(
        a.e1() * b.e2() - a.e2() * b.e1(),
        a.e2() * b.e0() - a.e0() * b.e2(),
        a.e0() * b.e1() - a.e1() * b.e0(),
    )
}

/// Outer product: Line ∧ Line = Point (meet)
#[inline]
pub fn outer_line_line<T: Float>(a: &Line<T>, b: &Line<T>) -> Point<T> {
    Point::new(
        a.e01() * b.e12() - a.e12() * b.e01(),
        a.e12() * b.e20() - a.e20() * b.e12(),
        a.e20() * b.e01() - a.e01() * b.e20(),
    )
}

/// Sandwich product: Motor * Point * reverse(Motor)
#[inline]
pub fn sandwich_motor_point<T: Float>(m: &Motor<T>, p: &Point<T>) -> Point<T> {
    let s = m.s();
    let d = m.e12();
    let tx = m.e20();
    let ty = m.e01();

    let px = p.e1();
    let py = p.e2();
    let pw = p.e0();

    // Rotation matrix elements
    let cos_full = s * s - d * d;
    let sin_full = T::TWO * s * d;

    // Translation
    let trans_x = T::TWO * (s * tx + d * ty);
    let trans_y = T::TWO * (s * ty - d * tx);

    Point::new(
        cos_full * px + sin_full * py + trans_x * pw,
        -sin_full * px + cos_full * py + trans_y * pw,
        pw,
    )
}

/// Sandwich product: Motor * Line * reverse(Motor)
#[inline]
pub fn sandwich_motor_line<T: Float>(m: &Motor<T>, l: &Line<T>) -> Line<T> {
    let s = m.s();
    let d = m.e12();
    let tx = m.e20();
    let ty = m.e01();

    let le12 = l.e12();
    let le20 = l.e20();
    let le01 = l.e01();

    // Rotation
    let cos_full = s * s - d * d;
    let sin_full = T::TWO * s * d;

    // Transformed normal
    let new_e20 = cos_full * le20 + sin_full * le01;
    let new_e01 = -sin_full * le20 + cos_full * le01;

    // Transformed distance (includes translation effect)
    let new_e12 = le12 - T::TWO * (tx * new_e20 + ty * new_e01);

    Line::new(new_e12, new_e20, new_e01)
}
```

## Phase 3: Create Extensions Module

### File: `src/specialized/projective/dim2/extensions.rs`

```rust
//! Domain-specific extensions for 2D PGA types.

use super::generated::products;
use super::generated::types::{Line, Motor, Point};
use crate::scalar::Float;

// ============================================================================
// Point extensions
// ============================================================================

impl<T: Float> Point<T> {
    /// Point at the origin.
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::one())
    }

    /// Ideal point (direction at infinity).
    #[inline]
    pub fn ideal(dx: T, dy: T) -> Self {
        Self::new(dx, dy, T::zero())
    }

    /// Create from homogeneous coordinates.
    #[inline]
    pub fn from_homogeneous(e1: T, e2: T, e0: T) -> Self {
        Self::new(e1, e2, e0)
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

    /// Cartesian x-coordinate.
    #[inline]
    pub fn x(&self) -> T {
        self.e1() / self.e0()
    }

    /// Cartesian y-coordinate.
    #[inline]
    pub fn y(&self) -> T {
        self.e2() / self.e0()
    }

    /// Homogeneous weight.
    #[inline]
    pub fn w(&self) -> T {
        self.e0()
    }

    /// Convert to Cartesian coordinates.
    pub fn to_cartesian(&self) -> Option<(T, T)> {
        if self.is_ideal() {
            None
        } else {
            Some((self.x(), self.y()))
        }
    }

    /// Direction of ideal point (attitude).
    pub fn attitude(&self) -> Option<(T, T)> {
        let norm = self.bulk_norm();
        if norm < T::epsilon() {
            None
        } else {
            Some((self.e1() / norm, self.e2() / norm))
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
        dx * dx + dy * dy
    }

    /// Midpoint between two finite points.
    pub fn midpoint(&self, other: &Point<T>) -> Point<T> {
        let two = T::TWO;
        Point::new(
            (self.e1() * other.e0() + other.e1() * self.e0()) / two,
            (self.e2() * other.e0() + other.e2() * self.e0()) / two,
            self.e0() * other.e0(),
        )
    }

    /// Normalize to unit weight.
    pub fn normalize(&self) -> Option<Self> {
        let w = self.e0();
        if w.abs() < T::epsilon() {
            None
        } else {
            Some(Self::new(self.e1() / w, self.e2() / w, T::one()))
        }
    }

    /// Geometric norm (for distance computations).
    pub fn geometric_norm(&self) -> T {
        (self.bulk_norm_squared() / (self.e0() * self.e0())).sqrt()
    }
}

// ============================================================================
// Line extensions
// ============================================================================

impl<T: Float> Line<T> {
    /// Create line from implicit equation ax + by + c = 0.
    #[inline]
    pub fn from_implicit(a: T, b: T, c: T) -> Self {
        Self::new(c, a, b)
    }

    /// X-axis (y = 0).
    #[inline]
    pub fn x_axis() -> Self {
        Self::from_implicit(T::zero(), T::one(), T::zero())
    }

    /// Y-axis (x = 0).
    #[inline]
    pub fn y_axis() -> Self {
        Self::from_implicit(T::one(), T::zero(), T::zero())
    }

    /// Normal vector (a, b) from ax + by + c = 0.
    #[inline]
    pub fn normal(&self) -> (T, T) {
        (self.e20(), self.e01())
    }

    /// Direction vector (perpendicular to normal).
    #[inline]
    pub fn direction(&self) -> (T, T) {
        (-self.e01(), self.e20())
    }

    /// Attitude (normalized normal).
    pub fn attitude(&self) -> Option<(T, T)> {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            None
        } else {
            Some((self.e20() / wn, self.e01() / wn))
        }
    }

    /// Signed distance from origin.
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            T::zero()
        } else {
            self.e12() / wn
        }
    }

    /// Meet of two lines (intersection point).
    #[inline]
    pub fn meet(&self, other: &Line<T>) -> Point<T> {
        products::outer_line_line(self, other)
    }

    /// Signed distance from a point to this line.
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        let (a, b) = self.normal();
        let c = self.e12();
        let wn = self.weight_norm();

        if wn < T::epsilon() || p.is_ideal() {
            return T::zero();
        }

        // Distance = (ax + by + c) / sqrt(a² + b²)
        let (px, py) = (p.x(), p.y());
        (a * px + b * py + c) / wn
    }

    /// Project a point onto this line.
    pub fn project(&self, p: &Point<T>) -> Point<T> {
        let dist = self.distance_to_point(p);
        let (nx, ny) = self.attitude().unwrap_or((T::zero(), T::zero()));

        if p.is_ideal() {
            return *p;
        }

        Point::new(
            p.e1() - dist * nx * p.e0(),
            p.e2() - dist * ny * p.e0(),
            p.e0(),
        )
    }

    /// Reflect a point across this line.
    pub fn reflect(&self, p: &Point<T>) -> Point<T> {
        let dist = self.distance_to_point(p);
        let (nx, ny) = self.attitude().unwrap_or((T::zero(), T::zero()));

        if p.is_ideal() {
            return *p;
        }

        let two = T::TWO;
        Point::new(
            p.e1() - two * dist * nx * p.e0(),
            p.e2() - two * dist * ny * p.e0(),
            p.e0(),
        )
    }

    /// Angle between two lines (unsigned).
    pub fn angle(&self, other: &Line<T>) -> T {
        let (a1, b1) = self.normal();
        let (a2, b2) = other.normal();

        let dot = a1 * a2 + b1 * b2;
        let n1 = self.weight_norm();
        let n2 = other.weight_norm();

        if n1 < T::epsilon() || n2 < T::epsilon() {
            return T::zero();
        }

        let cos_angle = (dot / (n1 * n2)).min(T::one()).max(-T::one());
        cos_angle.acos()
    }

    /// Unitize to unit normal.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new(self.e12() / wn, self.e20() / wn, self.e01() / wn)
    }

    /// Check if line is degenerate (zero normal).
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.weight_norm() < T::epsilon()
    }

    /// Reverse.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(-self.e12(), -self.e20(), -self.e01())
    }

    /// Geometric norm.
    pub fn geometric_norm(&self) -> T {
        let wn = self.weight_norm();
        let bn = self.bulk_norm();
        (wn * wn + bn * bn).sqrt()
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Identity motor (no transformation).
    #[inline]
    pub fn identity() -> Self {
        Self::new(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Pure rotation around origin.
    pub fn from_rotation(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new(half.cos(), half.sin(), T::zero(), T::zero())
    }

    /// Pure translation.
    pub fn from_translation(dx: T, dy: T) -> Self {
        let half = T::one() / T::TWO;
        Self::new(T::one(), T::zero(), dx * half, dy * half)
    }

    /// Rotation around an arbitrary point.
    pub fn from_rotation_around(center: &Point<T>, angle: T) -> Self {
        if center.is_ideal() {
            return Self::from_rotation(angle);
        }

        let (cx, cy) = (center.x(), center.y());

        // Translate to origin, rotate, translate back
        let to_origin = Self::from_translation(-cx, -cy);
        let rotate = Self::from_rotation(angle);
        let from_origin = Self::from_translation(cx, cy);

        from_origin.compose(&rotate).compose(&to_origin)
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

    /// Compose motors: self then other.
    #[inline]
    pub fn compose(&self, other: &Motor<T>) -> Motor<T> {
        products::geometric_motor_motor(other, self)
    }

    /// Reverse.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(self.s(), -self.e12(), -self.e20(), -self.e01())
    }

    /// Inverse motor.
    pub fn inverse(&self) -> Self {
        let wn_sq = self.weight_norm_squared();
        let rev = self.reverse();
        Self::new(
            rev.s() / wn_sq,
            rev.e12() / wn_sq,
            rev.e20() / wn_sq,
            rev.e01() / wn_sq,
        )
    }

    /// Check if motor is unitized.
    #[inline]
    pub fn is_unitized(&self, tolerance: T) -> bool {
        (self.weight_norm_squared() - T::one()).abs() < tolerance
    }

    /// Extract rotation angle.
    pub fn rotation_angle(&self) -> T {
        self.e12().atan2(self.s()) * T::TWO
    }

    /// Extract translation vector.
    pub fn translation(&self) -> (T, T) {
        let two = T::TWO;
        (
            two * (self.s() * self.e20() + self.e12() * self.e01()),
            two * (self.s() * self.e01() - self.e12() * self.e20()),
        )
    }

    /// Linear interpolation.
    pub fn lerp(&self, other: &Motor<T>, t: T) -> Self {
        let one_minus_t = T::one() - t;
        Self::new(
            self.s() * one_minus_t + other.s() * t,
            self.e12() * one_minus_t + other.e12() * t,
            self.e20() * one_minus_t + other.e20() * t,
            self.e01() * one_minus_t + other.e01() * t,
        )
    }

    /// Access components as tuple.
    #[inline]
    pub fn components(&self) -> (T, T, T, T) {
        (self.s(), self.e12(), self.e20(), self.e01())
    }
}
```

## Phase 4: Update Module Structure

### File: `src/specialized/projective/dim2/mod.rs`

```rust
//! 2D Projective Geometric Algebra.
//!
//! This module provides optimized types for 2D PGA (Cl(2,0,1)):
//!
//! - [`Point`] - Grade-1 elements (points in homogeneous coordinates)
//! - [`Line`] - Grade-2 elements (lines in implicit form)
//! - [`Motor`] - Even subalgebra elements (rigid transformations)
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim2::{Point, Line, Motor};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // Create points
//! let p1 = Point::new(1.0, 0.0);
//! let p2 = Point::new(0.0, 1.0);
//!
//! // Join points to get line through them
//! let line = p1.join(&p2);
//!
//! // Create a rotation around origin
//! let motor = Motor::from_rotation(FRAC_PI_2);
//! let rotated = motor.transform_point(&p1);
//! ```

mod generated;
mod extensions;
mod ops;
mod conversions;

pub use generated::types::{Line, Motor, Point};
pub use conversions::CONVERSION_TOLERANCE;

#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra;

#[cfg(feature = "rerun")]
mod rerun;

#[cfg(test)]
mod tests;
```

## Phase 5: Update Operators

### File: `src/specialized/projective/dim2/ops.rs`

```rust
//! Operator overloads for 2D PGA types.

use super::generated::products;
use super::generated::types::{Line, Motor, Point};
use crate::scalar::Float;
use std::ops::{Add, Mul, Neg, Sub};

// ============================================================================
// Point operators
// ============================================================================

impl<T: Float> Neg for Point<T> {
    type Output = Self;
    fn neg(self) -> Self { Self::new(-self.e1(), -self.e2(), -self.e0()) }
}

impl<T: Float> Add for Point<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.e1() + rhs.e1(), self.e2() + rhs.e2(), self.e0() + rhs.e0())
    }
}

impl<T: Float> Sub for Point<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.e1() - rhs.e1(), self.e2() - rhs.e2(), self.e0() - rhs.e0())
    }
}

impl<T: Float> Mul<T> for Point<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self::new(self.e1() * rhs, self.e2() * rhs, self.e0() * rhs)
    }
}

// ============================================================================
// Line operators
// ============================================================================

impl<T: Float> Neg for Line<T> {
    type Output = Self;
    fn neg(self) -> Self { Self::new(-self.e12(), -self.e20(), -self.e01()) }
}

impl<T: Float> Add for Line<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.e12() + rhs.e12(), self.e20() + rhs.e20(), self.e01() + rhs.e01())
    }
}

impl<T: Float> Sub for Line<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.e12() - rhs.e12(), self.e20() - rhs.e20(), self.e01() - rhs.e01())
    }
}

impl<T: Float> Mul<T> for Line<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self::new(self.e12() * rhs, self.e20() * rhs, self.e01() * rhs)
    }
}

// ============================================================================
// Motor operators
// ============================================================================

impl<T: Float> Neg for Motor<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.s(), -self.e12(), -self.e20(), -self.e01())
    }
}

impl<T: Float> Add for Motor<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(
            self.s() + rhs.s(),
            self.e12() + rhs.e12(),
            self.e20() + rhs.e20(),
            self.e01() + rhs.e01(),
        )
    }
}

impl<T: Float> Sub for Motor<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(
            self.s() - rhs.s(),
            self.e12() - rhs.e12(),
            self.e20() - rhs.e20(),
            self.e01() - rhs.e01(),
        )
    }
}

impl<T: Float> Mul<T> for Motor<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self::new(
            self.s() * rhs,
            self.e12() * rhs,
            self.e20() * rhs,
            self.e01() * rhs,
        )
    }
}

impl<T: Float> Mul<Motor<T>> for Motor<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        products::geometric_motor_motor(&self, &rhs)
    }
}
```

## Phase 6: Tests

### File: `src/specialized/projective/dim2/tests.rs`

```rust
//! Tests for 2D PGA codegen migration.

use super::*;
use crate::multivector::Multivector;
use crate::signature::Projective2;
use crate::test_utils::ABS_DIFF_EQ_EPS;
use approx::abs_diff_eq;
use proptest::prelude::*;

proptest! {
    // ========================================================================
    // Join/Meet operations
    // ========================================================================

    /// Join of two points gives line through them.
    #[test]
    fn join_gives_line_through_points(
        x1 in -10.0f64..10.0, y1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0, y2 in -10.0f64..10.0,
    ) {
        let p1 = Point::new(x1, y1, 1.0);
        let p2 = Point::new(x2, y2, 1.0);
        let line = p1.join(&p2);

        // Both points should be on the line (distance ≈ 0)
        let d1 = line.distance_to_point(&p1);
        let d2 = line.distance_to_point(&p2);

        prop_assert!(abs_diff_eq!(d1, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(d2, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Meet of two non-parallel lines gives intersection point.
    #[test]
    fn meet_gives_intersection(
        a1 in -1.0f64..1.0, b1 in -1.0f64..1.0, c1 in -10.0f64..10.0,
        a2 in -1.0f64..1.0, b2 in -1.0f64..1.0, c2 in -10.0f64..10.0,
    ) {
        // Ensure lines are not parallel
        let det = a1 * b2 - a2 * b1;
        if det.abs() < 0.1 {
            return Ok(());
        }

        let line1 = Line::from_implicit(a1, b1, c1);
        let line2 = Line::from_implicit(a2, b2, c2);
        let intersection = line1.meet(&line2);

        if let Some((x, y)) = intersection.to_cartesian() {
            // Point should be on both lines
            let d1 = line1.distance_to_point(&intersection);
            let d2 = line2.distance_to_point(&intersection);

            prop_assert!(abs_diff_eq!(d1, 0.0, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(d2, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Motor operations
    // ========================================================================

    /// Motor sandwich matches generic Multivector.
    #[test]
    fn motor_transform_point_matches_multivector(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0, ty in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0,
    ) {
        let motor = Motor::from_translation(tx, ty).compose(&Motor::from_rotation(angle));
        let point = Point::new(px, py, 1.0);

        let spec = motor.transform_point(&point);

        let mv_motor = Multivector::<f64, Projective2>::from(motor);
        let mv_point = Multivector::<f64, Projective2>::from(point);
        let gen = mv_motor.sandwich(&mv_point);
        let gen_point = Point::from_multivector_unchecked(&gen);

        if let (Some((sx, sy)), Some((gx, gy))) = (spec.to_cartesian(), gen_point.to_cartesian()) {
            prop_assert!(abs_diff_eq!(sx, gx, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(sy, gy, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    /// Motor preserves distance.
    #[test]
    fn motor_preserves_distance(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0, ty in -10.0f64..10.0,
        p1x in -10.0f64..10.0, p1y in -10.0f64..10.0,
        p2x in -10.0f64..10.0, p2y in -10.0f64..10.0,
    ) {
        let motor = Motor::from_translation(tx, ty)
            .compose(&Motor::from_rotation(angle))
            .normalized();
        let p1 = Point::new(p1x, p1y, 1.0);
        let p2 = Point::new(p2x, p2y, 1.0);

        let d_before = p1.distance(&p2);
        let t1 = motor.transform_point(&p1);
        let t2 = motor.transform_point(&p2);
        let d_after = t1.distance(&t2);

        prop_assert!(abs_diff_eq!(d_before, d_after, epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Motor composition is associative.
    #[test]
    fn motor_compose_associative(
        a1 in -std::f64::consts::PI..std::f64::consts::PI,
        a2 in -std::f64::consts::PI..std::f64::consts::PI,
        a3 in -std::f64::consts::PI..std::f64::consts::PI,
    ) {
        let m1 = Motor::from_rotation(a1);
        let m2 = Motor::from_rotation(a2);
        let m3 = Motor::from_rotation(a3);

        let lhs = m1.compose(&m2).compose(&m3);
        let rhs = m1.compose(&m2.compose(&m3));

        prop_assert!(abs_diff_eq!(lhs.s(), rhs.s(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(lhs.e12(), rhs.e12(), epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Motor inverse undoes transformation.
    #[test]
    fn motor_inverse(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0, ty in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0,
    ) {
        let motor = Motor::from_translation(tx, ty)
            .compose(&Motor::from_rotation(angle));
        let point = Point::new(px, py, 1.0);

        let transformed = motor.transform_point(&point);
        let back = motor.inverse().transform_point(&transformed);

        if let (Some((ox, oy)), Some((bx, by))) = (point.to_cartesian(), back.to_cartesian()) {
            prop_assert!(abs_diff_eq!(ox, bx, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(oy, by, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Conversion round-trips
    // ========================================================================

    /// Point round-trip through Multivector.
    #[test]
    fn point_roundtrip(
        e1 in -100.0f64..100.0,
        e2 in -100.0f64..100.0,
        e0 in -100.0f64..100.0,
    ) {
        let point = Point::new(e1, e2, e0);
        let mv = Multivector::<f64, Projective2>::from(point);
        let back = Point::from_multivector_unchecked(&mv);

        prop_assert!(abs_diff_eq!(point.e1(), back.e1(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(point.e2(), back.e2(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(point.e0(), back.e0(), epsilon = ABS_DIFF_EQ_EPS));
    }

    /// Line round-trip through Multivector.
    #[test]
    fn line_roundtrip(
        e12 in -100.0f64..100.0,
        e20 in -100.0f64..100.0,
        e01 in -100.0f64..100.0,
    ) {
        let line = Line::new(e12, e20, e01);
        let mv = Multivector::<f64, Projective2>::from(line);
        let back = Line::from_multivector_unchecked(&mv);

        prop_assert!(abs_diff_eq!(line.e12(), back.e12(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(line.e20(), back.e20(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(line.e01(), back.e01(), epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Verification Checklist

- [ ] TOML specification created
- [ ] Code generation succeeds
- [ ] Extensions module complete
- [ ] Operators updated
- [ ] Join/meet operations verified
- [ ] Motor transformations verified
- [ ] nalgebra conversions work
- [ ] All existing tests pass
- [ ] `cargo test --all-features` passes

## Files Changed

| File | Action |
|------|--------|
| `algebras/projective2.toml` | Create |
| `src/specialized/projective/dim2/mod.rs` | Update |
| `src/specialized/projective/dim2/generated/` | Create |
| `src/specialized/projective/dim2/extensions.rs` | Create |
| `src/specialized/projective/dim2/ops.rs` | Update |
| `src/specialized/projective/dim2/types.rs` | Delete |
| `src/specialized/projective/dim2/tests.rs` | Create |

## Summary

PRD-15.4 completes the migration of all Euclidean and Projective algebras to use the codegen framework:

| PRD | Algebra | Types | Complexity |
|-----|---------|-------|------------|
| 15.1 | Euclidean 3D | 5 | High (validates pattern) |
| 15.2 | Euclidean 2D | 3 | Low (confirms pattern) |
| 15.3 | Projective 3D | 5 | High (constraints) |
| 15.4 | Projective 2D | 3 | Low (confirms PGA pattern) |

After completing all four sub-PRDs, the specialized modules will be fully migrated to generated code while maintaining API compatibility.
