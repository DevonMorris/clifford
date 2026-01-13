//! 3D Projective Geometric Algebra types.
//!
//! This module provides specialized types for 3D PGA `Cl(3,0,1)`.
//!
//! # Types by Grade
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 1 | [`Point`] | e₁, e₂, e₃, e₀ | Point (homogeneous) |
//! | 2 | [`Line`] | 6 components | Line (Plücker coords) |
//! | 3 | [`Plane`] | 4 components | Plane |
//! | 0+2+4 | [`Motor`] | 8 components | Rigid transformation (rotation + translation) |
//! | 1+3 | [`Flector`] | 8 components | Reflection (improper isometry) |
//!
//! # Points as Homogeneous Coordinates
//!
//! A PGA point is just a homogeneous coordinate:
//!
//! ```text
//! P = x·e₁ + y·e₂ + z·e₃ + w·e₀
//! ```
//!
//! | Type | Condition | Cartesian |
//! |------|-----------|-----------|
//! | Finite point | `w ≠ 0` | `(x/w, y/w, z/w)` |
//! | Point at infinity | `w = 0` | Direction `(x, y, z)` |
//!
//! Use `Point::from_cartesian(x, y, z)` for finite points (sets `w = 1`).
//!
//! # Lines and Plücker Coordinates
//!
//! If you've worked with 3D line representations, you may know **Plücker coordinates**:
//! a 6-component representation `(d, m)` where `d` is direction and `m` is moment.
//!
//! In PGA, a line is a grade-2 blade (bivector) with the same 6 components:
//!
//! | Plücker | PGA | Meaning |
//! |---------|-----|---------|
//! | `d` (direction) | `e₂₃, e₁₃, e₁₂` | "bulk" part (Euclidean bivector) |
//! | `m` (moment) | `e₀₁, e₀₂, e₀₃` | "weight" part (involves e₀) |
//!
//! **Plücker constraint**: Valid lines satisfy `d · m = 0`. In PGA this is automatic
//! when lines are constructed geometrically (via join of points or meet of planes).
//!
//! # Join and Meet: Spanning and Intersecting
//!
//! The two fundamental operations in PGA are **join** (∨) and **meet** (∧):
//!
//! | Operation | Symbol | Result | Meaning |
//! |-----------|--------|--------|---------|
//! | **Join** | `a ∨ b` | Higher grade | Span (smallest object containing both) |
//! | **Meet** | `a ∧ b` | Lower grade | Intersection |
//!
//! ## Join Examples (building up)
//!
//! | Inputs | Result |
//! |--------|--------|
//! | Point ∨ Point | Line through both |
//! | Point ∨ Line | Plane containing both |
//! | Line ∨ Point | Plane containing both |
//!
//! ## Meet Examples (cutting down)
//!
//! | Inputs | Result |
//! |--------|--------|
//! | Plane ∧ Plane | Line of intersection |
//! | Plane ∧ Line | Point of intersection |
//! | Line ∧ Plane | Point of intersection |
//!
//! **Note**: In degenerate cases (parallel lines, point on plane, etc.),
//! the result is an "ideal" (infinite) element.
//!
//! # Motors: Better Than 4×4 Matrices
//!
//! A **motor** represents a rigid body transform (rotation + translation):
//!
//! | 4×4 Matrix | Motor |
//! |------------|-------|
//! | 16 floats | 8 floats |
//! | Orthogonality drift | Cheap normalization |
//! | Complex interpolation | Natural SLERP |
//! | `M₂ × M₁` | `M₂ * M₁` (same!) |
//!
//! Motors transform objects via the **sandwich product**: `X' = M X M̃`
//!
//! ## Motor Composition
//!
//! ```text
//! // Transform by M₁ then M₂
//! let combined = m2.compose(&m1);  // M₂ * M₁
//!
//! // Same as applying m1 then m2:
//! let p1 = m1.transform_point(&p);
//! let p2 = m2.transform_point(&p1);
//! // p2 == combined.transform_point(&p)
//! ```
//!
//! ## Motor from Components
//!
//! - `Motor::from_rotation_x(angle)`, `_y`, `_z` — pure rotation
//! - `Motor::from_translation(dx, dy, dz)` — pure translation
//! - `Motor::from_axis_angle(&axis, angle)` — rotation around arbitrary axis
//! - Compose with `.compose()` for combined transforms
//!
//! # Flectors: Reflections
//!
//! A **flector** is an improper isometry (includes a reflection). While motors
//! have `det = +1`, flectors have `det = -1`. Use for mirror transformations.
//!
//! # nalgebra Integration
//!
//! With `nalgebra-0_33` (or `0_32`/`0_34`), conversions are provided:
//!
//! | clifford | nalgebra |
//! |----------|----------|
//! | `Point` | `Point3<T>` |
//! | `Motor` | `Isometry3<T>` |
//! | `Line` | (no direct equivalent) |
//! | `Plane` | (no direct equivalent) |
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim3::{Point, Line, Motor};
//! use std::f64::consts::FRAC_PI_2;
//! use approx::abs_diff_eq;
//!
//! // Create points
//! let p1 = Point::from_cartesian(0.0, 0.0, 0.0);
//! let p2 = Point::from_cartesian(1.0, 0.0, 0.0);
//!
//! // Join two points to get a line
//! let line = p1.join(&p2);
//!
//! // Create a rotation motor
//! let rotor = Motor::from_rotation_z(FRAC_PI_2);
//!
//! // Create a translation motor
//! let translation = Motor::from_translation(1.0, 2.0, 3.0);
//!
//! // Compose motors
//! let combined = rotor.compose(&translation);
//! ```

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions
mod extensions;

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra;

#[cfg(feature = "rerun-0_28")]
mod rerun;

// Re-export generated types and products
pub use generated::products;
pub use generated::types::{Flector, Line, Motor, Plane, Point, Quadvector, Scalar};

// Re-export wrapper type aliases
pub use generated::types::{
    BulkFlector, BulkMotor, IdealLine, IdealPlane, IdealPoint, UnitizedLine, UnitizedPlane,
    UnitizedPoint,
};

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
pub use nalgebra::NalgebraConversionError;
