//! 3D Euclidean Geometric Algebra types.
//!
//! This module provides specialized types for 3D Euclidean space `Cl(3,0,0)`,
//! the most commonly used geometric algebra for graphics, physics, and robotics.
//!
//! # Types by Grade
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 0 | Scalar | 1 | Magnitude/weight |
//! | 1 | [`Vector`] | e₁, e₂, e₃ | Direction/position |
//! | 2 | [`Bivector`] | e₁₂, e₁₃, e₂₃ | Oriented plane/rotation axis |
//! | 3 | [`Trivector`] | e₁₂₃ | Oriented volume |
//!
//! # Rotors and Quaternions
//!
//! If you know quaternions, you already know rotors! A 3D rotor **is** a quaternion,
//! just written in a different basis:
//!
//! | Quaternion | Rotor | Meaning |
//! |------------|-------|---------|
//! | `w` | `s` (scalar) | cos(θ/2) |
//! | `x·i` | `e₂₃` | rotation in yz-plane |
//! | `y·j` | `e₃₁` (= -e₁₃) | rotation in zx-plane |
//! | `z·k` | `e₁₂` | rotation in xy-plane |
//!
//! The quaternion formula `q = w + xi + yj + zk` becomes:
//!
//! ```text
//! R = s + e₂₃·x + e₃₁·y + e₁₂·z
//!   = cos(θ/2) + sin(θ/2)·B̂
//! ```
//!
//! where `B̂` is the unit bivector of the rotation plane.
//!
//! ## Why Bivectors Beat Axis Vectors
//!
//! Traditional: "Rotate 90° around the z-axis"
//! GA: "Rotate 90° in the xy-plane"
//!
//! These describe the same rotation, but the GA version is more fundamental:
//!
//! - **Rotation happens IN a plane**, not AROUND an axis
//! - The "axis" is just the line perpendicular to that plane (only works in 3D!)
//! - In 2D there's no axis, but there IS a plane (the whole 2D space)
//! - In 4D+ there are multiple perpendicular directions—bivectors still work
//!
//! The axis vector is the **dual** of the rotation bivector: `axis = *B`
//!
//! # The Sandwich Product
//!
//! Rotors rotate vectors using the **sandwich product**:
//!
//! ```text
//! v' = R v R̃    (R̃ is the reverse of R)
//! ```
//!
//! This is identical to the quaternion formula `v' = q v q*` (where `q*` is conjugate).
//!
//! Why a sandwich? The geometric product `R v` would change the grade of `v`,
//! but `R v R̃` guarantees the result is still a vector.
//!
//! # Cross Product via GA
//!
//! The cross product is the dual of the wedge product:
//!
//! ```text
//! a × b = *(a ∧ b)
//! ```
//!
//! Where `*` is the Hodge dual (multiply by inverse pseudoscalar e₁₂₃⁻¹ = e₃₂₁).
//!
//! | Operation | Result | Type |
//! |-----------|--------|------|
//! | `a ∧ b` | Oriented plane containing a, b | Bivector |
//! | `*(a ∧ b)` | Vector perpendicular to that plane | Vector |
//! | `a × b` | Same as above | Vector |
//!
//! # nalgebra Integration
//!
//! With the `nalgebra-0_33` feature (or `0_32`/`0_34`), conversions are provided:
//!
//! | clifford | nalgebra |
//! |----------|----------|
//! | `Vector` | `Vector3<T>` |
//! | `Rotor` | `UnitQuaternion<T>` |
//! | `Bivector` | (no direct equivalent) |
//!
//! ```ignore
//! use clifford::specialized::euclidean::dim3::Rotor;
//! use nalgebra::UnitQuaternion;
//!
//! let rotor = Rotor::from_angle_z(std::f64::consts::FRAC_PI_2);
//! let quat: UnitQuaternion<f64> = rotor.into();
//! ```
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::euclidean::dim3::{Vector, Bivector, Rotor};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // Rotate 90° around the z-axis (in the xy-plane)
//! let rotor = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
//! let v = Vector::new(1.0, 0.0, 0.0);
//! let rotated = rotor.rotate(v);
//!
//! // x-axis rotated 90° becomes y-axis
//! assert!((rotated.y() - 1.0).abs() < 1e-10);
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

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
pub use nalgebra::NalgebraConversionError;

#[cfg(feature = "rerun-0_28")]
mod rerun;

// Re-export generated types and wrapper aliases
pub use generated::types::*;

// Re-export products module for direct access to algebraic products
pub use generated::products;

// Re-export Even as alias from extensions
pub use extensions::Even;
