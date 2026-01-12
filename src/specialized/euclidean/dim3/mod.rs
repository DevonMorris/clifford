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
//! # Rotors
//!
//! [`Rotor`] represents rotations as `scalar + bivector`. This is equivalent
//! to unit quaternions but expressed in the geometric algebra framework:
//!
//! ```
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

// Re-export generated types
pub use generated::types::{Bivector, Rotor, Scalar, Trivector, Vector};

// Re-export products module for direct access to algebraic products
pub use generated::products;

// Re-export Even as alias from extensions
pub use extensions::Even;
