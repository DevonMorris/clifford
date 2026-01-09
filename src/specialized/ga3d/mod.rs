//! 3D Euclidean Geometric Algebra types.
//!
//! This module provides specialized types for 3D Euclidean space `Cl(3,0,0)`,
//! the most commonly used geometric algebra for graphics, physics, and robotics.
//!
//! # Types by Grade
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 0 | [`Scalar`] | 1 | Magnitude/weight |
//! | 1 | [`Vec3`] | e₁, e₂, e₃ | Direction/position |
//! | 2 | [`Bivec3`] | e₁₂, e₁₃, e₂₃ | Oriented plane/rotation axis |
//! | 3 | [`Trivec3`] | e₁₂₃ | Oriented volume |
//!
//! # Rotors
//!
//! [`Rotor3`] represents rotations as `scalar + bivector`. This is equivalent
//! to unit quaternions but expressed in the geometric algebra framework:
//!
//! ```
//! use clifford::specialized::ga3d::{Vec3, Bivec3, Rotor3};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // Rotate 90° around the z-axis (in the xy-plane)
//! let rotor = Rotor3::from_angle_plane(FRAC_PI_2, Bivec3::unit_xy());
//! let v = Vec3::new(1.0, 0.0, 0.0);
//! let rotated = rotor.rotate(v);
//!
//! // x-axis rotated 90° becomes y-axis
//! assert!((rotated.y - 1.0).abs() < 1e-10);
//! ```

mod conversions;
mod ops;
mod types;

#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;

pub use conversions::{ConversionError, Specialized3};
pub use types::*;
