//! 3D Projective Geometric Algebra types.
//!
//! This module provides specialized types for 3D PGA `Cl(3,0,1)`.
//!
//! # Point-Based Formulation
//!
//! In the point-based formulation:
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 1 | [`Point`] | e₁, e₂, e₃, e₀ | Point (homogeneous) |
//! | 2 | [`Line`] | 6 components | Line (Plücker coords) |
//! | 3 | [`Plane`] | 4 components | Plane |
//! | 0+2+4 | [`Motor`] | 8 components | Rigid transformation (proper isometry) |
//! | 1+3 | [`Flector`] | 8 components | Reflection (improper isometry) |
//!
//! # Homogeneous Coordinates
//!
//! A point `P = x·e₁ + y·e₂ + z·e₃ + w·e₀` represents:
//! - Finite point `(x/w, y/w, z/w)` when `w ≠ 0`
//! - Ideal point (at infinity) in direction `(x, y, z)` when `w = 0`
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim3::{Point, Motor};
//! use std::f64::consts::FRAC_PI_2;
//! use approx::abs_diff_eq;
//!
//! // Create a point at (1, 0, 0)
//! let p = Point::from_cartesian(1.0, 0.0, 0.0);
//!
//! // 90° rotation around Z axis
//! let rotor = Motor::from_rotation_z(FRAC_PI_2);
//! let rotated = rotor.transform_point(&p);
//!
//! // (1,0,0) rotated 90° around Z becomes (0,1,0)
//! assert!(abs_diff_eq!(rotated.x(), 0.0, epsilon = 1e-10));
//! assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
//! assert!(abs_diff_eq!(rotated.z(), 0.0, epsilon = 1e-10));
//!
//! // Translation
//! let translation = Motor::from_translation(1.0, 2.0, 3.0);
//! let translated = translation.transform_point(&p);
//!
//! assert!(abs_diff_eq!(translated.x(), 2.0, epsilon = 1e-10));
//! assert!(abs_diff_eq!(translated.y(), 2.0, epsilon = 1e-10));
//! assert!(abs_diff_eq!(translated.z(), 3.0, epsilon = 1e-10));
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
pub use nalgebra::{FlectorConversionError, PointConversionError, Reflection3};

#[cfg(feature = "rerun-0_28")]
mod rerun;

#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;

// Re-export generated types
pub use generated::types::{Flector, Line, Motor, Plane, Point, Quadvector, Scalar};
