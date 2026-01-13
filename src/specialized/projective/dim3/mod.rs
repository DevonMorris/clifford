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

#[cfg(feature = "rerun-0_28")]
mod rerun;

// Re-export generated types and products
pub use generated::products;
pub use generated::types::{Flector, Line, Motor, Plane, Point, Quadvector, Scalar};

// Re-export wrapper type aliases
pub use generated::types::{BulkFlector, BulkMotor, IdealLine, IdealPlane, IdealPoint};
