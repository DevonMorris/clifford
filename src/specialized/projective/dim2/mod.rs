//! 2D Projective Geometric Algebra types.
//!
//! This module provides specialized types for 2D PGA `Cl(2,0,1)`.
//!
//! # Point-Based Formulation
//!
//! In the point-based formulation:
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 1 | [`Point`] | e₁, e₂, e₀ | Point (homogeneous) |
//! | 2 | [`Line`] | e₁₂, e₀₁, e₀₂ | Line |
//! | 0+2 | [`Motor`] | s, e₁₂, e₀₁, e₀₂ | Rigid transformation |
//!
//! # Homogeneous Coordinates
//!
//! A point `P = x·e₁ + y·e₂ + w·e₀` represents:
//! - Finite point `(x/w, y/w)` when `w ≠ 0`
//! - Ideal point (at infinity) in direction `(x, y)` when `w = 0`
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim2::{Point, Line, Motor};
//! use std::f64::consts::FRAC_PI_2;
//! use approx::abs_diff_eq;
//!
//! // Points
//! let origin: Point<f64> = Point::origin();
//! let p = Point::from_cartesian(3.0, 4.0);
//!
//! // Line meet - intersection of two lines
//! let x_axis: Line<f64> = Line::x_axis();
//! let y_axis: Line<f64> = Line::y_axis();
//! let intersection = x_axis.meet(&y_axis);
//!
//! // 90° rotation around origin
//! let rotor = Motor::from_rotation(FRAC_PI_2);
//!
//! // Translation
//! let translation = Motor::from_translation(1.0, 2.0);
//!
//! // Compose motors
//! let combined = rotor.compose(&translation);
//! ```

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions
mod extensions;

// TODO: Add nalgebra integration
// #[cfg(any(
//     feature = "nalgebra-0_32",
//     feature = "nalgebra-0_33",
//     feature = "nalgebra-0_34"
// ))]
// mod nalgebra;

// TODO: Add rerun visualization support
// #[cfg(feature = "rerun-0_28")]
// mod rerun;

// Re-export generated types and products
pub use generated::products;
pub use generated::types::{Line, Motor, Point, Scalar, Trivector};

// Re-export wrapper type aliases
pub use generated::types::{BulkMotor, IdealLine, IdealPoint, UnitizedLine, UnitizedPoint};
