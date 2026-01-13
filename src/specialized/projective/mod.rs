//! Specialized types for Projective Geometric Algebra (PGA).
//!
//! This module provides optimized, ergonomic types for point-based PGA organized
//! by dimension. Currently supported:
//!
//! - **[`dim2`]**: 2D PGA `Cl(2,0,1)`: Point, Line, Motor
//! - **[`dim3`]**: 3D PGA `Cl(3,0,1)`: Point, Line, Plane, Motor, Flector
//!
//! # Point-Based Formulation
//!
//! This library uses the point-based formulation of PGA where:
//! - **Grade 1 (vectors)**: Points in homogeneous coordinates
//! - **Grade 2 (bivectors)**: Lines
//! - **Grade 3 (trivectors)**: Planes (in 3D only)
//!
//! This is consistent with Euclidean GA conventions and nalgebra's `Point` types.
//!
//! # Type Organization
//!
//! ```text
//! projective/
//!   dim2/   - 2D PGA: Point, Line, Motor
//!   dim3/   - 3D PGA: Point, Line, Plane, Motor, Flector
//! ```
//!
//! # Motors
//!
//! Motors represent rigid body transformations (rotation + translation).
//! A motor `M` transforms geometric objects via the sandwich product:
//! `X' = M X M̃`
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim2::{Point, Line, Motor};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // Create points using from_cartesian (sets w=1)
//! let p1: Point<f64> = Point::from_cartesian(0.0, 0.0);  // Origin
//! let p2: Point<f64> = Point::from_cartesian(1.0, 0.0);  // Point at (1, 0)
//!
//! // Line meet - intersection of two lines
//! let x_axis: Line<f64> = Line::x_axis();
//! let y_axis: Line<f64> = Line::y_axis();
//! let intersection = x_axis.meet(&y_axis);
//!
//! // Motor: 90° rotation then translate by (1, 2)
//! let motor = Motor::from_translation(1.0, 2.0)
//!     .compose(&Motor::from_rotation(FRAC_PI_2));
//! ```

pub mod dim2;
pub mod dim3;
mod errors;

pub use errors::{LineConstraintError, MotorConstraintError};

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
pub use errors::PointConversionError;
