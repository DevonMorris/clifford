//! Specialized types for Projective Geometric Algebra (PGA).
//!
//! This module provides optimized, ergonomic types for point-based PGA organized
//! by dimension. Currently supported:
//!
//! - **[`dim2`]**: 2D PGA `Cl(2,0,1)`: Point, Line, Motor
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
//!   dim3/   - 3D PGA: Point, Line, Plane, Motor (future)
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
//! // Create points
//! let p1 = Point::new(0.0, 0.0);  // Origin
//! let p2 = Point::new(1.0, 0.0);  // Point at (1, 0)
//!
//! // Line through two points
//! let line = p1.join(&p2);
//!
//! // Motor: 90° rotation then translate by (1, 2)
//! let motor = Motor::from_translation(1.0, 2.0)
//!     .compose(&Motor::from_rotation(FRAC_PI_2));
//!
//! // Transform a point
//! let transformed = motor.transform_point(&p2);
//! ```

pub mod dim2;
