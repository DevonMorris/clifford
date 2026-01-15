//! Elliptic Projective Geometric Algebra.
//!
//! This module provides types for elliptic geometry, which is the geometry
//! of the projective plane with an elliptic metric.
//!
//! # Available Dimensions
//!
//! - **[`dim2`]**: 2D elliptic plane (RP²) using Cl(3,0,0)
//!
//! # What is Elliptic Geometry?
//!
//! Elliptic geometry is a non-Euclidean geometry where:
//! - There are **no parallel lines** (any two lines intersect)
//! - The sum of angles in a triangle **exceeds 180°**
//! - Space is **finite but unbounded** (like a sphere)
//!
//! The elliptic plane can be visualized as the unit sphere S² with
//! antipodal points identified (RP²), or equivalently as lines through
//! the origin in R³.
//!
//! # Algebra Structure
//!
//! The 2D elliptic plane uses Cl(3,0,0):
//!
//! | Grade | Elements | Geometric Meaning |
//! |-------|----------|-------------------|
//! | 0 | 1 | Scalar |
//! | 1 | e₁, e₂, e₃ | Points (homogeneous coords) |
//! | 2 | e₁₂, e₁₃, e₂₃ | Lines (great circles) |
//! | 3 | e₁₂₃ | Pseudoscalar |
//!
//! # Point-Line Duality
//!
//! Elliptic geometry exhibits perfect duality:
//! - Both points and lines have 3 components
//! - Operations on points mirror operations on lines
//! - Join of two points = line through them
//! - Meet of two lines = point of intersection
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::elliptic::dim2::{Point, Line, Rotor};
//! use clifford::ops::Wedge;
//!
//! // Two points on the elliptic plane
//! let p1 = Point::new(1.0, 0.0, 0.0);
//! let p2 = Point::new(0.0, 1.0, 0.0);
//!
//! // The line through both points
//! let line = p1.wedge(&p2);
//! ```

pub mod dim2;
