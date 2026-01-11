//! Conformal Geometric Algebra (CGA) types.
//!
//! This module provides specialized types for Conformal Geometric Algebra,
//! which extends Euclidean space with two extra dimensions to enable elegant
//! representations of circles, spheres, and transformations.
//!
//! # Currently Supported
//!
//! - **[`dim3`]**: 3D CGA (`Cl(4,1,0)`) with 32 basis blades
//!
//! # Conformal Model
//!
//! CGA embeds Euclidean points as null vectors in a higher-dimensional space.
//! A 3D point `(x, y, z)` becomes:
//!
//! ```text
//! P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
//! ```
//!
//! This embedding has remarkable properties:
//! - Points are null vectors: `P · P = 0`
//! - Distance from inner product: `d² = -2(P₁ · P₂)`
//! - Spheres are grade-4 blades
//! - Circles are grade-3 blades
//! - Lines and planes are "flat" versions of circles and spheres
//!
//! # References
//!
//! - <https://conformalgeometricalgebra.org/wiki/>
//!
//! # Example
//!
//! ```
//! use clifford::specialized::conformal::dim3::Point;
//!
//! // Create points at (0,0,0) and (3,4,0)
//! let p1 = Point::<f64>::origin();
//! let p2 = Point::new(3.0, 4.0, 0.0);
//!
//! // Compute distance using CGA inner product
//! let d = p1.distance(&p2);
//! assert!((d - 5.0).abs() < 1e-10);  // 3-4-5 triangle
//! ```

pub mod dim3;
