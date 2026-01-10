//! Metric signature definitions for Clifford algebras.
//!
//! This module provides the [`Signature`] trait and concrete signature types
//! that define the metric structure of different geometric algebras.
//!
//! # What is a Metric Signature?
//!
//! In a Clifford algebra, each basis vector `e_i` squares to a scalar:
//! `e_i² = +1`, `-1`, or `0`. The signature `(P, Q, R)` counts how many
//! basis vectors fall into each category:
//!
//! - `P`: Positive, `e_i² = +1` (space-like in physics)
//! - `Q`: Negative, `e_i² = -1` (time-like in physics)
//! - `R`: Zero, `e_i² = 0` (null/degenerate)
//!
//! # Available Signatures
//!
//! ## Euclidean (all positive)
//! - [`Euclidean2`]: 2D plane, `Cl(2,0,0)`
//! - [`Euclidean3`]: 3D space, `Cl(3,0,0)`
//! - [`Euclidean4`]: 4D space, `Cl(4,0,0)`
//!
//! ## Projective (Euclidean + null)
//! - [`Projective2`]: 2D PGA, `Cl(2,0,1)` - point-based formulation
//! - [`Projective3`]: 3D PGA, `Cl(3,0,1)` - point-based formulation
//!
//! ## Coming Soon
//! - CGA signatures for conformal geometry
//! - Minkowski signature for special relativity
//!
//! # Example
//!
//! ```
//! use clifford::prelude::*;
//!
//! // Query the algebra structure
//! assert_eq!(Euclidean3::DIM, 3);
//! assert_eq!(Euclidean3::num_blades(), 8);
//!
//! // All basis vectors in Euclidean space square to +1
//! assert_eq!(Euclidean3::metric(0), 1);
//! assert_eq!(Euclidean3::metric(1), 1);
//! assert_eq!(Euclidean3::metric(2), 1);
//!
//! // 3D PGA has one null vector
//! assert_eq!(Projective3::DIM, 4);
//! assert_eq!(Projective3::metric(3), 0);  // e₀² = 0
//! ```

mod euclidean;
mod metric;
mod projective;

pub use euclidean::{Euclidean2, Euclidean3, Euclidean4};
pub use metric::Signature;
pub use projective::{Projective2, Projective3};
