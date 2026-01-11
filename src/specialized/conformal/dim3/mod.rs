//! 3D Conformal Geometric Algebra types.
//!
//! This module provides specialized types for 3D CGA `Cl(4,1,0)`, which embeds
//! 3D Euclidean space into a 5D conformal space.
//!
//! # Types
//!
//! | Type | Grade | Description |
//! |------|-------|-------------|
//! | [`Point`] | 1 | Round point (null vector) |
//! | [`Sphere`] | 4 | Sphere (or imaginary sphere) |
//!
//! # Conformal Embedding
//!
//! A Euclidean point `(x, y, z)` is embedded as:
//!
//! ```text
//! P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
//! ```
//!
//! Where:
//! - `e₁, e₂, e₃` are the Euclidean basis vectors
//! - `e₀ = (e₋ - e₊)/2` is the origin
//! - `e∞ = e₋ + e₊` is the point at infinity
//!
//! # Key Properties
//!
//! - All conformal points are null: `P · P = 0`
//! - Distance from inner product: `d² = -2(P₁ · P₂)`
//!
//! # Example
//!
//! ```
//! use clifford::specialized::conformal::dim3::Point;
//!
//! let p = Point::<f64>::new(1.0, 2.0, 3.0);
//!
//! // Verify null constraint
//! assert!(p.is_null(1e-10));
//!
//! // Coordinates round-trip
//! assert!((p.x() - 1.0).abs() < 1e-10);
//! assert!((p.y() - 2.0).abs() < 1e-10);
//! assert!((p.z() - 3.0).abs() < 1e-10);
//! ```
//!
//! # References
//!
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Round_point>
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Sphere>

mod types;

#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra;

pub use types::*;
