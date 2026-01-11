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
//! | [`Dipole`] | 2 | Point pair (two points as bivector) |
//! | [`Circle`] | 3 | Circle (intersection of two spheres) |
//! | [`Line`] | 3 | Flat circle (circle through infinity) |
//! | [`Sphere`] | 4 | Sphere (or imaginary sphere) |
//! | [`Plane`] | 4 | Flat sphere (sphere through infinity) |
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
//! use clifford::specialized::conformal::dim3::{Point, Plane, Line};
//!
//! // Create a point
//! let p = Point::<f64>::new(1.0, 2.0, 3.0);
//!
//! // Verify null constraint
//! assert!(p.is_null(1e-10));
//!
//! // Create the xy-plane
//! let plane = Plane::<f64>::xy();
//! assert!(plane.contains(&Point::new(1.0, 2.0, 0.0), 1e-10));
//!
//! // Create the x-axis
//! let x_axis = Line::<f64>::x_axis();
//! assert!(x_axis.contains(&Point::new(5.0, 0.0, 0.0), 1e-10));
//! ```
//!
//! # References
//!
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Round_point>
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Sphere>
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Plane>
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Circle>
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Line>
//! - <https://conformalgeometricalgebra.org/wiki/index.php?title=Dipole>

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
