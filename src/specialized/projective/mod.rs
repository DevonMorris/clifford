//! Specialized types for Projective Geometric Algebra (PGA).
//!
//! This module provides optimized, ergonomic types for point-based PGA organized
//! by dimension. Currently supported:
//!
//! - **[`dim2`]**: 2D PGA `Cl(2,0,1)`: Point, Line, Motor
//! - **[`dim3`]**: 3D PGA `Cl(3,0,1)`: Point, Line, Plane, Motor, Flector
//!
//! # From Homogeneous Coordinates to PGA
//!
//! If you know homogeneous coordinates from graphics, PGA will feel familiar:
//!
//! | Concept | Homogeneous | PGA |
//! |---------|-------------|-----|
//! | 3D point | `(x, y, z, w)` | Grade-1 blade (vector) |
//! | Finite point | `w ≠ 0` | `e₁x + e₂y + e₃z + e₀w` |
//! | Point at infinity | `w = 0` | `e₁x + e₂y + e₃z` (no e₀) |
//! | Normalization | Divide by w | Weight norm = 1 |
//!
//! **Key insight**: In homogeneous coords, `w` is just a coordinate. In PGA,
//! it comes from a **degenerate basis vector** `e₀` where `e₀² = 0`.
//!
//! # The Degenerate Basis Vector
//!
//! PGA uses signature `Cl(n, 0, 1)`: n Euclidean basis vectors plus one
//! **degenerate** basis vector `e₀` that squares to zero:
//!
//! ```text
//! e₁² = +1   (Euclidean)
//! e₂² = +1   (Euclidean)
//! e₃² = +1   (Euclidean)
//! e₀² = 0    (Degenerate!)
//! ```
//!
//! This `e₀² = 0` is what makes projective geometry work:
//! - **Finite points**: have an `e₀` component (weight ≠ 0)
//! - **Ideal points** (at infinity): no `e₀` component (weight = 0)
//! - **Parallel lines**: meet at an ideal point
//!
//! # Geometric Objects as Blades
//!
//! In PGA, geometric objects are represented by blades of different grades:
//!
//! | Grade | 2D PGA | 3D PGA |
//! |-------|--------|--------|
//! | 1 | Point | Point |
//! | 2 | Line | Line |
//! | 3 | (pseudoscalar) | Plane |
//! | 4 | — | (pseudoscalar) |
//!
//! Higher-grade objects are built from lower-grade ones using the **wedge product**
//! (see Join and Meet in [`dim3`]).
//!
//! # Motors: 4×4 Matrices → 8 Numbers
//!
//! A **motor** is a PGA element that represents rigid body transforms
//! (rotation + translation). Motors unify what you'd do with 4×4 matrices:
//!
//! | Linear Algebra | PGA |
//! |----------------|-----|
//! | 4×4 homogeneous matrix | Motor (8 components in 3D) |
//! | Matrix multiply `M₂ × M₁` | Geometric product `M₂ * M₁` |
//! | Apply: `M × v` | Sandwich: `M v M̃` |
//! | Interpolation: complex | SLERP: natural |
//! | Drift: orthogonalize | Normalize: `M / |M|` |
//!
//! Motors are more compact (8 vs 16 numbers), don't drift as easily,
//! and interpolate cleanly.
//!
//! # Point-Based Formulation
//!
//! This library uses the **point-based** formulation of PGA where:
//! - **Grade 1 (vectors)**: Points in homogeneous coordinates
//! - **Grade 2 (bivectors)**: Lines
//! - **Grade 3 (trivectors)**: Planes (in 3D only)
//!
//! This is consistent with Euclidean GA conventions and nalgebra's `Point` types.
//!
//! **Note**: Some resources use a "plane-based" formulation where planes are grade 1.
//! The math is isomorphic but formulas differ.
//!
//! # Type Organization
//!
//! ```text
//! projective/
//!   dim2/   - 2D PGA: Point, Line, Motor
//!   dim3/   - 3D PGA: Point, Line, Plane, Motor, Flector
//! ```
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
//! // Motor: 90° rotation
//! let motor = Motor::from_rotation(FRAC_PI_2);
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
