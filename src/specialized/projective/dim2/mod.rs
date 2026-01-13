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
//! | 2 | [`Line`] | e₁₂, e₂₀, e₀₁ | Line |
//! | 0+2 | [`Motor`] | s, e₁₂, e₂₀, e₀₁ | Rigid transformation |
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
//!
//! // Points
//! let origin = Point::origin();
//! let p = Point::new(3.0, 4.0);
//!
//! // Line through two points
//! let line = origin.join(&p);
//!
//! // 90° rotation around origin
//! let rotor = Motor::from_rotation(FRAC_PI_2);
//! let rotated = rotor.transform_point(&p);
//!
//! // Translation
//! let translation = Motor::from_translation(1.0, 2.0);
//! let translated = translation.transform_point(&p);
//! ```

mod conversions;
mod ops;
mod types;

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra;

#[cfg(feature = "rerun-0_28")]
mod rerun;

pub use types::{Line, Motor, Point};
