//! 2D Euclidean Geometric Algebra types.
//!
//! This module provides specialized types for 2D Euclidean space `Cl(2,0,0)`.
//!
//! # Types by Grade
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 0 | scalar | 1 | Magnitude/weight |
//! | 1 | [`Vector`] | e₁, e₂ | Direction/position |
//! | 2 | [`Bivector`] | e₁₂ | Oriented area (pseudoscalar) |
//!
//! # Complex Number Analogy
//!
//! The 2D even subalgebra (scalar + bivector) is isomorphic to complex numbers:
//! - Scalar part → real part
//! - Bivector `e₁₂` → imaginary unit `i` (since `e₁₂² = -1`)
//!
//! A [`Rotor`] is essentially a complex number of unit magnitude.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim2::{Vector, Rotor};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // 90° rotation
//! let rotor = Rotor::from_angle(FRAC_PI_2);
//! let v = Vector::new(1.0, 0.0);
//! let rotated = rotor.rotate(v);
//!
//! assert!((rotated.y() - 1.0).abs() < 1e-10);
//! ```

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions
mod extensions;

#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra;

#[cfg(feature = "rerun-0_28")]
mod rerun;

// Re-export generated types
pub use generated::types::{Bivector, Rotor, Scalar, Vector};

// Re-export products module for direct access to algebraic products
pub use generated::products;

// Re-export Even alias from extensions
pub use extensions::Even;
