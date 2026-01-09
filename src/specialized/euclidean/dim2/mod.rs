//! 2D Euclidean Geometric Algebra types.
//!
//! This module provides specialized types for 2D Euclidean space `Cl(2,0,0)`.
//!
//! # Types by Grade
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 0 | scalar | 1 | Magnitude/weight |
//! | 1 | [`Vec2`] | e₁, e₂ | Direction/position |
//! | 2 | [`Bivec2`] | e₁₂ | Oriented area (pseudoscalar) |
//!
//! # Complex Number Analogy
//!
//! The 2D even subalgebra (scalar + bivector) is isomorphic to complex numbers:
//! - Scalar part → real part
//! - Bivector `e₁₂` → imaginary unit `i` (since `e₁₂² = -1`)
//!
//! A [`Rotor2`] is essentially a complex number of unit magnitude.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim2::{Vec2, Rotor2};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // 90° rotation
//! let rotor = Rotor2::from_angle(FRAC_PI_2);
//! let v = Vec2::new(1.0, 0.0);
//! let rotated = rotor.rotate(v);
//!
//! assert!((rotated.y - 1.0).abs() < 1e-10);
//! ```

mod conversions;
mod ops;
mod types;

#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;

pub use conversions::{ConversionError, Specialized2};
pub use types::*;
