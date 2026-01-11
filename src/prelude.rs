//! Convenient re-exports for common use.
//!
//! This module provides a "prelude" that re-exports the most commonly used
//! types and traits, allowing users to import everything they need with a
//! single `use` statement.
//!
//! # Example
//!
//! ```
//! use clifford::prelude::*;
//!
//! // Now you have access to common types and traits
//! let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
//! ```

pub use crate::algebra::Multivector;
pub use crate::basis::Blade;
pub use crate::scalar::Float;
pub use crate::signature::{
    Conformal2, Conformal3, Euclidean2, Euclidean3, Euclidean4, Projective2, Projective3, Signature,
};

// Re-export approx traits for approximate comparisons
pub use approx::{AbsDiffEq, RelativeEq, UlpsEq};
pub use approx::{abs_diff_eq, abs_diff_ne, relative_eq, relative_ne, ulps_eq, ulps_ne};
