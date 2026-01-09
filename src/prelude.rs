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
pub use crate::signature::{Euclidean2, Euclidean3, Euclidean4, Signature};

// Re-export Unsigned trait from typenum so users can access NumBlades::USIZE
// without importing typenum directly.
pub use typenum::Unsigned;
