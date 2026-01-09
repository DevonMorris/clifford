//! Clifford - A Rust library for Geometric Algebra (Clifford Algebra).
//!
//! This library provides types and operations for working with geometric algebra,
//! a mathematical framework that unifies and extends linear algebra, complex numbers,
//! quaternions, and more.
//!
//! # Overview
//!
//! Geometric Algebra (GA), also known as Clifford Algebra, provides a unified
//! mathematical language for geometry. It extends traditional vector algebra with
//! the geometric product, which combines the dot product and wedge product into
//! a single operation.
//!
//! # Module Structure
//!
//! - [`algebra`] - Core algebraic types ([`Multivector`](algebra::Multivector))
//! - [`scalar`] - Floating-point scalar type abstraction
//! - [`signature`] - Metric signatures defining the algebra
//! - [`basis`] - Basis blade representation and utilities
//!
//! # Getting Started
//!
//! ```
//! use clifford::algebra::Multivector;
//! use clifford::signature::Euclidean3;
//!
//! // Create two vectors
//! let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 0.0]);
//! let b: Multivector<f64, Euclidean3> = Multivector::vector(&[0.0, 1.0, 0.0]);
//!
//! // Geometric product combines dot and wedge products
//! let ab = &a * &b;
//!
//! // Dot product is the scalar part: a·b = 1*0 + 2*1 + 0*0 = 2
//! assert!((ab.scalar_part() - 2.0).abs() < 1e-10);
//! ```

pub mod algebra;
pub mod basis;
pub mod prelude;
pub mod scalar;
pub mod signature;
pub mod specialized;

/// Test utilities available only during testing.
#[cfg(test)]
pub(crate) mod test_utils {
    /// Standard epsilon for absolute difference comparisons in tests.
    ///
    /// Use this constant instead of magic numbers like `1e-10` or `1e-9`.
    /// This value is chosen to be strict enough for most operations while
    /// allowing for reasonable floating-point accumulation in compound operations.
    pub const ABS_DIFF_EQ_EPS: f64 = 1e-10;
}
