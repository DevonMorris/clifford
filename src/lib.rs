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
//! - [`scalar`] - Floating-point scalar type abstraction
//! - [`signature`] - Metric signatures defining the algebra
//! - [`basis`] - Basis blade indexing and grade utilities
//!
//! # Getting Started
//!
//! ```
//! use clifford::scalar::Float;
//!
//! // The Float trait provides generic floating-point operations
//! fn magnitude<T: Float>(x: T, y: T) -> T {
//!     (x * x + y * y).sqrt()
//! }
//!
//! let m: f64 = magnitude(3.0, 4.0);
//! assert!((m - 5.0).abs() < f64::EPSILON);
//! ```

pub mod basis;
pub mod scalar;
pub mod signature;
