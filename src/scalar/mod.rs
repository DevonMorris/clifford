//! Scalar type abstractions for geometric algebra.
//!
//! This module provides the [`Float`] trait which abstracts over floating-point
//! scalar types, enabling the library to be generic over numeric precision.
//!
//! # Overview
//!
//! In geometric algebra, scalars form the grade-0 elements of the algebra.
//! All higher-grade elements (vectors, bivectors, pseudoscalars, etc.) have
//! coefficients that are scalars. This module defines the requirements for
//! those scalar types.
//!
//! # Supported Types
//!
//! The [`Float`] trait is implemented for:
//! - [`f32`] - Single precision (32-bit)
//! - [`f64`] - Double precision (64-bit)
//!
//! # Example
//!
//! ```
//! use clifford::scalar::Float;
//!
//! fn dot_product<T: Float>(a: &[T], b: &[T]) -> T {
//!     a.iter()
//!         .zip(b.iter())
//!         .map(|(&x, &y)| x * y)
//!         .fold(T::ZERO, |acc, x| acc + x)
//! }
//!
//! let a = [1.0_f64, 2.0, 3.0];
//! let b = [4.0_f64, 5.0, 6.0];
//! assert!((dot_product(&a, &b) - 32.0).abs() < f64::EPSILON);
//! ```

mod float;

pub use float::Float;
