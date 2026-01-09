//! Core algebraic types for geometric algebra.
//!
//! This module provides the [`Multivector`] type, the fundamental element of
//! geometric algebra that can represent scalars, vectors, bivectors, and
//! arbitrary combinations thereof.
//!
//! # The Multivector
//!
//! A multivector is a linear combination of basis blades:
//!
//! ```text
//! M = s + v₁e₁ + v₂e₂ + v₃e₃ + b₁e₁₂ + b₂e₂₃ + b₃e₁₃ + pe₁₂₃
//!     ─   ─────────────────   ───────────────────────   ──────
//!     │         │                      │                  │
//!   scalar   vector              bivector          pseudoscalar
//!   (grade 0) (grade 1)          (grade 2)          (grade 3)
//! ```
//!
//! # Operations
//!
//! The main operations on multivectors are:
//!
//! - **Geometric product** (`*`): The fundamental product combining dot and wedge
//! - **Addition/Subtraction** (`+`, `-`): Component-wise
//! - **Scalar multiplication** (`* scalar`): Scale all components
//! - **Reverse** (`.reverse()`): Reverses blade orientation
//! - **Inverse** (`.inverse()`): Multiplicative inverse (when it exists)
//!
//! # Example
//!
//! ```
//! use clifford::algebra::Multivector;
//! use clifford::signature::Euclidean3;
//!
//! // Create two vectors
//! let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 0.0, 0.0]);
//! let b: Multivector<f64, Euclidean3> = Multivector::vector(&[0.0, 1.0, 0.0]);
//!
//! // Geometric product: ab = a·b + a∧b = 0 + e₁₂
//! let ab = &a * &b;
//!
//! // The result is a pure bivector (oriented plane)
//! assert!(ab.scalar_part().abs() < 1e-10); // No scalar part (perpendicular)
//! ```

mod multivector;

pub use multivector::{MAX_BLADES, Multivector};
