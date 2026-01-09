//! Basis blades and grade utilities for geometric algebra.
//!
//! This module provides the [`Blade`] type and utilities for working with
//! basis blades, the fundamental building blocks of geometric algebra.
//!
//! # What is a Blade?
//!
//! A **blade** is a geometric object that represents an oriented subspace.
//! Every multivector (the general elements of geometric algebra) can be
//! written as a sum of blades.
//!
//! | Grade | Name | Represents | Example |
//! |-------|------|------------|---------|
//! | 0 | Scalar | Pure magnitude | `1`, `2.5` |
//! | 1 | Vector | Directed line | `e₁`, `3e₂` |
//! | 2 | Bivector | Oriented plane | `e₁₂`, `e₁∧e₂` |
//! | 3 | Trivector | Oriented volume | `e₁₂₃` |
//! | k | k-vector | k-dimensional subspace | ... |
//!
//! # The Geometric Product
//!
//! Blades multiply using the **geometric product**, which unifies:
//!
//! - **Dot product** (contraction): Measures alignment, reduces grade
//! - **Wedge product** (extension): Measures perpendicularity, increases grade
//!
//! For two vectors `a` and `b`:
//! ```text
//! ab = a·b + a∧b
//!    = |a||b|cos(θ) + |a||b|sin(θ)·B
//! ```
//! where `B` is the unit bivector of the plane containing `a` and `b`.
//!
//! # Key Insight: Why This Matters
//!
//! The geometric product encodes both "how parallel" (dot) and "how perpendicular"
//! (wedge) two vectors are in a single algebraic operation. This is why geometric
//! algebra is so powerful for geometry—rotations, reflections, and projections
//! all emerge naturally from this single product.
//!
//! # Representation
//!
//! We represent basis blades using a bitmask index. For a blade index `k`,
//! bit `i` is set if and only if basis vector `eᵢ` is present in the blade.
//!
//! ```text
//! Index  Binary  Blade   Grade   Meaning
//! ─────  ──────  ──────  ─────   ────────────────
//!   0    000     1       0       Scalar
//!   1    001     e₁      1       x-direction
//!   2    010     e₂      1       y-direction
//!   3    011     e₁₂     2       xy-plane
//!   4    100     e₃      1       z-direction
//!   5    101     e₁₃     2       xz-plane
//!   6    110     e₂₃     2       yz-plane
//!   7    111     e₁₂₃    3       Pseudoscalar (volume)
//! ```
//!
//! # Example
//!
//! ```
//! use clifford::basis::Blade;
//!
//! // Create basis vectors
//! let e1 = Blade::basis_vector(0);
//! let e2 = Blade::basis_vector(1);
//!
//! // Euclidean metric: all basis vectors square to +1
//! let euclidean = |_: usize| 1i8;
//!
//! // Vector squared gives a scalar: e₁² = 1
//! let (sign, result) = e1.product(&e1, euclidean);
//! assert_eq!(sign, 1);
//! assert_eq!(result, Blade::scalar());
//!
//! // Different vectors give a bivector: e₁e₂ = e₁₂
//! let (sign, result) = e1.product(&e2, euclidean);
//! assert_eq!(sign, 1);
//! assert_eq!(result.grade(), 2);
//!
//! // Order matters! e₂e₁ = -e₁₂
//! let (sign, result) = e2.product(&e1, euclidean);
//! assert_eq!(sign, -1);
//! println!("e₂e₁ = -{}", result); // prints: e₂e₁ = -e₁₂
//! ```

mod blade;
mod index;

pub use blade::{anticommutes, basis_product, Blade};
pub use index::{binomial, blades_of_grade, grade_of_blade, grade_start_index};
