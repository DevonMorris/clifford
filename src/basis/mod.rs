//! Basis blade representation and operations.
//!
//! This module provides utilities for working with basis blades, the fundamental
//! building blocks of geometric algebra. Every multivector can be expressed as
//! a linear combination of basis blades.
//!
//! # Blade Indexing
//!
//! Basis blades are indexed using a bitmask representation. For a blade index `k`,
//! bit `i` is set if and only if basis vector `e_i` is part of the blade.
//!
//! For example, in 3D:
//! - `0b000 = 0` → scalar `1`
//! - `0b001 = 1` → `e₁`
//! - `0b010 = 2` → `e₂`
//! - `0b011 = 3` → `e₁₂`
//! - `0b100 = 4` → `e₃`
//! - `0b101 = 5` → `e₁₃`
//! - `0b110 = 6` → `e₂₃`
//! - `0b111 = 7` → `e₁₂₃`
//!
//! # Grade
//!
//! The grade of a blade is the number of basis vectors it contains. This equals
//! the population count (number of set bits) of its index.
//!
//! | Grade | Name | Example (3D) |
//! |-------|------|--------------|
//! | 0 | Scalar | `1` |
//! | 1 | Vector | `e₁`, `e₂`, `e₃` |
//! | 2 | Bivector | `e₁₂`, `e₂₃`, `e₁₃` |
//! | 3 | Trivector | `e₁₂₃` |
//!
//! # Blade Products
//!
//! The [`basis_product`] function computes the geometric product of two basis blades,
//! returning both the result index and the sign. The sign accounts for:
//!
//! 1. **Swaps**: Reordering basis vectors to canonical order introduces `-1` for each swap
//! 2. **Metric**: Contracting matching vectors applies the metric coefficient
//!
//! # Example
//!
//! ```
//! use clifford::basis::{grade_of_blade, basis_product, blades_of_grade};
//!
//! // Grade is population count
//! assert_eq!(grade_of_blade(0b111), 3); // e₁₂₃ has grade 3
//!
//! // Blade count by grade in 3D
//! assert_eq!(blades_of_grade(3, 2), 3); // 3 bivectors
//!
//! // Geometric product with Euclidean metric
//! let euclidean = |_| 1i8;
//! let (sign, result) = basis_product(0b001, 0b010, euclidean);
//! assert_eq!(sign, 1);   // e₁ * e₂ = +e₁₂
//! assert_eq!(result, 0b011);
//! ```

mod blade;
mod index;

pub use blade::{anticommutes, basis_product};
pub use index::{binomial, blades_of_grade, grade_of_blade, grade_start_index};
