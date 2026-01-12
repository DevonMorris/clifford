//! Core blade algebra computations.
//!
//! This module provides the mathematical foundation for code generation:
//!
//! - [`Blade`]: Represents a basis blade with canonical ordering
//! - [`Algebra`]: Encapsulates a geometric algebra with its metric signature
//! - [`ProductTable`]: Precomputed product lookup for efficient code generation
//!
//! # Canonical Ordering
//!
//! Blades are represented as bitmasks where bit `i` indicates the presence
//! of basis vector `eᵢ`. This representation is inherently canonical:
//!
//! | Blade | Binary | Index |
//! |-------|--------|-------|
//! | 1     | `0000` | 0     |
//! | e₁    | `0001` | 1     |
//! | e₂    | `0010` | 2     |
//! | e₁₂   | `0011` | 3     |
//! | e₃    | `0100` | 4     |
//!
//! # Sign Computation
//!
//! When computing `eₐ * eᵦ`, the sign comes from:
//! 1. Swaps needed to reorder into canonical form
//! 2. Metric contributions when basis vectors square
//!
//! # Example
//!
//! ```
//! use clifford_codegen::algebra::{Algebra, Blade, ProductTable};
//!
//! let algebra = Algebra::euclidean(3);
//! let table = ProductTable::new(&algebra);
//!
//! // e1 * e2 = e12 with sign +1
//! let (sign, result) = table.geometric(1, 2);
//! assert_eq!(sign, 1);
//! assert_eq!(result, 3);
//!
//! // e2 * e1 = -e12 (anticommutative)
//! let (sign, result) = table.geometric(2, 1);
//! assert_eq!(sign, -1);
//! assert_eq!(result, 3);
//! ```

mod blade;
mod constraints;
mod grade;
#[cfg(test)]
mod proptest;
mod sign;
mod signature;
mod table;
mod versor;
#[cfg(test)]
mod verification;

pub use blade::Blade;
pub use constraints::{
    satisfies_all_constraints, satisfies_antiproduct_constraint, satisfies_geometric_constraint,
};
pub use grade::{
    antireverse_sign, binomial, blades_of_grade, blades_of_grades, geometric_grades, grade,
    inner_grade, left_contraction_grade, outer_grade, reverse_sign,
};
pub use sign::basis_product;
pub use signature::Algebra;
pub use table::ProductTable;
pub use versor::{VersorInfo, VersorParity, versor_parity};
