//! Metric signature trait for Clifford algebras.
//!
//! This module provides the [`Signature`] trait which encodes the metric signature
//! of a geometric algebra as the triple `(P, Q, R)`:
//!
//! - `P`: number of basis vectors squaring to `+1`
//! - `Q`: number of basis vectors squaring to `-1`
//! - `R`: number of basis vectors squaring to `0` (degenerate/null)
//!
//! # Notation
//!
//! A Clifford algebra with signature `(P, Q, R)` is commonly written as:
//! - `Cl(P, Q)` when `R = 0`
//! - `Cl(P, Q, R)` for the general case
//!
//! The total dimension is `N = P + Q + R`, and the algebra has `2^N` basis blades.
//!
//! # Common Signatures
//!
//! | Algebra | Signature | Description |
//! |---------|-----------|-------------|
//! | Euclidean 2D | `Cl(2,0,0)` | Standard 2D plane |
//! | Euclidean 3D | `Cl(3,0,0)` | Standard 3D space |
//! | Minkowski | `Cl(1,3,0)` | Special relativity |
//! | PGA 3D | `Cl(3,0,1)` | Projective geometry |
//! | CGA 3D | `Cl(4,1,0)` | Conformal geometry |

use generic_array::ArrayLength;
use typenum::Unsigned;

/// Trait defining the metric signature of a Clifford algebra.
///
/// The signature determines how basis vectors square under the geometric product:
/// - Positive basis vectors: `e_i² = +1` (for `i < P`)
/// - Negative basis vectors: `e_i² = -1` (for `P ≤ i < P + Q`)
/// - Null basis vectors: `e_i² = 0` (for `P + Q ≤ i < P + Q + R`)
///
/// # Type-Level Constants
///
/// The signature is encoded as compile-time constants, enabling the compiler
/// to optimize based on the specific algebra being used. The [`NumBlades`](Self::NumBlades)
/// associated type enables zero-overhead generic array sizing.
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
/// struct Minkowski;
///
/// impl Signature for Minkowski {
///     type NumBlades = typenum::U16; // 2^4 = 16 blades
///
///     const P: usize = 1;
///     const Q: usize = 3;
///     const R: usize = 0;
///
///     fn metric(i: usize) -> i8 {
///         if i == 0 { 1 } else { -1 }
///     }
/// }
///
/// assert_eq!(Minkowski::DIM, 4);
/// assert_eq!(Minkowski::num_blades(), 16);
/// assert_eq!(Minkowski::metric(0), 1);  // time-like
/// assert_eq!(Minkowski::metric(1), -1); // space-like
/// ```
pub trait Signature: Copy + Clone + Default + 'static {
    /// Type-level number of basis blades for generic array sizing.
    ///
    /// This should be `typenum::U{2^DIM}`, e.g., `U4` for 2D, `U8` for 3D, `U16` for 4D.
    /// This enables [`Multivector`](crate::algebra::Multivector) to allocate exactly
    /// the right amount of space at compile time.
    type NumBlades: ArrayLength + Unsigned;

    /// Number of basis vectors squaring to `+1`.
    const P: usize;

    /// Number of basis vectors squaring to `-1`.
    const Q: usize;

    /// Number of basis vectors squaring to `0` (degenerate/null).
    const R: usize;

    /// Total dimension of the vector space: `P + Q + R`.
    const DIM: usize = Self::P + Self::Q + Self::R;

    /// Returns the number of basis blades in the algebra: `2^DIM`.
    ///
    /// This is a convenience method that returns `NumBlades::USIZE`.
    #[inline]
    fn num_blades() -> usize {
        Self::NumBlades::USIZE
    }

    /// Returns the metric coefficient for basis vector `e_i`.
    ///
    /// # Arguments
    ///
    /// * `i` - Index of the basis vector (0-indexed)
    ///
    /// # Returns
    ///
    /// - `+1` if `e_i² = +1` (positive/space-like)
    /// - `-1` if `e_i² = -1` (negative/time-like)
    /// - `0` if `e_i² = 0` (null/degenerate)
    ///
    /// # Panics
    ///
    /// Implementations should panic if `i >= DIM`.
    fn metric(i: usize) -> i8;
}
