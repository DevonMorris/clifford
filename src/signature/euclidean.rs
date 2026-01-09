//! Euclidean signature types for standard geometric algebras.
//!
//! This module provides signature types for Euclidean spaces where all basis
//! vectors square to `+1`. These are the most commonly used signatures in
//! practical geometric algebra applications.
//!
//! # Available Signatures
//!
//! - [`Euclidean2`]: 2D plane, `Cl(2,0,0)` with 4 basis blades
//! - [`Euclidean3`]: 3D space, `Cl(3,0,0)` with 8 basis blades
//! - [`Euclidean4`]: 4D space, `Cl(4,0,0)` with 16 basis blades
//!
//! # Basis Blade Structure
//!
//! For Euclidean 3D (`Cl(3,0,0)`), the 8 basis blades are:
//!
//! | Index | Binary | Blade | Grade | Name |
//! |-------|--------|-------|-------|------|
//! | 0 | `000` | `1` | 0 | Scalar |
//! | 1 | `001` | `e₁` | 1 | Vector |
//! | 2 | `010` | `e₂` | 1 | Vector |
//! | 3 | `011` | `e₁₂` | 2 | Bivector |
//! | 4 | `100` | `e₃` | 1 | Vector |
//! | 5 | `101` | `e₁₃` | 2 | Bivector |
//! | 6 | `110` | `e₂₃` | 2 | Bivector |
//! | 7 | `111` | `e₁₂₃` | 3 | Pseudoscalar |

use super::Signature;

/// 2D Euclidean signature: `Cl(2,0,0)`.
///
/// Represents the standard 2D plane with basis vectors `e₁` and `e₂`.
///
/// # Basis Blades (4 total)
///
/// | Index | Blade | Grade | Description |
/// |-------|-------|-------|-------------|
/// | 0 | `1` | 0 | Scalar |
/// | 1 | `e₁` | 1 | x-direction vector |
/// | 2 | `e₂` | 1 | y-direction vector |
/// | 3 | `e₁₂` | 2 | Pseudoscalar (oriented area) |
///
/// # Properties
///
/// - `e₁² = e₂² = +1`
/// - `e₁₂² = e₁e₂e₁e₂ = -e₁e₁e₂e₂ = -1`
/// - The pseudoscalar `e₁₂` behaves like the imaginary unit `i`
///
/// # Example
///
/// ```
/// use clifford::signature::{Euclidean2, Signature};
///
/// assert_eq!(Euclidean2::DIM, 2);
/// assert_eq!(Euclidean2::NUM_BLADES, 4);
/// assert_eq!(Euclidean2::metric(0), 1);
/// assert_eq!(Euclidean2::metric(1), 1);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Euclidean2;

impl Signature for Euclidean2 {
    const P: usize = 2;
    const Q: usize = 0;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Euclidean2");
        1
    }
}

/// 3D Euclidean signature: `Cl(3,0,0)`.
///
/// Represents standard 3D space with basis vectors `e₁`, `e₂`, and `e₃`.
/// This is the most commonly used geometric algebra for 3D graphics,
/// physics simulations, and robotics.
///
/// # Basis Blades (8 total)
///
/// | Index | Blade | Grade | Description |
/// |-------|-------|-------|-------------|
/// | 0 | `1` | 0 | Scalar |
/// | 1 | `e₁` | 1 | x-direction vector |
/// | 2 | `e₂` | 1 | y-direction vector |
/// | 3 | `e₁₂` | 2 | xy-plane bivector |
/// | 4 | `e₃` | 1 | z-direction vector |
/// | 5 | `e₁₃` | 2 | xz-plane bivector |
/// | 6 | `e₂₃` | 2 | yz-plane bivector |
/// | 7 | `e₁₂₃` | 3 | Pseudoscalar (oriented volume) |
///
/// # Bivectors and Rotation
///
/// The three bivectors (`e₁₂`, `e₁₃`, `e₂₃`) represent oriented planes.
/// They are dual to the three axes and are used to construct rotors
/// (rotation operators). Each bivector squares to `-1`:
///
/// - `e₁₂² = -1` (rotation in xy-plane, around z-axis)
/// - `e₁₃² = -1` (rotation in xz-plane, around y-axis)
/// - `e₂₃² = -1` (rotation in yz-plane, around x-axis)
///
/// # Example
///
/// ```
/// use clifford::signature::{Euclidean3, Signature};
///
/// assert_eq!(Euclidean3::DIM, 3);
/// assert_eq!(Euclidean3::NUM_BLADES, 8);
///
/// // All basis vectors square to +1
/// for i in 0..3 {
///     assert_eq!(Euclidean3::metric(i), 1);
/// }
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Euclidean3;

impl Signature for Euclidean3 {
    const P: usize = 3;
    const Q: usize = 0;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Euclidean3");
        1
    }
}

/// 4D Euclidean signature: `Cl(4,0,0)`.
///
/// Represents 4D Euclidean space with basis vectors `e₁`, `e₂`, `e₃`, and `e₄`.
///
/// # Basis Blades (16 total)
///
/// | Grade | Count | Description |
/// |-------|-------|-------------|
/// | 0 | 1 | Scalar |
/// | 1 | 4 | Vectors |
/// | 2 | 6 | Bivectors |
/// | 3 | 4 | Trivectors |
/// | 4 | 1 | Pseudoscalar (4-volume) |
///
/// # Applications
///
/// While less common than 3D, 4D Euclidean GA is useful for:
/// - Conformal model preprocessing
/// - Quaternion-like operations with additional structure
/// - Higher-dimensional geometry research
///
/// # Example
///
/// ```
/// use clifford::signature::{Euclidean4, Signature};
///
/// assert_eq!(Euclidean4::DIM, 4);
/// assert_eq!(Euclidean4::NUM_BLADES, 16);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Euclidean4;

impl Signature for Euclidean4 {
    const P: usize = 4;
    const Q: usize = 0;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Euclidean4");
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn euclidean2_metric_always_positive(i in 0usize..2) {
            prop_assert_eq!(Euclidean2::metric(i), 1);
        }

        #[test]
        fn euclidean3_metric_always_positive(i in 0usize..3) {
            prop_assert_eq!(Euclidean3::metric(i), 1);
        }

        #[test]
        fn euclidean4_metric_always_positive(i in 0usize..4) {
            prop_assert_eq!(Euclidean4::metric(i), 1);
        }
    }

    #[test]
    fn euclidean_dimensions() {
        assert_eq!(Euclidean2::DIM, 2);
        assert_eq!(Euclidean2::NUM_BLADES, 4);

        assert_eq!(Euclidean3::DIM, 3);
        assert_eq!(Euclidean3::NUM_BLADES, 8);

        assert_eq!(Euclidean4::DIM, 4);
        assert_eq!(Euclidean4::NUM_BLADES, 16);
    }

    #[test]
    fn euclidean_pqr_consistency() {
        assert_eq!(Euclidean2::P + Euclidean2::Q + Euclidean2::R, Euclidean2::DIM);
        assert_eq!(Euclidean3::P + Euclidean3::Q + Euclidean3::R, Euclidean3::DIM);
        assert_eq!(Euclidean4::P + Euclidean4::Q + Euclidean4::R, Euclidean4::DIM);
    }
}
