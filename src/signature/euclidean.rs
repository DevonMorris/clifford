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
use typenum::{U4, U8, U16};

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
/// use clifford::prelude::*;
///
/// assert_eq!(Euclidean2::DIM, 2);
/// assert_eq!(Euclidean2::num_blades(), 4);
/// assert_eq!(Euclidean2::metric(0), 1);
/// assert_eq!(Euclidean2::metric(1), 1);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Euclidean2;

impl Signature for Euclidean2 {
    type NumBlades = U4; // 2^2 = 4

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
/// use clifford::prelude::*;
///
/// assert_eq!(Euclidean3::DIM, 3);
/// assert_eq!(Euclidean3::num_blades(), 8);
///
/// // All basis vectors square to +1
/// for i in 0..3 {
///     assert_eq!(Euclidean3::metric(i), 1);
/// }
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Euclidean3;

impl Signature for Euclidean3 {
    type NumBlades = U8; // 2^3 = 8

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
/// use clifford::prelude::*;
///
/// assert_eq!(Euclidean4::DIM, 4);
/// assert_eq!(Euclidean4::num_blades(), 16);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Euclidean4;

impl Signature for Euclidean4 {
    type NumBlades = U16; // 2^4 = 16

    const P: usize = 4;
    const Q: usize = 0;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Euclidean4");
        1
    }
}

/// 1D Hyperbolic signature: `Cl(1,0,0)`.
///
/// This is used for hyperbolic numbers (split-complex numbers).
/// With one positive basis vector `e₁` where `e₁² = +1`, this
/// algebra produces hyperbolic numbers `a + bj` where `j² = +1`.
///
/// # Basis Blades (2 total)
///
/// | Index | Blade | Grade | Description |
/// |-------|-------|-------|-------------|
/// | 0 | `1` | 0 | Scalar (real part) |
/// | 1 | `e₁` | 1 | Hyperbolic unit (j) |
///
/// # Properties
///
/// - `e₁² = +1` (unlike complex numbers where `i² = -1`)
/// - Hyperbolic numbers have zero divisors: `(1+j)(1-j) = 0`
/// - The norm uses grade involution: `z * involute(z) = a² - b²`
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Cl1_0_0::DIM, 1);
/// assert_eq!(Cl1_0_0::num_blades(), 2);
/// assert_eq!(Cl1_0_0::metric(0), 1);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Cl1_0_0;

impl Signature for Cl1_0_0 {
    type NumBlades = typenum::U2; // 2^1 = 2

    const P: usize = 1;
    const Q: usize = 0;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Cl1_0_0");
        1
    }
}

/// Type alias for the hyperbolic numbers signature.
pub type Hyperbolic1 = Cl1_0_0;

/// 1D Anti-Euclidean signature: `Cl(0,1,0)`.
///
/// This is used for complex numbers.
/// With one negative basis vector `e₁` where `e₁² = -1`, this
/// algebra produces complex numbers `a + bi` where `i² = -1`.
///
/// # Basis Blades (2 total)
///
/// | Index | Blade | Grade | Description |
/// |-------|-------|-------|-------------|
/// | 0 | `1` | 0 | Scalar (real part) |
/// | 1 | `e₁` | 1 | Imaginary unit (i) |
///
/// # Properties
///
/// - `e₁² = -1` (like complex numbers)
/// - Complex numbers have no zero divisors
/// - The norm uses Clifford conjugate: `z * conjugate(z) = a² + b²`
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Cl0_1_0::DIM, 1);
/// assert_eq!(Cl0_1_0::num_blades(), 2);
/// assert_eq!(Cl0_1_0::metric(0), -1);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Cl0_1_0;

impl Signature for Cl0_1_0 {
    type NumBlades = typenum::U2; // 2^1 = 2

    const P: usize = 0;
    const Q: usize = 1;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Cl0_1_0");
        -1
    }
}

/// Type alias for the complex numbers signature.
pub type Complex1 = Cl0_1_0;

/// 2D Negative-definite signature: `Cl(0,2,0)`.
///
/// This is the quaternion algebra. With two basis vectors `e₁, e₂`
/// where `e₁² = e₂² = -1`, this gives the standard quaternion algebra
/// with basis mapping: `i → e₁`, `j → e₂`, `k → e₁e₂`.
///
/// # Basis Blades (4 total)
///
/// | Index | Blade | Grade | Description |
/// |-------|-------|-------|-------------|
/// | 0 | `1` | 0 | Scalar (real part, w) |
/// | 1 | `e₁` | 1 | Imaginary i |
/// | 2 | `e₂` | 1 | Imaginary j |
/// | 3 | `e₁₂` | 2 | Imaginary k = ij |
///
/// # Properties
///
/// - `e₁² = e₂² = -1`
/// - `e₁e₂ = -e₂e₁` (anti-commutative)
/// - Quaternions have no zero divisors (division algebra)
/// - The norm uses Clifford conjugate: `q * conjugate(q) = w² + x² + y² + z²`
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Cl0_2_0::DIM, 2);
/// assert_eq!(Cl0_2_0::num_blades(), 4);
/// assert_eq!(Cl0_2_0::metric(0), -1);
/// assert_eq!(Cl0_2_0::metric(1), -1);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Cl0_2_0;

impl Signature for Cl0_2_0 {
    type NumBlades = typenum::U4; // 2^2 = 4

    const P: usize = 0;
    const Q: usize = 2;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Cl0_2_0");
        -1
    }
}

/// Type alias for the quaternion signature.
pub type Quaternion2 = Cl0_2_0;

/// 2D Indefinite signature: `Cl(1,1,0)`.
///
/// This is the Minkowski plane (2D spacetime algebra). With one
/// spacelike basis vector `e₁` (e₁² = +1) and one timelike basis
/// vector `e₂` (e₂² = -1), this gives a 2D relativistic algebra.
///
/// # Basis Blades (4 total)
///
/// | Index | Blade | Grade | Description |
/// |-------|-------|-------|-------------|
/// | 0 | `1` | 0 | Scalar |
/// | 1 | `e₁` | 1 | Spacelike vector |
/// | 2 | `e₂` | 1 | Timelike vector |
/// | 3 | `e₁₂` | 2 | Bivector (pseudoscalar) |
///
/// # Properties
///
/// - `e₁² = +1` (spacelike)
/// - `e₂² = -1` (timelike)
/// - Indefinite norm: v² can be positive, negative, or zero
/// - Null vectors exist: e.g., `e₁ + e₂` has v² = 0
/// - Even subalgebra is isomorphic to hyperbolic numbers
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Cl1_1_0::DIM, 2);
/// assert_eq!(Cl1_1_0::num_blades(), 4);
/// assert_eq!(Cl1_1_0::metric(0), 1);   // e1² = +1
/// assert_eq!(Cl1_1_0::metric(1), -1);  // e2² = -1
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Cl1_1_0;

impl Signature for Cl1_1_0 {
    type NumBlades = typenum::U4; // 2^2 = 4

    const P: usize = 1;
    const Q: usize = 1;
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Cl1_1_0");
        match i {
            0 => 1,  // e1² = +1 (spacelike)
            1 => -1, // e2² = -1 (timelike)
            _ => unreachable!(),
        }
    }
}

/// Type alias for the Minkowski plane signature.
pub type Minkowski2 = Cl1_1_0;

/// 1D Degenerate signature: `Cl(0,0,1)`.
///
/// This is used for dual numbers (automatic differentiation).
/// With one null basis vector `e₁` where `e₁² = 0`, this
/// algebra produces dual numbers `a + bε` where `ε² = 0`.
///
/// # Basis Blades (2 total)
///
/// | Index | Blade | Grade | Description |
/// |-------|-------|-------|-------------|
/// | 0 | `1` | 0 | Scalar (real part) |
/// | 1 | `e₁` | 1 | Dual unit (ε, nilpotent) |
///
/// # Properties
///
/// - `e₁² = 0` (nilpotent)
/// - Dual numbers have zero divisors: all `bε` have zero norm
/// - The norm is degenerate: `z * involute(z) = a²`
/// - Foundation for forward-mode automatic differentiation
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Cl0_0_1::DIM, 1);
/// assert_eq!(Cl0_0_1::num_blades(), 2);
/// assert_eq!(Cl0_0_1::metric(0), 0);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Cl0_0_1;

impl Signature for Cl0_0_1 {
    type NumBlades = typenum::U2; // 2^1 = 2

    const P: usize = 0;
    const Q: usize = 0;
    const R: usize = 1;

    #[inline]
    fn metric(i: usize) -> i8 {
        debug_assert!(i < Self::DIM, "basis index {i} out of range for Cl0_0_1");
        0
    }
}

/// Type alias for the dual numbers signature.
pub type Dual1 = Cl0_0_1;

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
        assert_eq!(Euclidean2::num_blades(), 4);

        assert_eq!(Euclidean3::DIM, 3);
        assert_eq!(Euclidean3::num_blades(), 8);

        assert_eq!(Euclidean4::DIM, 4);
        assert_eq!(Euclidean4::num_blades(), 16);
    }

    #[test]
    fn euclidean_pqr_consistency() {
        assert_eq!(
            Euclidean2::P + Euclidean2::Q + Euclidean2::R,
            Euclidean2::DIM
        );
        assert_eq!(
            Euclidean3::P + Euclidean3::Q + Euclidean3::R,
            Euclidean3::DIM
        );
        assert_eq!(
            Euclidean4::P + Euclidean4::Q + Euclidean4::R,
            Euclidean4::DIM
        );
    }
}
