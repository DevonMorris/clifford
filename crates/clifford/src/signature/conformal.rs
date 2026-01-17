//! Conformal signature types for Conformal Geometric Algebra (CGA).
//!
//! This module provides signature types for Conformal Geometric Algebras,
//! which extend Euclidean space by embedding it into a higher-dimensional space
//! with one extra positive basis vector (`e₊²=+1`) and one negative basis vector
//! (`e₋²=-1`).
//!
//! # Null Basis Convention
//!
//! CGA uses a "null basis" where two special vectors are constructed:
//! - `e∞ = e₋ + e₊` (point at infinity)
//! - `e₀ = (e₋ - e₊) / 2` (origin)
//!
//! These have the key properties:
//! - `e∞² = 0`, `e₀² = 0` (both are null vectors)
//! - `e∞ · e₀ = -1`
//!
//! # Available Signatures
//!
//! - [`Conformal2`]: 2D CGA, `Cl(3,1,0)` with 16 basis blades
//! - [`Conformal3`]: 3D CGA, `Cl(4,1,0)` with 32 basis blades
//!
//! # Basis Vector Ordering
//!
//! Basis vectors are ordered with Euclidean vectors first, then conformal vectors:
//! - 2D CGA: `e₁, e₂, e₊, e₋` (indices 0, 1, 2, 3)
//! - 3D CGA: `e₁, e₂, e₃, e₊, e₋` (indices 0, 1, 2, 3, 4)
//!
//! # Reference
//!
//! Primary resource: <https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page>
//!
//! # Example
//!
//! ```
//! use clifford::prelude::*;
//!
//! // 3D Conformal GA
//! assert_eq!(Conformal3::DIM, 5);
//! assert_eq!(Conformal3::num_blades(), 32);
//!
//! // Euclidean basis vectors square to +1
//! assert_eq!(Conformal3::metric(0), 1);  // e₁² = 1
//! assert_eq!(Conformal3::metric(1), 1);  // e₂² = 1
//! assert_eq!(Conformal3::metric(2), 1);  // e₃² = 1
//!
//! // Conformal basis vectors
//! assert_eq!(Conformal3::metric(3), 1);  // e₊² = +1
//! assert_eq!(Conformal3::metric(4), -1); // e₋² = -1
//! ```

use super::Signature;
use crate::algebra::Multivector;
use crate::basis::Blade;
use crate::scalar::Float;
use typenum::{U16, U32};

/// 2D Conformal signature: `Cl(3,1,0)`.
///
/// Embeds 2D Euclidean space into a 4D conformal space with basis vectors
/// `e₁`, `e₂`, `e₊`, and `e₋`, where `e₊² = +1` and `e₋² = -1`.
///
/// # Null Basis
///
/// The null basis vectors are:
/// - `e∞ = e₋ + e₊` (point at infinity)
/// - `e₀ = (e₋ - e₊) / 2` (origin)
///
/// Use [`Conformal2::e_infinity()`] and [`Conformal2::e_origin()`] to create these.
///
/// # Conformal Embedding
///
/// A 2D Euclidean point `(x, y)` is embedded as:
/// ```text
/// P = x·e₁ + y·e₂ + e₀ + ½(x² + y²)·e∞
/// ```
///
/// # Basis Blades (16 total)
///
/// | Grade | Count | Description |
/// |-------|-------|-------------|
/// | 0 | 1 | Scalar |
/// | 1 | 4 | Vectors (including conformal) |
/// | 2 | 6 | Bivectors |
/// | 3 | 4 | Trivectors |
/// | 4 | 1 | Pseudoscalar |
///
/// # Properties
///
/// - `e₁² = e₂² = +1` (Euclidean)
/// - `e₊² = +1`, `e₋² = -1` (conformal)
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Conformal2::DIM, 4);
/// assert_eq!(Conformal2::num_blades(), 16);
///
/// // Euclidean basis vectors
/// assert_eq!(Conformal2::metric(0), 1);  // e₁² = 1
/// assert_eq!(Conformal2::metric(1), 1);  // e₂² = 1
///
/// // Conformal basis vectors
/// assert_eq!(Conformal2::metric(2), 1);  // e₊² = +1
/// assert_eq!(Conformal2::metric(3), -1); // e₋² = -1
/// ```
///
/// Reference: <https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page>
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Conformal2;

impl Signature for Conformal2 {
    type NumBlades = U16; // 2^4 = 16

    const P: usize = 3; // e₁, e₂, e₊
    const Q: usize = 1; // e₋
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        match i {
            0..=2 => 1, // e₁² = e₂² = e₊² = +1
            3 => -1,    // e₋² = -1
            _ => panic!("basis index {i} out of range for Conformal2 (DIM=4)"),
        }
    }
}

impl Conformal2 {
    /// Index of the positive conformal basis vector `e₊`.
    pub const E_PLUS: usize = 2;

    /// Index of the negative conformal basis vector `e₋`.
    pub const E_MINUS: usize = 3;

    /// Creates `e∞ = e₋ + e₊` (point at infinity).
    ///
    /// This is a null vector (`e∞² = 0`) representing the point at infinity.
    ///
    /// # Key Properties
    ///
    /// - `e∞ · e∞ = 0`
    /// - `e∞ · e₀ = -1`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::prelude::*;
    ///
    /// let e_inf = Conformal2::e_infinity::<f64>();
    /// let sq = &e_inf * &e_inf;
    /// assert!(sq.is_zero(1e-10));
    /// ```
    pub fn e_infinity<T: Float>() -> Multivector<T, Self> {
        let mut mv = Multivector::zero();
        mv.set(Blade::basis_vector(Self::E_PLUS), T::one());
        mv.set(Blade::basis_vector(Self::E_MINUS), T::one());
        mv
    }

    /// Creates `e₀ = (e₋ - e₊) / 2` (origin).
    ///
    /// This is a null vector (`e₀² = 0`) representing the origin.
    ///
    /// # Key Properties
    ///
    /// - `e₀ · e₀ = 0`
    /// - `e∞ · e₀ = -1`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::prelude::*;
    ///
    /// let e_o = Conformal2::e_origin::<f64>();
    /// let sq = &e_o * &e_o;
    /// assert!(sq.is_zero(1e-10));
    /// ```
    pub fn e_origin<T: Float>() -> Multivector<T, Self> {
        let mut mv = Multivector::zero();
        let half = T::one() / T::TWO;
        mv.set(Blade::basis_vector(Self::E_MINUS), half);
        mv.set(Blade::basis_vector(Self::E_PLUS), -half);
        mv
    }
}

/// 3D Conformal signature: `Cl(4,1,0)`.
///
/// Embeds 3D Euclidean space into a 5D conformal space with basis vectors
/// `e₁`, `e₂`, `e₃`, `e₊`, and `e₋`, where `e₊² = +1` and `e₋² = -1`.
///
/// # Null Basis
///
/// The null basis vectors are:
/// - `e∞ = e₋ + e₊` (point at infinity)
/// - `e₀ = (e₋ - e₊) / 2` (origin)
///
/// Use [`Conformal3::e_infinity()`] and [`Conformal3::e_origin()`] to create these.
///
/// # Conformal Embedding
///
/// A 3D Euclidean point `(x, y, z)` is embedded as:
/// ```text
/// P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
/// ```
///
/// This embedding has the property that the inner product of two conformal
/// points gives (negative) twice the squared Euclidean distance:
/// ```text
/// P₁ · P₂ = -½ |p₁ - p₂|²
/// ```
///
/// # Basis Blades (32 total)
///
/// | Grade | Count | Description |
/// |-------|-------|-------------|
/// | 0 | 1 | Scalar |
/// | 1 | 5 | Vectors (including conformal) |
/// | 2 | 10 | Bivectors (dipoles, lines) |
/// | 3 | 10 | Trivectors (circles, planes) |
/// | 4 | 5 | 4-vectors (spheres) |
/// | 5 | 1 | Pseudoscalar |
///
/// # Geometric Objects
///
/// CGA represents geometric objects as blades:
/// - **Grade 1**: Round points (null vectors)
/// - **Grade 2**: Dipoles (point pairs), flat points
/// - **Grade 3**: Circles, lines
/// - **Grade 4**: Spheres, planes
///
/// # Versors
///
/// CGA supports various transformations via versors:
/// - **Rotors**: Rotations
/// - **Translators**: Translations
/// - **Motors**: Combined rotation + translation
/// - **Dilators**: Uniform scaling
/// - **Inversors**: Sphere inversion
///
/// # Properties
///
/// - `e₁² = e₂² = e₃² = +1` (Euclidean)
/// - `e₊² = +1`, `e₋² = -1` (conformal)
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Conformal3::DIM, 5);
/// assert_eq!(Conformal3::num_blades(), 32);
///
/// // Verify null basis properties
/// let e_inf = Conformal3::e_infinity::<f64>();
/// let e_o = Conformal3::e_origin::<f64>();
///
/// // e∞² = 0
/// let inf_sq = &e_inf * &e_inf;
/// assert!(inf_sq.is_zero(1e-10));
///
/// // e₀² = 0
/// let o_sq = &e_o * &e_o;
/// assert!(o_sq.is_zero(1e-10));
///
/// // e∞ · e₀ = -1
/// let dot = e_inf.inner(&e_o).scalar_part();
/// assert!((dot - (-1.0)).abs() < 1e-10);
/// ```
///
/// Reference: <https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page>
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Conformal3;

impl Signature for Conformal3 {
    type NumBlades = U32; // 2^5 = 32

    const P: usize = 4; // e₁, e₂, e₃, e₊
    const Q: usize = 1; // e₋
    const R: usize = 0;

    #[inline]
    fn metric(i: usize) -> i8 {
        match i {
            0..=3 => 1, // e₁² = e₂² = e₃² = e₊² = +1
            4 => -1,    // e₋² = -1
            _ => panic!("basis index {i} out of range for Conformal3 (DIM=5)"),
        }
    }
}

/// Type alias for the 3D Conformal GA signature.
///
/// `Cl4_1_0` follows the naming convention `Cl{P}_{Q}_{R}` for the signature `(P, Q, R)`.
/// This is an alias for [`Conformal3`].
pub type Cl4_1_0 = Conformal3;

impl Conformal3 {
    /// Index of the positive conformal basis vector `e₊`.
    pub const E_PLUS: usize = 3;

    /// Index of the negative conformal basis vector `e₋`.
    pub const E_MINUS: usize = 4;

    /// Creates `e∞ = e₋ + e₊` (point at infinity).
    ///
    /// This is a null vector (`e∞² = 0`) representing the point at infinity.
    ///
    /// # Key Properties
    ///
    /// - `e∞ · e∞ = 0`
    /// - `e∞ · e₀ = -1`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::prelude::*;
    ///
    /// let e_inf = Conformal3::e_infinity::<f64>();
    /// let sq = &e_inf * &e_inf;
    /// assert!(sq.is_zero(1e-10));
    /// ```
    pub fn e_infinity<T: Float>() -> Multivector<T, Self> {
        let mut mv = Multivector::zero();
        mv.set(Blade::basis_vector(Self::E_PLUS), T::one());
        mv.set(Blade::basis_vector(Self::E_MINUS), T::one());
        mv
    }

    /// Creates `e₀ = (e₋ - e₊) / 2` (origin).
    ///
    /// This is a null vector (`e₀² = 0`) representing the origin.
    ///
    /// # Key Properties
    ///
    /// - `e₀ · e₀ = 0`
    /// - `e∞ · e₀ = -1`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::prelude::*;
    ///
    /// let e_o = Conformal3::e_origin::<f64>();
    /// let sq = &e_o * &e_o;
    /// assert!(sq.is_zero(1e-10));
    /// ```
    pub fn e_origin<T: Float>() -> Multivector<T, Self> {
        let mut mv = Multivector::zero();
        let half = T::one() / T::TWO;
        mv.set(Blade::basis_vector(Self::E_MINUS), half);
        mv.set(Blade::basis_vector(Self::E_PLUS), -half);
        mv
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;
    use proptest::prelude::*;

    // ========================================================================
    // Basic signature properties
    // ========================================================================

    #[test]
    fn conformal2_dimensions() {
        assert_eq!(Conformal2::DIM, 4);
        assert_eq!(Conformal2::num_blades(), 16);
    }

    #[test]
    fn conformal3_dimensions() {
        assert_eq!(Conformal3::DIM, 5);
        assert_eq!(Conformal3::num_blades(), 32);
    }

    #[test]
    fn conformal2_pqr_consistency() {
        assert_eq!(
            Conformal2::P + Conformal2::Q + Conformal2::R,
            Conformal2::DIM
        );
    }

    #[test]
    fn conformal3_pqr_consistency() {
        assert_eq!(
            Conformal3::P + Conformal3::Q + Conformal3::R,
            Conformal3::DIM
        );
    }

    proptest! {
        // ====================================================================
        // Metric properties
        // ====================================================================

        #[test]
        fn conformal2_euclidean_metric_positive(i in 0usize..2) {
            prop_assert_eq!(Conformal2::metric(i), 1);
        }

        #[test]
        fn conformal3_euclidean_metric_positive(i in 0usize..3) {
            prop_assert_eq!(Conformal3::metric(i), 1);
        }

        // ====================================================================
        // Null basis properties for Conformal2
        // ====================================================================

        #[test]
        fn conformal2_e_infinity_squares_to_zero(s in -100.0f64..100.0) {
            // e∞ = e₋ + e₊ should square to zero
            let e_inf = Conformal2::e_infinity::<f64>() * s;
            let sq = e_inf * e_inf;
            prop_assert!(sq.is_zero(RELATIVE_EQ_EPS));
        }

        #[test]
        fn conformal2_e_origin_squares_to_zero(s in -100.0f64..100.0) {
            // e₀ = (e₋ - e₊) / 2 should square to zero
            let e_o = Conformal2::e_origin::<f64>() * s;
            let sq = e_o * e_o;
            prop_assert!(sq.is_zero(RELATIVE_EQ_EPS));
        }

        #[test]
        fn conformal2_e_inf_dot_e_o_equals_minus_one(_dummy in 0..1i32) {
            // e∞ · e₀ = -1
            let e_inf = Conformal2::e_infinity::<f64>();
            let e_o = Conformal2::e_origin::<f64>();
            let dot = e_inf.inner(&e_o).scalar_part();
            prop_assert!(relative_eq!(dot, -1.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // Null basis properties for Conformal3
        // ====================================================================

        #[test]
        fn conformal3_e_infinity_squares_to_zero(s in -100.0f64..100.0) {
            // e∞ = e₋ + e₊ should square to zero
            let e_inf = Conformal3::e_infinity::<f64>() * s;
            let sq = e_inf * e_inf;
            prop_assert!(sq.is_zero(RELATIVE_EQ_EPS));
        }

        #[test]
        fn conformal3_e_origin_squares_to_zero(s in -100.0f64..100.0) {
            // e₀ = (e₋ - e₊) / 2 should square to zero
            let e_o = Conformal3::e_origin::<f64>() * s;
            let sq = e_o * e_o;
            prop_assert!(sq.is_zero(RELATIVE_EQ_EPS));
        }

        #[test]
        fn conformal3_e_inf_dot_e_o_equals_minus_one(_dummy in 0..1i32) {
            // e∞ · e₀ = -1
            let e_inf = Conformal3::e_infinity::<f64>();
            let e_o = Conformal3::e_origin::<f64>();
            let dot = e_inf.inner(&e_o).scalar_part();
            prop_assert!(relative_eq!(dot, -1.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // Algebraic properties for Conformal3
        // ====================================================================

        #[test]
        fn conformal3_geometric_product_associative(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
            c in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = (a * b) * c;
            let rhs = a * (b * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn conformal3_geometric_product_distributive(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
            c in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = a * (b + c);
            let rhs = (a * b) + (a * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn conformal3_reverse_involutory(a in any::<Multivector<f64, Conformal3>>()) {
            prop_assert!(relative_eq!(
                a.reverse().reverse(),
                a,
                max_relative = RELATIVE_EQ_EPS
            ));
        }

        #[test]
        fn conformal3_reverse_antimorphism(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = (a * b).reverse();
            let rhs = b.reverse() * a.reverse();
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn conformal3_outer_associative(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
            c in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = a.exterior(&b).exterior(&c);
            let rhs = a.exterior(&b.exterior(&c));
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // Euclidean subspace
        // ====================================================================

        #[test]
        fn conformal3_euclidean_subspace_inner_product(
            ax in -10.0f64..10.0, ay in -10.0f64..10.0, az in -10.0f64..10.0,
            bx in -10.0f64..10.0, by in -10.0f64..10.0, bz in -10.0f64..10.0,
        ) {
            // Pure Euclidean vectors (no e₊, e₋ components)
            let cga_a: Multivector<f64, Conformal3> =
                Multivector::vector(&[ax, ay, az, 0.0, 0.0]);
            let cga_b: Multivector<f64, Conformal3> =
                Multivector::vector(&[bx, by, bz, 0.0, 0.0]);

            // Inner product should match Euclidean dot product
            let cga_dot = cga_a.inner(&cga_b).scalar_part();
            let expected_dot = ax * bx + ay * by + az * bz;
            prop_assert!(relative_eq!(cga_dot, expected_dot, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }
    }

    // ========================================================================
    // Specific metric value tests
    // ========================================================================

    #[test]
    fn conformal2_metric_values() {
        assert_eq!(Conformal2::metric(0), 1); // e₁² = +1
        assert_eq!(Conformal2::metric(1), 1); // e₂² = +1
        assert_eq!(Conformal2::metric(2), 1); // e₊² = +1
        assert_eq!(Conformal2::metric(3), -1); // e₋² = -1
    }

    #[test]
    fn conformal3_metric_values() {
        assert_eq!(Conformal3::metric(0), 1); // e₁² = +1
        assert_eq!(Conformal3::metric(1), 1); // e₂² = +1
        assert_eq!(Conformal3::metric(2), 1); // e₃² = +1
        assert_eq!(Conformal3::metric(3), 1); // e₊² = +1
        assert_eq!(Conformal3::metric(4), -1); // e₋² = -1
    }

    #[test]
    fn conformal3_basis_vector_indices() {
        assert_eq!(Conformal3::E_PLUS, 3);
        assert_eq!(Conformal3::E_MINUS, 4);
    }

    #[test]
    fn conformal2_basis_vector_indices() {
        assert_eq!(Conformal2::E_PLUS, 2);
        assert_eq!(Conformal2::E_MINUS, 3);
    }
}
