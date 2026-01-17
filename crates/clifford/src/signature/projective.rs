//! Projective signature types for Projective Geometric Algebra (PGA).
//!
//! This module provides signature types for Projective Geometric Algebras,
//! which extend Euclidean space with a degenerate (null) basis vector `e₀`
//! that squares to zero (`e₀² = 0`).
//!
//! # Point-Based Formulation
//!
//! This library uses the **point-based** formulation of PGA where:
//! - Grade 1 (vectors): Points in homogeneous coordinates
//! - Grade 2 (bivectors): Lines
//! - Grade 3 (trivectors): Planes (in 3D)
//!
//! This is consistent with Euclidean GA conventions and nalgebra's `Point` types.
//!
//! # Available Signatures
//!
//! - [`Projective2`]: 2D PGA, `Cl(2,0,1)` with 8 basis blades
//! - [`Projective3`]: 3D PGA, `Cl(3,0,1)` with 16 basis blades
//!
//! # Basis Vector Ordering
//!
//! Basis vectors are ordered with Euclidean vectors first, then the null vector:
//! - 2D PGA: `e₁, e₂, e₀` (indices 0, 1, 2)
//! - 3D PGA: `e₁, e₂, e₃, e₀` (indices 0, 1, 2, 3)
//!
//! # Example
//!
//! ```
//! use clifford::prelude::*;
//!
//! // 3D Projective GA
//! assert_eq!(Projective3::DIM, 4);
//! assert_eq!(Projective3::num_blades(), 16);
//!
//! // Euclidean basis vectors square to +1
//! assert_eq!(Projective3::metric(0), 1);  // e₁² = 1
//! assert_eq!(Projective3::metric(1), 1);  // e₂² = 1
//! assert_eq!(Projective3::metric(2), 1);  // e₃² = 1
//!
//! // Null basis vector squares to 0
//! assert_eq!(Projective3::metric(3), 0);  // e₀² = 0
//! ```

use super::Signature;
use typenum::{U8, U16};

/// 2D Projective signature: `Cl(2,0,1)`.
///
/// Represents 2D projective space with basis vectors `e₁`, `e₂`, and `e₀`,
/// where `e₀² = 0` is the null (degenerate) basis vector.
///
/// # Point-Based Interpretation
///
/// In the point-based formulation:
/// - **Grade 1 (vectors)**: Points `P = x·e₁ + y·e₂ + w·e₀`
/// - **Grade 2 (bivectors)**: Lines
/// - **Grade 3 (pseudoscalar)**: Oriented area
///
/// A point with `w = 0` represents an ideal point (point at infinity).
///
/// # Basis Blades (8 total)
///
/// | Index | Binary | Blade | Grade | Description |
/// |-------|--------|-------|-------|-------------|
/// | 0 | `000` | `1` | 0 | Scalar |
/// | 1 | `001` | `e₁` | 1 | x-component of point |
/// | 2 | `010` | `e₂` | 1 | y-component of point |
/// | 3 | `011` | `e₁₂` | 2 | Line component |
/// | 4 | `100` | `e₀` | 1 | Homogeneous weight |
/// | 5 | `101` | `e₀₁` | 2 | Line component |
/// | 6 | `110` | `e₀₂` | 2 | Line component |
/// | 7 | `111` | `e₀₁₂` | 3 | Pseudoscalar |
///
/// # Properties
///
/// - `e₁² = e₂² = +1` (Euclidean)
/// - `e₀² = 0` (null/degenerate)
/// - `e₁₂² = -1` (rotation generator)
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Projective2::DIM, 3);
/// assert_eq!(Projective2::num_blades(), 8);
///
/// // Euclidean basis vectors
/// assert_eq!(Projective2::metric(0), 1);  // e₁² = 1
/// assert_eq!(Projective2::metric(1), 1);  // e₂² = 1
///
/// // Null basis vector
/// assert_eq!(Projective2::metric(2), 0);  // e₀² = 0
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Projective2;

impl Signature for Projective2 {
    type NumBlades = U8; // 2^3 = 8

    const P: usize = 2;
    const Q: usize = 0;
    const R: usize = 1;

    #[inline]
    fn metric(i: usize) -> i8 {
        match i {
            0 | 1 => 1, // e₁² = e₂² = 1
            2 => 0,     // e₀² = 0
            _ => panic!("basis index {i} out of range for Projective2 (DIM=3)"),
        }
    }
}

/// 3D Projective signature: `Cl(3,0,1)`.
///
/// Represents 3D projective space with basis vectors `e₁`, `e₂`, `e₃`, and `e₀`,
/// where `e₀² = 0` is the null (degenerate) basis vector.
///
/// # Point-Based Interpretation
///
/// In the point-based formulation:
/// - **Grade 1 (vectors)**: Points `P = x·e₁ + y·e₂ + z·e₃ + w·e₀`
/// - **Grade 2 (bivectors)**: Lines (Plücker coordinates)
/// - **Grade 3 (trivectors)**: Planes
/// - **Grade 4 (pseudoscalar)**: Oriented volume
///
/// A point with `w = 0` represents an ideal point (point at infinity).
///
/// # Basis Blades (16 total)
///
/// | Grade | Count | Blades | Description |
/// |-------|-------|--------|-------------|
/// | 0 | 1 | `1` | Scalar |
/// | 1 | 4 | `e₁, e₂, e₃, e₀` | Points |
/// | 2 | 6 | `e₁₂, e₁₃, e₂₃, e₀₁, e₀₂, e₀₃` | Lines |
/// | 3 | 4 | `e₁₂₃, e₀₁₂, e₀₁₃, e₀₂₃` | Planes |
/// | 4 | 1 | `e₀₁₂₃` | Pseudoscalar |
///
/// # Motors
///
/// The even subalgebra (grades 0, 2, 4) contains **motors**, which represent
/// rigid body transformations (rotation + translation). A motor `M` transforms
/// a point `P` via the sandwich product: `P' = M P M̃`.
///
/// # Properties
///
/// - `e₁² = e₂² = e₃² = +1` (Euclidean)
/// - `e₀² = 0` (null/degenerate)
/// - Bivectors `e₁₂, e₁₃, e₂₃` each square to `-1` (rotation generators)
/// - Bivectors `e₀₁, e₀₂, e₀₃` square to `0` (translation generators)
///
/// # Example
///
/// ```
/// use clifford::prelude::*;
///
/// assert_eq!(Projective3::DIM, 4);
/// assert_eq!(Projective3::num_blades(), 16);
///
/// // Euclidean basis vectors
/// assert_eq!(Projective3::metric(0), 1);  // e₁² = 1
/// assert_eq!(Projective3::metric(1), 1);  // e₂² = 1
/// assert_eq!(Projective3::metric(2), 1);  // e₃² = 1
///
/// // Null basis vector
/// assert_eq!(Projective3::metric(3), 0);  // e₀² = 0
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Projective3;

impl Signature for Projective3 {
    type NumBlades = U16; // 2^4 = 16

    const P: usize = 3;
    const Q: usize = 0;
    const R: usize = 1;

    #[inline]
    fn metric(i: usize) -> i8 {
        match i {
            0..=2 => 1, // e₁² = e₂² = e₃² = 1
            3 => 0,     // e₀² = 0
            _ => panic!("basis index {i} out of range for Projective3 (DIM=4)"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::Multivector;
    use crate::basis::Blade;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;
    use proptest::prelude::*;

    // ========================================================================
    // Basic signature tests
    // ========================================================================

    proptest! {
        #[test]
        fn projective2_euclidean_metric_positive(i in 0usize..2) {
            prop_assert_eq!(Projective2::metric(i), 1);
        }

        #[test]
        fn projective3_euclidean_metric_positive(i in 0usize..3) {
            prop_assert_eq!(Projective3::metric(i), 1);
        }
    }

    #[test]
    fn projective2_null_metric() {
        assert_eq!(Projective2::metric(2), 0);
    }

    #[test]
    fn projective3_null_metric() {
        assert_eq!(Projective3::metric(3), 0);
    }

    #[test]
    fn projective_dimensions() {
        assert_eq!(Projective2::DIM, 3);
        assert_eq!(Projective2::num_blades(), 8);

        assert_eq!(Projective3::DIM, 4);
        assert_eq!(Projective3::num_blades(), 16);
    }

    #[test]
    fn projective_pqr_consistency() {
        assert_eq!(
            Projective2::P + Projective2::Q + Projective2::R,
            Projective2::DIM
        );
        assert_eq!(
            Projective3::P + Projective3::Q + Projective3::R,
            Projective3::DIM
        );
    }

    // ========================================================================
    // Null vector behavior tests
    // ========================================================================

    #[test]
    fn pga2_null_vector_squares_to_zero() {
        // e₀ is at index 2 (basis vector index), blade index is 1 << 2 = 4
        let e0: Multivector<f64, Projective2> = Multivector::basis_vector(2);
        let e0_sq = e0 * e0;
        assert!(e0_sq.is_zero(RELATIVE_EQ_EPS));
    }

    #[test]
    fn pga3_null_vector_squares_to_zero() {
        // e₀ is at index 3 (basis vector index), blade index is 1 << 3 = 8
        let e0: Multivector<f64, Projective3> = Multivector::basis_vector(3);
        let e0_sq = e0 * e0;
        assert!(e0_sq.is_zero(RELATIVE_EQ_EPS));
    }

    #[test]
    fn pga3_euclidean_vectors_square_to_one() {
        let e1: Multivector<f64, Projective3> = Multivector::basis_vector(0);
        let e2: Multivector<f64, Projective3> = Multivector::basis_vector(1);
        let e3: Multivector<f64, Projective3> = Multivector::basis_vector(2);

        assert!(relative_eq!(
            (e1 * e1).scalar_part(),
            1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            (e2 * e2).scalar_part(),
            1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            (e3 * e3).scalar_part(),
            1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn pga3_euclidean_bivector_squares_to_minus_one() {
        // e₁₂ should square to -1
        let e1: Multivector<f64, Projective3> = Multivector::basis_vector(0);
        let e2: Multivector<f64, Projective3> = Multivector::basis_vector(1);
        let e12 = e1 * e2;
        let e12_sq = e12 * e12;

        assert!(relative_eq!(
            e12_sq.scalar_part(),
            -1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn pga3_null_bivector_squares_to_zero() {
        // e₀₁ should square to 0 (contains null vector)
        let e0: Multivector<f64, Projective3> = Multivector::basis_vector(3);
        let e1: Multivector<f64, Projective3> = Multivector::basis_vector(0);
        let e01 = e0 * e1;
        let e01_sq = e01 * e01;

        assert!(e01_sq.is_zero(RELATIVE_EQ_EPS));
    }

    // ========================================================================
    // Point representation tests
    // ========================================================================

    #[test]
    fn pga3_point_at_origin() {
        // Origin is represented as e₀ (just the homogeneous weight)
        let origin: Multivector<f64, Projective3> = Multivector::basis_vector(3);

        // Should be grade 1
        assert_eq!(origin.grade(RELATIVE_EQ_EPS), Some(1));
    }

    #[test]
    fn pga3_finite_point() {
        // Point at (1, 2, 3) is e₁ + 2e₂ + 3e₃ + e₀
        let mut point: Multivector<f64, Projective3> = Multivector::zero();
        point.set(Blade::basis_vector(0), 1.0); // e₁
        point.set(Blade::basis_vector(1), 2.0); // e₂
        point.set(Blade::basis_vector(2), 3.0); // e₃
        point.set(Blade::basis_vector(3), 1.0); // e₀ (weight = 1)

        // Should be grade 1
        assert_eq!(point.grade(RELATIVE_EQ_EPS), Some(1));
    }

    #[test]
    fn pga3_ideal_point() {
        // Ideal point (at infinity) in direction (1, 0, 0) has no e₀ component
        let mut ideal: Multivector<f64, Projective3> = Multivector::zero();
        ideal.set(Blade::basis_vector(0), 1.0); // e₁ only

        // Should be grade 1
        assert_eq!(ideal.grade(RELATIVE_EQ_EPS), Some(1));

        // No e₀ component (index 3 -> blade index 8)
        assert!(relative_eq!(
            ideal.get(Blade::basis_vector(3)),
            0.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    // ========================================================================
    // Property-based algebraic tests for PGA2
    // ========================================================================

    proptest! {
        #[test]
        fn pga2_geometric_product_associative(
            a in any::<Multivector<f64, Projective2>>(),
            b in any::<Multivector<f64, Projective2>>(),
            c in any::<Multivector<f64, Projective2>>(),
        ) {
            let lhs = (a * b) * c;
            let rhs = a * (b * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_geometric_product_distributive(
            a in any::<Multivector<f64, Projective2>>(),
            b in any::<Multivector<f64, Projective2>>(),
            c in any::<Multivector<f64, Projective2>>(),
        ) {
            let lhs = a * (b + c);
            let rhs = (a * b) + (a * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_scaled_null_vector_squares_to_zero(s in -100.0f64..100.0) {
            // Any scalar multiple of e₀ should square to zero
            let e0: Multivector<f64, Projective2> = Multivector::basis_vector(2);
            let scaled = e0 * s;
            let sq = scaled * scaled;
            prop_assert!(sq.is_zero(RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_reverse_involutory(a in any::<Multivector<f64, Projective2>>()) {
            prop_assert!(relative_eq!(a.reverse().reverse(), a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_reverse_antimorphism(
            a in any::<Multivector<f64, Projective2>>(),
            b in any::<Multivector<f64, Projective2>>(),
        ) {
            let lhs = (a * b).reverse();
            let rhs = b.reverse() * a.reverse();
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_involute_involutory(a in any::<Multivector<f64, Projective2>>()) {
            prop_assert!(relative_eq!(a.involute().involute(), a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_outer_associative(
            a in any::<Multivector<f64, Projective2>>(),
            b in any::<Multivector<f64, Projective2>>(),
            c in any::<Multivector<f64, Projective2>>(),
        ) {
            let lhs = a.exterior(&b).exterior(&c);
            let rhs = a.exterior(&b.exterior(&c));
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_one_is_multiplicative_identity(a in any::<Multivector<f64, Projective2>>()) {
            let one = Multivector::<f64, Projective2>::one();
            prop_assert!(relative_eq!(a * one, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(one * a, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_zero_is_additive_identity(a in any::<Multivector<f64, Projective2>>()) {
            let zero = Multivector::<f64, Projective2>::zero();
            prop_assert!(relative_eq!(a + zero, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_add_commutative(
            a in any::<Multivector<f64, Projective2>>(),
            b in any::<Multivector<f64, Projective2>>()
        ) {
            prop_assert!(relative_eq!(a + b, b + a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga2_grade_decomposition_complete(a in any::<Multivector<f64, Projective2>>()) {
            // Sum of all grade projections equals original (grades 0, 1, 2, 3)
            let sum = (0..=3)
                .map(|k| a.grade_select(k))
                .fold(Multivector::<f64, Projective2>::zero(), |acc, x| acc + x);
            prop_assert!(relative_eq!(sum, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }
    }

    // ========================================================================
    // Property-based algebraic tests for PGA3
    // ========================================================================

    proptest! {
        #[test]
        fn pga3_geometric_product_associative(
            a in any::<Multivector<f64, Projective3>>(),
            b in any::<Multivector<f64, Projective3>>(),
            c in any::<Multivector<f64, Projective3>>(),
        ) {
            let lhs = (a * b) * c;
            let rhs = a * (b * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_geometric_product_distributive(
            a in any::<Multivector<f64, Projective3>>(),
            b in any::<Multivector<f64, Projective3>>(),
            c in any::<Multivector<f64, Projective3>>(),
        ) {
            let lhs = a * (b + c);
            let rhs = (a * b) + (a * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_scaled_null_vector_squares_to_zero(s in -100.0f64..100.0) {
            // Any scalar multiple of e₀ should square to zero
            let e0: Multivector<f64, Projective3> = Multivector::basis_vector(3);
            let scaled = e0 * s;
            let sq = scaled * scaled;
            prop_assert!(sq.is_zero(RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_reverse_involutory(a in any::<Multivector<f64, Projective3>>()) {
            prop_assert!(relative_eq!(a.reverse().reverse(), a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_reverse_antimorphism(
            a in any::<Multivector<f64, Projective3>>(),
            b in any::<Multivector<f64, Projective3>>(),
        ) {
            let lhs = (a * b).reverse();
            let rhs = b.reverse() * a.reverse();
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_involute_involutory(a in any::<Multivector<f64, Projective3>>()) {
            prop_assert!(relative_eq!(a.involute().involute(), a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_outer_associative(
            a in any::<Multivector<f64, Projective3>>(),
            b in any::<Multivector<f64, Projective3>>(),
            c in any::<Multivector<f64, Projective3>>(),
        ) {
            let lhs = a.exterior(&b).exterior(&c);
            let rhs = a.exterior(&b.exterior(&c));
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_one_is_multiplicative_identity(a in any::<Multivector<f64, Projective3>>()) {
            let one = Multivector::<f64, Projective3>::one();
            prop_assert!(relative_eq!(a * one, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(one * a, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_zero_is_additive_identity(a in any::<Multivector<f64, Projective3>>()) {
            let zero = Multivector::<f64, Projective3>::zero();
            prop_assert!(relative_eq!(a + zero, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_add_commutative(
            a in any::<Multivector<f64, Projective3>>(),
            b in any::<Multivector<f64, Projective3>>()
        ) {
            prop_assert!(relative_eq!(a + b, b + a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_grade_decomposition_complete(a in any::<Multivector<f64, Projective3>>()) {
            // Sum of all grade projections equals original (grades 0, 1, 2, 3, 4)
            let sum = (0..=4)
                .map(|k| a.grade_select(k))
                .fold(Multivector::<f64, Projective3>::zero(), |acc, x| acc + x);
            prop_assert!(relative_eq!(sum, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn pga3_even_plus_odd_equals_original(a in any::<Multivector<f64, Projective3>>()) {
            let reconstructed = a.even() + a.odd();
            prop_assert!(relative_eq!(reconstructed, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // Note: dual/undual tests are NOT included because in PGA the pseudoscalar
        // contains the null vector e₀ and is therefore not invertible (I² = 0).
        // The dual operation requires an invertible pseudoscalar.
    }

    // ========================================================================
    // nalgebra consistency tests
    // ========================================================================

    #[cfg(feature = "nalgebra-0_33")]
    mod nalgebra_consistency_0_33 {
        use super::*;
        use nalgebra_0_33 as na;

        proptest! {
            /// Euclidean subspace dot product matches nalgebra
            #[test]
            fn pga3_euclidean_dot_matches_nalgebra(
                ax in -100.0f64..100.0, ay in -100.0f64..100.0, az in -100.0f64..100.0,
                bx in -100.0f64..100.0, by in -100.0f64..100.0, bz in -100.0f64..100.0,
            ) {
                // Pure Euclidean vectors (no e₀ component) - components for e₁, e₂, e₃, e₀
                let pga_a: Multivector<f64, Projective3> =
                    Multivector::vector(&[ax, ay, az, 0.0]);
                let pga_b: Multivector<f64, Projective3> =
                    Multivector::vector(&[bx, by, bz, 0.0]);

                let na_a = na::Vector3::new(ax, ay, az);
                let na_b = na::Vector3::new(bx, by, bz);

                let pga_dot = pga_a.inner(&pga_b).scalar_part();
                let na_dot = na_a.dot(&na_b);

                prop_assert!(relative_eq!(pga_dot, na_dot, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            }

            /// Euclidean subspace geometric product gives correct scalar + bivector
            #[test]
            fn pga3_euclidean_vectors_product(
                ax in -10.0f64..10.0, ay in -10.0f64..10.0, az in -10.0f64..10.0,
                bx in -10.0f64..10.0, by in -10.0f64..10.0, bz in -10.0f64..10.0,
            ) {
                let pga_a: Multivector<f64, Projective3> =
                    Multivector::vector(&[ax, ay, az, 0.0]);
                let pga_b: Multivector<f64, Projective3> =
                    Multivector::vector(&[bx, by, bz, 0.0]);

                let na_a = na::Vector3::new(ax, ay, az);
                let na_b = na::Vector3::new(bx, by, bz);

                // Geometric product of vectors: ab = a·b + a∧b
                let product = pga_a * pga_b;

                // Scalar part should be dot product
                let expected_scalar = na_a.dot(&na_b);
                prop_assert!(relative_eq!(
                    product.scalar_part(),
                    expected_scalar,
                    max_relative = RELATIVE_EQ_EPS
                ));

                // Should have no grade-1 or grade-3 parts (only scalar and bivector)
                prop_assert!(product.grade_select(1).is_zero(RELATIVE_EQ_EPS));
                prop_assert!(product.grade_select(3).is_zero(RELATIVE_EQ_EPS));
            }

            /// Cross product can be extracted from wedge product in 3D
            #[test]
            fn pga3_wedge_matches_cross_product(
                ax in -10.0f64..10.0, ay in -10.0f64..10.0, az in -10.0f64..10.0,
                bx in -10.0f64..10.0, by in -10.0f64..10.0, bz in -10.0f64..10.0,
            ) {
                let pga_a: Multivector<f64, Projective3> =
                    Multivector::vector(&[ax, ay, az, 0.0]);
                let pga_b: Multivector<f64, Projective3> =
                    Multivector::vector(&[bx, by, bz, 0.0]);

                let na_a = na::Vector3::new(ax, ay, az);
                let na_b = na::Vector3::new(bx, by, bz);
                let na_cross = na_a.cross(&na_b);

                // The wedge product a∧b gives a bivector
                let wedge = pga_a.exterior(&pga_b);

                // In 3D Euclidean subspace, the bivector components are related to cross product:
                // e₁₂ component = ax*by - ay*bx = (a×b)_z
                // e₁₃ component = ax*bz - az*bx = -(a×b)_y
                // e₂₃ component = ay*bz - az*by = (a×b)_x
                //
                // Blade indices: e₁₂ = 0b0011 = 3, e₁₃ = 0b0101 = 5, e₂₃ = 0b0110 = 6
                let biv_12 = wedge.get(Blade::from_index(0b0011)); // e₁₂
                let biv_13 = wedge.get(Blade::from_index(0b0101)); // e₁₃
                let biv_23 = wedge.get(Blade::from_index(0b0110)); // e₂₃

                prop_assert!(relative_eq!(biv_23, na_cross.x, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
                prop_assert!(relative_eq!(-biv_13, na_cross.y, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
                prop_assert!(relative_eq!(biv_12, na_cross.z, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            }
        }
    }

    #[cfg(feature = "nalgebra-0_34")]
    mod nalgebra_consistency_0_34 {
        use super::*;
        use nalgebra_0_34 as na;

        proptest! {
            /// Euclidean subspace dot product matches nalgebra
            #[test]
            fn pga3_euclidean_dot_matches_nalgebra(
                ax in -100.0f64..100.0, ay in -100.0f64..100.0, az in -100.0f64..100.0,
                bx in -100.0f64..100.0, by in -100.0f64..100.0, bz in -100.0f64..100.0,
            ) {
                // Pure Euclidean vectors (no e₀ component) - components for e₁, e₂, e₃, e₀
                let pga_a: Multivector<f64, Projective3> =
                    Multivector::vector(&[ax, ay, az, 0.0]);
                let pga_b: Multivector<f64, Projective3> =
                    Multivector::vector(&[bx, by, bz, 0.0]);

                let na_a = na::Vector3::new(ax, ay, az);
                let na_b = na::Vector3::new(bx, by, bz);

                let pga_dot = pga_a.inner(&pga_b).scalar_part();
                let na_dot = na_a.dot(&na_b);

                prop_assert!(relative_eq!(pga_dot, na_dot, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            }

            /// Euclidean subspace geometric product gives correct scalar + bivector
            #[test]
            fn pga3_euclidean_vectors_product(
                ax in -10.0f64..10.0, ay in -10.0f64..10.0, az in -10.0f64..10.0,
                bx in -10.0f64..10.0, by in -10.0f64..10.0, bz in -10.0f64..10.0,
            ) {
                let pga_a: Multivector<f64, Projective3> =
                    Multivector::vector(&[ax, ay, az, 0.0]);
                let pga_b: Multivector<f64, Projective3> =
                    Multivector::vector(&[bx, by, bz, 0.0]);

                let na_a = na::Vector3::new(ax, ay, az);
                let na_b = na::Vector3::new(bx, by, bz);

                // Geometric product of vectors: ab = a·b + a∧b
                let product = &pga_a * &pga_b;

                // Scalar part should be dot product
                let expected_scalar = na_a.dot(&na_b);
                prop_assert!(relative_eq!(
                    product.scalar_part(),
                    expected_scalar,
                    max_relative = RELATIVE_EQ_EPS
                ));

                // Should have no grade-1 or grade-3 parts (only scalar and bivector)
                prop_assert!(product.grade_select(1).is_zero(RELATIVE_EQ_EPS));
                prop_assert!(product.grade_select(3).is_zero(RELATIVE_EQ_EPS));
            }

            /// Cross product can be extracted from wedge product in 3D
            #[test]
            fn pga3_wedge_matches_cross_product(
                ax in -10.0f64..10.0, ay in -10.0f64..10.0, az in -10.0f64..10.0,
                bx in -10.0f64..10.0, by in -10.0f64..10.0, bz in -10.0f64..10.0,
            ) {
                let pga_a: Multivector<f64, Projective3> =
                    Multivector::vector(&[ax, ay, az, 0.0]);
                let pga_b: Multivector<f64, Projective3> =
                    Multivector::vector(&[bx, by, bz, 0.0]);

                let na_a = na::Vector3::new(ax, ay, az);
                let na_b = na::Vector3::new(bx, by, bz);
                let na_cross = na_a.cross(&na_b);

                // The wedge product a∧b gives a bivector
                let wedge = pga_a.exterior(&pga_b);

                // In 3D Euclidean subspace, the bivector components are related to cross product:
                // e₁₂ component = ax*by - ay*bx = (a×b)_z
                // e₁₃ component = ax*bz - az*bx = -(a×b)_y
                // e₂₃ component = ay*bz - az*by = (a×b)_x
                //
                // Blade indices: e₁₂ = 0b0011 = 3, e₁₃ = 0b0101 = 5, e₂₃ = 0b0110 = 6
                let biv_12 = wedge.get(Blade::from_index(0b0011)); // e₁₂
                let biv_13 = wedge.get(Blade::from_index(0b0101)); // e₁₃
                let biv_23 = wedge.get(Blade::from_index(0b0110)); // e₂₃

                prop_assert!(relative_eq!(biv_23, na_cross.x, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
                prop_assert!(relative_eq!(-biv_13, na_cross.y, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
                prop_assert!(relative_eq!(biv_12, na_cross.z, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            }
        }
    }
}
