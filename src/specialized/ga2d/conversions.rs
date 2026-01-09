//! Conversions between 2D specialized types and generic `Multivector`.
//!
//! This module provides bidirectional conversions:
//! - `From<SpecializedType> for Multivector<T, Euclidean2>` (infallible)
//! - `TryFrom<Multivector<T, Euclidean2>> for SpecializedType` (fallible, validates structure)
//! - `from_multivector_unchecked` methods (infallible, no validation)
//!
//! # Choosing a Conversion Method
//!
//! | Method | Validates | Branches | Use Case |
//! |--------|-----------|----------|----------|
//! | `TryFrom` | Yes | Yes | Runtime validation of computed results |
//! | `from_multivector_unchecked` | No | No | AD/dual numbers, performance-critical code |
//!
//! # Automatic Differentiation (AD) Compatibility
//!
//! The `TryFrom` implementations use branching to validate that coefficients
//! in unexpected grades are below a tolerance threshold. This branching is
//! **not compatible with automatic differentiation** using dual numbers, as:
//!
//! - Comparisons like `x.abs() > threshold` lose derivative information
//! - The branch creates a non-differentiable discontinuity
//!
//! For AD-compatible code paths, use the `from_multivector_unchecked` methods
//! which perform direct field extraction without any branching or validation.
//!
//! # Tolerance
//!
//! The `TryFrom` implementations use [`CONVERSION_TOLERANCE`] (1e-10) rather than
//! machine epsilon to account for accumulated floating-point errors in computed
//! results. This is intentionally more lenient than `T::EPSILON`.

use core::fmt;

use crate::algebra::Multivector;
use crate::basis::Blade;
use crate::scalar::Float;
use crate::signature::Euclidean2;

use super::{Bivec2, Rotor2, Vec2};

/// Error type for conversion from `Multivector` to specialized types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConversionError {
    /// The multivector has non-zero coefficients in unexpected grades.
    InvalidGrade,
}

impl fmt::Display for ConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGrade => write!(
                f,
                "multivector has non-zero coefficients in unexpected grades"
            ),
        }
    }
}

impl std::error::Error for ConversionError {}

// ============================================================================
// Blade indices for Euclidean2
// ============================================================================

/// Index for scalar (grade 0).
const SCALAR_IDX: usize = 0;
/// Index for e₁ (grade 1).
const E1_IDX: usize = 1;
/// Index for e₂ (grade 1).
const E2_IDX: usize = 2;
/// Index for e₁₂ (grade 2).
const E12_IDX: usize = 3;

/// Tolerance for `TryFrom` grade validation.
///
/// This is intentionally larger than machine epsilon to accommodate
/// accumulated floating-point errors from arithmetic operations.
/// Using `T::EPSILON` (~2.2e-16 for f64) would be too strict for
/// practical use with computed results.
///
/// For exact conversions without tolerance checking, use the
/// `from_multivector_unchecked` methods.
pub const CONVERSION_TOLERANCE: f64 = 1e-10;

// ============================================================================
// Vec2 conversions
// ============================================================================

impl<T: Float> From<Vec2<T>> for Multivector<T, Euclidean2> {
    /// Converts a 2D vector to a generic multivector.
    ///
    /// The vector components map to grade-1 blades:
    /// - `x` → `e₁` (index 1)
    /// - `y` → `e₂` (index 2)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Vec2;
    ///
    /// let v = Vec2::new(3.0_f64, 4.0);
    /// let mv: Multivector<f64, Euclidean2> = v.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(1)); // Pure vector
    /// ```
    fn from(v: Vec2<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E1_IDX), v.x);
        mv.set(Blade::from_index(E2_IDX), v.y);
        mv
    }
}

impl<T: Float> Vec2<T> {
    /// Converts a generic multivector to a 2D vector without validation.
    ///
    /// This method extracts the grade-1 components directly without checking
    /// that other grades are zero. It is branch-free and suitable for use
    /// with automatic differentiation (dual numbers).
    ///
    /// # Safety
    ///
    /// This is not unsafe in the Rust sense, but it may produce unexpected
    /// results if the input multivector has non-zero components in grades
    /// other than 1. Those components will be silently ignored.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Vec2;
    ///
    /// let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[3.0, 4.0]);
    /// let v = Vec2::from_multivector_unchecked(&mv);
    /// assert_eq!(v.x, 3.0);
    /// assert_eq!(v.y, 4.0);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean2>) -> Self {
        Self::new(
            mv.get(Blade::from_index(E1_IDX)),
            mv.get(Blade::from_index(E2_IDX)),
        )
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Vec2<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 2D vector.
    ///
    /// Succeeds only if the multivector is a pure vector (grade 1 only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Vec2::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// scalar or bivector components (above [`CONVERSION_TOLERANCE`]).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Vec2;
    ///
    /// let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[3.0, 4.0]);
    /// let v = Vec2::try_from(mv).unwrap();
    /// assert_eq!(v.x, 3.0);
    /// assert_eq!(v.y, 4.0);
    /// ```
    fn try_from(mv: Multivector<T, Euclidean2>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only grade-1 components are non-zero
        if mv.get(Blade::from_index(SCALAR_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }
        if mv.get(Blade::from_index(E12_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

// ============================================================================
// Bivec2 conversions
// ============================================================================

impl<T: Float> From<Bivec2<T>> for Multivector<T, Euclidean2> {
    /// Converts a 2D bivector to a generic multivector.
    ///
    /// The bivector coefficient maps to the grade-2 blade:
    /// - `value` → `e₁₂` (index 3)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Bivec2;
    ///
    /// let b = Bivec2::new(5.0_f64);
    /// let mv: Multivector<f64, Euclidean2> = b.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(2)); // Pure bivector
    /// ```
    fn from(b: Bivec2<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E12_IDX), b.0);
        mv
    }
}

impl<T: Float> Bivec2<T> {
    /// Converts a generic multivector to a 2D bivector without validation.
    ///
    /// This method extracts the grade-2 component directly without checking
    /// that other grades are zero. It is branch-free and suitable for use
    /// with automatic differentiation (dual numbers).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Bivec2;
    ///
    /// let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
    /// mv.set(Blade::from_index(3), 5.0); // e₁₂
    /// let b = Bivec2::from_multivector_unchecked(&mv);
    /// assert_eq!(b.value(), 5.0);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean2>) -> Self {
        Self::new(mv.get(Blade::from_index(E12_IDX)))
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Bivec2<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 2D bivector.
    ///
    /// Succeeds only if the multivector is a pure bivector (grade 2 only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Bivec2::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// scalar or vector components (above [`CONVERSION_TOLERANCE`]).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Bivec2;
    ///
    /// let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
    /// mv.set(Blade::from_index(3), 5.0); // e₁₂
    /// let b = Bivec2::try_from(mv).unwrap();
    /// assert_eq!(b.value(), 5.0);
    /// ```
    fn try_from(mv: Multivector<T, Euclidean2>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only grade-2 component is non-zero
        if mv.get(Blade::from_index(SCALAR_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }
        if mv.get(Blade::from_index(E1_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }
        if mv.get(Blade::from_index(E2_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

// ============================================================================
// Rotor2 conversions
// ============================================================================

impl<T: Float> From<Rotor2<T>> for Multivector<T, Euclidean2> {
    /// Converts a 2D rotor to a generic multivector.
    ///
    /// The rotor components map to:
    /// - `s` → scalar (index 0)
    /// - `xy` → `e₁₂` (index 3)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Rotor2;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let r = Rotor2::from_angle(FRAC_PI_4);
    /// let mv: Multivector<f64, Euclidean2> = r.into();
    ///
    /// // Rotor has scalar and bivector parts only
    /// assert!(mv.grade_select(1).is_zero(1e-10));
    /// ```
    fn from(r: Rotor2<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(SCALAR_IDX), r.s);
        mv.set(Blade::from_index(E12_IDX), r.xy);
        mv
    }
}

impl<T: Float> Rotor2<T> {
    /// Converts a generic multivector to a 2D rotor without validation.
    ///
    /// This method extracts the even-grade (scalar + bivector) components
    /// directly without checking that odd grades are zero. It is branch-free
    /// and suitable for use with automatic differentiation (dual numbers).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Rotor2;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let original = Rotor2::from_angle(FRAC_PI_4);
    /// let mv: Multivector<f64, Euclidean2> = original.into();
    /// let recovered = Rotor2::from_multivector_unchecked(&mv);
    ///
    /// assert!((original.s - recovered.s).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean2>) -> Self {
        Self::new(
            mv.get(Blade::from_index(SCALAR_IDX)),
            mv.get(Blade::from_index(E12_IDX)),
        )
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Rotor2<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 2D rotor.
    ///
    /// Succeeds only if the multivector is an even element (scalar + bivector only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Rotor2::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// vector components (above [`CONVERSION_TOLERANCE`]).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::ga2d::Rotor2;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let original = Rotor2::from_angle(FRAC_PI_4);
    /// let mv: Multivector<f64, Euclidean2> = original.into();
    /// let recovered = Rotor2::try_from(mv).unwrap();
    ///
    /// assert!((original.s - recovered.s).abs() < 1e-10);
    /// assert!((original.xy - recovered.xy).abs() < 1e-10);
    /// ```
    fn try_from(mv: Multivector<T, Euclidean2>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only even-grade components (scalar, bivector) are non-zero
        if mv.get(Blade::from_index(E1_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }
        if mv.get(Blade::from_index(E2_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;
    use std::f64::consts::FRAC_PI_4;

    use super::super::arbitrary::{UnitRotor2, UnitVec2};

    // ========================================================================
    // Vec2 tests
    // ========================================================================

    proptest! {
        #[test]
        fn vec2_roundtrip(v in any::<Vec2<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = v.into();
            let back = Vec2::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(v, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec2_to_multivector_is_grade_1(v in any::<Vec2<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = v.into();
            // Should be pure vector (grade 1) or zero
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(2).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn vec2_conversion_fails_with_scalar() {
        let mv: Multivector<f64, Euclidean2> = Multivector::scalar(1.0);
        assert!(Vec2::try_from(mv).is_err());
    }

    // ========================================================================
    // Bivec2 tests
    // ========================================================================

    proptest! {
        #[test]
        fn bivec2_roundtrip(val in -100.0f64..100.0) {
            let b = Bivec2::new(val);
            let mv: Multivector<f64, Euclidean2> = b.into();
            let back = Bivec2::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(b, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivec2_to_multivector_is_grade_2(val in -100.0f64..100.0) {
            let b = Bivec2::new(val);
            let mv: Multivector<f64, Euclidean2> = b.into();
            // Should be pure bivector (grade 2) or zero
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn bivec2_conversion_fails_with_vector() {
        let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[1.0, 0.0]);
        assert!(Bivec2::try_from(mv).is_err());
    }

    // ========================================================================
    // Rotor2 tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor2_roundtrip(r in any::<UnitRotor2<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = (*r).into();
            let back = Rotor2::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(*r, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor2_to_multivector_is_even(r in any::<UnitRotor2<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = (*r).into();
            // Should have no vector (odd) part
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn rotor2_conversion_fails_with_vector() {
        let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[1.0, 0.0]);
        assert!(Rotor2::try_from(mv).is_err());
    }

    // ========================================================================
    // Operation consistency tests
    // ========================================================================

    proptest! {
        #[test]
        fn dot_consistency(
            a in any::<Vec2<f64>>(),
            b in any::<Vec2<f64>>(),
        ) {
            let spec_result = a.dot(b);

            let gen_a: Multivector<f64, Euclidean2> = a.into();
            let gen_b: Multivector<f64, Euclidean2> = b.into();
            let gen_result = gen_a.inner(&gen_b).scalar_part();

            prop_assert!(abs_diff_eq!(spec_result, gen_result, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn wedge_consistency(
            a in any::<Vec2<f64>>(),
            b in any::<Vec2<f64>>(),
        ) {
            let spec_result = a.wedge(b);

            let gen_a: Multivector<f64, Euclidean2> = a.into();
            let gen_b: Multivector<f64, Euclidean2> = b.into();
            let gen_result = gen_a.outer(&gen_b);
            let gen_as_bivec = Bivec2::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_bivec, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_rotation_consistency(
            r in any::<UnitRotor2<f64>>(),
            v in any::<UnitVec2<f64>>(),
        ) {
            // Specialized rotation
            let spec_result = r.rotate(*v);

            // Generic sandwich product: R̃ v R
            let gen_r: Multivector<f64, Euclidean2> = (*r).into();
            let gen_v: Multivector<f64, Euclidean2> = (*v).into();
            let gen_r_rev = gen_r.reverse();
            let gen_result = &(&gen_r_rev * &gen_v) * &gen_r;
            let gen_as_vec = Vec2::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_vec, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn rotor_from_angle_consistency() {
        let r = Rotor2::from_angle(FRAC_PI_4);
        let mv: Multivector<f64, Euclidean2> = r.into();

        // cos(π/8) and sin(π/8)
        let expected_s = (FRAC_PI_4 / 2.0).cos();
        let expected_xy = (FRAC_PI_4 / 2.0).sin();

        assert!(abs_diff_eq!(
            mv.scalar_part(),
            expected_s,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E12_IDX)),
            expected_xy,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // ========================================================================
    // Tests starting from generic Multivector (TryFrom consistency)
    // ========================================================================

    proptest! {
        /// Test that a generic vector Multivector converts to Vec2 correctly.
        #[test]
        fn generic_vector_to_vec2(x in -100.0f64..100.0, y in -100.0f64..100.0) {
            let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[x, y]);
            let v = Vec2::try_from(mv).expect("pure vector should convert");
            prop_assert!(abs_diff_eq!(v.x, x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(v.y, y, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic bivector Multivector converts to Bivec2 correctly.
        #[test]
        fn generic_bivector_to_bivec2(xy in -100.0f64..100.0) {
            let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
            mv.set(Blade::from_index(E12_IDX), xy);
            let b = Bivec2::try_from(mv).expect("pure bivector should convert");
            prop_assert!(abs_diff_eq!(b.value(), xy, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic even Multivector converts to Rotor2 correctly.
        #[test]
        fn generic_even_to_rotor2(s in -100.0f64..100.0, xy in -100.0f64..100.0) {
            let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
            mv.set(Blade::from_index(SCALAR_IDX), s);
            mv.set(Blade::from_index(E12_IDX), xy);
            let r = Rotor2::try_from(mv).expect("even element should convert");
            prop_assert!(abs_diff_eq!(r.s, s, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(r.xy, xy, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that full generic Multivector fails TryFrom for specialized types.
        #[test]
        fn generic_mixed_fails_tryfrom(mv in any::<Multivector<f64, Euclidean2>>()) {
            // A random multivector almost certainly has mixed grades
            // Only check if it has more than one grade
            let has_scalar = mv.get(Blade::from_index(SCALAR_IDX)).abs() > ABS_DIFF_EQ_EPS;
            let has_vec = mv.get(Blade::from_index(E1_IDX)).abs() > ABS_DIFF_EQ_EPS
                || mv.get(Blade::from_index(E2_IDX)).abs() > ABS_DIFF_EQ_EPS;
            let has_bivec = mv.get(Blade::from_index(E12_IDX)).abs() > ABS_DIFF_EQ_EPS;

            if has_vec && (has_scalar || has_bivec) {
                // Has vector and something else - should fail vec conversion
                prop_assert!(Vec2::try_from(mv.clone()).is_err());
            }
            if has_bivec && (has_scalar || has_vec) {
                // Has bivector and scalar or vector - should fail bivec conversion
                prop_assert!(Bivec2::try_from(mv).is_err());
            }
        }

        /// Test that vector wedge product (generic) converts to Bivec2.
        #[test]
        fn generic_wedge_to_bivec2(
            ax in -10.0f64..10.0, ay in -10.0f64..10.0,
            bx in -10.0f64..10.0, by in -10.0f64..10.0,
        ) {
            let a: Multivector<f64, Euclidean2> = Multivector::vector(&[ax, ay]);
            let b: Multivector<f64, Euclidean2> = Multivector::vector(&[bx, by]);
            let wedge = a.outer(&b);

            // Should be a pure bivector, convertible to Bivec2
            let bivec = Bivec2::try_from(wedge).expect("wedge of vectors should be bivector");

            // Check value matches ax*by - ay*bx
            let expected = ax * by - ay * bx;
            prop_assert!(abs_diff_eq!(bivec.value(), expected, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that vector geometric product (generic) converts to Rotor2.
        #[test]
        fn generic_geometric_to_rotor2(
            ax in -10.0f64..10.0, ay in -10.0f64..10.0,
            bx in -10.0f64..10.0, by in -10.0f64..10.0,
        ) {
            let a: Multivector<f64, Euclidean2> = Multivector::vector(&[ax, ay]);
            let b: Multivector<f64, Euclidean2> = Multivector::vector(&[bx, by]);
            let product = &a * &b;

            // Should be scalar + bivector (even), convertible to Rotor2
            let rotor = Rotor2::try_from(product).expect("geometric product should be even");

            // Check components
            let expected_s = ax * bx + ay * by; // dot product
            let expected_xy = ax * by - ay * bx; // wedge product
            prop_assert!(abs_diff_eq!(rotor.s, expected_s, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotor.xy, expected_xy, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
