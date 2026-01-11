//! Conversions between 3D CGA specialized types and generic `Multivector`.
//!
//! This module provides bidirectional conversions:
//! - `From<SpecializedType> for Multivector<T, Conformal3>` (infallible)
//! - `TryFrom<Multivector<T, Conformal3>> for SpecializedType` (fallible, validates structure)
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
//! **not compatible with automatic differentiation** using dual numbers.
//!
//! For AD-compatible code paths, use the `from_multivector_unchecked` methods
//! which perform direct field extraction without any branching or validation.
//!
//! # CGA Basis Structure
//!
//! Conformal 3D uses basis vectors: `e₁, e₂, e₃, e₊, e₋`
//! - Euclidean: `e₁² = e₂² = e₃² = +1`
//! - Conformal: `e₊² = +1`, `e₋² = -1`
//!
//! Null basis:
//! - `e∞ = e₋ + e₊` (point at infinity)
//! - `e₀ = (e₋ - e₊) / 2` (origin)

use core::fmt;

use crate::algebra::Multivector;
use crate::basis::Blade;
use crate::scalar::Float;
use crate::signature::Conformal3;

use super::{FlatPoint, Point};

/// Error type for conversion from `Multivector` to specialized CGA types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConversionError {
    /// The multivector has non-zero coefficients in unexpected grades.
    InvalidGrade,
    /// The null constraint is violated (for Point).
    NullConstraintViolated,
}

impl fmt::Display for ConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGrade => write!(
                f,
                "multivector has non-zero coefficients in unexpected grades"
            ),
            Self::NullConstraintViolated => write!(f, "point does not satisfy null constraint"),
        }
    }
}

impl std::error::Error for ConversionError {}

// ============================================================================
// Blade indices for Conformal3 (Cl(4,1,0))
//
// Binary encoding: bit i = 1 means basis vector e_i is present
// Basis order: e₁ (bit 0), e₂ (bit 1), e₃ (bit 2), e₊ (bit 3), e₋ (bit 4)
// ============================================================================

/// Index for scalar (grade 0).
const SCALAR_IDX: usize = 0b00000; // 0

// Grade 1 basis vectors
/// Index for e₁ (grade 1).
const E1_IDX: usize = 0b00001; // 1
/// Index for e₂ (grade 1).
const E2_IDX: usize = 0b00010; // 2
/// Index for e₃ (grade 1).
const E3_IDX: usize = 0b00100; // 4
/// Index for e₊ (grade 1).
const EP_IDX: usize = 0b01000; // 8
/// Index for e₋ (grade 1).
const EM_IDX: usize = 0b10000; // 16

// Grade 2 bivectors (10 total)
/// Index for e₁₂ (grade 2).
const E12_IDX: usize = 0b00011; // 3
/// Index for e₁₃ (grade 2).
const E13_IDX: usize = 0b00101; // 5
/// Index for e₂₃ (grade 2).
const E23_IDX: usize = 0b00110; // 6
/// Index for e₁₊ (grade 2).
const E1P_IDX: usize = 0b01001; // 9
/// Index for e₂₊ (grade 2).
const E2P_IDX: usize = 0b01010; // 10
/// Index for e₃₊ (grade 2).
const E3P_IDX: usize = 0b01100; // 12
/// Index for e₁₋ (grade 2).
const E1M_IDX: usize = 0b10001; // 17
/// Index for e₂₋ (grade 2).
const E2M_IDX: usize = 0b10010; // 18
/// Index for e₃₋ (grade 2).
const E3M_IDX: usize = 0b10100; // 20
/// Index for e₊₋ (grade 2).
const EPM_IDX: usize = 0b11000; // 24

/// Tolerance for `TryFrom` grade validation.
///
/// This is intentionally larger than machine epsilon to accommodate
/// accumulated floating-point errors from arithmetic operations.
pub const CONVERSION_TOLERANCE: f64 = 1e-10;

// ============================================================================
// Point conversions (grade 1)
// ============================================================================

impl<T: Float> From<Point<T>> for Multivector<T, Conformal3> {
    /// Converts a CGA point to a generic multivector.
    ///
    /// The point components map to grade-1 blades:
    /// - `e1` → `e₁` (index 1)
    /// - `e2` → `e₂` (index 2)
    /// - `e3` → `e₃` (index 4)
    /// - `ep` → `e₊` (index 8)
    /// - `em` → `e₋` (index 16)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Conformal3;
    /// use clifford::specialized::conformal::dim3::Point;
    ///
    /// let p = Point::<f64>::new(1.0, 2.0, 3.0);
    /// let mv: Multivector<f64, Conformal3> = p.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(1)); // Pure vector
    /// ```
    fn from(p: Point<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E1_IDX), p.e1());
        mv.set(Blade::from_index(E2_IDX), p.e2());
        mv.set(Blade::from_index(E3_IDX), p.e3());
        mv.set(Blade::from_index(EP_IDX), p.ep());
        mv.set(Blade::from_index(EM_IDX), p.em());
        mv
    }
}

impl<T: Float> Point<T> {
    /// Converts a generic multivector to a CGA point without validation.
    ///
    /// This method extracts the grade-1 components directly without checking
    /// that other grades are zero or that the null constraint holds. It is
    /// branch-free and suitable for use with automatic differentiation.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Conformal3;
    /// use clifford::specialized::conformal::dim3::Point;
    ///
    /// let p = Point::<f64>::new(1.0, 2.0, 3.0);
    /// let mv: Multivector<f64, Conformal3> = p.into();
    /// let back = Point::from_multivector_unchecked(&mv);
    ///
    /// assert!((p.x() - back.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Conformal3>) -> Self {
        Self::from_conformal_unchecked(
            mv.get(Blade::from_index(E1_IDX)),
            mv.get(Blade::from_index(E2_IDX)),
            mv.get(Blade::from_index(E3_IDX)),
            mv.get(Blade::from_index(EP_IDX)),
            mv.get(Blade::from_index(EM_IDX)),
        )
    }
}

impl<T: Float> TryFrom<Multivector<T, Conformal3>> for Point<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a CGA point.
    ///
    /// Succeeds only if:
    /// 1. The multivector is a pure vector (grade 1 only)
    /// 2. The null constraint is satisfied (P · P ≈ 0)
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Point::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// - `ConversionError::InvalidGrade` if non-grade-1 components are present
    /// - `ConversionError::NullConstraintViolated` if P · P ≠ 0
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Conformal3;
    /// use clifford::specialized::conformal::dim3::Point;
    ///
    /// let p = Point::<f64>::new(1.0, 2.0, 3.0);
    /// let mv: Multivector<f64, Conformal3> = p.into();
    /// let back = Point::try_from(mv).unwrap();
    ///
    /// assert!((p.x() - back.x()).abs() < 1e-10);
    /// ```
    fn try_from(mv: Multivector<T, Conformal3>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only grade-1 components are non-zero
        // We check a representative sample of other grades
        if mv.get(Blade::from_index(SCALAR_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        // Check grade 2 (some representative bivectors)
        if mv.get(Blade::from_index(E12_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E13_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E23_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(EPM_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }

        let point = Self::from_multivector_unchecked(&mv);

        // Verify null constraint
        if !point.is_null(tolerance) {
            return Err(ConversionError::NullConstraintViolated);
        }

        Ok(point)
    }
}

// ============================================================================
// FlatPoint conversions (grade 2 - bivector with e∞ factor)
// ============================================================================

impl<T: Float> From<FlatPoint<T>> for Multivector<T, Conformal3> {
    /// Converts a flat point to a generic multivector.
    ///
    /// FlatPoint components are bivectors formed by wedging with e∞ = e₋ + e₊:
    /// - `e1i` (e₁∞) → `e₁₋ + e₁₊` (indices 17 and 9)
    /// - `e2i` (e₂∞) → `e₂₋ + e₂₊` (indices 18 and 10)
    /// - `e3i` (e₃∞) → `e₃₋ + e₃₊` (indices 20 and 12)
    /// - `e0i` (e₀∞) → `-e₊₋` (index 24, negated)
    ///
    /// # Sign Convention
    ///
    /// The e₀∞ component maps to `-e₊₋` because:
    /// ```text
    /// e₀∞ = e₀ ∧ e∞ = ((e₋ - e₊)/2) ∧ (e₋ + e₊)
    ///     = (e₋ ∧ e₊) = -e₊₋
    /// ```
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Conformal3;
    /// use clifford::specialized::conformal::dim3::FlatPoint;
    ///
    /// let fp = FlatPoint::<f64>::new(1.0, 2.0, 3.0);
    /// let mv: Multivector<f64, Conformal3> = fp.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(2)); // Pure bivector
    /// ```
    fn from(fp: FlatPoint<T>) -> Self {
        let mut mv = Multivector::zero();

        // e₁∞ = e₁ ∧ (e₋ + e₊) = e₁₋ + e₁₊
        mv.set(Blade::from_index(E1M_IDX), fp.e1i());
        mv.set(Blade::from_index(E1P_IDX), fp.e1i());

        // e₂∞ = e₂ ∧ (e₋ + e₊) = e₂₋ + e₂₊
        mv.set(Blade::from_index(E2M_IDX), fp.e2i());
        mv.set(Blade::from_index(E2P_IDX), fp.e2i());

        // e₃∞ = e₃ ∧ (e₋ + e₊) = e₃₋ + e₃₊
        mv.set(Blade::from_index(E3M_IDX), fp.e3i());
        mv.set(Blade::from_index(E3P_IDX), fp.e3i());

        // e₀∞ = e₀ ∧ e∞ = ((e₋ - e₊)/2) ∧ (e₋ + e₊) = e₋ ∧ e₊ = -e₊₋
        // So coefficient of e₊₋ is -e0i
        mv.set(Blade::from_index(EPM_IDX), -fp.e0i());

        mv
    }
}

impl<T: Float> FlatPoint<T> {
    /// Converts a generic multivector to a flat point without validation.
    ///
    /// This method extracts the grade-2 components that correspond to
    /// the flat point structure. It is branch-free and suitable for AD.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Conformal3;
    /// use clifford::specialized::conformal::dim3::FlatPoint;
    ///
    /// let fp = FlatPoint::<f64>::new(1.0, 2.0, 3.0);
    /// let mv: Multivector<f64, Conformal3> = fp.into();
    /// let back = FlatPoint::from_multivector_unchecked(&mv);
    ///
    /// assert!((fp.x() - back.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Conformal3>) -> Self {
        // Extract e₁∞ from e₁₋ (or e₁₊, they should be equal for valid flat points)
        let e1i = mv.get(Blade::from_index(E1M_IDX));
        let e2i = mv.get(Blade::from_index(E2M_IDX));
        let e3i = mv.get(Blade::from_index(E3M_IDX));
        // e₀∞ = -e₊₋, so negate to get the e0i component
        let e0i = -mv.get(Blade::from_index(EPM_IDX));

        Self::from_components_unchecked(e1i, e2i, e3i, e0i)
    }
}

impl<T: Float> TryFrom<Multivector<T, Conformal3>> for FlatPoint<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a flat point.
    ///
    /// Succeeds only if the multivector has the structure of a flat point
    /// (specific grade-2 pattern with e∞ factor).
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the structure doesn't match.
    fn try_from(mv: Multivector<T, Conformal3>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that scalar is zero
        if mv.get(Blade::from_index(SCALAR_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        // Check that grade-1 is zero
        if mv.get(Blade::from_index(E1_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E2_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E3_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(EP_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(EM_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }

        // Check Euclidean bivectors are zero (not part of flat point structure)
        if mv.get(Blade::from_index(E12_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E13_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E23_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }

        // Verify flat point structure: e₁₊ = e₁₋, e₂₊ = e₂₋, e₃₊ = e₃₋
        let e1p = mv.get(Blade::from_index(E1P_IDX));
        let e1m = mv.get(Blade::from_index(E1M_IDX));
        if (e1p - e1m).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        let e2p = mv.get(Blade::from_index(E2P_IDX));
        let e2m = mv.get(Blade::from_index(E2M_IDX));
        if (e2p - e2m).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        let e3p = mv.get(Blade::from_index(E3P_IDX));
        let e3m = mv.get(Blade::from_index(E3M_IDX));
        if (e3p - e3m).abs() > tolerance {
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

    // ========================================================================
    // Point tests
    // ========================================================================

    proptest! {
        /// Point roundtrip through Multivector preserves coordinates.
        #[test]
        fn point_roundtrip(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let p = Point::new(x, y, z);
            let mv: Multivector<f64, Conformal3> = p.into();
            let back = Point::try_from(mv).unwrap();

            prop_assert!(abs_diff_eq!(p.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Point converts to pure grade-1 multivector.
        #[test]
        fn point_to_multivector_is_grade_1(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let p = Point::new(x, y, z);
            let mv: Multivector<f64, Conformal3> = p.into();

            // Should be pure vector (grade 1)
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(2).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(3).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(4).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(5).is_zero(ABS_DIFF_EQ_EPS));
        }

        /// Point null constraint verified via Multivector inner product.
        #[test]
        fn point_null_via_multivector_inner(x in -10.0f64..10.0, y in -10.0f64..10.0, z in -10.0f64..10.0) {
            let p = Point::new(x, y, z);
            let mv: Multivector<f64, Conformal3> = p.into();

            // P · P = 0 for conformal points
            let inner = mv.inner(&mv);
            let norm_sq = inner.scalar_part();

            // Use relative tolerance for large coordinates
            let scale = (1.0 + x * x + y * y + z * z).max(1.0);
            prop_assert!(norm_sq.abs() < ABS_DIFF_EQ_EPS * scale,
                "null constraint violated: {} for point ({}, {}, {})", norm_sq, x, y, z);
        }

        /// Distance formula works via Multivector inner product.
        #[test]
        fn point_distance_via_multivector_inner(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);

            let mv1: Multivector<f64, Conformal3> = p1.into();
            let mv2: Multivector<f64, Conformal3> = p2.into();

            // For unit-weight conformal points: P₁ · P₂ = -½|p₁ - p₂|²
            let inner = mv1.inner(&mv2);
            let inner_scalar = inner.scalar_part();

            let euclidean_dist_sq = (x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2);
            let expected_inner = -euclidean_dist_sq / 2.0;

            prop_assert!(abs_diff_eq!(inner_scalar, expected_inner, epsilon = ABS_DIFF_EQ_EPS * 10.0),
                "inner product mismatch: {} vs expected {} for points ({},{},{}) and ({},{},{})",
                inner_scalar, expected_inner, x1, y1, z1, x2, y2, z2);
        }

        /// Point from_multivector_unchecked matches roundtrip.
        #[test]
        fn point_unchecked_roundtrip(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let p = Point::new(x, y, z);
            let mv: Multivector<f64, Conformal3> = p.into();
            let back = Point::from_multivector_unchecked(&mv);

            prop_assert!(abs_diff_eq!(p.e1(), back.e1(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.e2(), back.e2(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.e3(), back.e3(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.ep(), back.ep(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.em(), back.em(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn point_conversion_fails_with_scalar() {
        let mv: Multivector<f64, Conformal3> = Multivector::scalar(1.0);
        assert!(Point::try_from(mv).is_err());
    }

    #[test]
    fn point_origin_is_null() {
        let p = Point::<f64>::origin();
        let mv: Multivector<f64, Conformal3> = p.into();
        let inner = mv.inner(&mv);
        assert!(inner.scalar_part().abs() < ABS_DIFF_EQ_EPS);
    }

    // ========================================================================
    // FlatPoint tests
    // ========================================================================

    proptest! {
        /// FlatPoint roundtrip through Multivector preserves coordinates.
        #[test]
        fn flatpoint_roundtrip(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let fp = FlatPoint::new(x, y, z);
            let mv: Multivector<f64, Conformal3> = fp.into();
            let back = FlatPoint::try_from(mv).unwrap();

            prop_assert!(abs_diff_eq!(fp.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.z(), back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// FlatPoint converts to pure grade-2 multivector.
        #[test]
        fn flatpoint_to_multivector_is_grade_2(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let fp = FlatPoint::new(x, y, z);
            let mv: Multivector<f64, Conformal3> = fp.into();

            // Should be pure bivector (grade 2)
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(3).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(4).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(5).is_zero(ABS_DIFF_EQ_EPS));
        }

        /// FlatPoint unchecked roundtrip matches.
        #[test]
        fn flatpoint_unchecked_roundtrip(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let fp = FlatPoint::new(x, y, z);
            let mv: Multivector<f64, Conformal3> = fp.into();
            let back = FlatPoint::from_multivector_unchecked(&mv);

            prop_assert!(abs_diff_eq!(fp.e1i(), back.e1i(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.e2i(), back.e2i(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.e3i(), back.e3i(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.e0i(), back.e0i(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn flatpoint_conversion_fails_with_vector() {
        let mv: Multivector<f64, Conformal3> = Multivector::vector(&[1.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(FlatPoint::try_from(mv).is_err());
    }

    // ========================================================================
    // Cross-type consistency tests
    // ========================================================================

    proptest! {
        /// Point and FlatPoint relationship via outer product with e∞.
        ///
        /// FlatPoint = Point ∧ e∞ should hold algebraically.
        #[test]
        fn point_wedge_einf_equals_flatpoint(x in -10.0f64..10.0, y in -10.0f64..10.0, z in -10.0f64..10.0) {
            let p = Point::new(x, y, z);
            let fp = FlatPoint::from_round(&p);

            let mv_point: Multivector<f64, Conformal3> = p.into();
            let e_inf = Conformal3::e_infinity::<f64>();

            // Compute P ∧ e∞
            let wedge = mv_point.outer(&e_inf);

            let mv_fp: Multivector<f64, Conformal3> = fp.into();

            // They should be equal
            prop_assert!(abs_diff_eq!(wedge, mv_fp, epsilon = ABS_DIFF_EQ_EPS * 10.0),
                "Point ∧ e∞ ≠ FlatPoint for ({}, {}, {})", x, y, z);
        }

        /// Distance between specialized Point matches generic inner product.
        #[test]
        fn specialized_distance_matches_generic(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);

            // Specialized distance
            let specialized_dist = p1.distance(&p2);

            // Generic distance via inner product
            let mv1: Multivector<f64, Conformal3> = p1.into();
            let mv2: Multivector<f64, Conformal3> = p2.into();
            let inner = mv1.inner(&mv2).scalar_part();
            let generic_dist = (-2.0 * inner).max(0.0).sqrt();

            prop_assert!(abs_diff_eq!(specialized_dist, generic_dist, epsilon = ABS_DIFF_EQ_EPS * 10.0));
        }
    }
}
