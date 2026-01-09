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
//! The `TryFrom` implementations use `CONVERSION_TOLERANCE` (1e-10) rather than
//! machine epsilon to account for accumulated floating-point errors in computed
//! results. This is intentionally more lenient than `T::EPSILON`.

use core::fmt;

use crate::algebra::Multivector;
use crate::basis::Blade;
use crate::scalar::Float;
use crate::signature::Euclidean2;

use super::{Bivector, Rotor, Vector};

/// Tuple type for decomposed 2D multivector components.
///
/// Contains `(scalar, vector, bivector)` as Options. Each component is `Some`
/// if that grade has non-zero coefficients (above tolerance), `None` otherwise.
pub type Specialized<T> = (Option<T>, Option<Vector<T>>, Option<Bivector<T>>);

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
// Vector conversions
// ============================================================================

impl<T: Float> From<Vector<T>> for Multivector<T, Euclidean2> {
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
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let v = Vector::new(3.0_f64, 4.0);
    /// let mv: Multivector<f64, Euclidean2> = v.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(1)); // Pure vector
    /// ```
    fn from(v: Vector<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E1_IDX), v.x);
        mv.set(Blade::from_index(E2_IDX), v.y);
        mv
    }
}

impl<T: Float> Vector<T> {
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
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[3.0, 4.0]);
    /// let v = Vector::from_multivector_unchecked(&mv);
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

impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Vector<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 2D vector.
    ///
    /// Succeeds only if the multivector is a pure vector (grade 1 only),
    /// using `CONVERSION_TOLERANCE` to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Vector::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// scalar or bivector components (above `CONVERSION_TOLERANCE`).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[3.0, 4.0]);
    /// let v = Vector::try_from(mv).unwrap();
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
// Bivector conversions
// ============================================================================

impl<T: Float> From<Bivector<T>> for Multivector<T, Euclidean2> {
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
    /// use clifford::specialized::euclidean::dim2::Bivector;
    ///
    /// let b = Bivector::new(5.0_f64);
    /// let mv: Multivector<f64, Euclidean2> = b.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(2)); // Pure bivector
    /// ```
    fn from(b: Bivector<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E12_IDX), b.0);
        mv
    }
}

impl<T: Float> Bivector<T> {
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
    /// use clifford::specialized::euclidean::dim2::Bivector;
    ///
    /// let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
    /// mv.set(Blade::from_index(3), 5.0); // e₁₂
    /// let b = Bivector::from_multivector_unchecked(&mv);
    /// assert_eq!(b.value(), 5.0);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean2>) -> Self {
        Self::new(mv.get(Blade::from_index(E12_IDX)))
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Bivector<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 2D bivector.
    ///
    /// Succeeds only if the multivector is a pure bivector (grade 2 only),
    /// using `CONVERSION_TOLERANCE` to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Bivector::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// scalar or vector components (above `CONVERSION_TOLERANCE`).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::euclidean::dim2::Bivector;
    ///
    /// let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
    /// mv.set(Blade::from_index(3), 5.0); // e₁₂
    /// let b = Bivector::try_from(mv).unwrap();
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
// Rotor conversions
// ============================================================================

impl<T: Float> From<Rotor<T>> for Multivector<T, Euclidean2> {
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
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let r = Rotor::from_angle(FRAC_PI_4);
    /// let mv: Multivector<f64, Euclidean2> = r.into();
    ///
    /// // Rotor has scalar and bivector parts only
    /// assert!(mv.grade_select(1).is_zero(1e-10));
    /// ```
    fn from(r: Rotor<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(SCALAR_IDX), r.s);
        mv.set(Blade::from_index(E12_IDX), r.xy);
        mv
    }
}

impl<T: Float> Rotor<T> {
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
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let original = Rotor::from_angle(FRAC_PI_4);
    /// let mv: Multivector<f64, Euclidean2> = original.into();
    /// let recovered = Rotor::from_multivector_unchecked(&mv);
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

impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Rotor<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 2D rotor.
    ///
    /// Succeeds only if the multivector is an even element (scalar + bivector only),
    /// using `CONVERSION_TOLERANCE` to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Rotor::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// vector components (above `CONVERSION_TOLERANCE`).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let original = Rotor::from_angle(FRAC_PI_4);
    /// let mv: Multivector<f64, Euclidean2> = original.into();
    /// let recovered = Rotor::try_from(mv).unwrap();
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

// ============================================================================
// Multivector to_specialized decomposition
// ============================================================================

impl<T: Float> Multivector<T, Euclidean2> {
    /// Decomposes a 2D multivector into its specialized grade components.
    ///
    /// Returns a tuple `(scalar, vector, bivector)` where each component is `Some`
    /// if that grade has non-zero coefficients (above tolerance), `None` otherwise.
    ///
    /// This method is branch-free in extraction (uses `from_multivector_unchecked`)
    /// but branches on the tolerance check to determine `Some` vs `None`.
    ///
    /// # Arguments
    ///
    /// * `tolerance` - Threshold below which coefficients are considered zero
    ///
    /// # Returns
    ///
    /// A [`Specialized`] tuple containing:
    /// - `Option<T>`: Scalar (grade 0) if non-zero
    /// - `Option<Vector<T>>`: Vector (grade 1) if non-zero
    /// - `Option<Bivector<T>>`: Bivector (grade 2) if non-zero
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// // Create a pure vector multivector
    /// let v = Vector::new(3.0_f64, 4.0);
    /// let mv: Multivector<f64, Euclidean2> = v.into();
    ///
    /// let (scalar, vector, bivector) = mv.to_specialized(1e-10);
    ///
    /// assert!(scalar.is_none()); // No scalar part
    /// assert!(vector.is_some()); // Has vector part
    /// assert!(bivector.is_none()); // No bivector part
    ///
    /// let vec = vector.unwrap();
    /// assert!((vec.x - 3.0).abs() < 1e-10);
    /// assert!((vec.y - 4.0).abs() < 1e-10);
    /// ```
    ///
    /// # Example: Mixed multivector
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean2;
    /// use clifford::specialized::euclidean::dim2::{Vector, Bivector};
    ///
    /// // Create a mixed multivector: 2 + 3e₁ + 4e₂ + 5e₁₂
    /// let mut mv: Multivector<f64, Euclidean2> = Multivector::scalar(2.0);
    /// mv = &mv + &Multivector::from(Vector::new(3.0, 4.0));
    /// mv = &mv + &Multivector::from(Bivector::new(5.0));
    ///
    /// let (scalar, vector, bivector) = mv.to_specialized(1e-10);
    ///
    /// assert_eq!(scalar, Some(2.0));
    /// assert_eq!(vector.map(|v| (v.x, v.y)), Some((3.0, 4.0)));
    /// assert_eq!(bivector.map(|b| b.value()), Some(5.0));
    /// ```
    pub fn to_specialized(&self, tolerance: T) -> Specialized<T> {
        // Extract scalar (grade 0)
        let scalar_val = self.get(Blade::from_index(SCALAR_IDX));
        let scalar = if scalar_val.abs() > tolerance {
            Some(scalar_val)
        } else {
            None
        };

        // Extract vector (grade 1)
        let vec = Vector::from_multivector_unchecked(self);
        let vector = if vec.x.abs() > tolerance || vec.y.abs() > tolerance {
            Some(vec)
        } else {
            None
        };

        // Extract bivector (grade 2)
        let biv = Bivector::from_multivector_unchecked(self);
        let bivector = if biv.0.abs() > tolerance {
            Some(biv)
        } else {
            None
        };

        (scalar, vector, bivector)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;
    use std::f64::consts::FRAC_PI_4;

    use super::super::arbitrary::{UnitRotor, UnitVector};

    // ========================================================================
    // Vector tests
    // ========================================================================

    proptest! {
        #[test]
        fn vec2_roundtrip(v in any::<Vector<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = v.into();
            let back = Vector::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(v, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec2_to_multivector_is_grade_1(v in any::<Vector<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = v.into();
            // Should be pure vector (grade 1) or zero
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(2).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn vec2_conversion_fails_with_scalar() {
        let mv: Multivector<f64, Euclidean2> = Multivector::scalar(1.0);
        assert!(Vector::try_from(mv).is_err());
    }

    // ========================================================================
    // Bivector tests
    // ========================================================================

    proptest! {
        #[test]
        fn bivec2_roundtrip(val in -100.0f64..100.0) {
            let b = Bivector::new(val);
            let mv: Multivector<f64, Euclidean2> = b.into();
            let back = Bivector::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(b, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivec2_to_multivector_is_grade_2(val in -100.0f64..100.0) {
            let b = Bivector::new(val);
            let mv: Multivector<f64, Euclidean2> = b.into();
            // Should be pure bivector (grade 2) or zero
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn bivec2_conversion_fails_with_vector() {
        let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[1.0, 0.0]);
        assert!(Bivector::try_from(mv).is_err());
    }

    // ========================================================================
    // Rotor tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor2_roundtrip(r in any::<UnitRotor<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = (*r).into();
            let back = Rotor::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(*r, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor2_to_multivector_is_even(r in any::<UnitRotor<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = (*r).into();
            // Should have no vector (odd) part
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn rotor2_conversion_fails_with_vector() {
        let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[1.0, 0.0]);
        assert!(Rotor::try_from(mv).is_err());
    }

    // ========================================================================
    // Operation consistency tests
    // ========================================================================

    proptest! {
        #[test]
        fn dot_consistency(
            a in any::<Vector<f64>>(),
            b in any::<Vector<f64>>(),
        ) {
            let spec_result = a.dot(b);

            let gen_a: Multivector<f64, Euclidean2> = a.into();
            let gen_b: Multivector<f64, Euclidean2> = b.into();
            let gen_result = gen_a.inner(&gen_b).scalar_part();

            prop_assert!(abs_diff_eq!(spec_result, gen_result, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn wedge_consistency(
            a in any::<Vector<f64>>(),
            b in any::<Vector<f64>>(),
        ) {
            let spec_result = a.wedge(b);

            let gen_a: Multivector<f64, Euclidean2> = a.into();
            let gen_b: Multivector<f64, Euclidean2> = b.into();
            let gen_result = gen_a.outer(&gen_b);
            let gen_as_bivec = Bivector::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_bivec, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_rotation_consistency(
            r in any::<UnitRotor<f64>>(),
            v in any::<UnitVector<f64>>(),
        ) {
            // Specialized rotation
            let spec_result = r.rotate(*v);

            // Generic sandwich product: R̃ v R
            let gen_r: Multivector<f64, Euclidean2> = (*r).into();
            let gen_v: Multivector<f64, Euclidean2> = (*v).into();
            let gen_r_rev = gen_r.reverse();
            let gen_result = &(&gen_r_rev * &gen_v) * &gen_r;
            let gen_as_vec = Vector::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_vec, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn rotor_from_angle_consistency() {
        let r = Rotor::from_angle(FRAC_PI_4);
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
        /// Test that a generic vector Multivector converts to Vector correctly.
        #[test]
        fn generic_vector_to_vec2(x in -100.0f64..100.0, y in -100.0f64..100.0) {
            let mv: Multivector<f64, Euclidean2> = Multivector::vector(&[x, y]);
            let v = Vector::try_from(mv).expect("pure vector should convert");
            prop_assert!(abs_diff_eq!(v.x, x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(v.y, y, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic bivector Multivector converts to Bivector correctly.
        #[test]
        fn generic_bivector_to_bivec2(xy in -100.0f64..100.0) {
            let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
            mv.set(Blade::from_index(E12_IDX), xy);
            let b = Bivector::try_from(mv).expect("pure bivector should convert");
            prop_assert!(abs_diff_eq!(b.value(), xy, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic even Multivector converts to Rotor correctly.
        #[test]
        fn generic_even_to_rotor2(s in -100.0f64..100.0, xy in -100.0f64..100.0) {
            let mut mv: Multivector<f64, Euclidean2> = Multivector::zero();
            mv.set(Blade::from_index(SCALAR_IDX), s);
            mv.set(Blade::from_index(E12_IDX), xy);
            let r = Rotor::try_from(mv).expect("even element should convert");
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
                prop_assert!(Vector::try_from(mv.clone()).is_err());
            }
            if has_bivec && (has_scalar || has_vec) {
                // Has bivector and scalar or vector - should fail bivec conversion
                prop_assert!(Bivector::try_from(mv).is_err());
            }
        }

        /// Test that vector wedge product (generic) converts to Bivector.
        #[test]
        fn generic_wedge_to_bivec2(
            ax in -10.0f64..10.0, ay in -10.0f64..10.0,
            bx in -10.0f64..10.0, by in -10.0f64..10.0,
        ) {
            let a: Multivector<f64, Euclidean2> = Multivector::vector(&[ax, ay]);
            let b: Multivector<f64, Euclidean2> = Multivector::vector(&[bx, by]);
            let wedge = a.outer(&b);

            // Should be a pure bivector, convertible to Bivector
            let bivec = Bivector::try_from(wedge).expect("wedge of vectors should be bivector");

            // Check value matches ax*by - ay*bx
            let expected = ax * by - ay * bx;
            prop_assert!(abs_diff_eq!(bivec.value(), expected, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that vector geometric product (generic) converts to Rotor.
        #[test]
        fn generic_geometric_to_rotor2(
            ax in -10.0f64..10.0, ay in -10.0f64..10.0,
            bx in -10.0f64..10.0, by in -10.0f64..10.0,
        ) {
            let a: Multivector<f64, Euclidean2> = Multivector::vector(&[ax, ay]);
            let b: Multivector<f64, Euclidean2> = Multivector::vector(&[bx, by]);
            let product = &a * &b;

            // Should be scalar + bivector (even), convertible to Rotor
            let rotor = Rotor::try_from(product).expect("geometric product should be even");

            // Check components
            let expected_s = ax * bx + ay * by; // dot product
            let expected_xy = ax * by - ay * bx; // wedge product
            prop_assert!(abs_diff_eq!(rotor.s, expected_s, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotor.xy, expected_xy, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // to_specialized tests
    // ========================================================================

    proptest! {
        /// Test that to_specialized correctly identifies pure vectors.
        #[test]
        fn to_specialized_pure_vector(v in any::<Vector<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = v.into();
            let (scalar, vector, bivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            prop_assert!(scalar.is_none());
            prop_assert!(bivector.is_none());

            if v.x.abs().max(v.y.abs()) > ABS_DIFF_EQ_EPS {
                prop_assert!(vector.is_some());
                let vec = vector.unwrap();
                prop_assert!(abs_diff_eq!(vec, v, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized correctly identifies pure bivectors.
        #[test]
        fn to_specialized_pure_bivector(xy in -100.0f64..100.0) {
            let b = Bivector::new(xy);
            let mv: Multivector<f64, Euclidean2> = b.into();
            let (scalar, vector, bivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            prop_assert!(scalar.is_none());
            prop_assert!(vector.is_none());

            if xy.abs() > ABS_DIFF_EQ_EPS {
                prop_assert!(bivector.is_some());
                let biv = bivector.unwrap();
                prop_assert!(abs_diff_eq!(biv, b, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized correctly identifies pure scalars.
        #[test]
        fn to_specialized_pure_scalar(s in -100.0f64..100.0) {
            let mv: Multivector<f64, Euclidean2> = Multivector::scalar(s);
            let (scalar, vector, bivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            prop_assert!(vector.is_none());
            prop_assert!(bivector.is_none());

            if s.abs() > ABS_DIFF_EQ_EPS {
                prop_assert!(scalar.is_some());
                prop_assert!(abs_diff_eq!(scalar.unwrap(), s, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized correctly identifies rotors (even multivectors).
        #[test]
        fn to_specialized_rotor(r in any::<UnitRotor<f64>>()) {
            let mv: Multivector<f64, Euclidean2> = (*r).into();
            let (scalar, vector, bivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            // Rotor has scalar and bivector, no vector
            prop_assert!(vector.is_none());

            // Check scalar part
            if r.s.abs() > ABS_DIFF_EQ_EPS {
                prop_assert!(scalar.is_some());
                prop_assert!(abs_diff_eq!(scalar.unwrap(), r.s, epsilon = ABS_DIFF_EQ_EPS));
            }

            // Check bivector part
            if r.xy.abs() > ABS_DIFF_EQ_EPS {
                prop_assert!(bivector.is_some());
                prop_assert!(abs_diff_eq!(bivector.unwrap().value(), r.xy, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized handles arbitrary multivectors correctly.
        #[test]
        fn to_specialized_arbitrary(mv in any::<Multivector<f64, Euclidean2>>()) {
            let (scalar, vector, bivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            // Reconstruct from components and verify
            let mut reconstructed: Multivector<f64, Euclidean2> = Multivector::zero();

            if let Some(s) = scalar {
                reconstructed = &reconstructed + &Multivector::scalar(s);
            }
            if let Some(v) = vector {
                reconstructed = &reconstructed + &Multivector::from(v);
            }
            if let Some(b) = bivector {
                reconstructed = &reconstructed + &Multivector::from(b);
            }

            // Reconstructed should be approximately equal to original (within tolerance)
            // Note: small values below tolerance are zeroed out, so we check equivalence
            for i in 0..4 {
                let orig = mv.get(Blade::from_index(i));
                let recon = reconstructed.get(Blade::from_index(i));
                if orig.abs() > ABS_DIFF_EQ_EPS {
                    prop_assert!(abs_diff_eq!(orig, recon, epsilon = ABS_DIFF_EQ_EPS));
                } else {
                    // Original was below tolerance, reconstructed should be zero
                    prop_assert!(abs_diff_eq!(recon, 0.0, epsilon = ABS_DIFF_EQ_EPS));
                }
            }
        }
    }

    #[test]
    fn to_specialized_zero_multivector() {
        let mv: Multivector<f64, Euclidean2> = Multivector::zero();
        let (scalar, vector, bivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

        assert!(scalar.is_none());
        assert!(vector.is_none());
        assert!(bivector.is_none());
    }

    #[test]
    fn to_specialized_full_multivector() {
        // Create a multivector with all grades non-zero
        let mut mv: Multivector<f64, Euclidean2> = Multivector::scalar(2.0);
        mv = &mv + &Multivector::from(Vector::new(3.0, 4.0));
        mv = &mv + &Multivector::from(Bivector::new(5.0));

        let (scalar, vector, bivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

        assert_eq!(scalar, Some(2.0));
        assert!(vector.is_some());
        let v = vector.unwrap();
        assert!(abs_diff_eq!(v.x, 3.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(v.y, 4.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(bivector.is_some());
        assert!(abs_diff_eq!(
            bivector.unwrap().value(),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }
}
