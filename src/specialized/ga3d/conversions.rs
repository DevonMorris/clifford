//! Conversions between 3D specialized types and generic `Multivector`.
//!
//! This module provides bidirectional conversions:
//! - `From<SpecializedType> for Multivector<T, Euclidean3>` (infallible)
//! - `TryFrom<Multivector<T, Euclidean3>> for SpecializedType` (fallible, validates structure)
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
use crate::signature::Euclidean3;

use super::{Bivec3, Even3, Rotor3, Trivec3, Vec3};

/// Tuple type for decomposed 3D multivector components.
///
/// Contains `(scalar, vector, bivector, trivector)` as Options. Each component is `Some`
/// if that grade has non-zero coefficients (above tolerance), `None` otherwise.
pub type Specialized3<T> = (
    Option<T>,
    Option<Vec3<T>>,
    Option<Bivec3<T>>,
    Option<Trivec3<T>>,
);

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
// Blade indices for Euclidean3
// ============================================================================

/// Index for scalar (grade 0).
const SCALAR_IDX: usize = 0;
/// Index for e₁ (grade 1).
const E1_IDX: usize = 1;
/// Index for e₂ (grade 1).
const E2_IDX: usize = 2;
/// Index for e₁₂ (grade 2).
const E12_IDX: usize = 3;
/// Index for e₃ (grade 1).
const E3_IDX: usize = 4;
/// Index for e₁₃ (grade 2).
const E13_IDX: usize = 5;
/// Index for e₂₃ (grade 2).
const E23_IDX: usize = 6;
/// Index for e₁₂₃ (grade 3).
const E123_IDX: usize = 7;

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
// Vec3 conversions
// ============================================================================

impl<T: Float> From<Vec3<T>> for Multivector<T, Euclidean3> {
    /// Converts a 3D vector to a generic multivector.
    ///
    /// The vector components map to grade-1 blades:
    /// - `x` → `e₁` (index 1)
    /// - `y` → `e₂` (index 2)
    /// - `z` → `e₃` (index 4)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Vec3;
    ///
    /// let v = Vec3::new(1.0_f64, 2.0, 3.0);
    /// let mv: Multivector<f64, Euclidean3> = v.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(1)); // Pure vector
    /// ```
    fn from(v: Vec3<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E1_IDX), v.x);
        mv.set(Blade::from_index(E2_IDX), v.y);
        mv.set(Blade::from_index(E3_IDX), v.z);
        mv
    }
}

impl<T: Float> Vec3<T> {
    /// Converts a generic multivector to a 3D vector without validation.
    ///
    /// This method extracts the grade-1 components directly without checking
    /// that other grades are zero. It is branch-free and suitable for use
    /// with automatic differentiation (dual numbers).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Vec3;
    ///
    /// let mv: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    /// let v = Vec3::from_multivector_unchecked(&mv);
    /// assert_eq!(v.x, 1.0);
    /// assert_eq!(v.y, 2.0);
    /// assert_eq!(v.z, 3.0);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean3>) -> Self {
        Self::new(
            mv.get(Blade::from_index(E1_IDX)),
            mv.get(Blade::from_index(E2_IDX)),
            mv.get(Blade::from_index(E3_IDX)),
        )
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Vec3<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 3D vector.
    ///
    /// Succeeds only if the multivector is a pure vector (grade 1 only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Vec3::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// scalar, bivector, or trivector components (above [`CONVERSION_TOLERANCE`]).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Vec3;
    ///
    /// let mv: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    /// let v = Vec3::try_from(mv).unwrap();
    /// assert_eq!(v.x, 1.0);
    /// assert_eq!(v.y, 2.0);
    /// assert_eq!(v.z, 3.0);
    /// ```
    fn try_from(mv: Multivector<T, Euclidean3>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only grade-1 components are non-zero
        // Grade 0
        if mv.get(Blade::from_index(SCALAR_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 2
        if mv.get(Blade::from_index(E12_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E13_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E23_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 3
        if mv.get(Blade::from_index(E123_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

// ============================================================================
// Bivec3 conversions
// ============================================================================

impl<T: Float> From<Bivec3<T>> for Multivector<T, Euclidean3> {
    /// Converts a 3D bivector to a generic multivector.
    ///
    /// The bivector components map to grade-2 blades:
    /// - `xy` → `e₁₂` (index 3)
    /// - `xz` → `e₁₃` (index 5)
    /// - `yz` → `e₂₃` (index 6)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Bivec3;
    ///
    /// let b = Bivec3::new(1.0_f64, 2.0, 3.0);
    /// let mv: Multivector<f64, Euclidean3> = b.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(2)); // Pure bivector
    /// ```
    fn from(b: Bivec3<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E12_IDX), b.xy);
        mv.set(Blade::from_index(E13_IDX), b.xz);
        mv.set(Blade::from_index(E23_IDX), b.yz);
        mv
    }
}

impl<T: Float> Bivec3<T> {
    /// Converts a generic multivector to a 3D bivector without validation.
    ///
    /// This method extracts the grade-2 components directly without checking
    /// that other grades are zero. It is branch-free and suitable for use
    /// with automatic differentiation (dual numbers).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Bivec3;
    ///
    /// let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
    /// mv.set(Blade::from_index(3), 1.0); // e12
    /// mv.set(Blade::from_index(5), 2.0); // e13
    /// mv.set(Blade::from_index(6), 3.0); // e23
    /// let b = Bivec3::from_multivector_unchecked(&mv);
    /// assert_eq!(b.xy, 1.0);
    /// assert_eq!(b.xz, 2.0);
    /// assert_eq!(b.yz, 3.0);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean3>) -> Self {
        Self::new(
            mv.get(Blade::from_index(E12_IDX)),
            mv.get(Blade::from_index(E13_IDX)),
            mv.get(Blade::from_index(E23_IDX)),
        )
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Bivec3<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 3D bivector.
    ///
    /// Succeeds only if the multivector is a pure bivector (grade 2 only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Bivec3::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// scalar, vector, or trivector components (above [`CONVERSION_TOLERANCE`]).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Bivec3;
    ///
    /// let b = Bivec3::new(1.0_f64, 2.0, 3.0);
    /// let mv: Multivector<f64, Euclidean3> = b.into();
    /// let back = Bivec3::try_from(mv).unwrap();
    /// assert_eq!(back.xy, 1.0);
    /// assert_eq!(back.xz, 2.0);
    /// assert_eq!(back.yz, 3.0);
    /// ```
    fn try_from(mv: Multivector<T, Euclidean3>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only grade-2 components are non-zero
        // Grade 0
        if mv.get(Blade::from_index(SCALAR_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 1
        if mv.get(Blade::from_index(E1_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E2_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E3_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 3
        if mv.get(Blade::from_index(E123_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

// ============================================================================
// Trivec3 conversions
// ============================================================================

impl<T: Float> From<Trivec3<T>> for Multivector<T, Euclidean3> {
    /// Converts a 3D trivector to a generic multivector.
    ///
    /// The trivector coefficient maps to the grade-3 blade:
    /// - `value` → `e₁₂₃` (index 7)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Trivec3;
    ///
    /// let t = Trivec3::new(5.0_f64);
    /// let mv: Multivector<f64, Euclidean3> = t.into();
    ///
    /// assert_eq!(mv.grade(1e-10), Some(3)); // Pure trivector
    /// ```
    fn from(t: Trivec3<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E123_IDX), t.0);
        mv
    }
}

impl<T: Float> Trivec3<T> {
    /// Converts a generic multivector to a 3D trivector without validation.
    ///
    /// This method extracts the grade-3 component directly without checking
    /// that other grades are zero. It is branch-free and suitable for use
    /// with automatic differentiation (dual numbers).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Trivec3;
    ///
    /// let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
    /// mv.set(Blade::from_index(7), 5.0); // e123
    /// let t = Trivec3::from_multivector_unchecked(&mv);
    /// assert_eq!(t.value(), 5.0);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean3>) -> Self {
        Self::new(mv.get(Blade::from_index(E123_IDX)))
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Trivec3<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 3D trivector.
    ///
    /// Succeeds only if the multivector is a pure trivector (grade 3 only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Trivec3::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// scalar, vector, or bivector components (above [`CONVERSION_TOLERANCE`]).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Trivec3;
    ///
    /// let t = Trivec3::new(5.0_f64);
    /// let mv: Multivector<f64, Euclidean3> = t.into();
    /// let back = Trivec3::try_from(mv).unwrap();
    /// assert_eq!(back.value(), 5.0);
    /// ```
    fn try_from(mv: Multivector<T, Euclidean3>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only grade-3 component is non-zero
        // Grade 0
        if mv.get(Blade::from_index(SCALAR_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 1
        if mv.get(Blade::from_index(E1_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E2_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E3_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 2
        if mv.get(Blade::from_index(E12_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E13_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E23_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

// ============================================================================
// Rotor3 conversions
// ============================================================================

impl<T: Float> From<Rotor3<T>> for Multivector<T, Euclidean3> {
    /// Converts a 3D rotor to a generic multivector.
    ///
    /// The rotor components map to:
    /// - `s` → scalar (index 0)
    /// - `b.xy` → `e₁₂` (index 3)
    /// - `b.xz` → `e₁₃` (index 5)
    /// - `b.yz` → `e₂₃` (index 6)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::{Bivec3, Rotor3};
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let r = Rotor3::from_angle_plane(FRAC_PI_4, Bivec3::unit_xy());
    /// let mv: Multivector<f64, Euclidean3> = r.into();
    ///
    /// // Rotor has scalar and bivector parts only (even grades)
    /// assert!(mv.grade_select(1).is_zero(1e-10));
    /// assert!(mv.grade_select(3).is_zero(1e-10));
    /// ```
    fn from(r: Rotor3<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(SCALAR_IDX), r.s);
        mv.set(Blade::from_index(E12_IDX), r.b.xy);
        mv.set(Blade::from_index(E13_IDX), r.b.xz);
        mv.set(Blade::from_index(E23_IDX), r.b.yz);
        mv
    }
}

impl<T: Float> Rotor3<T> {
    /// Converts a generic multivector to a 3D rotor without validation.
    ///
    /// This method extracts the even-grade components directly without checking
    /// that odd grades are zero. It is branch-free and suitable for use
    /// with automatic differentiation (dual numbers).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Rotor3;
    ///
    /// let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
    /// mv.set(Blade::from_index(0), 0.9238); // scalar (cos(π/8))
    /// mv.set(Blade::from_index(3), 0.3827); // e12 (sin(π/8))
    /// let r = Rotor3::from_multivector_unchecked(&mv);
    /// assert!((r.s - 0.9238).abs() < 1e-4);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean3>) -> Self {
        Self::new(
            mv.get(Blade::from_index(SCALAR_IDX)),
            Bivec3::new(
                mv.get(Blade::from_index(E12_IDX)),
                mv.get(Blade::from_index(E13_IDX)),
                mv.get(Blade::from_index(E23_IDX)),
            ),
        )
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Rotor3<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 3D rotor.
    ///
    /// Succeeds only if the multivector is an even element (scalar + bivector only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Rotor3::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// vector or trivector components (above [`CONVERSION_TOLERANCE`]).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::{Bivec3, Rotor3};
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let original = Rotor3::from_angle_plane(FRAC_PI_4, Bivec3::unit_xy());
    /// let mv: Multivector<f64, Euclidean3> = original.into();
    /// let recovered = Rotor3::try_from(mv).unwrap();
    ///
    /// assert!((original.s - recovered.s).abs() < 1e-10);
    /// ```
    fn try_from(mv: Multivector<T, Euclidean3>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only even-grade components are non-zero
        // Grade 1
        if mv.get(Blade::from_index(E1_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E2_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E3_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 3
        if mv.get(Blade::from_index(E123_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

// ============================================================================
// Even3 conversions
// ============================================================================

impl<T: Float> From<Even3<T>> for Multivector<T, Euclidean3> {
    /// Converts a 3D even multivector to a generic multivector.
    ///
    /// The components map identically to `Rotor3`:
    /// - `s` → scalar (index 0)
    /// - `b.xy` → `e₁₂` (index 3)
    /// - `b.xz` → `e₁₃` (index 5)
    /// - `b.yz` → `e₂₃` (index 6)
    fn from(e: Even3<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(SCALAR_IDX), e.s);
        mv.set(Blade::from_index(E12_IDX), e.b.xy);
        mv.set(Blade::from_index(E13_IDX), e.b.xz);
        mv.set(Blade::from_index(E23_IDX), e.b.yz);
        mv
    }
}

impl<T: Float> Even3<T> {
    /// Converts a generic multivector to a 3D even multivector without validation.
    ///
    /// This method extracts the even-grade components directly without checking
    /// that odd grades are zero. It is branch-free and suitable for use
    /// with automatic differentiation (dual numbers).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Even3;
    ///
    /// let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
    /// mv.set(Blade::from_index(0), 1.0);  // scalar
    /// mv.set(Blade::from_index(3), 0.5); // e12
    /// let e = Even3::from_multivector_unchecked(&mv);
    /// assert_eq!(e.s, 1.0);
    /// assert_eq!(e.b.xy, 0.5);
    /// ```
    #[inline]
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean3>) -> Self {
        Self::new(
            mv.get(Blade::from_index(SCALAR_IDX)),
            Bivec3::new(
                mv.get(Blade::from_index(E12_IDX)),
                mv.get(Blade::from_index(E13_IDX)),
                mv.get(Blade::from_index(E23_IDX)),
            ),
        )
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Even3<T> {
    type Error = ConversionError;

    /// Attempts to convert a generic multivector to a 3D even multivector.
    ///
    /// Succeeds only if the multivector is an even element (scalar + bivector only),
    /// using [`CONVERSION_TOLERANCE`] to check for near-zero components.
    ///
    /// # Note on Automatic Differentiation
    ///
    /// This method uses branching for validation and is **not suitable for AD**.
    /// Use [`Even3::from_multivector_unchecked`] for AD-compatible conversions.
    ///
    /// # Errors
    ///
    /// Returns `ConversionError::InvalidGrade` if the multivector has non-zero
    /// vector or trivector components (above [`CONVERSION_TOLERANCE`]).
    fn try_from(mv: Multivector<T, Euclidean3>) -> Result<Self, Self::Error> {
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check that only even-grade components are non-zero
        // Grade 1
        if mv.get(Blade::from_index(E1_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E2_IDX)).abs() > tolerance
            || mv.get(Blade::from_index(E3_IDX)).abs() > tolerance
        {
            return Err(ConversionError::InvalidGrade);
        }
        // Grade 3
        if mv.get(Blade::from_index(E123_IDX)).abs() > tolerance {
            return Err(ConversionError::InvalidGrade);
        }

        Ok(Self::from_multivector_unchecked(&mv))
    }
}

// ============================================================================
// Multivector to_specialized decomposition
// ============================================================================

impl<T: Float> Multivector<T, Euclidean3> {
    /// Decomposes a 3D multivector into its specialized grade components.
    ///
    /// Returns a tuple `(scalar, vector, bivector, trivector)` where each component
    /// is `Some` if that grade has non-zero coefficients (above tolerance), `None`
    /// otherwise.
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
    /// A [`Specialized3`] tuple containing:
    /// - `Option<T>`: Scalar (grade 0) if non-zero
    /// - `Option<Vec3<T>>`: Vector (grade 1) if non-zero
    /// - `Option<Bivec3<T>>`: Bivector (grade 2) if non-zero
    /// - `Option<Trivec3<T>>`: Trivector (grade 3) if non-zero
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::Vec3;
    ///
    /// // Create a pure vector multivector
    /// let v = Vec3::new(1.0_f64, 2.0, 3.0);
    /// let mv: Multivector<f64, Euclidean3> = v.into();
    ///
    /// let (scalar, vector, bivector, trivector) = mv.to_specialized(1e-10);
    ///
    /// assert!(scalar.is_none()); // No scalar part
    /// assert!(vector.is_some()); // Has vector part
    /// assert!(bivector.is_none()); // No bivector part
    /// assert!(trivector.is_none()); // No trivector part
    ///
    /// let vec = vector.unwrap();
    /// assert!((vec.x - 1.0).abs() < 1e-10);
    /// assert!((vec.y - 2.0).abs() < 1e-10);
    /// assert!((vec.z - 3.0).abs() < 1e-10);
    /// ```
    ///
    /// # Example: Mixed multivector
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use clifford::specialized::ga3d::{Vec3, Bivec3, Trivec3};
    ///
    /// // Create a mixed multivector: 2 + 1e₁ + 2e₂ + 3e₃ + e₁₂ + 5e₁₂₃
    /// let mut mv: Multivector<f64, Euclidean3> = Multivector::scalar(2.0);
    /// mv = &mv + &Multivector::from(Vec3::new(1.0, 2.0, 3.0));
    /// mv = &mv + &Multivector::from(Bivec3::new(1.0, 0.0, 0.0));
    /// mv = &mv + &Multivector::from(Trivec3::new(5.0));
    ///
    /// let (scalar, vector, bivector, trivector) = mv.to_specialized(1e-10);
    ///
    /// assert_eq!(scalar, Some(2.0));
    /// assert!(vector.is_some());
    /// assert!(bivector.is_some());
    /// assert_eq!(trivector.map(|t| t.value()), Some(5.0));
    /// ```
    pub fn to_specialized(&self, tolerance: T) -> Specialized3<T> {
        // Extract scalar (grade 0)
        let scalar_val = self.get(Blade::from_index(SCALAR_IDX));
        let scalar = if scalar_val.abs() > tolerance {
            Some(scalar_val)
        } else {
            None
        };

        // Extract vector (grade 1)
        let vec = Vec3::from_multivector_unchecked(self);
        let vector =
            if vec.x.abs() > tolerance || vec.y.abs() > tolerance || vec.z.abs() > tolerance {
                Some(vec)
            } else {
                None
            };

        // Extract bivector (grade 2)
        let biv = Bivec3::from_multivector_unchecked(self);
        let bivector =
            if biv.xy.abs() > tolerance || biv.xz.abs() > tolerance || biv.yz.abs() > tolerance {
                Some(biv)
            } else {
                None
            };

        // Extract trivector (grade 3)
        let tri = Trivec3::from_multivector_unchecked(self);
        let trivector = if tri.0.abs() > tolerance {
            Some(tri)
        } else {
            None
        };

        (scalar, vector, bivector, trivector)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;
    use std::f64::consts::FRAC_PI_4;

    use super::super::arbitrary::{UnitBivec3, UnitRotor3, UnitVec3};

    // ========================================================================
    // Vec3 tests
    // ========================================================================

    proptest! {
        #[test]
        fn vec3_roundtrip(v in any::<Vec3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = v.into();
            let back = Vec3::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(v, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec3_to_multivector_is_grade_1(v in any::<Vec3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = v.into();
            // Should be pure vector (grade 1) or zero
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(2).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(3).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn vec3_conversion_fails_with_scalar() {
        let mv: Multivector<f64, Euclidean3> = Multivector::scalar(1.0);
        assert!(Vec3::try_from(mv).is_err());
    }

    // ========================================================================
    // Bivec3 tests
    // ========================================================================

    proptest! {
        #[test]
        fn bivec3_roundtrip(b in any::<Bivec3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = b.into();
            let back = Bivec3::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(b, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivec3_to_multivector_is_grade_2(b in any::<Bivec3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = b.into();
            // Should be pure bivector (grade 2) or zero
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(3).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn bivec3_conversion_fails_with_vector() {
        let mv: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 0.0, 0.0]);
        assert!(Bivec3::try_from(mv).is_err());
    }

    // ========================================================================
    // Trivec3 tests
    // ========================================================================

    proptest! {
        #[test]
        fn trivec3_roundtrip(val in -100.0f64..100.0) {
            let t = Trivec3::new(val);
            let mv: Multivector<f64, Euclidean3> = t.into();
            let back = Trivec3::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(t, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn trivec3_to_multivector_is_grade_3(val in -100.0f64..100.0) {
            let t = Trivec3::new(val);
            let mv: Multivector<f64, Euclidean3> = t.into();
            // Should be pure trivector (grade 3) or zero
            prop_assert!(mv.grade_select(0).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(2).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn trivec3_conversion_fails_with_vector() {
        let mv: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 0.0, 0.0]);
        assert!(Trivec3::try_from(mv).is_err());
    }

    // ========================================================================
    // Rotor3 tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor3_roundtrip(r in any::<UnitRotor3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = (*r).into();
            let back = Rotor3::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(*r, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor3_to_multivector_is_even(r in any::<UnitRotor3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = (*r).into();
            // Should have no odd (vector, trivector) parts
            prop_assert!(mv.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(mv.grade_select(3).is_zero(ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn rotor3_conversion_fails_with_vector() {
        let mv: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 0.0, 0.0]);
        assert!(Rotor3::try_from(mv).is_err());
    }

    // ========================================================================
    // Even3 tests
    // ========================================================================

    proptest! {
        #[test]
        fn even3_roundtrip(
            s in -100.0f64..100.0,
            xy in -100.0f64..100.0,
            xz in -100.0f64..100.0,
            yz in -100.0f64..100.0,
        ) {
            let e = Even3::new(s, Bivec3::new(xy, xz, yz));
            let mv: Multivector<f64, Euclidean3> = e.into();
            let back = Even3::try_from(mv).unwrap();
            prop_assert!(abs_diff_eq!(e, back, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Operation consistency tests
    // ========================================================================

    proptest! {
        #[test]
        fn dot_consistency(
            a in any::<Vec3<f64>>(),
            b in any::<Vec3<f64>>(),
        ) {
            let spec_result = a.dot(b);

            let gen_a: Multivector<f64, Euclidean3> = a.into();
            let gen_b: Multivector<f64, Euclidean3> = b.into();
            let gen_result = gen_a.inner(&gen_b).scalar_part();

            prop_assert!(abs_diff_eq!(spec_result, gen_result, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn wedge_consistency(
            a in any::<Vec3<f64>>(),
            b in any::<Vec3<f64>>(),
        ) {
            let spec_result = a.wedge(b);

            let gen_a: Multivector<f64, Euclidean3> = a.into();
            let gen_b: Multivector<f64, Euclidean3> = b.into();
            let gen_result = gen_a.outer(&gen_b);
            let gen_as_bivec = Bivec3::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_bivec, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn geometric_product_consistency(
            a in any::<UnitVec3<f64>>(),
            b in any::<UnitVec3<f64>>(),
        ) {
            let spec_result = a.geometric(*b);

            let gen_a: Multivector<f64, Euclidean3> = (*a).into();
            let gen_b: Multivector<f64, Euclidean3> = (*b).into();
            let gen_result = &gen_a * &gen_b;
            let gen_as_even = Even3::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_even, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_rotation_consistency(
            r in any::<UnitRotor3<f64>>(),
            v in any::<UnitVec3<f64>>(),
        ) {
            // Specialized rotation
            let spec_result = r.rotate(*v);

            // Generic sandwich product: R̃ v R
            let gen_r: Multivector<f64, Euclidean3> = (*r).into();
            let gen_v: Multivector<f64, Euclidean3> = (*v).into();
            let gen_r_rev = gen_r.reverse();
            let gen_result = &(&gen_r_rev * &gen_v) * &gen_r;
            let gen_as_vec = Vec3::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_vec, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivec3_dual_consistency(
            b in any::<UnitBivec3<f64>>(),
        ) {
            // Specialized dual
            let spec_result = b.dual();

            // Generic dual
            let gen_b: Multivector<f64, Euclidean3> = (*b).into();
            let gen_result = gen_b.dual();
            let gen_as_vec = Vec3::try_from(gen_result).unwrap();

            prop_assert!(abs_diff_eq!(spec_result, gen_as_vec, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn rotor_from_angle_plane_consistency() {
        let plane = Bivec3::unit_xy();
        let r = Rotor3::from_angle_plane(FRAC_PI_4, plane);
        let mv: Multivector<f64, Euclidean3> = r.into();

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
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E13_IDX)),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E23_IDX)),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // ========================================================================
    // Tests starting from generic Multivector (TryFrom consistency)
    // ========================================================================

    // Import generic arbitrary types
    use crate::algebra::arbitrary::{NonZeroVectorE3, VectorE3};

    proptest! {
        /// Test that a generic VectorE3 converts to Vec3 correctly.
        #[test]
        fn generic_vectore3_to_vec3(v in any::<VectorE3>()) {
            // Get components before moving the value
            let x = v.get(Blade::from_index(E1_IDX));
            let y = v.get(Blade::from_index(E2_IDX));
            let z = v.get(Blade::from_index(E3_IDX));

            let vec3 = Vec3::try_from(v.into_inner()).expect("VectorE3 should convert to Vec3");

            // Verify components match
            prop_assert!(abs_diff_eq!(vec3.x, x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(vec3.y, y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(vec3.z, z, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic vector Multivector converts to Vec3 correctly.
        #[test]
        fn generic_vector_to_vec3(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let mv: Multivector<f64, Euclidean3> = Multivector::vector(&[x, y, z]);
            let v = Vec3::try_from(mv).expect("pure vector should convert");
            prop_assert!(abs_diff_eq!(v.x, x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(v.y, y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(v.z, z, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic bivector Multivector converts to Bivec3 correctly.
        #[test]
        fn generic_bivector_to_bivec3(xy in -100.0f64..100.0, xz in -100.0f64..100.0, yz in -100.0f64..100.0) {
            let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
            mv.set(Blade::from_index(E12_IDX), xy);
            mv.set(Blade::from_index(E13_IDX), xz);
            mv.set(Blade::from_index(E23_IDX), yz);
            let b = Bivec3::try_from(mv).expect("pure bivector should convert");
            prop_assert!(abs_diff_eq!(b.xy, xy, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(b.xz, xz, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(b.yz, yz, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic trivector Multivector converts to Trivec3 correctly.
        #[test]
        fn generic_trivector_to_trivec3(xyz in -100.0f64..100.0) {
            let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
            mv.set(Blade::from_index(E123_IDX), xyz);
            let t = Trivec3::try_from(mv).expect("pure trivector should convert");
            prop_assert!(abs_diff_eq!(t.value(), xyz, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that a generic even Multivector converts to Rotor3 correctly.
        #[test]
        fn generic_even_to_rotor3(
            s in -100.0f64..100.0,
            xy in -100.0f64..100.0,
            xz in -100.0f64..100.0,
            yz in -100.0f64..100.0,
        ) {
            let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
            mv.set(Blade::from_index(SCALAR_IDX), s);
            mv.set(Blade::from_index(E12_IDX), xy);
            mv.set(Blade::from_index(E13_IDX), xz);
            mv.set(Blade::from_index(E23_IDX), yz);
            let r = Rotor3::try_from(mv).expect("even element should convert");
            prop_assert!(abs_diff_eq!(r.s, s, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(r.b.xy, xy, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(r.b.xz, xz, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(r.b.yz, yz, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that full generic Multivector fails TryFrom for specialized types.
        #[test]
        fn generic_mixed_fails_vec3(mv in any::<Multivector<f64, Euclidean3>>()) {
            let has_scalar = mv.get(Blade::from_index(SCALAR_IDX)).abs() > ABS_DIFF_EQ_EPS;
            let has_bivec = mv.get(Blade::from_index(E12_IDX)).abs() > ABS_DIFF_EQ_EPS
                || mv.get(Blade::from_index(E13_IDX)).abs() > ABS_DIFF_EQ_EPS
                || mv.get(Blade::from_index(E23_IDX)).abs() > ABS_DIFF_EQ_EPS;
            let has_trivec = mv.get(Blade::from_index(E123_IDX)).abs() > ABS_DIFF_EQ_EPS;
            let has_vec = mv.get(Blade::from_index(E1_IDX)).abs() > ABS_DIFF_EQ_EPS
                || mv.get(Blade::from_index(E2_IDX)).abs() > ABS_DIFF_EQ_EPS
                || mv.get(Blade::from_index(E3_IDX)).abs() > ABS_DIFF_EQ_EPS;

            if has_vec && (has_scalar || has_bivec || has_trivec) {
                // Has vector and something else - should fail vec conversion
                prop_assert!(Vec3::try_from(mv).is_err());
            }
        }

        /// Test that vector wedge product (generic) converts to Bivec3.
        #[test]
        fn generic_wedge_to_bivec3(a in any::<VectorE3>(), b in any::<VectorE3>()) {
            let wedge = a.outer(&*b);

            // Should be a pure bivector, convertible to Bivec3
            let bivec = Bivec3::try_from(wedge).expect("wedge of vectors should be bivector");

            // Check values match component-wise wedge formula
            let ax = a.get(Blade::from_index(E1_IDX));
            let ay = a.get(Blade::from_index(E2_IDX));
            let az = a.get(Blade::from_index(E3_IDX));
            let bx = b.get(Blade::from_index(E1_IDX));
            let by = b.get(Blade::from_index(E2_IDX));
            let bz = b.get(Blade::from_index(E3_IDX));

            let expected_xy = ax * by - ay * bx;
            let expected_xz = ax * bz - az * bx;
            let expected_yz = ay * bz - az * by;

            prop_assert!(abs_diff_eq!(bivec.xy, expected_xy, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(bivec.xz, expected_xz, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(bivec.yz, expected_yz, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test that vector geometric product (generic) converts to Even3/Rotor3.
        #[test]
        fn generic_geometric_to_even3(a in any::<VectorE3>(), b in any::<VectorE3>()) {
            let product = &*a * &*b;

            // Should be scalar + bivector (even), convertible to Even3
            let even = Even3::try_from(product).expect("geometric product should be even");

            // Check scalar = dot product
            let ax = a.get(Blade::from_index(E1_IDX));
            let ay = a.get(Blade::from_index(E2_IDX));
            let az = a.get(Blade::from_index(E3_IDX));
            let bx = b.get(Blade::from_index(E1_IDX));
            let by = b.get(Blade::from_index(E2_IDX));
            let bz = b.get(Blade::from_index(E3_IDX));

            let expected_s = ax * bx + ay * by + az * bz;
            prop_assert!(abs_diff_eq!(even.s, expected_s, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Test sandwich product via generic Multivector produces result with vector grade.
        #[test]
        fn generic_sandwich_produces_vector_grade(
            n in any::<NonZeroVectorE3>(),
            v in any::<VectorE3>(),
        ) {
            // Reflection: n v n (with normalized n)
            let n_normalized = n.normalize().unwrap();
            let reflected = n_normalized.sandwich(&*v);

            // The result should be predominantly grade-1
            // Due to floating-point arithmetic, tiny components in other grades may appear
            // So we verify the grade-1 extraction rather than strict conversion
            let grade1 = reflected.grade_select(1);

            // The non-grade-1 parts should be negligible (within numerical tolerance)
            let grade0 = reflected.grade_select(0);
            let grade2 = reflected.grade_select(2);
            let grade3 = reflected.grade_select(3);

            prop_assert!(grade0.is_zero(1e-9), "grade 0 should be negligible");
            prop_assert!(grade2.is_zero(1e-9), "grade 2 should be negligible");
            prop_assert!(grade3.is_zero(1e-9), "grade 3 should be negligible");

            // The grade-1 part should be extractable to Vec3
            let result = Vec3::try_from(grade1).expect("grade-1 part should be vector");

            // Verify finite values
            prop_assert!(result.x.is_finite());
            prop_assert!(result.y.is_finite());
            prop_assert!(result.z.is_finite());
        }

        /// Test dual of generic bivector produces valid vector.
        #[test]
        fn generic_bivector_dual_to_vec3(
            xy in -10.0f64..10.0,
            xz in -10.0f64..10.0,
            yz in -10.0f64..10.0,
        ) {
            let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
            mv.set(Blade::from_index(E12_IDX), xy);
            mv.set(Blade::from_index(E13_IDX), xz);
            mv.set(Blade::from_index(E23_IDX), yz);

            let dual = mv.dual();

            // Dual of bivector should be a vector
            let vec = Vec3::try_from(dual).expect("dual of bivector should be vector");

            // The specialized Bivec3::dual formula: (yz, -xz, xy)
            let bivec = Bivec3::new(xy, xz, yz);
            let spec_dual = bivec.dual();

            prop_assert!(abs_diff_eq!(vec, spec_dual, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // to_specialized tests
    // ========================================================================

    proptest! {
        /// Test that to_specialized correctly identifies pure vectors.
        #[test]
        fn to_specialized_pure_vector(v in any::<Vec3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = v.into();
            let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            prop_assert!(scalar.is_none());
            prop_assert!(bivector.is_none());
            prop_assert!(trivector.is_none());

            if v.x.abs().max(v.y.abs()).max(v.z.abs()) > ABS_DIFF_EQ_EPS {
                prop_assert!(vector.is_some());
                let vec = vector.unwrap();
                prop_assert!(abs_diff_eq!(vec, v, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized correctly identifies pure bivectors.
        #[test]
        fn to_specialized_pure_bivector(b in any::<Bivec3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = b.into();
            let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            prop_assert!(scalar.is_none());
            prop_assert!(vector.is_none());
            prop_assert!(trivector.is_none());

            if b.xy.abs().max(b.xz.abs()).max(b.yz.abs()) > ABS_DIFF_EQ_EPS {
                prop_assert!(bivector.is_some());
                let biv = bivector.unwrap();
                prop_assert!(abs_diff_eq!(biv, b, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized correctly identifies pure trivectors.
        #[test]
        fn to_specialized_pure_trivector(xyz in -100.0f64..100.0) {
            let t = Trivec3::new(xyz);
            let mv: Multivector<f64, Euclidean3> = t.into();
            let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            prop_assert!(scalar.is_none());
            prop_assert!(vector.is_none());
            prop_assert!(bivector.is_none());

            if xyz.abs() > ABS_DIFF_EQ_EPS {
                prop_assert!(trivector.is_some());
                let tri = trivector.unwrap();
                prop_assert!(abs_diff_eq!(tri, t, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized correctly identifies pure scalars.
        #[test]
        fn to_specialized_pure_scalar(s in -100.0f64..100.0) {
            let mv: Multivector<f64, Euclidean3> = Multivector::scalar(s);
            let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            prop_assert!(vector.is_none());
            prop_assert!(bivector.is_none());
            prop_assert!(trivector.is_none());

            if s.abs() > ABS_DIFF_EQ_EPS {
                prop_assert!(scalar.is_some());
                prop_assert!(abs_diff_eq!(scalar.unwrap(), s, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized correctly identifies rotors (even multivectors).
        #[test]
        fn to_specialized_rotor(r in any::<UnitRotor3<f64>>()) {
            let mv: Multivector<f64, Euclidean3> = (*r).into();
            let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            // Rotor has scalar and bivector, no vector or trivector
            prop_assert!(vector.is_none());
            prop_assert!(trivector.is_none());

            // Check scalar part
            if r.s.abs() > ABS_DIFF_EQ_EPS {
                prop_assert!(scalar.is_some());
                prop_assert!(abs_diff_eq!(scalar.unwrap(), r.s, epsilon = ABS_DIFF_EQ_EPS));
            }

            // Check bivector part
            if r.b.xy.abs().max(r.b.xz.abs()).max(r.b.yz.abs()) > ABS_DIFF_EQ_EPS {
                prop_assert!(bivector.is_some());
                prop_assert!(abs_diff_eq!(bivector.unwrap(), r.b, epsilon = ABS_DIFF_EQ_EPS));
            }
        }

        /// Test that to_specialized handles arbitrary multivectors correctly.
        #[test]
        fn to_specialized_arbitrary(mv in any::<Multivector<f64, Euclidean3>>()) {
            let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

            // Reconstruct from components and verify
            let mut reconstructed: Multivector<f64, Euclidean3> = Multivector::zero();

            if let Some(s) = scalar {
                reconstructed = &reconstructed + &Multivector::scalar(s);
            }
            if let Some(v) = vector {
                reconstructed = &reconstructed + &Multivector::from(v);
            }
            if let Some(b) = bivector {
                reconstructed = &reconstructed + &Multivector::from(b);
            }
            if let Some(t) = trivector {
                reconstructed = &reconstructed + &Multivector::from(t);
            }

            // Reconstructed should be approximately equal to original (within tolerance)
            // Note: small values below tolerance are zeroed out, so we check equivalence
            for i in 0..8 {
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
        let mv: Multivector<f64, Euclidean3> = Multivector::zero();
        let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

        assert!(scalar.is_none());
        assert!(vector.is_none());
        assert!(bivector.is_none());
        assert!(trivector.is_none());
    }

    #[test]
    fn to_specialized_full_multivector() {
        // Create a multivector with all grades non-zero
        let mut mv: Multivector<f64, Euclidean3> = Multivector::scalar(2.0);
        mv = &mv + &Multivector::from(Vec3::new(1.0, 2.0, 3.0));
        mv = &mv + &Multivector::from(Bivec3::new(4.0, 5.0, 6.0));
        mv = &mv + &Multivector::from(Trivec3::new(7.0));

        let (scalar, vector, bivector, trivector) = mv.to_specialized(ABS_DIFF_EQ_EPS);

        assert_eq!(scalar, Some(2.0));

        assert!(vector.is_some());
        let v = vector.unwrap();
        assert!(abs_diff_eq!(v.x, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(v.y, 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(v.z, 3.0, epsilon = ABS_DIFF_EQ_EPS));

        assert!(bivector.is_some());
        let b = bivector.unwrap();
        assert!(abs_diff_eq!(b.xy, 4.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(b.xz, 5.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(b.yz, 6.0, epsilon = ABS_DIFF_EQ_EPS));

        assert!(trivector.is_some());
        assert!(abs_diff_eq!(
            trivector.unwrap().value(),
            7.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }
}
