//! Conversions between specialized 2D PGA types and other representations.
//!
//! This module provides:
//!
//! - Bidirectional conversions between specialized types ([`Point`], [`Line`],
//!   [`Motor`]) and generic [`Multivector<T, Projective2>`]
//! - Conversions from Euclidean [`Vector`](crate::specialized::euclidean::dim2::Vector)
//!   to projective [`Point`] (embedding with w=1)
//! - Extraction from projective [`Point`] back to Euclidean
//!   [`Vector`](crate::specialized::euclidean::dim2::Vector)

use crate::algebra::Multivector;
use crate::basis::Blade;
use crate::scalar::Float;
use crate::signature::Projective2;
use crate::specialized::euclidean::dim2::Vector as EuclideanVector;

use super::types::{Line, Motor, Point};

// ============================================================================
// Blade index constants for Projective2
// ============================================================================
//
// These constants map blade indices (bitmask representation) to their
// corresponding basis elements in the Projective2 signature.

/// Grade-0 scalar blade index.
const SCALAR: usize = 0b000;

/// Grade-1 e₁ basis vector index.
const E1: usize = 0b001;
/// Grade-1 e₂ basis vector index.
const E2: usize = 0b010;
/// Grade-1 e₀ (null) basis vector index.
const E0: usize = 0b100;

/// Grade-2 e₁₂ bivector index.
const E12: usize = 0b011;
/// Grade-2 e₂₀ = e₂ ∧ e₀ bivector index.
const E20: usize = 0b110;
/// Grade-2 e₀₁ = e₀ ∧ e₁ bivector index.
const E01: usize = 0b101;

/// Grade-3 e₀₁₂ pseudoscalar index.
const E012: usize = 0b111;

// ============================================================================
// Point conversions
// ============================================================================

impl<T: Float> From<Point<T>> for Multivector<T, Projective2> {
    /// Converts a [`Point`] to a [`Multivector`].
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Projective2;
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let point = Point::new(1.0, 2.0);
    /// let mv: Multivector<f64, Projective2> = point.into();
    /// ```
    fn from(p: Point<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E1), p.e1);
        mv.set(Blade::from_index(E2), p.e2);
        mv.set(Blade::from_index(E0), p.e0);
        mv
    }
}

impl<T: Float> TryFrom<Multivector<T, Projective2>> for Point<T> {
    type Error = ConversionError;

    /// Tries to convert a [`Multivector`] to a [`Point`].
    ///
    /// Returns an error if the multivector has non-zero components outside grade 1.
    fn try_from(mv: Multivector<T, Projective2>) -> Result<Self, Self::Error> {
        // Check that only grade-1 components are non-zero
        let eps = T::epsilon();
        if mv.get(Blade::from_index(SCALAR)).abs() > eps {
            return Err(ConversionError::NonZeroScalar);
        }
        if mv.get(Blade::from_index(E12)).abs() > eps
            || mv.get(Blade::from_index(E20)).abs() > eps
            || mv.get(Blade::from_index(E01)).abs() > eps
        {
            return Err(ConversionError::NonZeroBivector);
        }
        if mv.get(Blade::from_index(E012)).abs() > eps {
            return Err(ConversionError::NonZeroPseudoscalar);
        }

        Ok(Point {
            e1: mv.get(Blade::from_index(E1)),
            e2: mv.get(Blade::from_index(E2)),
            e0: mv.get(Blade::from_index(E0)),
        })
    }
}

// ============================================================================
// Line conversions
// ============================================================================

impl<T: Float> From<Line<T>> for Multivector<T, Projective2> {
    /// Converts a [`Line`] to a [`Multivector`].
    fn from(l: Line<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E12), l.e12);
        mv.set(Blade::from_index(E20), l.e20);
        mv.set(Blade::from_index(E01), l.e01);
        mv
    }
}

impl<T: Float> TryFrom<Multivector<T, Projective2>> for Line<T> {
    type Error = ConversionError;

    /// Tries to convert a [`Multivector`] to a [`Line`].
    ///
    /// Returns an error if the multivector has non-zero components outside grade 2.
    fn try_from(mv: Multivector<T, Projective2>) -> Result<Self, Self::Error> {
        let eps = T::epsilon();
        if mv.get(Blade::from_index(SCALAR)).abs() > eps {
            return Err(ConversionError::NonZeroScalar);
        }
        if mv.get(Blade::from_index(E1)).abs() > eps
            || mv.get(Blade::from_index(E2)).abs() > eps
            || mv.get(Blade::from_index(E0)).abs() > eps
        {
            return Err(ConversionError::NonZeroVector);
        }
        if mv.get(Blade::from_index(E012)).abs() > eps {
            return Err(ConversionError::NonZeroPseudoscalar);
        }

        Ok(Line {
            e12: mv.get(Blade::from_index(E12)),
            e20: mv.get(Blade::from_index(E20)),
            e01: mv.get(Blade::from_index(E01)),
        })
    }
}

// ============================================================================
// Motor conversions
// ============================================================================

impl<T: Float> From<Motor<T>> for Multivector<T, Projective2> {
    /// Converts a [`Motor`] to a [`Multivector`].
    fn from(m: Motor<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(SCALAR), m.s);
        mv.set(Blade::from_index(E12), m.e12);
        mv.set(Blade::from_index(E20), m.e20);
        mv.set(Blade::from_index(E01), m.e01);
        mv
    }
}

impl<T: Float> TryFrom<Multivector<T, Projective2>> for Motor<T> {
    type Error = ConversionError;

    /// Tries to convert a [`Multivector`] to a [`Motor`].
    ///
    /// Returns an error if the multivector has non-zero odd-grade components.
    fn try_from(mv: Multivector<T, Projective2>) -> Result<Self, Self::Error> {
        let eps = T::epsilon();
        // Motors are even-grade: scalar + bivector
        // Check that vector and pseudoscalar parts are zero
        if mv.get(Blade::from_index(E1)).abs() > eps
            || mv.get(Blade::from_index(E2)).abs() > eps
            || mv.get(Blade::from_index(E0)).abs() > eps
        {
            return Err(ConversionError::NonZeroVector);
        }
        if mv.get(Blade::from_index(E012)).abs() > eps {
            return Err(ConversionError::NonZeroPseudoscalar);
        }

        Ok(Motor {
            s: mv.get(Blade::from_index(SCALAR)),
            e12: mv.get(Blade::from_index(E12)),
            e20: mv.get(Blade::from_index(E20)),
            e01: mv.get(Blade::from_index(E01)),
        })
    }
}

// ============================================================================
// Euclidean conversions
// ============================================================================

impl<T: Float> From<EuclideanVector<T>> for Point<T> {
    /// Embeds a Euclidean 2D vector as a projective point with weight w=1.
    ///
    /// This creates a finite point at the Cartesian coordinates (x, y).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector as EucVec;
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let v = EucVec::new(3.0, 4.0);
    /// let p: Point<f64> = v.into();
    ///
    /// assert_eq!(p.x(), 3.0);
    /// assert_eq!(p.y(), 4.0);
    /// assert_eq!(p.w(), 1.0);
    /// ```
    #[inline]
    fn from(v: EuclideanVector<T>) -> Self {
        Point::new(v.x, v.y)
    }
}

impl<T: Float> TryFrom<Point<T>> for EuclideanVector<T> {
    type Error = ConversionError;

    /// Extracts a Euclidean 2D vector from a projective point.
    ///
    /// Returns an error if the point is ideal (w ≈ 0).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector as EucVec;
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let p = Point::new(3.0, 4.0);
    /// let v: EucVec<f64> = p.try_into().unwrap();
    ///
    /// assert_eq!(v.x, 3.0);
    /// assert_eq!(v.y, 4.0);
    /// ```
    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0.abs() < T::epsilon() {
            return Err(ConversionError::IdealPoint);
        }
        Ok(EuclideanVector::new(p.e1 / p.e0, p.e2 / p.e0))
    }
}

impl<T: Float> Motor<T> {
    /// Transforms a Euclidean 2D vector using this motor.
    ///
    /// This is a convenience method that:
    /// 1. Embeds the vector as a projective point (w=1)
    /// 2. Applies the motor transformation
    /// 3. Extracts the Euclidean coordinates
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector as EucVec;
    /// use clifford::specialized::projective::dim2::Motor;
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let v = EucVec::new(1.0, 0.0);
    ///
    /// // Rotate 90° then translate
    /// let rotation = Motor::from_rotation(FRAC_PI_2);
    /// let translation = Motor::from_translation(1.0, 2.0);
    /// let motor = rotation.compose(&translation);
    ///
    /// let result = motor.transform_euclidean(&v);
    /// assert!(abs_diff_eq!(result.x, 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.y, 3.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn transform_euclidean(&self, v: &EuclideanVector<T>) -> EuclideanVector<T> {
        let p = Point::new(v.x, v.y);
        let transformed = self.transform_point(&p);
        // Safe because motor preserves weight of finite points
        EuclideanVector::new(
            transformed.e1 / transformed.e0,
            transformed.e2 / transformed.e0,
        )
    }
}

// ============================================================================
// ConversionError
// ============================================================================

/// Error type for conversions from [`Multivector`] or specialized types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)]
pub enum ConversionError {
    /// Multivector has non-zero scalar component where not expected.
    NonZeroScalar,
    /// Multivector has non-zero vector (grade 1) components where not expected.
    NonZeroVector,
    /// Multivector has non-zero bivector (grade 2) components where not expected.
    NonZeroBivector,
    /// Multivector has non-zero pseudoscalar component where not expected.
    NonZeroPseudoscalar,
    /// Projective point is ideal (at infinity) and cannot be converted to Euclidean.
    IdealPoint,
}

impl core::fmt::Display for ConversionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonZeroScalar => write!(f, "multivector has non-zero scalar component"),
            Self::NonZeroVector => write!(f, "multivector has non-zero vector components"),
            Self::NonZeroBivector => write!(f, "multivector has non-zero bivector components"),
            Self::NonZeroPseudoscalar => {
                write!(f, "multivector has non-zero pseudoscalar component")
            }
            Self::IdealPoint => {
                write!(f, "projective point is ideal (at infinity)")
            }
        }
    }
}

impl std::error::Error for ConversionError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;

    #[test]
    fn point_roundtrip() {
        let point = Point::new(3.0, 4.0);
        let mv: Multivector<f64, Projective2> = point.into();
        let back = Point::try_from(mv).unwrap();

        assert!(abs_diff_eq!(point.e1, back.e1, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e2, back.e2, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e0, back.e0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_roundtrip() {
        let line = Line::new(1.0, 2.0, 3.0);
        let mv: Multivector<f64, Projective2> = line.into();
        let back = Line::try_from(mv).unwrap();

        assert!(abs_diff_eq!(line.e12, back.e12, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(line.e20, back.e20, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(line.e01, back.e01, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_roundtrip() {
        let motor = Motor::new(0.5, 0.1, 0.2, 0.3);
        let mv: Multivector<f64, Projective2> = motor.into();
        let back = Motor::try_from(mv).unwrap();

        assert!(abs_diff_eq!(motor.s, back.s, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e12, back.e12, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e20, back.e20, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e01, back.e01, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn conversion_error_for_wrong_grade() {
        // Try to convert a motor (even) to a point (grade 1)
        let motor = Motor::from_rotation(0.5);
        let mv: Multivector<f64, Projective2> = motor.into();
        let result = Point::try_from(mv);
        assert!(result.is_err());
    }
}
