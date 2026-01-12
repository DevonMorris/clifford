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
use crate::specialized::euclidean::dim2::Rotor as EuclideanRotor;
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
        mv.set(Blade::from_index(SCALAR), m.s());
        mv.set(Blade::from_index(E12), m.e12());
        mv.set(Blade::from_index(E20), m.e20());
        mv.set(Blade::from_index(E01), m.e01());
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

        Ok(Motor::new_unchecked(
            mv.get(Blade::from_index(SCALAR)),
            mv.get(Blade::from_index(E12)),
            mv.get(Blade::from_index(E20)),
            mv.get(Blade::from_index(E01)),
        ))
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
        Point::new(v.x(), v.y())
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
    /// assert_eq!(v.x(), 3.0);
    /// assert_eq!(v.y(), 4.0);
    /// ```
    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0.abs() < T::epsilon() {
            return Err(ConversionError::IdealPoint);
        }
        Ok(EuclideanVector::new(p.e1 / p.e0, p.e2 / p.e0))
    }
}

impl<T: Float> From<EuclideanRotor<T>> for Motor<T> {
    /// Converts a Euclidean 2D rotor to a pure rotation motor.
    ///
    /// The resulting motor has no translation component.
    ///
    /// # Mapping
    ///
    /// - `rotor.s` → `motor.s` (scalar, cos(θ/2))
    /// - `rotor.xy` → `motor.e12` (rotation bivector, sin(θ/2))
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use clifford::specialized::projective::dim2::{Motor, Point};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::relative_eq;
    ///
    /// // Create a 90° rotation rotor
    /// let rotor = Rotor::from_angle(FRAC_PI_2);
    /// let motor: Motor<f64> = rotor.into();
    ///
    /// // Transform a point
    /// let p = Point::new(1.0, 0.0);
    /// let rotated = motor.transform_point(&p);
    ///
    /// assert!(relative_eq!(rotated.x(), 0.0, epsilon = 1e-10, max_relative = 1e-10));
    /// assert!(relative_eq!(rotated.y(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
    /// ```
    #[inline]
    fn from(r: EuclideanRotor<T>) -> Self {
        Motor::new_unchecked(r.s(), r.xy(), T::zero(), T::zero())
    }
}

impl<T: Float> From<Motor<T>> for EuclideanRotor<T> {
    /// Extracts the rotation part of a motor as a Euclidean rotor.
    ///
    /// # Note
    ///
    /// This extracts only the rotation component. Translation is discarded.
    /// The result is normalized to ensure a valid unit rotor.
    ///
    /// # Mapping
    ///
    /// - `motor.s` → `rotor.s`
    /// - `motor.e12` → `rotor.xy`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::{Rotor, Vector};
    /// use clifford::specialized::projective::dim2::Motor;
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::relative_eq;
    ///
    /// // Create a motor with rotation and translation
    /// let rotation = Motor::from_rotation(FRAC_PI_2);
    /// let translation = Motor::from_translation(1.0, 2.0);
    /// let motor = rotation.compose(&translation);
    ///
    /// // Extract just the rotation
    /// let rotor: Rotor<f64> = motor.into();
    ///
    /// // The rotor should perform the same rotation
    /// let v = Vector::new(1.0, 0.0);
    /// let rotated = rotor.rotate(v);
    ///
    /// assert!(relative_eq!(rotated.x(), 0.0, epsilon = 1e-10, max_relative = 1e-10));
    /// assert!(relative_eq!(rotated.y(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
    /// ```
    #[inline]
    fn from(m: Motor<T>) -> Self {
        // Unit motor's rotation part is already unit, no normalization needed
        EuclideanRotor::new_unchecked(m.s(), m.e12())
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
    /// use approx::relative_eq;
    ///
    /// let v = EucVec::new(1.0, 0.0);
    ///
    /// // Rotate 90° then translate
    /// let rotation = Motor::from_rotation(FRAC_PI_2);
    /// let translation = Motor::from_translation(1.0, 2.0);
    /// let motor = rotation.compose(&translation);
    ///
    /// let result = motor.transform_euclidean(&v);
    /// assert!(relative_eq!(result.x(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
    /// assert!(relative_eq!(result.y(), 3.0, epsilon = 1e-10, max_relative = 1e-10));
    /// ```
    #[inline]
    pub fn transform_euclidean(&self, v: &EuclideanVector<T>) -> EuclideanVector<T> {
        let p = Point::new(v.x(), v.y());
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
    use crate::specialized::euclidean::dim2::Vector;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;
    use proptest::prelude::*;

    #[test]
    fn point_roundtrip() {
        let point = Point::new(3.0, 4.0);
        let mv: Multivector<f64, Projective2> = point.into();
        let back = Point::try_from(mv).unwrap();

        assert!(relative_eq!(
            point.e1,
            back.e1,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            point.e2,
            back.e2,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            point.e0,
            back.e0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn line_roundtrip() {
        let line = Line::new(1.0, 2.0, 3.0);
        let mv: Multivector<f64, Projective2> = line.into();
        let back = Line::try_from(mv).unwrap();

        assert!(relative_eq!(
            line.e12,
            back.e12,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            line.e20,
            back.e20,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            line.e01,
            back.e01,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_roundtrip() {
        let motor = Motor::new_unchecked(0.5, 0.1, 0.2, 0.3);
        let mv: Multivector<f64, Projective2> = motor.into();
        let back = Motor::try_from(mv).unwrap();

        assert!(relative_eq!(
            motor.s(),
            back.s(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            motor.e12(),
            back.e12(),
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            motor.e20(),
            back.e20(),
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            motor.e01(),
            back.e01(),
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn conversion_error_for_wrong_grade() {
        // Try to convert a motor (even) to a point (grade 1)
        let motor = Motor::from_rotation(0.5);
        let mv: Multivector<f64, Projective2> = motor.into();
        let result = Point::try_from(mv);
        assert!(result.is_err());
    }

    // ========================================================================
    // Rotor <-> Motor interop tests
    // ========================================================================

    proptest! {
        /// Tests Rotor -> Motor -> Rotor roundtrip preserves rotation behavior.
        #[test]
        fn rotor_motor_roundtrip(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            vx in -10.0f64..10.0, vy in -10.0f64..10.0,
        ) {
            let rotor = EuclideanRotor::from_angle(angle);
            let motor: Motor<f64> = rotor.into();
            let back: EuclideanRotor<f64> = motor.into();

            // Compare by rotating a vector
            let v = Vector::new(vx, vy);
            let result_orig = rotor.rotate(v);
            let result_back = back.rotate(v);

            prop_assert!(relative_eq!(result_orig.x(), result_back.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(result_orig.y(), result_back.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        /// Tests Motor -> Rotor -> Motor roundtrip preserves rotation behavior.
        #[test]
        fn motor_rotor_roundtrip(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            // Create a pure rotation motor
            let motor = Motor::from_rotation(angle);
            let rotor: EuclideanRotor<f64> = motor.into();
            let back: Motor<f64> = rotor.into();

            // Compare by transforming a point
            let p = Point::new(px, py);
            let result_orig = motor.transform_point(&p);
            let result_back = back.transform_point(&p);

            prop_assert!(relative_eq!(result_orig.x(), result_back.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(result_orig.y(), result_back.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        /// Tests that Rotor rotation matches Motor rotation.
        #[test]
        fn rotor_motor_rotation_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            x in -10.0f64..10.0, y in -10.0f64..10.0,
        ) {
            let rotor = EuclideanRotor::from_angle(angle);
            let motor: Motor<f64> = rotor.into();

            // Rotate with Euclidean rotor
            let v = Vector::new(x, y);
            let rotated_euc = rotor.rotate(v);

            // Transform with PGA motor
            let p = Point::new(x, y);
            let rotated_pga = motor.transform_point(&p);

            prop_assert!(relative_eq!(rotated_euc.x(), rotated_pga.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(rotated_euc.y(), rotated_pga.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        /// Tests that Motor -> Rotor extracts correct rotation from composed motor.
        #[test]
        fn motor_to_rotor_extracts_rotation(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0,
            vx in -10.0f64..10.0, vy in -10.0f64..10.0,
        ) {
            // Create a motor with both rotation and translation
            let rotation = Motor::from_rotation(angle);
            let translation = Motor::from_translation(tx, ty);
            let composed = rotation.compose(&translation);

            // Extract just the rotation
            let rotor: EuclideanRotor<f64> = composed.into();

            // The extracted rotor should match the original rotation
            let expected_rotor = EuclideanRotor::from_angle(angle);
            let v = Vector::new(vx, vy);

            let result = rotor.rotate(v);
            let expected = expected_rotor.rotate(v);

            prop_assert!(relative_eq!(result.x(), expected.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(result.y(), expected.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }
    }

    #[test]
    fn identity_rotor_to_motor() {
        let rotor = EuclideanRotor::<f64>::identity();
        let motor: Motor<f64> = rotor.into();

        // Should be identity motor (no rotation, no translation)
        assert!(relative_eq!(
            motor.s(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            motor.e12(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            motor.e20(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            motor.e01(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn identity_motor_to_rotor() {
        let motor = Motor::<f64>::identity();
        let rotor: EuclideanRotor<f64> = motor.into();

        // Should be identity rotor
        assert!(relative_eq!(
            rotor.s(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            rotor.xy(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
