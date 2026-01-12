//! Conversions between specialized 3D PGA types and other representations.
//!
//! This module provides:
//!
//! - Bidirectional conversions between specialized types ([`Point`], [`Motor`])
//!   and generic [`Multivector<T, Projective3>`]
//! - Conversions from Euclidean [`Vector`](crate::specialized::euclidean::dim3::Vector)
//!   to projective [`Point`] (embedding with w=1)
//! - Extraction from projective [`Point`] back to Euclidean
//!   [`Vector`](crate::specialized::euclidean::dim3::Vector)

use crate::algebra::Multivector;
use crate::basis::Blade;
use crate::scalar::Float;
use crate::signature::Projective3;
use crate::specialized::euclidean::dim3::Rotor as EuclideanRotor;
use crate::specialized::euclidean::dim3::Vector as EuclideanVector;

use super::types::{Motor, Point};

// ============================================================================
// Blade index constants for Projective3
// ============================================================================
//
// These constants map blade indices (bitmask representation) to their
// corresponding basis elements in the Projective3 signature.
// Bit ordering: bit 0 = e₁, bit 1 = e₂, bit 2 = e₃, bit 3 = e₀

/// Grade-0 scalar blade index.
const SCALAR: usize = 0b0000;

/// Grade-1 e₁ basis vector index.
const E1: usize = 0b0001;
/// Grade-1 e₂ basis vector index.
const E2: usize = 0b0010;
/// Grade-1 e₃ basis vector index.
const E3: usize = 0b0100;
/// Grade-1 e₀ (null) basis vector index.
const E0: usize = 0b1000;

/// Grade-2 e₁₂ bivector index.
const E12: usize = 0b0011;
/// Grade-2 e₂₃ bivector index.
const E23: usize = 0b0110;
/// Grade-2 e₃₁ bivector index.
const E31: usize = 0b0101;
/// Grade-2 e₀₁ bivector index.
const E01: usize = 0b1001;
/// Grade-2 e₀₂ bivector index.
const E02: usize = 0b1010;
/// Grade-2 e₀₃ bivector index.
const E03: usize = 0b1100;

/// Grade-3 e₁₂₃ trivector index.
const E123: usize = 0b0111;
/// Grade-3 e₀₁₂ trivector index.
const E012: usize = 0b1011;
/// Grade-3 e₀₁₃ trivector index.
const E013: usize = 0b1101;
/// Grade-3 e₀₂₃ trivector index.
const E023: usize = 0b1110;

/// Grade-4 e₀₁₂₃ pseudoscalar index.
const E0123: usize = 0b1111;

// ============================================================================
// Point conversions
// ============================================================================

impl<T: Float> From<Point<T>> for Multivector<T, Projective3> {
    /// Converts a [`Point`] to a [`Multivector`].
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Projective3;
    /// use clifford::specialized::projective::dim3::Point;
    ///
    /// let point = Point::new(1.0, 2.0, 3.0);
    /// let mv: Multivector<f64, Projective3> = point.into();
    /// ```
    fn from(p: Point<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(E1), p.e1);
        mv.set(Blade::from_index(E2), p.e2);
        mv.set(Blade::from_index(E3), p.e3);
        mv.set(Blade::from_index(E0), p.e0);
        mv
    }
}

impl<T: Float> TryFrom<Multivector<T, Projective3>> for Point<T> {
    type Error = ConversionError;

    /// Tries to convert a [`Multivector`] to a [`Point`].
    ///
    /// Returns an error if the multivector has non-zero components outside grade 1.
    fn try_from(mv: Multivector<T, Projective3>) -> Result<Self, Self::Error> {
        // Check that only grade-1 components are non-zero
        let eps = T::epsilon();

        if mv.get(Blade::from_index(SCALAR)).abs() > eps {
            return Err(ConversionError::NonZeroScalar);
        }

        // Check bivectors
        if mv.get(Blade::from_index(E12)).abs() > eps
            || mv.get(Blade::from_index(E23)).abs() > eps
            || mv.get(Blade::from_index(E31)).abs() > eps
            || mv.get(Blade::from_index(E01)).abs() > eps
            || mv.get(Blade::from_index(E02)).abs() > eps
            || mv.get(Blade::from_index(E03)).abs() > eps
        {
            return Err(ConversionError::NonZeroBivector);
        }

        // Check trivectors
        if mv.get(Blade::from_index(E123)).abs() > eps
            || mv.get(Blade::from_index(E012)).abs() > eps
            || mv.get(Blade::from_index(E013)).abs() > eps
            || mv.get(Blade::from_index(E023)).abs() > eps
        {
            return Err(ConversionError::NonZeroTrivector);
        }

        // Check pseudoscalar
        if mv.get(Blade::from_index(E0123)).abs() > eps {
            return Err(ConversionError::NonZeroPseudoscalar);
        }

        Ok(Point {
            e1: mv.get(Blade::from_index(E1)),
            e2: mv.get(Blade::from_index(E2)),
            e3: mv.get(Blade::from_index(E3)),
            e0: mv.get(Blade::from_index(E0)),
        })
    }
}

// ============================================================================
// Motor conversions
// ============================================================================

impl<T: Float> From<Motor<T>> for Multivector<T, Projective3> {
    /// Converts a [`Motor`] to a [`Multivector`].
    ///
    /// # Sign Convention
    ///
    /// The Motor struct uses the convention where `e01`, `e02`, `e03` represent
    /// e₀∧e₁, e₀∧e₂, e₀∧e₃ respectively. However, blade indices use ascending
    /// order: blade E01 (0b1001) represents e₁∧e₀ = -e₀∧e₁.
    ///
    /// Therefore, translation bivectors are negated during conversion:
    /// - Motor.e01 (representing e₀∧e₁) → -coefficient at blade index E01 (e₁∧e₀)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Projective3;
    /// use clifford::specialized::projective::dim3::Motor;
    ///
    /// let motor = Motor::from_translation(1.0, 2.0, 3.0);
    /// let mv: Multivector<f64, Projective3> = motor.into();
    /// ```
    fn from(m: Motor<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_index(SCALAR), m.s());
        mv.set(Blade::from_index(E23), m.e23());
        mv.set(Blade::from_index(E31), m.e31());
        mv.set(Blade::from_index(E12), m.e12());
        // Translation bivectors: e₀ᵢ = -eᵢ₀, so negate
        mv.set(Blade::from_index(E01), -m.e01());
        mv.set(Blade::from_index(E02), -m.e02());
        mv.set(Blade::from_index(E03), -m.e03());
        mv.set(Blade::from_index(E0123), m.e0123());
        mv
    }
}

impl<T: Float> TryFrom<Multivector<T, Projective3>> for Motor<T> {
    type Error = ConversionError;

    /// Tries to convert a [`Multivector`] to a [`Motor`].
    ///
    /// Returns an error if the multivector has non-zero odd-grade components
    /// (vectors or trivectors).
    ///
    /// # Sign Convention
    ///
    /// Translation bivectors are negated during extraction (inverse of From).
    /// See [`From<Motor<T>> for Multivector`] for details.
    fn try_from(mv: Multivector<T, Projective3>) -> Result<Self, Self::Error> {
        let eps = T::epsilon();

        // Motors are even-grade: scalar + bivector + pseudoscalar
        // Check that vectors are zero
        if mv.get(Blade::from_index(E1)).abs() > eps
            || mv.get(Blade::from_index(E2)).abs() > eps
            || mv.get(Blade::from_index(E3)).abs() > eps
            || mv.get(Blade::from_index(E0)).abs() > eps
        {
            return Err(ConversionError::NonZeroVector);
        }

        // Check that trivectors are zero
        if mv.get(Blade::from_index(E123)).abs() > eps
            || mv.get(Blade::from_index(E012)).abs() > eps
            || mv.get(Blade::from_index(E013)).abs() > eps
            || mv.get(Blade::from_index(E023)).abs() > eps
        {
            return Err(ConversionError::NonZeroTrivector);
        }

        Ok(Motor::new_unchecked(
            mv.get(Blade::from_index(SCALAR)),
            mv.get(Blade::from_index(E23)),
            mv.get(Blade::from_index(E31)),
            mv.get(Blade::from_index(E12)),
            // Translation bivectors: negate (inverse of From conversion)
            -mv.get(Blade::from_index(E01)),
            -mv.get(Blade::from_index(E02)),
            -mv.get(Blade::from_index(E03)),
            mv.get(Blade::from_index(E0123)),
        ))
    }
}

// ============================================================================
// Euclidean conversions
// ============================================================================

impl<T: Float> From<EuclideanVector<T>> for Point<T> {
    /// Embeds a Euclidean 3D vector as a projective point with weight w=1.
    ///
    /// This creates a finite point at the Cartesian coordinates (x, y, z).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector as EucVec;
    /// use clifford::specialized::projective::dim3::Point;
    ///
    /// let v = EucVec::new(3.0, 4.0, 5.0);
    /// let p: Point<f64> = v.into();
    ///
    /// assert_eq!(p.x(), 3.0);
    /// assert_eq!(p.y(), 4.0);
    /// assert_eq!(p.z(), 5.0);
    /// assert_eq!(p.w(), 1.0);
    /// ```
    #[inline]
    fn from(v: EuclideanVector<T>) -> Self {
        Point::new(v.x(), v.y(), v.z())
    }
}

impl<T: Float> TryFrom<Point<T>> for EuclideanVector<T> {
    type Error = ConversionError;

    /// Extracts a Euclidean 3D vector from a projective point.
    ///
    /// Returns an error if the point is ideal (w ≈ 0).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector as EucVec;
    /// use clifford::specialized::projective::dim3::Point;
    ///
    /// let p = Point::new(3.0, 4.0, 5.0);
    /// let v: EucVec<f64> = p.try_into().unwrap();
    ///
    /// assert_eq!(v.x(), 3.0);
    /// assert_eq!(v.y(), 4.0);
    /// assert_eq!(v.z(), 5.0);
    /// ```
    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0.abs() < T::epsilon() {
            return Err(ConversionError::IdealPoint);
        }
        Ok(EuclideanVector::new(p.e1 / p.e0, p.e2 / p.e0, p.e3 / p.e0))
    }
}

impl<T: Float> From<EuclideanRotor<T>> for Motor<T> {
    /// Converts a Euclidean 3D rotor to a pure rotation motor.
    ///
    /// The resulting motor has no translation component.
    ///
    /// # Mapping
    ///
    /// The Euclidean bivector components map to PGA motor bivectors:
    /// - `rotor.s` → `motor.s` (scalar)
    /// - `rotor.b.xy` → `motor.e12` (xy-plane rotation)
    /// - `rotor.b.xz` → `-motor.e31` (xz-plane, sign due to e₁₃ = -e₃₁)
    /// - `rotor.b.yz` → `motor.e23` (yz-plane rotation)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Bivector, Rotor};
    /// use clifford::specialized::projective::dim3::{Motor, Point};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// // Create a 90° rotation around Z axis
    /// let rotor = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let motor: Motor<f64> = rotor.into();
    ///
    /// // Transform a point
    /// let p = Point::new(1.0, 0.0, 0.0);
    /// let rotated = motor.transform_point(&p);
    ///
    /// assert!(abs_diff_eq!(rotated.x(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.z(), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    fn from(r: EuclideanRotor<T>) -> Self {
        // Mapping: xy -> e12, xz -> -e31 (sign!), yz -> e23
        Motor::new_unchecked(
            r.s(),
            r.b().yz(),  // e23: yz-plane rotation
            -r.b().xz(), // e31: note sign flip (e13 = -e31)
            r.b().xy(),  // e12: xy-plane rotation
            T::zero(),   // e01: no translation
            T::zero(),   // e02: no translation
            T::zero(),   // e03: no translation
            T::zero(),   // e0123: no pseudoscalar
        )
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
    /// - `motor.e12` → `rotor.b.xy`
    /// - `-motor.e31` → `rotor.b.xz` (sign due to e₃₁ = -e₁₃)
    /// - `motor.e23` → `rotor.b.yz`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Vector};
    /// use clifford::specialized::projective::dim3::Motor;
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// // Create a motor with rotation and translation
    /// let rotation = Motor::from_rotation_z(FRAC_PI_2);
    /// let translation = Motor::from_translation(1.0, 2.0, 3.0);
    /// let motor = rotation.compose(&translation);
    ///
    /// // Extract just the rotation
    /// let rotor: Rotor<f64> = motor.into();
    ///
    /// // The rotor should perform the same rotation
    /// let v = Vector::new(1.0, 0.0, 0.0);
    /// let rotated = rotor.rotate(v);
    ///
    /// assert!(abs_diff_eq!(rotated.x(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.z(), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    fn from(m: Motor<T>) -> Self {
        // Inverse mapping: e12 -> xy, -e31 -> xz, e23 -> yz
        EuclideanRotor::new(m.s(), m.e12(), -m.e31(), m.e23()).normalized()
    }
}

impl<T: Float> Motor<T> {
    /// Transforms a Euclidean 3D vector using this motor.
    ///
    /// This is a convenience method that:
    /// 1. Embeds the vector as a projective point (w=1)
    /// 2. Applies the motor transformation
    /// 3. Extracts the Euclidean coordinates
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector as EucVec;
    /// use clifford::specialized::projective::dim3::Motor;
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let v = EucVec::new(1.0, 0.0, 0.0);
    ///
    /// // Rotate 90° around Z then translate
    /// let rotation = Motor::from_rotation_z(FRAC_PI_2);
    /// let translation = Motor::from_translation(1.0, 2.0, 3.0);
    /// let motor = rotation.compose(&translation);
    ///
    /// let result = motor.transform_euclidean(&v);
    /// assert!(abs_diff_eq!(result.x(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.y(), 3.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.z(), 3.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn transform_euclidean(&self, v: &EuclideanVector<T>) -> EuclideanVector<T> {
        let p = Point::new(v.x(), v.y(), v.z());
        let transformed = self.transform_point(&p);
        // Safe because motor preserves weight of finite points
        EuclideanVector::new(
            transformed.e1 / transformed.e0,
            transformed.e2 / transformed.e0,
            transformed.e3 / transformed.e0,
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
    /// Multivector has non-zero trivector (grade 3) components where not expected.
    NonZeroTrivector,
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
            Self::NonZeroTrivector => write!(f, "multivector has non-zero trivector components"),
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
    use crate::specialized::euclidean::dim3::Bivector as EucBivector;
    use crate::specialized::euclidean::dim3::Rotor as EucRotor;
    use crate::specialized::euclidean::dim3::Vector;
    use crate::specialized::projective::dim3::arbitrary::UnitMotor;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    #[test]
    fn point_roundtrip() {
        let point = Point::new(3.0, 4.0, 5.0);
        let mv: Multivector<f64, Projective3> = point.into();
        let back = Point::try_from(mv).unwrap();

        assert!(abs_diff_eq!(point.e1, back.e1, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e2, back.e2, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e3, back.e3, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e0, back.e0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn ideal_point_roundtrip() {
        let point = Point::ideal(1.0, 2.0, 3.0);
        let mv: Multivector<f64, Projective3> = point.into();
        let back = Point::try_from(mv).unwrap();

        assert!(abs_diff_eq!(point.e1, back.e1, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e2, back.e2, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e3, back.e3, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(point.e0, back.e0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_roundtrip() {
        let motor = Motor::new_unchecked(0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.05);
        let mv: Multivector<f64, Projective3> = motor.into();
        let back = Motor::try_from(mv).unwrap();

        assert!(abs_diff_eq!(motor.s(), back.s(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(
            motor.e23(),
            back.e23(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e31(),
            back.e31(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e12(),
            back.e12(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e01(),
            back.e01(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e02(),
            back.e02(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e03(),
            back.e03(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e0123(),
            back.e0123(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn rotation_motor_roundtrip() {
        let motor = Motor::from_rotation_z(std::f64::consts::FRAC_PI_4);
        let mv: Multivector<f64, Projective3> = motor.into();
        let back = Motor::try_from(mv).unwrap();

        assert!(abs_diff_eq!(motor.s(), back.s(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(
            motor.e12(),
            back.e12(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn translation_motor_roundtrip() {
        let motor = Motor::from_translation(1.0, 2.0, 3.0);
        let mv: Multivector<f64, Projective3> = motor.into();
        let back = Motor::try_from(mv).unwrap();

        assert!(abs_diff_eq!(motor.s(), back.s(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(
            motor.e01(),
            back.e01(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e02(),
            back.e02(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            motor.e03(),
            back.e03(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn conversion_error_motor_to_point() {
        // Try to convert a motor (even) to a point (grade 1)
        let motor = Motor::from_rotation_z(0.5);
        let mv: Multivector<f64, Projective3> = motor.into();
        let result = Point::try_from(mv);
        assert!(result.is_err());
    }

    #[test]
    fn conversion_error_point_to_motor() {
        // Try to convert a point (grade 1) to a motor (even)
        let point = Point::new(1.0, 2.0, 3.0);
        let mv: Multivector<f64, Projective3> = point.into();
        let result = Motor::try_from(mv);
        assert!(result.is_err());
    }

    #[test]
    fn multivector_grade_extraction() {
        // Verify that conversion correctly extracts only the relevant grades
        let motor = Motor::from_translation(1.0, 2.0, 3.0);
        let mv: Multivector<f64, Projective3> = motor.into();

        // Grade 1 should be zero
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E1)),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E2)),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E3)),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Grade 0 should be 1 (identity rotation)
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(SCALAR)),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Translation bivectors: Motor stores d/2, and conversion negates
        // Motor has e01=0.5, e02=1.0, e03=1.5 (half of 1, 2, 3)
        // In multivector, these are negated (e₀ᵢ = -eᵢ₀ convention)
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E01)),
            -0.5,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E02)),
            -1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            mv.get(Blade::from_index(E03)),
            -1.5,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // Property-based tests
    proptest! {
        #[test]
        fn point_multivector_roundtrip(p in any::<Point<f64>>()) {
            let mv: Multivector<f64, Projective3> = p.into();
            let back = Point::try_from(mv).unwrap();

            prop_assert!(abs_diff_eq!(p.e1, back.e1, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.e2, back.e2, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.e3, back.e3, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.e0, back.e0, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn motor_multivector_roundtrip(m in any::<UnitMotor<f64>>()) {
            let mv: Multivector<f64, Projective3> = (*m).into();
            let back = Motor::try_from(mv).unwrap();

            prop_assert!(abs_diff_eq!(m.s(), back.s(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(m.e23(), back.e23(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(m.e31(), back.e31(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(m.e12(), back.e12(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(m.e01(), back.e01(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(m.e02(), back.e02(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(m.e03(), back.e03(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(m.e0123(), back.e0123(), epsilon = ABS_DIFF_EQ_EPS));
        }

        // ====================================================================
        // Rotor <-> Motor interop tests
        // ====================================================================

        /// Tests Rotor -> Motor -> Rotor roundtrip preserves rotation behavior.
        #[test]
        fn rotor_motor_roundtrip(
            r in any::<EucRotor<f64>>(),
            vx in -10.0f64..10.0, vy in -10.0f64..10.0, vz in -10.0f64..10.0,
        ) {
            let r = r.normalized();
            let motor: Motor<f64> = r.into();
            let back: EuclideanRotor<f64> = motor.into();

            // Compare by rotating a vector
            let v = Vector::new(vx, vy, vz);
            let result_orig = r.rotate(v);
            let result_back = back.rotate(v);

            prop_assert!(abs_diff_eq!(result_orig.x(), result_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.y(), result_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.z(), result_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests Motor -> Rotor -> Motor roundtrip preserves rotation behavior.
        #[test]
        fn motor_rotor_roundtrip(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let axis = Vector::new(ax / len, ay / len, az / len);

            // Create a pure rotation motor
            let motor = Motor::from_axis_angle(&axis, angle);
            let rotor: EuclideanRotor<f64> = motor.into();
            let back: Motor<f64> = rotor.into();

            // Compare by transforming a point
            let p = Point::new(px, py, pz);
            let result_orig = motor.transform_point(&p);
            let result_back = back.transform_point(&p);

            prop_assert!(abs_diff_eq!(result_orig.x(), result_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.y(), result_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.z(), result_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests that Rotor rotation matches Motor rotation.
        #[test]
        fn rotor_motor_rotation_equivalence(
            r in any::<EucRotor<f64>>(),
            x in -10.0f64..10.0, y in -10.0f64..10.0, z in -10.0f64..10.0,
        ) {
            let r = r.normalized();
            let motor: Motor<f64> = r.into();

            // Rotate with Euclidean rotor
            let v = Vector::new(x, y, z);
            let rotated_euc = r.rotate(v);

            // Transform with PGA motor
            let p = Point::new(x, y, z);
            let rotated_pga = motor.transform_point(&p);

            prop_assert!(abs_diff_eq!(rotated_euc.x(), rotated_pga.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotated_euc.y(), rotated_pga.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotated_euc.z(), rotated_pga.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests that Motor -> Rotor extracts correct rotation from composed motor.
        #[test]
        fn motor_to_rotor_extracts_rotation(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0, tz in -10.0f64..10.0,
            vx in -10.0f64..10.0, vy in -10.0f64..10.0, vz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let axis = Vector::new(ax / len, ay / len, az / len);

            // Create a motor with both rotation and translation
            let rotation = Motor::from_axis_angle(&axis, angle);
            let translation = Motor::from_translation(tx, ty, tz);
            let composed = rotation.compose(&translation);

            // Extract just the rotation
            let rotor: EuclideanRotor<f64> = composed.into();

            // The extracted rotor should match the original rotation
            let expected_rotor = EuclideanRotor::from_angle_plane(
                angle,
                EucBivector::new(axis.z(), -axis.y(), axis.x())  // Map axis to bivector
            );
            let v = Vector::new(vx, vy, vz);

            let result = rotor.rotate(v);
            let expected = expected_rotor.rotate(v);

            prop_assert!(abs_diff_eq!(result.x(), expected.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result.y(), expected.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result.z(), expected.z(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn identity_rotor_to_motor() {
        let rotor = EuclideanRotor::<f64>::identity();
        let motor: Motor<f64> = rotor.into();

        // Should be identity motor (no rotation, no translation)
        assert!(abs_diff_eq!(motor.s(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e23(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e31(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e12(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e01(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e02(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e03(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn identity_motor_to_rotor() {
        let motor = Motor::<f64>::identity();
        let rotor: EuclideanRotor<f64> = motor.into();

        // Should be identity rotor
        assert!(abs_diff_eq!(rotor.s(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rotor.b().xy(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rotor.b().xz(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rotor.b().yz(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}

// Note: A test comparing `transform_point` with generic sandwich product
// `M P M̃` is not included because the generic geometric product doesn't
// correctly implement PGA point transformation. In PGA Cl(3,0,1), the null
// vector e₀² = 0 causes certain products to vanish, making the naive
// sandwich product insufficient for translation. The specialized
// `transform_point` uses an optimized formula that correctly handles the
// PGA algebra.
