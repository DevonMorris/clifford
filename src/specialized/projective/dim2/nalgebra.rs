//! nalgebra interoperability for 2D Projective GA types.
//!
//! This module provides bidirectional conversions between clifford's 2D PGA types
//! and nalgebra's equivalent types.
//!
//! Enable with feature `nalgebra-0_32`, `nalgebra-0_33`, or `nalgebra-0_34`.
//!
//! # Conversions
//!
//! | clifford | nalgebra | Notes |
//! |----------|----------|-------|
//! | [`Point<T>`] | [`na::Point2<T>`] | Homogeneous ↔ Cartesian coordinates |
//! | [`Motor<T>`] | [`na::Isometry2<T>`] | Rigid transformation (rotation + translation) |
//!
//! # Motor ↔ Isometry Correspondence
//!
//! Both PGA motors and nalgebra isometries represent rigid transformations
//! (rotation + translation), but with different internal representations.
//!
//! **Motor representation** (2D PGA with dual form, grades 1+3):
//! - Components `(e0, e012)` encode the rotation: e0 = -sin(θ/2), e012 = cos(θ/2)
//! - Components `(e1, e2)` encode the translation: e1 = dy/2, e2 = -dx/2
//!
//! This dual form matches 3D PGA where motors include the pseudoscalar,
//! allowing the antisandwich product to work directly for transformations.
//!
//! **Isometry representation**:
//! - `UnitComplex` for rotation
//! - `Translation2` for translation

use core::fmt;

#[cfg(feature = "nalgebra-0_32")]
use nalgebra_0_32 as na;
#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use crate::scalar::Float;

use super::{Motor, Point};

// ============================================================================
// Error type
// ============================================================================

/// Error when converting from clifford types to nalgebra types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NalgebraConversionError {
    /// Point is ideal (at infinity), cannot convert to finite coordinates.
    IdealPoint,
}

impl fmt::Display for NalgebraConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IdealPoint => write!(f, "point is ideal (at infinity)"),
        }
    }
}

impl std::error::Error for NalgebraConversionError {}

// ============================================================================
// Point <-> Point2
// ============================================================================

impl<T: Float + na::Scalar> From<na::Point2<T>> for Point<T> {
    /// Converts a nalgebra 2D point to a PGA point.
    ///
    /// Creates a finite point with homogeneous weight `w = 1`.
    #[inline]
    fn from(p: na::Point2<T>) -> Self {
        Point::from_cartesian(p.x, p.y)
    }
}

impl<T: Float + na::Scalar> TryFrom<Point<T>> for na::Point2<T> {
    type Error = NalgebraConversionError;

    /// Attempts to convert a PGA point to a nalgebra 2D point.
    ///
    /// Returns an error if the point is ideal (at infinity).
    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0().abs() < T::epsilon() {
            return Err(NalgebraConversionError::IdealPoint);
        }
        Ok(na::Point2::new(p.x(), p.y()))
    }
}

// ============================================================================
// Motor <-> Isometry2
// ============================================================================

impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry2<T> {
    /// Converts a PGA motor to a nalgebra 2D isometry.
    ///
    /// # Mathematical Correspondence
    ///
    /// The motor's internal representation is "translate then rotate" (translation
    /// applied first, then rotation). nalgebra's Isometry2 uses "rotate then translate".
    ///
    /// For these to produce the same transformation:
    /// - Motor effect: Rot(θ) × (p + t') = Rot(θ)×p + Rot(θ)×t'
    /// - nalgebra effect: Rot(θ)×p + t
    ///
    /// So: t = Rot(θ) × t' (we post-rotate the decoded translation)
    fn from(motor: Motor<T>) -> Self {
        // Extract rotation from (e0, e012) = (-sin(θ/2), cos(θ/2))
        let rotation = na::UnitComplex::new(motor.rotation_angle());

        // Decode motor's internal translation (t')
        let cos_half = motor.e012();
        let sin_half = -motor.e0();
        let e1_t = motor.e1() * cos_half + motor.e2() * sin_half;
        let e2_t = motor.e2() * cos_half - motor.e1() * sin_half;

        // Internal translation t' = (-2*e2_t, 2*e1_t)
        let t_prime_x = -e2_t * T::TWO;
        let t_prime_y = e1_t * T::TWO;

        // Post-rotate by θ to get nalgebra translation: t = Rot(θ) × t'
        // cos(θ) = cos²(θ/2) - sin²(θ/2), sin(θ) = 2*sin(θ/2)*cos(θ/2)
        let cos_theta = cos_half * cos_half - sin_half * sin_half;
        let sin_theta = T::TWO * sin_half * cos_half;

        let tx = t_prime_x * cos_theta - t_prime_y * sin_theta;
        let ty = t_prime_x * sin_theta + t_prime_y * cos_theta;

        let translation = na::Translation2::new(tx, ty);

        na::Isometry2::from_parts(translation, rotation)
    }
}

impl<T: Float + na::RealField> From<na::Isometry2<T>> for Motor<T> {
    /// Converts a nalgebra 2D isometry to a PGA motor.
    ///
    /// # Mathematical Correspondence
    ///
    /// nalgebra's Isometry2 uses "rotate then translate": Rot(θ)×p + t
    /// The motor's internal representation is "translate then rotate": Rot(θ)×(p + t')
    ///
    /// For these to produce the same transformation:
    /// - t = Rot(θ) × t'
    /// - So: t' = Rot(-θ) × t (we pre-rotate the translation)
    fn from(iso: na::Isometry2<T>) -> Self {
        let angle = iso.rotation.angle();
        let t = iso.translation.vector;
        let half = angle / T::TWO;
        let cos_half = num_traits::Float::cos(half);
        let sin_half = num_traits::Float::sin(half);

        // Pre-rotate translation by -θ to get internal translation t'
        // Rot(-θ) = [[cos, sin], [-sin, cos]]
        let cos_theta = cos_half * cos_half - sin_half * sin_half;
        let sin_theta = T::TWO * sin_half * cos_half;

        let t_prime_x = t.x * cos_theta + t.y * sin_theta;
        let t_prime_y = -t.x * sin_theta + t.y * cos_theta;

        // Encode t' using standard motor encoding
        // e1 = t'_y/2 * cos(θ/2) + t'_x/2 * sin(θ/2)
        // e2 = -t'_x/2 * cos(θ/2) + t'_y/2 * sin(θ/2)
        Motor::new(
            t_prime_y / T::TWO * cos_half + t_prime_x / T::TWO * sin_half,
            -t_prime_x / T::TWO * cos_half + t_prime_y / T::TWO * sin_half,
            -sin_half,
            cos_half,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::Transform;
    use crate::specialized::projective::dim2::UnitizedPoint;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;
    use proptest::prelude::*;

    /// Epsilon for floating-point comparisons in tests.
    const EPS: f64 = RELATIVE_EQ_EPS;

    // ========================================================================
    // Point round-trip tests
    // ========================================================================

    proptest! {
        #[test]
        fn point_roundtrip(p in any::<UnitizedPoint<f64>>()) {
            let na_p: na::Point2<f64> = (*p).try_into().unwrap();
            let back: Point<f64> = na_p.into();

            // Compare Cartesian coordinates
            prop_assert!(relative_eq!(p.x(), back.x(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(p.y(), back.y(), epsilon = EPS, max_relative = EPS));
        }

        #[test]
        fn point_try_from_finite(p in any::<UnitizedPoint<f64>>()) {
            let result: Result<na::Point2<f64>, _> = (*p).try_into();
            prop_assert!(result.is_ok());

            let na_p = result.unwrap();
            prop_assert!(relative_eq!(p.x(), na_p.x, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(p.y(), na_p.y, epsilon = EPS, max_relative = EPS));
        }
    }

    #[test]
    fn point_try_from_ideal_fails() {
        let ideal = Point::<f64>::ideal(1.0, 0.0);
        let result: Result<na::Point2<f64>, _> = ideal.try_into();
        assert!(matches!(result, Err(NalgebraConversionError::IdealPoint)));
    }

    // ========================================================================
    // Motor conversion tests
    // ========================================================================

    #[test]
    fn identity_motor_to_isometry() {
        let motor = Motor::<f64>::identity();
        let iso: na::Isometry2<f64> = motor.into();

        // Identity isometry should have zero rotation and zero translation
        assert!(relative_eq!(
            iso.rotation.angle(),
            0.0,
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            iso.translation.vector.x,
            0.0,
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            iso.translation.vector.y,
            0.0,
            epsilon = EPS,
            max_relative = EPS
        ));
    }

    #[test]
    fn isometry_identity_to_motor() {
        let iso = na::Isometry2::<f64>::identity();
        let motor: Motor<f64> = iso.into();

        // Identity motor in dual form: e1=0, e2=0, e0=0, e012=1
        assert!(relative_eq!(
            motor.e1(),
            0.0,
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            motor.e2(),
            0.0,
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            motor.e0(),
            0.0,
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            motor.e012(),
            1.0,
            epsilon = EPS,
            max_relative = EPS
        ));
    }

    proptest! {
        /// Test that translation factories produce equivalent results.
        #[test]
        fn translation_factory_equivalence(
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
        ) {
            let motor = Motor::from_translation(tx, ty);
            let iso: na::Isometry2<f64> = motor.into();

            // Translation should match
            prop_assert!(
                relative_eq!(iso.translation.vector.x, tx, epsilon = EPS, max_relative = EPS),
                "tx mismatch: {} vs {}", iso.translation.vector.x, tx
            );
            prop_assert!(
                relative_eq!(iso.translation.vector.y, ty, epsilon = EPS, max_relative = EPS),
                "ty mismatch: {} vs {}", iso.translation.vector.y, ty
            );
            // Rotation should be identity
            prop_assert!(
                relative_eq!(iso.rotation.angle(), 0.0, epsilon = EPS, max_relative = EPS),
                "rotation should be 0, got {}", iso.rotation.angle()
            );
        }

        /// Test that rotation factories produce equivalent rotation angles.
        #[test]
        fn rotation_factory_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
        ) {
            let motor = Motor::from_rotation(angle);
            let iso: na::Isometry2<f64> = motor.into();

            // Rotation angle should match
            prop_assert!(
                relative_eq!(iso.rotation.angle(), angle, epsilon = EPS, max_relative = EPS),
                "angle mismatch: {} vs {}", iso.rotation.angle(), angle
            );
            // Translation should be zero
            prop_assert!(
                relative_eq!(iso.translation.vector.x, 0.0, epsilon = EPS, max_relative = EPS),
                "tx should be 0, got {}", iso.translation.vector.x
            );
            prop_assert!(
                relative_eq!(iso.translation.vector.y, 0.0, epsilon = EPS, max_relative = EPS),
                "ty should be 0, got {}", iso.translation.vector.y
            );
        }

        /// Test Motor → Isometry → Motor roundtrip for pure translations.
        #[test]
        fn translation_roundtrip(
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
        ) {
            let motor = Motor::from_translation(tx, ty);
            let iso: na::Isometry2<f64> = motor.into();
            let back: Motor<f64> = iso.into();

            // Motor components should match (dual form: e1, e2, e0, e012)
            prop_assert!(relative_eq!(motor.e1(), back.e1(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(motor.e2(), back.e2(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(motor.e0(), back.e0(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(motor.e012(), back.e012(), epsilon = EPS, max_relative = EPS));
        }

        /// Test Motor → Isometry → Motor roundtrip for pure rotations.
        #[test]
        fn rotation_roundtrip(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
        ) {
            let motor = Motor::from_rotation(angle);
            let iso: na::Isometry2<f64> = motor.into();
            let back: Motor<f64> = iso.into();

            // Motor components should match (dual form: e1, e2, e0, e012)
            prop_assert!(relative_eq!(motor.e1(), back.e1(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(motor.e2(), back.e2(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(motor.e0(), back.e0(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(motor.e012(), back.e012(), epsilon = EPS, max_relative = EPS));
        }
    }

    // ========================================================================
    // Operational equivalence tests
    // ========================================================================

    proptest! {
        /// Test that transform_point gives equivalent results via clifford and nalgebra.
        #[test]
        fn transform_point_equivalence(
            px in -100.0f64..100.0,
            py in -100.0f64..100.0,
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
        ) {
            // Build isometry in nalgebra: rotation then translation
            let iso = na::Isometry2::new(na::Vector2::new(tx, ty), angle);

            // Convert to motor
            let motor: Motor<f64> = iso.into();

            // Transform a point using both
            let p_cliff = Point::from_cartesian(px, py);
            let p_na = na::Point2::new(px, py);

            let result_cliff = Transform::transform(&motor, &p_cliff);
            let result_na = iso.transform_point(&p_na);

            prop_assert!(
                relative_eq!(result_cliff.x(), result_na.x, epsilon = EPS, max_relative = EPS),
                "x mismatch: clifford {} vs nalgebra {}", result_cliff.x(), result_na.x
            );
            prop_assert!(
                relative_eq!(result_cliff.y(), result_na.y, epsilon = EPS, max_relative = EPS),
                "y mismatch: clifford {} vs nalgebra {}", result_cliff.y(), result_na.y
            );
        }

        /// Test that motor inverse matches isometry inverse.
        #[test]
        fn inverse_equivalence(
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
        ) {
            let iso = na::Isometry2::new(na::Vector2::new(tx, ty), angle);
            let motor: Motor<f64> = iso.into();

            // Transform a point, then apply inverse
            let test_pt = na::Point2::new(3.0, 4.0);
            let test_cliff = Point::from_cartesian(3.0, 4.0);

            let transformed_na = iso.transform_point(&test_pt);
            let transformed_cliff = Transform::transform(&motor, &test_cliff);

            // Apply inverse
            let back_na = iso.inverse().transform_point(&transformed_na);
            let back_cliff = motor.inverse().transform(&transformed_cliff);

            prop_assert!(
                relative_eq!(back_cliff.x(), back_na.x, epsilon = EPS, max_relative = EPS),
                "inverse x mismatch: clifford {} vs nalgebra {}", back_cliff.x(), back_na.x
            );
            prop_assert!(
                relative_eq!(back_cliff.y(), back_na.y, epsilon = EPS, max_relative = EPS),
                "inverse y mismatch: clifford {} vs nalgebra {}", back_cliff.y(), back_na.y
            );

            // Also verify round-trip back to original point
            prop_assert!(
                relative_eq!(back_cliff.x(), test_cliff.x(), epsilon = EPS, max_relative = EPS),
                "round-trip x mismatch: {} vs {}", back_cliff.x(), test_cliff.x()
            );
            prop_assert!(
                relative_eq!(back_cliff.y(), test_cliff.y(), epsilon = EPS, max_relative = EPS),
                "round-trip y mismatch: {} vs {}", back_cliff.y(), test_cliff.y()
            );
        }
    }
}
