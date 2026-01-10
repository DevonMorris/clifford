//! nalgebra interoperability for 2D PGA types.
//!
//! This module provides conversions between 2D PGA types and nalgebra types:
//! - [`Point`] ↔ [`nalgebra::Point2`]
//! - [`Motor`] ↔ [`nalgebra::Isometry2`]

use crate::scalar::Float;

use super::types::{Motor, Point};

#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

// ============================================================================
// Point <-> Point2
// ============================================================================

impl<T> From<na::Point2<T>> for Point<T>
where
    T: Float + na::Scalar,
{
    /// Converts a nalgebra [`Point2`](na::Point2) to a PGA [`Point`].
    ///
    /// The resulting point has homogeneous weight `w = 1`.
    fn from(p: na::Point2<T>) -> Self {
        Point::new(p.x, p.y)
    }
}

impl<T> TryFrom<Point<T>> for na::Point2<T>
where
    T: Float + na::Scalar,
{
    type Error = PointConversionError;

    /// Tries to convert a PGA [`Point`] to a nalgebra [`Point2`](na::Point2).
    ///
    /// Returns an error if the point is at infinity (weight ≈ 0).
    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0.abs() < T::epsilon() {
            return Err(PointConversionError::PointAtInfinity);
        }
        Ok(na::Point2::new(p.x(), p.y()))
    }
}

// ============================================================================
// Motor <-> Isometry2
// ============================================================================

impl<T> From<na::Isometry2<T>> for Motor<T>
where
    T: Float + na::RealField + Copy,
{
    /// Converts a nalgebra [`Isometry2`](na::Isometry2) to a PGA [`Motor`].
    fn from(iso: na::Isometry2<T>) -> Self {
        let angle = iso.rotation.angle();
        let translation = iso.translation.vector;

        let rotation = Motor::from_rotation(angle);
        let trans = Motor::from_translation(translation.x, translation.y);

        // nalgebra applies rotation first, then translation
        // In PGA compose, leftmost motor acts first: rotation.compose(&trans) = R * T
        rotation.compose(&trans)
    }
}

impl<T> From<Motor<T>> for na::Isometry2<T>
where
    T: Float + na::RealField + Copy,
{
    /// Converts a PGA [`Motor`] to a nalgebra [`Isometry2`](na::Isometry2).
    ///
    /// Note: This extracts the rotation angle and translation from the motor.
    /// For non-unit motors, the result may not be exact.
    fn from(m: Motor<T>) -> Self {
        let angle = m.rotation_angle();
        let (tx, ty) = m.translation();

        // Note: for a composed motor, the translation extraction is approximate
        // For accurate extraction, we'd need to decompose the motor properly
        na::Isometry2::new(na::Vector2::new(tx, ty), angle)
    }
}

// ============================================================================
// Error types
// ============================================================================

/// Error type for point conversions to nalgebra types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointConversionError {
    /// The point is at infinity (homogeneous weight ≈ 0).
    PointAtInfinity,
}

impl core::fmt::Display for PointConversionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::PointAtInfinity => write!(f, "point is at infinity (w ≈ 0)"),
        }
    }
}

impl std::error::Error for PointConversionError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn point_roundtrip(x in -100.0f64..100.0, y in -100.0f64..100.0) {
            let pga_point = Point::new(x, y);
            let na_point: na::Point2<f64> = pga_point.try_into().unwrap();
            let back: Point<f64> = na_point.into();

            prop_assert!(abs_diff_eq!(pga_point.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_point.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn motor_transform_matches_isometry(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            // Create motor and isometry representing the same transformation
            // nalgebra: first rotate, then translate
            // PGA compose: leftmost acts first, so rotation.compose(&translation)
            let rotation = Motor::from_rotation(angle);
            let translation = Motor::from_translation(tx, ty);
            let motor = rotation.compose(&translation);

            let iso = na::Isometry2::new(na::Vector2::new(tx, ty), angle);

            // Transform a point using both
            let pga_point = Point::new(px, py);
            let pga_result = motor.transform_point(&pga_point);

            let na_point = na::Point2::new(px, py);
            let na_result = iso * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn isometry_to_motor_roundtrip(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            // Create an isometry and convert to motor
            let iso = na::Isometry2::new(na::Vector2::new(tx, ty), angle);
            let motor: Motor<f64> = iso.into();

            // Transform a point with both and compare
            let na_point = na::Point2::new(px, py);
            let na_result = iso * na_point;

            let pga_point = Point::new(px, py);
            let pga_result = motor.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn ideal_point_conversion_fails() {
        let ideal = Point::<f64>::ideal(1.0, 0.0);
        let result: Result<na::Point2<f64>, _> = ideal.try_into();
        assert!(result.is_err());
    }

    // ========================================================================
    // Composition order consistency tests
    // ========================================================================
    //
    // In nalgebra: (iso1 * iso2) * p = iso1 * (iso2 * p)
    //   - iso1 * iso2 applies iso2 first, then iso1
    //
    // In PGA: m1.compose(&m2) applies m1 first, then m2
    //   - So nalgebra's iso1 * iso2 corresponds to our m2.compose(&m1)
    //   - Or equivalently, our m1.compose(&m2) corresponds to nalgebra's m2 * m1

    proptest! {
        /// Verifies that Motor composition matches nalgebra Isometry multiplication
        /// with the appropriate order reversal.
        #[test]
        fn composition_order_matches_nalgebra(
            angle1 in -std::f64::consts::PI..std::f64::consts::PI,
            tx1 in -10.0f64..10.0, ty1 in -10.0f64..10.0,
            angle2 in -std::f64::consts::PI..std::f64::consts::PI,
            tx2 in -10.0f64..10.0, ty2 in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            // Create two isometries and motors
            let iso1 = na::Isometry2::new(na::Vector2::new(tx1, ty1), angle1);
            let iso2 = na::Isometry2::new(na::Vector2::new(tx2, ty2), angle2);
            let motor1: Motor<f64> = iso1.into();
            let motor2: Motor<f64> = iso2.into();

            // nalgebra: iso1 * iso2 applies iso2 first, then iso1
            let na_composed = iso1 * iso2;
            let na_point = na::Point2::new(px, py);
            let na_result = na_composed * na_point;

            // PGA: motor2.compose(&motor1) applies motor2 first, then motor1
            // This should match nalgebra's iso1 * iso2
            let pga_composed = motor2.compose(&motor1);
            let pga_point = Point::new(px, py);
            let pga_result = pga_composed.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Verifies that Motor inverse matches nalgebra Isometry inverse.
        #[test]
        fn inverse_matches_nalgebra(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            let iso = na::Isometry2::new(na::Vector2::new(tx, ty), angle);
            let motor: Motor<f64> = iso.into();

            // Transform a point, then apply inverse
            let na_point = na::Point2::new(px, py);
            let na_transformed = iso * na_point;
            let na_back = iso.inverse() * na_transformed;

            let pga_point = Point::new(px, py);
            let pga_transformed = motor.transform_point(&pga_point);
            let pga_back = motor.inverse().transform_point(&pga_transformed);

            // Both should return to the original point
            prop_assert!(abs_diff_eq!(pga_back.x(), na_back.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_back.y(), na_back.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_back.x(), px, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_back.y(), py, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Verifies that pure rotations match between PGA and nalgebra.
        #[test]
        fn pure_rotation_matches_nalgebra(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            let na_rot = na::Isometry2::rotation(angle);
            let pga_rot = Motor::from_rotation(angle);

            let na_point = na::Point2::new(px, py);
            let na_result = na_rot * na_point;

            let pga_point = Point::new(px, py);
            let pga_result = pga_rot.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Verifies that pure translations match between PGA and nalgebra.
        #[test]
        fn pure_translation_matches_nalgebra(
            tx in -10.0f64..10.0, ty in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            let na_trans = na::Isometry2::translation(tx, ty);
            let pga_trans = Motor::from_translation(tx, ty);

            let na_point = na::Point2::new(px, py);
            let na_result = na_trans * na_point;

            let pga_point = Point::new(px, py);
            let pga_result = pga_trans.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Verifies triple composition consistency.
        #[test]
        fn triple_composition_matches_nalgebra(
            angle1 in -std::f64::consts::PI..std::f64::consts::PI,
            tx1 in -5.0f64..5.0, ty1 in -5.0f64..5.0,
            angle2 in -std::f64::consts::PI..std::f64::consts::PI,
            tx2 in -5.0f64..5.0, ty2 in -5.0f64..5.0,
            angle3 in -std::f64::consts::PI..std::f64::consts::PI,
            tx3 in -5.0f64..5.0, ty3 in -5.0f64..5.0,
            px in -5.0f64..5.0, py in -5.0f64..5.0,
        ) {
            let iso1 = na::Isometry2::new(na::Vector2::new(tx1, ty1), angle1);
            let iso2 = na::Isometry2::new(na::Vector2::new(tx2, ty2), angle2);
            let iso3 = na::Isometry2::new(na::Vector2::new(tx3, ty3), angle3);

            let motor1: Motor<f64> = iso1.into();
            let motor2: Motor<f64> = iso2.into();
            let motor3: Motor<f64> = iso3.into();

            // nalgebra: iso1 * iso2 * iso3 applies iso3 first, then iso2, then iso1
            let na_composed = iso1 * iso2 * iso3;
            let na_point = na::Point2::new(px, py);
            let na_result = na_composed * na_point;

            // PGA: to match, we need motor3.compose(&motor2).compose(&motor1)
            let pga_composed = motor3.compose(&motor2).compose(&motor1);
            let pga_point = Point::new(px, py);
            let pga_result = pga_composed.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn identity_motor_matches_nalgebra_identity() {
        let na_identity = na::Isometry2::<f64>::identity();
        let pga_identity = Motor::<f64>::identity();

        // Both should leave points unchanged
        let points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (3.0, 4.0), (-2.5, 7.1)];
        for (px, py) in points {
            let na_point = na::Point2::new(px, py);
            let na_result = na_identity * na_point;

            let pga_point = Point::new(px, py);
            let pga_result = pga_identity.transform_point(&pga_point);

            assert!(abs_diff_eq!(
                pga_result.x(),
                na_result.x,
                epsilon = ABS_DIFF_EQ_EPS
            ));
            assert!(abs_diff_eq!(
                pga_result.y(),
                na_result.y,
                epsilon = ABS_DIFF_EQ_EPS
            ));
            assert!(abs_diff_eq!(pga_result.x(), px, epsilon = ABS_DIFF_EQ_EPS));
            assert!(abs_diff_eq!(pga_result.y(), py, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
