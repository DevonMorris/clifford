//! nalgebra interoperability for 3D PGA types.
//!
//! This module provides conversions between 3D PGA types and nalgebra types:
//! - [`Point`] ↔ [`nalgebra::Point3`]
//! - [`Motor`] ↔ [`nalgebra::Isometry3`]

use crate::scalar::Float;

use super::types::Point;

// Motor will be used once we add Motor <-> Isometry3 conversions
#[allow(unused_imports)]
use super::types::Motor;

#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

// ============================================================================
// Point <-> Point3
// ============================================================================

impl<T> From<na::Point3<T>> for Point<T>
where
    T: Float + na::Scalar,
{
    /// Converts a nalgebra [`Point3`](na::Point3) to a PGA [`Point`].
    ///
    /// The resulting point has homogeneous weight `w = 1`.
    fn from(p: na::Point3<T>) -> Self {
        Point::new(p.x, p.y, p.z)
    }
}

impl<T> TryFrom<Point<T>> for na::Point3<T>
where
    T: Float + na::Scalar,
{
    type Error = PointConversionError;

    /// Tries to convert a PGA [`Point`] to a nalgebra [`Point3`](na::Point3).
    ///
    /// Returns an error if the point is at infinity (weight ≈ 0).
    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0.abs() < T::epsilon() {
            return Err(PointConversionError::PointAtInfinity);
        }
        Ok(na::Point3::new(p.x(), p.y(), p.z()))
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
    use std::f64::consts::FRAC_PI_2;

    // ========================================================================
    // Point conversion tests
    // ========================================================================

    proptest! {
        #[test]
        fn point_roundtrip(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
            let pga_point = Point::new(x, y, z);
            let na_point: na::Point3<f64> = pga_point.try_into().unwrap();
            let back: Point<f64> = na_point.into();

            prop_assert!(abs_diff_eq!(pga_point.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_point.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_point.z(), back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn ideal_point_conversion_fails() {
        let ideal = Point::<f64>::ideal(1.0, 0.0, 0.0);
        let result: Result<na::Point3<f64>, _> = ideal.try_into();
        assert!(result.is_err());
    }

    // ========================================================================
    // Identity motor tests
    // ========================================================================

    #[test]
    fn identity_motor_preserves_points() {
        let motor = Motor::<f64>::identity();
        let points = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
            Point::new(3.0, 4.0, 5.0),
            Point::new(-2.5, 7.1, -1.3),
        ];

        for p in points {
            let result = motor.transform_point(&p);
            assert!(
                abs_diff_eq!(result.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS),
                "x mismatch for point {:?}",
                p
            );
            assert!(
                abs_diff_eq!(result.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS),
                "y mismatch for point {:?}",
                p
            );
            assert!(
                abs_diff_eq!(result.z(), p.z(), epsilon = ABS_DIFF_EQ_EPS),
                "z mismatch for point {:?}",
                p
            );
        }
    }

    // ========================================================================
    // Translation tests
    // ========================================================================

    #[test]
    fn pure_translation_concrete_cases() {
        // Translate (1, 2, 3) by (10, 20, 30) -> (11, 22, 33)
        let motor = Motor::from_translation(10.0, 20.0, 30.0);
        let p = Point::new(1.0, 2.0, 3.0);
        let result = motor.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), 11.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 22.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 33.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    proptest! {
        #[test]
        fn pure_translation_matches_nalgebra(
            dx in -10.0f64..10.0, dy in -10.0f64..10.0, dz in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let pga_trans = Motor::from_translation(dx, dy, dz);
            let na_trans = na::Translation3::new(dx, dy, dz);

            let pga_point = Point::new(px, py, pz);
            let pga_result = pga_trans.transform_point(&pga_point);

            let na_point = na::Point3::new(px, py, pz);
            let na_result = na_trans * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Rotation tests - concrete cases first
    // ========================================================================

    #[test]
    fn rotation_z_90_degrees() {
        // (1, 0, 0) rotated 90° around Z -> (0, 1, 0)
        let motor = Motor::from_rotation_z(FRAC_PI_2);
        let p = Point::new(1.0, 0.0, 0.0);
        let result = motor.transform_point(&p);

        assert!(
            abs_diff_eq!(result.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "x = {} expected 0",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 1.0, epsilon = ABS_DIFF_EQ_EPS),
            "y = {} expected 1",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "z = {} expected 0",
            result.z()
        );
    }

    #[test]
    fn rotation_x_90_degrees() {
        // (0, 1, 0) rotated 90° around X -> (0, 0, 1)
        let motor = Motor::from_rotation_x(FRAC_PI_2);
        let p = Point::new(0.0, 1.0, 0.0);
        let result = motor.transform_point(&p);

        assert!(
            abs_diff_eq!(result.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "x = {} expected 0",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "y = {} expected 0",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), 1.0, epsilon = ABS_DIFF_EQ_EPS),
            "z = {} expected 1",
            result.z()
        );
    }

    #[test]
    fn rotation_y_90_degrees() {
        // (0, 0, 1) rotated 90° around Y -> (1, 0, 0)
        let motor = Motor::from_rotation_y(FRAC_PI_2);
        let p = Point::new(0.0, 0.0, 1.0);
        let result = motor.transform_point(&p);

        assert!(
            abs_diff_eq!(result.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS),
            "x = {} expected 1",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "y = {} expected 0",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "z = {} expected 0",
            result.z()
        );
    }

    // ========================================================================
    // Rotation tests - property tests against nalgebra
    // ========================================================================

    proptest! {
        #[test]
        fn rotation_z_matches_nalgebra(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let pga_rot = Motor::from_rotation_z(angle);
            let na_rot = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), angle);

            let pga_point = Point::new(px, py, pz);
            let pga_result = pga_rot.transform_point(&pga_point);

            let na_point = na::Point3::new(px, py, pz);
            let na_result = na_rot * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotation_x_matches_nalgebra(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let pga_rot = Motor::from_rotation_x(angle);
            let na_rot = na::UnitQuaternion::from_axis_angle(&na::Vector3::x_axis(), angle);

            let pga_point = Point::new(px, py, pz);
            let pga_result = pga_rot.transform_point(&pga_point);

            let na_point = na::Point3::new(px, py, pz);
            let na_result = na_rot * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotation_y_matches_nalgebra(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let pga_rot = Motor::from_rotation_y(angle);
            let na_rot = na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), angle);

            let pga_point = Point::new(px, py, pz);
            let pga_result = pga_rot.transform_point(&pga_point);

            let na_point = na::Point3::new(px, py, pz);
            let na_result = na_rot * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn from_axis_angle_matches_nalgebra(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let axis = (ax / len, ay / len, az / len);

            let pga_rot = Motor::from_axis_angle(axis, angle);
            let na_axis = na::Unit::new_normalize(na::Vector3::new(axis.0, axis.1, axis.2));
            let na_rot = na::UnitQuaternion::from_axis_angle(&na_axis, angle);

            let pga_point = Point::new(px, py, pz);
            let pga_result = pga_rot.transform_point(&pga_point);

            let na_point = na::Point3::new(px, py, pz);
            let na_result = na_rot * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Composition tests - concrete cases first
    // ========================================================================

    #[test]
    fn compose_rotation_then_translation() {
        // (1,0,0) -> rotate 90° around Z to (0,1,0) -> translate by (1,2,3) to (1,3,3)
        let rotation = Motor::from_rotation_z(FRAC_PI_2);
        let translation = Motor::from_translation(1.0, 2.0, 3.0);
        let combined = rotation.compose(&translation);

        let p = Point::new(1.0, 0.0, 0.0);
        let result = combined.transform_point(&p);

        assert!(
            abs_diff_eq!(result.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS),
            "x = {} expected 1",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 3.0, epsilon = ABS_DIFF_EQ_EPS),
            "y = {} expected 3",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS),
            "z = {} expected 3",
            result.z()
        );
    }

    #[test]
    fn compose_translation_then_rotation() {
        // (1,0,0) -> translate by (1,2,3) to (2,2,3) -> rotate 90° around Z to (-2,2,3)
        let translation = Motor::from_translation(1.0, 2.0, 3.0);
        let rotation = Motor::from_rotation_z(FRAC_PI_2);
        let combined = translation.compose(&rotation);

        let p = Point::new(1.0, 0.0, 0.0);
        let result = combined.transform_point(&p);

        assert!(
            abs_diff_eq!(result.x(), -2.0, epsilon = ABS_DIFF_EQ_EPS),
            "x = {} expected -2",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS),
            "y = {} expected 2",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS),
            "z = {} expected 3",
            result.z()
        );
    }

    #[test]
    fn compose_two_rotations() {
        // Rotate 90° around X, then 90° around Y
        // (1,0,0) -> rotate X 90° -> (1,0,0) (unchanged, on axis)
        // Then rotate Y 90° -> (0,0,-1)
        let rot_x = Motor::from_rotation_x(FRAC_PI_2);
        let rot_y = Motor::from_rotation_y(FRAC_PI_2);
        let combined = rot_x.compose(&rot_y);

        let p = Point::new(1.0, 0.0, 0.0);
        let result = combined.transform_point(&p);

        assert!(
            abs_diff_eq!(result.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "x = {} expected 0",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS),
            "y = {} expected 0",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), -1.0, epsilon = ABS_DIFF_EQ_EPS),
            "z = {} expected -1",
            result.z()
        );
    }

    // ========================================================================
    // Composition tests - property tests against nalgebra
    // ========================================================================

    proptest! {
        /// Verifies that Motor composition produces the same result as applying
        /// transformations in sequence.
        #[test]
        fn compose_equals_sequential_application(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0, tz in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let rotation = Motor::from_rotation_z(angle);
            let translation = Motor::from_translation(tx, ty, tz);
            let composed = rotation.compose(&translation);

            let p = Point::new(px, py, pz);

            // Composed transformation
            let result_composed = composed.transform_point(&p);

            // Sequential application: rotate first, then translate
            let intermediate = rotation.transform_point(&p);
            let result_sequential = translation.transform_point(&intermediate);

            prop_assert!(abs_diff_eq!(result_composed.x(), result_sequential.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_composed.y(), result_sequential.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_composed.z(), result_sequential.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Verifies composition order mapping to nalgebra.
        /// nalgebra: iso1 * iso2 applies iso2 first, then iso1
        /// PGA: motor2.compose(&motor1) applies motor2 first, then motor1
        /// So motor2.compose(&motor1) should match iso1 * iso2
        #[test]
        fn composition_order_matches_nalgebra(
            angle1 in -std::f64::consts::PI..std::f64::consts::PI,
            tx1 in -5.0f64..5.0, ty1 in -5.0f64..5.0, tz1 in -5.0f64..5.0,
            angle2 in -std::f64::consts::PI..std::f64::consts::PI,
            tx2 in -5.0f64..5.0, ty2 in -5.0f64..5.0, tz2 in -5.0f64..5.0,
            px in -5.0f64..5.0, py in -5.0f64..5.0, pz in -5.0f64..5.0,
        ) {
            // Create isometries and motors
            let iso1 = na::Isometry3::new(
                na::Vector3::new(tx1, ty1, tz1),
                na::Vector3::new(0.0, 0.0, angle1),
            );
            let iso2 = na::Isometry3::new(
                na::Vector3::new(tx2, ty2, tz2),
                na::Vector3::new(0.0, 0.0, angle2),
            );

            // Convert to motors by building equivalent transformation
            // Isometry3::new applies rotation first, then translation
            let rot1 = Motor::from_rotation_z(angle1);
            let trans1 = Motor::from_translation(tx1, ty1, tz1);
            let motor1 = rot1.compose(&trans1);

            let rot2 = Motor::from_rotation_z(angle2);
            let trans2 = Motor::from_translation(tx2, ty2, tz2);
            let motor2 = rot2.compose(&trans2);

            // nalgebra: iso1 * iso2 applies iso2 first, then iso1
            let na_composed = iso1 * iso2;
            let na_point = na::Point3::new(px, py, pz);
            let na_result = na_composed * na_point;

            // PGA: to match iso1 * iso2, we need motor2.compose(&motor1)
            // because compose applies left motor first
            let pga_composed = motor2.compose(&motor1);
            let pga_point = Point::new(px, py, pz);
            let pga_result = pga_composed.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
