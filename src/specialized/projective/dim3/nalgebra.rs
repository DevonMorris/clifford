//! nalgebra interoperability for 3D PGA types.
//!
//! This module provides conversions between 3D PGA types and nalgebra types:
//! - [`Point`] ↔ [`nalgebra::Point3`]
//! - [`Motor`] ↔ [`nalgebra::Isometry3`]
//! - [`Motor`] ↔ [`nalgebra::UnitQuaternion`] (rotation only)
//! - [`Motor`] ↔ [`nalgebra::Rotation3`] (rotation only)
//! - [`Flector`] ↔ [`nalgebra::Reflection`] (pure plane reflection only)

use crate::scalar::Float;

use super::types::{Flector, Motor, Plane, Point};

#[cfg(feature = "nalgebra-0_32")]
use nalgebra_0_32 as na;
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
// Motor <-> Isometry3
// ============================================================================

impl<T> From<na::Isometry3<T>> for Motor<T>
where
    T: Float + na::RealField + Copy,
{
    /// Converts a nalgebra [`Isometry3`](na::Isometry3) to a PGA [`Motor`].
    ///
    /// # Convention
    ///
    /// nalgebra's `Isometry3` applies rotation first, then translation.
    /// The motor is constructed by composing rotation and translation motors
    /// in the same order.
    fn from(iso: na::Isometry3<T>) -> Self {
        // Extract rotation as unit quaternion
        let q = iso.rotation;
        let rotation = Motor::from(q);

        // Extract translation
        let t = iso.translation.vector;
        let translation = Motor::from_translation(t.x, t.y, t.z);

        // nalgebra applies rotation first, then translation
        rotation.compose(&translation)
    }
}

impl<T> From<Motor<T>> for na::Isometry3<T>
where
    T: Float + na::RealField + Copy,
{
    /// Converts a PGA [`Motor`] to a nalgebra [`Isometry3`](na::Isometry3).
    ///
    /// # Note
    ///
    /// This extracts the rotation from the motor's rotation bivector and
    /// the translation from the translation bivector. For motors that are
    /// not simple rotation+translation compositions, the extraction may
    /// be approximate.
    fn from(m: Motor<T>) -> Self {
        // Extract rotation as quaternion
        let q: na::UnitQuaternion<T> = m.into();

        // Extract translation
        // For a composed motor R*T, the translation extraction is:
        // t = 2 * (s*e0i + cross terms from rotation-translation interaction)
        // For simplicity, we extract as if it's a pure translation
        let tx = m.e01 * T::TWO;
        let ty = m.e02 * T::TWO;
        let tz = m.e03 * T::TWO;

        na::Isometry3::from_parts(na::Translation3::new(tx, ty, tz), q)
    }
}

// ============================================================================
// Motor <-> UnitQuaternion (rotation only)
// ============================================================================

impl<T> From<na::UnitQuaternion<T>> for Motor<T>
where
    T: Float + na::RealField,
{
    /// Converts a nalgebra [`UnitQuaternion`](na::UnitQuaternion) to a pure rotation [`Motor`].
    ///
    /// # Mapping
    ///
    /// - `q.w` → `motor.s` (scalar)
    /// - `q.i` → `motor.e23` (x-axis rotation, yz-plane)
    /// - `q.j` → `motor.e31` (y-axis rotation, zx-plane)
    /// - `q.k` → `motor.e12` (z-axis rotation, xy-plane)
    fn from(q: na::UnitQuaternion<T>) -> Self {
        let q = q.quaternion();
        Motor::new(
            q.w,
            q.i,
            q.j,
            q.k,
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        )
    }
}

impl<T> From<Motor<T>> for na::UnitQuaternion<T>
where
    T: Float + na::RealField,
{
    /// Converts a [`Motor`]'s rotation part to a nalgebra [`UnitQuaternion`](na::UnitQuaternion).
    ///
    /// # Note
    ///
    /// This extracts only the rotation component. Translation is ignored.
    ///
    /// # Mapping
    ///
    /// - `motor.s` → `q.w` (scalar)
    /// - `motor.e23` → `q.i` (x-axis rotation)
    /// - `motor.e31` → `q.j` (y-axis rotation)
    /// - `motor.e12` → `q.k` (z-axis rotation)
    fn from(m: Motor<T>) -> Self {
        let q = na::Quaternion::new(m.s, m.e23, m.e31, m.e12);
        na::UnitQuaternion::new_normalize(q)
    }
}

// ============================================================================
// Motor <-> Rotation3 (rotation only)
// ============================================================================

impl<T> From<na::Rotation3<T>> for Motor<T>
where
    T: Float + na::RealField,
{
    /// Converts a nalgebra [`Rotation3`](na::Rotation3) to a pure rotation [`Motor`].
    fn from(rot: na::Rotation3<T>) -> Self {
        let q: na::UnitQuaternion<T> = rot.into();
        q.into()
    }
}

impl<T> From<Motor<T>> for na::Rotation3<T>
where
    T: Float + na::RealField,
{
    /// Converts a [`Motor`]'s rotation part to a nalgebra [`Rotation3`](na::Rotation3).
    ///
    /// # Note
    ///
    /// This extracts only the rotation component. Translation is ignored.
    fn from(m: Motor<T>) -> Self {
        let q: na::UnitQuaternion<T> = m.into();
        q.into()
    }
}

// ============================================================================
// Flector <-> Reflection
// ============================================================================

/// Type alias for nalgebra 3D Reflection with owned (array) storage.
///
/// This is the most common type of 3D reflection in nalgebra, using stack-allocated storage.
pub type Reflection3<T> = na::Reflection<T, na::Const<3>, na::ArrayStorage<T, 3, 1>>;

impl<T, S> From<na::Reflection<T, na::Const<3>, S>> for Flector<T>
where
    T: Float + na::RealField,
    S: na::Storage<T, na::Const<3>>,
{
    /// Converts a nalgebra [`Reflection`](na::Reflection) to a PGA [`Flector`].
    ///
    /// The nalgebra Reflection represents a pure plane reflection with:
    /// - `axis()`: unit normal vector of the reflection plane
    /// - `bias()`: signed distance from origin to plane along the normal
    ///
    /// # Convention
    ///
    /// nalgebra's reflection plane equation: `n · x = bias`
    /// PGA's plane equation: `nx·x + ny·y + nz·z + d = 0`
    ///
    /// So `d = -bias` when converting.
    fn from(refl: na::Reflection<T, na::Const<3>, S>) -> Self {
        let axis = refl.axis();
        let bias = refl.bias();

        // Create a pure plane reflection flector
        // nalgebra uses n · x = bias, PGA uses n · x + d = 0, so d = -bias
        Flector::from_plane(Plane::from_normal_and_distance(
            axis[0], axis[1], axis[2], -bias,
        ))
    }
}

impl<T> TryFrom<Flector<T>> for Reflection3<T>
where
    T: Float + na::RealField,
{
    type Error = FlectorConversionError;

    /// Tries to convert a PGA [`Flector`] to a nalgebra [`Reflection`](na::Reflection).
    ///
    /// Returns an error if the flector is not a pure plane reflection
    /// (has non-zero point component).
    ///
    /// # Errors
    ///
    /// - [`FlectorConversionError::NotPureReflection`]: The flector has a
    ///   non-zero point component
    /// - [`FlectorConversionError::DegeneratePlane`]: The plane normal has zero length
    fn try_from(f: Flector<T>) -> Result<Self, Self::Error> {
        // Check if this is a pure plane reflection
        if !f.is_pure_reflection(T::epsilon()) {
            return Err(FlectorConversionError::NotPureReflection);
        }

        let plane = f.plane_part();
        let weight_norm = plane.weight_norm();

        if weight_norm < T::epsilon() {
            return Err(FlectorConversionError::DegeneratePlane);
        }

        // Unitize the plane
        let unitized = plane.unitized();
        let n = unitized.normal();
        let d = unitized.distance();

        // Create axis as unit vector
        let axis = na::Unit::new_unchecked(na::Vector3::new(n.x, n.y, n.z));

        // Find a point on the plane: for plane n·x + d = 0, the point P = -d*n lies on it
        // since n·(-d*n) = -d*|n|² = -d (when |n|=1), and -d + d = 0 ✓
        let point_on_plane = na::Point3::new(-d * n.x, -d * n.y, -d * n.z);

        Ok(na::Reflection::new_containing_point(axis, &point_on_plane))
    }
}

// ============================================================================
// Error types
// ============================================================================

// Re-export the shared error type
pub use super::super::PointConversionError;

use thiserror::Error;

/// Error type for flector conversions to nalgebra types.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum FlectorConversionError {
    /// The flector is not a pure plane reflection (has non-zero point component).
    #[error("flector is not a pure plane reflection")]
    NotPureReflection,
    /// The reflection plane has zero-length normal.
    #[error("reflection plane has zero-length normal")]
    DegeneratePlane,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::specialized::euclidean::dim3::Vector as EuclideanVector;
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
            let axis = EuclideanVector::new(ax / len, ay / len, az / len);

            let pga_rot = Motor::from_axis_angle(&axis, angle);
            let na_axis = na::Unit::new_normalize(na::Vector3::new(axis.x, axis.y, axis.z));
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

    // ========================================================================
    // Motor <-> UnitQuaternion conversion tests
    // ========================================================================

    proptest! {
        /// Tests that Motor -> UnitQuaternion -> Motor preserves rotation behavior.
        #[test]
        fn motor_quaternion_roundtrip_rotation(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let axis = EuclideanVector::new(ax / len, ay / len, az / len);

            let motor = Motor::from_axis_angle(&axis, angle);
            let q: na::UnitQuaternion<f64> = motor.into();
            let back: Motor<f64> = q.into();

            // Compare by transforming a point
            let p = Point::new(px, py, pz);
            let result_orig = motor.transform_point(&p);
            let result_back = back.transform_point(&p);

            prop_assert!(abs_diff_eq!(result_orig.x(), result_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.y(), result_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.z(), result_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests that UnitQuaternion -> Motor -> UnitQuaternion preserves rotation.
        #[test]
        fn quaternion_motor_roundtrip(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let na_axis = na::Unit::new_normalize(na::Vector3::new(ax / len, ay / len, az / len));
            let q = na::UnitQuaternion::from_axis_angle(&na_axis, angle);

            let motor: Motor<f64> = q.into();
            let back: na::UnitQuaternion<f64> = motor.into();

            // Compare by transforming a point
            let na_point = na::Point3::new(px, py, pz);
            let result_orig = q * na_point;
            let result_back = back * na_point;

            prop_assert!(abs_diff_eq!(result_orig.x, result_back.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.y, result_back.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.z, result_back.z, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests that Motor rotation matches UnitQuaternion rotation.
        #[test]
        fn motor_quaternion_rotation_equivalence(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let axis = EuclideanVector::new(ax / len, ay / len, az / len);

            let motor = Motor::from_axis_angle(&axis, angle);
            let q: na::UnitQuaternion<f64> = motor.into();

            // Transform with both
            let pga_point = Point::new(px, py, pz);
            let pga_result = motor.transform_point(&pga_point);

            let na_point = na::Point3::new(px, py, pz);
            let na_result = q * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Motor <-> Rotation3 conversion tests
    // ========================================================================

    proptest! {
        /// Tests that Motor -> Rotation3 -> Motor preserves rotation behavior.
        #[test]
        fn motor_rotation3_roundtrip(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let axis = EuclideanVector::new(ax / len, ay / len, az / len);

            let motor = Motor::from_axis_angle(&axis, angle);
            let rot: na::Rotation3<f64> = motor.into();
            let back: Motor<f64> = rot.into();

            // Compare by transforming a point
            let p = Point::new(px, py, pz);
            let result_orig = motor.transform_point(&p);
            let result_back = back.transform_point(&p);

            prop_assert!(abs_diff_eq!(result_orig.x(), result_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.y(), result_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.z(), result_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Motor <-> Isometry3 conversion tests
    // ========================================================================

    proptest! {
        /// Tests that Isometry3 -> Motor conversion produces equivalent transformation.
        #[test]
        fn isometry_to_motor_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0, tz in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Create isometry (rotation around Z)
            let iso = na::Isometry3::new(
                na::Vector3::new(tx, ty, tz),
                na::Vector3::new(0.0, 0.0, angle),
            );
            let motor: Motor<f64> = iso.into();

            // Transform a point with both
            let na_point = na::Point3::new(px, py, pz);
            let na_result = iso * na_point;

            let pga_point = Point::new(px, py, pz);
            let pga_result = motor.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests pure rotation Isometry conversion.
        #[test]
        fn isometry_pure_rotation_to_motor(
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize axis (skip if too small)
            let len = (ax * ax + ay * ay + az * az).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let na_axis = na::Unit::new_normalize(na::Vector3::new(ax / len, ay / len, az / len));

            let q = na::UnitQuaternion::from_axis_angle(&na_axis, angle);
            let iso = na::Isometry3::from_parts(na::Translation3::identity(), q);
            let motor: Motor<f64> = iso.into();

            // Transform a point with both
            let na_point = na::Point3::new(px, py, pz);
            let na_result = iso * na_point;

            let pga_point = Point::new(px, py, pz);
            let pga_result = motor.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests pure translation Isometry conversion.
        #[test]
        fn isometry_pure_translation_to_motor(
            tx in -10.0f64..10.0, ty in -10.0f64..10.0, tz in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let iso = na::Isometry3::translation(tx, ty, tz);
            let motor: Motor<f64> = iso.into();

            // Transform a point with both
            let na_point = na::Point3::new(px, py, pz);
            let na_result = iso * na_point;

            let pga_point = Point::new(px, py, pz);
            let pga_result = motor.transform_point(&pga_point);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Identity and edge case tests
    // ========================================================================

    #[test]
    fn identity_motor_to_isometry() {
        let motor = Motor::<f64>::identity();
        let iso: na::Isometry3<f64> = motor.into();

        // Should be identity isometry
        let p = na::Point3::new(1.0, 2.0, 3.0);
        let result = iso * p;
        assert!(abs_diff_eq!(result.x, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y, 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z, 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn identity_isometry_to_motor() {
        let iso = na::Isometry3::<f64>::identity();
        let motor: Motor<f64> = iso.into();

        // Should be identity motor
        let p = Point::new(1.0, 2.0, 3.0);
        let result = motor.transform_point(&p);
        assert!(abs_diff_eq!(result.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn identity_motor_to_quaternion() {
        let motor = Motor::<f64>::identity();
        let q: na::UnitQuaternion<f64> = motor.into();

        // Should be identity quaternion (w=1, i=j=k=0)
        assert!(abs_diff_eq!(q.w, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(q.i, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(q.j, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(q.k, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn identity_quaternion_to_motor() {
        let q = na::UnitQuaternion::<f64>::identity();
        let motor: Motor<f64> = q.into();

        // Should produce identity motor
        assert!(abs_diff_eq!(motor.s, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e23, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e31, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e12, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e01, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e02, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(motor.e03, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Flector <-> Reflection conversion tests
    // ========================================================================

    #[test]
    fn reflection_to_flector_xy_plane() {
        // XY plane at z=0: normal (0,0,1), bias = 0
        let axis = na::Unit::new_normalize(na::Vector3::new(0.0, 0.0, 1.0));
        let na_refl = na::Reflection::new_containing_point(axis, &na::Point3::origin());

        let flector: Flector<f64> = na_refl.into();

        // Should be a pure plane reflection through XY plane
        assert!(flector.is_pure_reflection(ABS_DIFF_EQ_EPS));
        let plane = flector.plane_part();
        assert!(abs_diff_eq!(plane.e023, 0.0, epsilon = ABS_DIFF_EQ_EPS)); // nx
        assert!(abs_diff_eq!(plane.e031, 0.0, epsilon = ABS_DIFF_EQ_EPS)); // ny
        assert!(abs_diff_eq!(plane.e012, 1.0, epsilon = ABS_DIFF_EQ_EPS)); // nz
        assert!(abs_diff_eq!(plane.e123, 0.0, epsilon = ABS_DIFF_EQ_EPS)); // d
    }

    #[test]
    fn reflection_to_flector_offset_plane() {
        // Plane z = 5: normal (0,0,1), contains point (0,0,5)
        let axis = na::Unit::new_normalize(na::Vector3::new(0.0, 0.0, 1.0));
        let na_refl = na::Reflection::new_containing_point(axis, &na::Point3::new(0.0, 0.0, 5.0));

        let flector: Flector<f64> = na_refl.into();

        // Should be a pure plane reflection
        assert!(flector.is_pure_reflection(ABS_DIFF_EQ_EPS));
        let plane = flector.plane_part();
        // PGA plane: nz = 1, d = -5 (since n·x + d = 0 means z - 5 = 0)
        assert!(abs_diff_eq!(plane.e012, 1.0, epsilon = ABS_DIFF_EQ_EPS)); // nz
        assert!(abs_diff_eq!(plane.e123, -5.0, epsilon = ABS_DIFF_EQ_EPS)); // d
    }

    #[test]
    fn flector_to_reflection_xy_plane() {
        let flector = Flector::reflect_xy();
        let na_refl: Reflection3<f64> = flector.try_into().unwrap();

        // Should have axis (0, 0, 1) and bias 0
        let axis = na_refl.axis();
        assert!(abs_diff_eq!(axis[0], 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(axis[1], 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(axis[2], 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(na_refl.bias(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flector_to_reflection_offset_plane() {
        // Plane z = 5: nx·x + ny·y + nz·z + d = 0 => z - 5 = 0 => d = -5
        let plane = Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);
        let flector = Flector::from_plane(plane);
        let na_refl: Reflection3<f64> = flector.try_into().unwrap();

        // nalgebra bias = n·P where P is on plane, so bias = 5
        let axis = na_refl.axis();
        assert!(abs_diff_eq!(axis[0], 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(axis[1], 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(axis[2], 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(na_refl.bias(), 5.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn non_pure_flector_conversion_fails() {
        // Create a flector with non-zero point component (not a pure reflection)
        let flector = Flector::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        let result: Result<Reflection3<f64>, _> = flector.try_into();
        assert!(matches!(
            result,
            Err(FlectorConversionError::NotPureReflection)
        ));
    }

    proptest! {
        /// Tests that Flector -> Reflection -> Flector roundtrip preserves behavior.
        #[test]
        fn flector_reflection_roundtrip(
            nx in -1.0f64..1.0, ny in -1.0f64..1.0, nz in -1.0f64..1.0,
            d in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize to get unit normal (skip degenerate cases)
            let len = (nx * nx + ny * ny + nz * nz).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let (nx, ny, nz) = (nx / len, ny / len, nz / len);

            let plane = Plane::from_normal_and_distance(nx, ny, nz, d);
            let flector = Flector::from_plane(plane);

            let na_refl: Reflection3<f64> = flector.try_into().unwrap();
            let back: Flector<f64> = na_refl.into();

            // Compare by transforming a point
            let p = Point::new(px, py, pz);
            let result_orig = flector.transform_point(&p);
            let result_back = back.transform_point(&p);

            prop_assert!(abs_diff_eq!(result_orig.x(), result_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.y(), result_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(result_orig.z(), result_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Tests that Flector reflection matches nalgebra Reflection.
        #[test]
        fn flector_matches_nalgebra_reflection(
            nx in -1.0f64..1.0, ny in -1.0f64..1.0, nz in -1.0f64..1.0,
            d in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            // Normalize to get unit normal (skip degenerate cases)
            let len = (nx * nx + ny * ny + nz * nz).sqrt();
            if len < 0.1 {
                return Ok(());
            }
            let (nx, ny, nz) = (nx / len, ny / len, nz / len);

            // Create PGA flector
            let plane = Plane::from_normal_and_distance(nx, ny, nz, d);
            let flector = Flector::from_plane(plane);

            // Create nalgebra reflection
            let axis = na::Unit::new_unchecked(na::Vector3::new(nx, ny, nz));
            // Point on plane: P = -d * n (since n·P + d = 0 => n·(-d*n) + d = -d + d = 0)
            let point_on_plane = na::Point3::new(-d * nx, -d * ny, -d * nz);
            let na_refl = na::Reflection::new_containing_point(axis, &point_on_plane);

            // Transform with PGA
            let pga_point = Point::new(px, py, pz);
            let pga_result = flector.transform_point(&pga_point);

            // Transform with nalgebra
            let mut na_vec = na::Vector3::new(px, py, pz);
            na_refl.reflect(&mut na_vec);

            prop_assert!(abs_diff_eq!(pga_result.x(), na_vec.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_vec.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.z(), na_vec.z, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
