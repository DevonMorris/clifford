//! nalgebra interoperability for 3D Projective GA types.
//!
//! This module provides bidirectional conversions between clifford's 3D PGA types
//! and nalgebra's equivalent types.
//!
//! Enable with feature `nalgebra-0_32`, `nalgebra-0_33`, or `nalgebra-0_34`.
//!
//! # Conversions
//!
//! | clifford | nalgebra | Notes |
//! |----------|----------|-------|
//! | [`Point<T>`] | [`na::Point3<T>`] | Homogeneous ↔ Cartesian coordinates |
//! | [`Motor<T>`] | [`na::Isometry3<T>`] | Rigid transformation (rotation + translation) |
//!
//! # Motor ↔ Isometry Correspondence
//!
//! Both PGA motors and nalgebra isometries represent rigid transformations
//! (rotation + translation), but with different internal representations.
//!
//! **Motor representation** (3D PGA with point-based formulation):
//! - Components `(e01, e02, e03, e0123)` encode the rotation (like a quaternion)
//! - Components `(s, e23, e31, e12)` encode the translation
//!
//! **Isometry representation**:
//! - `UnitQuaternion` for rotation
//! - `Translation3` for translation
//!
//! # Rotation Convention Difference
//!
//! **Important:** PGA rotation factories use the opposite direction convention from nalgebra.
//! - `Motor::from_rotation_z(θ)` rotates by `θ` in PGA convention
//! - This equals `na::UnitQuaternion::from_axis_angle(&z, -θ)` in nalgebra convention
//!
//! The Motor ↔ Isometry conversion handles this automatically: converting a motor to
//! an isometry and back preserves the transformation semantics perfectly. The convention
//! difference only matters when creating rotations directly from angle values.
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::projective::dim3::{Motor, Point};
//! use nalgebra as na;
//!
//! // Create a motor and convert to isometry
//! let motor = Motor::from_translation(1.0, 2.0, 3.0);
//! let iso: na::Isometry3<f64> = motor.into();
//!
//! // Both represent the same transformation
//! let p = Point::from_cartesian(0.0, 0.0, 0.0);
//! let na_p = na::Point3::new(0.0, 0.0, 0.0);
//!
//! let result_ga = motor.transform_point(&p);
//! let result_na = iso.transform_point(&na_p);
//!
//! // Results should match
//! assert!((result_ga.x() - result_na.x).abs() < 1e-10);
//! ```

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
// Point <-> Point3
// ============================================================================

impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T> {
    /// Converts a nalgebra 3D point to a PGA point.
    ///
    /// Creates a finite point with homogeneous weight `w = 1`.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::Point;
    /// use nalgebra::Point3;
    ///
    /// let na_p = Point3::new(1.0, 2.0, 3.0);
    /// let p: Point<f64> = na_p.into();
    /// assert_eq!(p.x(), 1.0);
    /// assert_eq!(p.y(), 2.0);
    /// assert_eq!(p.z(), 3.0);
    /// ```
    #[inline]
    fn from(p: na::Point3<T>) -> Self {
        Point::from_cartesian(p.x, p.y, p.z)
    }
}

impl<T: Float + na::Scalar> TryFrom<Point<T>> for na::Point3<T> {
    type Error = NalgebraConversionError;

    /// Attempts to convert a PGA point to a nalgebra 3D point.
    ///
    /// Returns an error if the point is ideal (at infinity).
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::Point;
    /// use nalgebra::Point3;
    ///
    /// let finite = Point::from_cartesian(1.0, 2.0, 3.0);
    /// let na_p: Point3<f64> = finite.try_into().unwrap();
    ///
    /// let ideal = Point::ideal(1.0, 0.0, 0.0);
    /// assert!(Point3::<f64>::try_from(ideal).is_err());
    /// ```
    fn try_from(p: Point<T>) -> Result<Self, Self::Error> {
        if p.e0().abs() < T::epsilon() {
            return Err(NalgebraConversionError::IdealPoint);
        }
        Ok(na::Point3::new(p.x(), p.y(), p.z()))
    }
}

// ============================================================================
// Motor <-> Isometry3
// ============================================================================

impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry3<T> {
    /// Converts a PGA motor to a nalgebra 3D isometry.
    ///
    /// # Mathematical Correspondence
    ///
    /// In the point-based PGA formulation used by this library:
    /// - Rotation is encoded in `(e01, e02, e03, e0123)` similar to a quaternion
    /// - Translation is encoded in `(s, e23, e31, e12)`
    ///
    /// The motor is unitized before conversion to ensure a valid isometry.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::Motor;
    /// use nalgebra::Isometry3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let motor = Motor::from_rotation_z(FRAC_PI_2);
    /// let iso: Isometry3<f64> = motor.into();
    ///
    /// // Both represent a 90° rotation around the z-axis
    /// ```
    fn from(motor: Motor<T>) -> Self {
        let m = motor.unitized();

        // Extract rotation quaternion from (e01, e02, e03, e0123)
        // PGA antisandwich uses antireverse which differs from quaternion conjugate.
        // Negate the bivector components to match nalgebra's rotation direction.
        // Mapping: (w, i, j, k) = (e0123, -e01, -e02, -e03)
        let rotation = na::UnitQuaternion::from_quaternion(na::Quaternion::new(
            m.e0123(),
            -m.e01(),
            -m.e02(),
            -m.e03(),
        ));

        // Extract translation by applying motor to origin and reading result
        // This handles the complex interaction between rotation and translation
        let origin = Point::origin();
        let transformed = m.transform_point(&origin);
        let translation = na::Translation3::new(transformed.x(), transformed.y(), transformed.z());

        na::Isometry3::from_parts(translation, rotation)
    }
}

impl<T: Float + na::RealField> From<na::Isometry3<T>> for Motor<T> {
    /// Converts a nalgebra 3D isometry to a PGA motor.
    ///
    /// # Mathematical Correspondence
    ///
    /// The motor is constructed by composing a rotation motor and translation motor.
    /// Since nalgebra applies rotation first then translation, we compose them as
    /// `rotation.compose(&translation)` which applies rotation first, then translation.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::{Motor, Point};
    /// use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};
    ///
    /// let iso = Isometry3::from_parts(
    ///     Translation3::new(1.0, 2.0, 3.0),
    ///     UnitQuaternion::identity(),
    /// );
    /// let motor: Motor<f64> = iso.into();
    ///
    /// // Both represent the same transformation
    /// ```
    fn from(iso: na::Isometry3<T>) -> Self {
        let q = iso.rotation.quaternion();
        let t = iso.translation.vector;

        // Build rotation motor from quaternion
        // PGA antisandwich uses antireverse which differs from quaternion conjugate.
        // Negate the bivector components to match nalgebra's rotation direction.
        // Mapping: (e01, e02, e03, e0123) = (-i, -j, -k, w)
        let rotation = Motor::new_unchecked(
            T::zero(), // s
            T::zero(), // e23
            T::zero(), // e31
            T::zero(), // e12
            -q.i,      // e01 (negated)
            -q.j,      // e02 (negated)
            -q.k,      // e03 (negated)
            q.w,       // e0123
        );

        // Build translation motor
        let translation = Motor::from_translation(t.x, t.y, t.z);

        // Compose: rotation first, then translation (matches nalgebra semantics)
        // In compose, self.compose(&other) applies other first, then self
        // So translation.compose(&rotation) applies rotation first, then translation
        translation.compose(&rotation)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::specialized::projective::dim3::{Flector, Plane, UnitizedPoint};
    use crate::test_utils::RELATIVE_EQ_EPS;
    use crate::wrappers::Unitized;
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
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let back: Point<f64> = na_p.into();

            // Compare Cartesian coordinates
            prop_assert!(relative_eq!(p.x(), back.x(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(p.y(), back.y(), epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(p.z(), back.z(), epsilon = EPS, max_relative = EPS));
        }

        #[test]
        fn point_try_from_finite(p in any::<UnitizedPoint<f64>>()) {
            let result: Result<na::Point3<f64>, _> = (*p).try_into();
            prop_assert!(result.is_ok());

            let na_p = result.unwrap();
            prop_assert!(relative_eq!(p.x(), na_p.x, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(p.y(), na_p.y, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(p.z(), na_p.z, epsilon = EPS, max_relative = EPS));
        }
    }

    #[test]
    fn point_try_from_ideal_fails() {
        let ideal = Point::<f64>::ideal(1.0, 0.0, 0.0);
        let result: Result<na::Point3<f64>, _> = ideal.try_into();
        assert!(matches!(result, Err(NalgebraConversionError::IdealPoint)));
    }

    // ========================================================================
    // Motor round-trip tests
    // ========================================================================

    proptest! {
        /// Motor round-trip through Isometry for arbitrary unitized motors.
        ///
        /// Tests that Motor → Isometry → Motor preserves the transformation.
        /// Any unitized motor represents a valid rigid transformation.
        #[test]
        fn motor_roundtrip(m in any::<Unitized<Motor<f64>>>()) {
            // Round-trip through Isometry
            let iso: na::Isometry3<f64> = (*m).into();
            let back: Motor<f64> = iso.into();

            // Test with several points to ensure transformation equivalence
            let test_points = [
                Point::from_cartesian(1.0, 0.0, 0.0),
                Point::from_cartesian(0.0, 1.0, 0.0),
                Point::from_cartesian(0.0, 0.0, 1.0),
                Point::from_cartesian(1.0, 2.0, 3.0),
            ];

            for test_p in test_points {
                let result1 = m.transform_point(&test_p);
                let result2 = back.transform_point(&test_p);

                prop_assert!(
                    relative_eq!(result1.x(), result2.x(), epsilon = EPS, max_relative = EPS),
                    "x mismatch: {} vs {}", result1.x(), result2.x()
                );
                prop_assert!(
                    relative_eq!(result1.y(), result2.y(), epsilon = EPS, max_relative = EPS),
                    "y mismatch: {} vs {}", result1.y(), result2.y()
                );
                prop_assert!(
                    relative_eq!(result1.z(), result2.z(), epsilon = EPS, max_relative = EPS),
                    "z mismatch: {} vs {}", result1.z(), result2.z()
                );
            }
        }
    }

    // ========================================================================
    // Operational equivalence tests
    // ========================================================================

    proptest! {
        /// Verify motor.transform_point() == isometry.transform_point()
        #[test]
        fn transform_point_equivalence(
            m in any::<Unitized<Motor<f64>>>(),
            p in any::<UnitizedPoint<f64>>(),
        ) {
            // Transform with clifford motor
            let result_ga = m.transform_point(&*p);

            // Transform with nalgebra isometry
            let iso: na::Isometry3<f64> = (*m).into();
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = iso.transform_point(&na_p);

            // Compare results
            prop_assert!(
                relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS),
                "x mismatch: GA={} vs NA={}", result_ga.x(), na_result.x
            );
            prop_assert!(
                relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS),
                "y mismatch: GA={} vs NA={}", result_ga.y(), na_result.y
            );
            prop_assert!(
                relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS),
                "z mismatch: GA={} vs NA={}", result_ga.z(), na_result.z
            );
        }

        /// Verify motor composition == isometry composition
        #[test]
        fn composition_equivalence(
            m1 in any::<Unitized<Motor<f64>>>(),
            m2 in any::<Unitized<Motor<f64>>>(),
            p in any::<UnitizedPoint<f64>>(),
        ) {
            // Compose with clifford
            let composed_ga = m1.compose(&*m2);
            let result_ga = composed_ga.transform_point(&*p);

            // Compose with nalgebra
            let iso1: na::Isometry3<f64> = (*m1).into();
            let iso2: na::Isometry3<f64> = (*m2).into();
            let composed_na = iso1 * iso2;
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = composed_na.transform_point(&na_p);

            // Compare results
            prop_assert!(
                relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS),
                "x mismatch: GA={} vs NA={}", result_ga.x(), na_result.x
            );
            prop_assert!(
                relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS),
                "y mismatch: GA={} vs NA={}", result_ga.y(), na_result.y
            );
            prop_assert!(
                relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS),
                "z mismatch: GA={} vs NA={}", result_ga.z(), na_result.z
            );
        }

        /// Verify motor inverse == isometry inverse
        #[test]
        fn inverse_equivalence(
            m in any::<Unitized<Motor<f64>>>(),
            p in any::<UnitizedPoint<f64>>(),
        ) {
            // Inverse with clifford
            let inv_ga = m.inverse();
            let result_ga = inv_ga.transform_point(&*p);

            // Inverse with nalgebra
            let iso: na::Isometry3<f64> = (*m).into();
            let inv_na = iso.inverse();
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = inv_na.transform_point(&na_p);

            // Compare results
            prop_assert!(
                relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS),
                "x mismatch: GA={} vs NA={}", result_ga.x(), na_result.x
            );
            prop_assert!(
                relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS),
                "y mismatch: GA={} vs NA={}", result_ga.y(), na_result.y
            );
            prop_assert!(
                relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS),
                "z mismatch: GA={} vs NA={}", result_ga.z(), na_result.z
            );
        }
    }

    // ========================================================================
    // Factory method equivalence tests
    // ========================================================================

    proptest! {
        /// Verify from_translation produces equivalent results
        #[test]
        fn translation_factory_equivalence(
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            tz in -100.0f64..100.0,
            p in any::<UnitizedPoint<f64>>(),
        ) {
            let motor = Motor::from_translation(tx, ty, tz);
            let iso = na::Isometry3::from_parts(
                na::Translation3::new(tx, ty, tz),
                na::UnitQuaternion::identity(),
            );

            let result_ga = motor.transform_point(&*p);
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = iso.transform_point(&na_p);

            prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
        }

        /// Verify from_rotation_x produces equivalent results.
        ///
        /// Note: PGA rotation convention is opposite to nalgebra quaternion convention.
        /// Motor::from_rotation_x(angle) equals nalgebra rotation by -angle.
        #[test]
        fn rotation_x_factory_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            p in any::<UnitizedPoint<f64>>(),
        ) {
            let motor = Motor::from_rotation_x(angle);
            // PGA rotates opposite to nalgebra convention, so use -angle
            let iso = na::Isometry3::from_parts(
                na::Translation3::identity(),
                na::UnitQuaternion::from_axis_angle(&na::Vector3::x_axis(), -angle),
            );

            let result_ga = motor.transform_point(&*p);
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = iso.transform_point(&na_p);

            prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
        }

        /// Verify from_rotation_y produces equivalent results.
        ///
        /// Note: PGA rotation convention is opposite to nalgebra quaternion convention.
        /// Motor::from_rotation_y(angle) equals nalgebra rotation by -angle.
        #[test]
        fn rotation_y_factory_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            p in any::<UnitizedPoint<f64>>(),
        ) {
            let motor = Motor::from_rotation_y(angle);
            // PGA rotates opposite to nalgebra convention, so use -angle
            let iso = na::Isometry3::from_parts(
                na::Translation3::identity(),
                na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), -angle),
            );

            let result_ga = motor.transform_point(&*p);
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = iso.transform_point(&na_p);

            prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
        }

        /// Verify from_rotation_z produces equivalent results.
        ///
        /// Note: PGA rotation convention is opposite to nalgebra quaternion convention.
        /// Motor::from_rotation_z(angle) equals nalgebra rotation by -angle.
        #[test]
        fn rotation_z_factory_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            p in any::<UnitizedPoint<f64>>(),
        ) {
            let motor = Motor::from_rotation_z(angle);
            // PGA rotates opposite to nalgebra convention, so use -angle
            let iso = na::Isometry3::from_parts(
                na::Translation3::identity(),
                na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), -angle),
            );

            let result_ga = motor.transform_point(&*p);
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = iso.transform_point(&na_p);

            prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
        }

        /// Verify from_axis_angle produces equivalent results.
        ///
        /// Note: PGA rotation convention is opposite to nalgebra quaternion convention.
        /// Motor::from_axis_angle(axis, angle) equals nalgebra rotation by -angle.
        #[test]
        fn axis_angle_factory_equivalence(
            ax in -1.0f64..1.0,
            ay in -1.0f64..1.0,
            az in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            p in any::<UnitizedPoint<f64>>(),
        ) {
            // Skip degenerate axes
            let axis_len = (ax*ax + ay*ay + az*az).sqrt();
            prop_assume!(axis_len > 0.1);

            let axis = crate::specialized::euclidean::dim3::Vector::new(ax, ay, az);
            let motor = Motor::from_axis_angle(&axis, angle);

            // PGA rotates opposite to nalgebra convention, so use -angle
            let na_axis = na::Unit::new_normalize(na::Vector3::new(ax, ay, az));
            let iso = na::Isometry3::from_parts(
                na::Translation3::identity(),
                na::UnitQuaternion::from_axis_angle(&na_axis, -angle),
            );

            let result_ga = motor.transform_point(&*p);
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_result = iso.transform_point(&na_p);

            prop_assert!(relative_eq!(result_ga.x(), na_result.x, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.y(), na_result.y, epsilon = EPS, max_relative = EPS));
            prop_assert!(relative_eq!(result_ga.z(), na_result.z, epsilon = EPS, max_relative = EPS));
        }
    }

    // ========================================================================
    // Identity and edge case tests
    // ========================================================================

    #[test]
    fn identity_motor_equivalence() {
        let motor = Motor::<f64>::identity();
        let iso: na::Isometry3<f64> = motor.into();

        // Identity isometry should have identity rotation and zero translation
        let p = na::Point3::new(1.0, 2.0, 3.0);
        let result = iso.transform_point(&p);
        assert!(relative_eq!(
            result.x,
            p.x,
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            result.y,
            p.y,
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            result.z,
            p.z,
            epsilon = EPS,
            max_relative = EPS
        ));
    }

    // ========================================================================
    // Flector (reflection) tests - nalgebra doesn't have a direct equivalent
    // ========================================================================

    /// Helper: compute reflection of a point through a plane through the origin.
    ///
    /// In PGA, the plane representation uses basis elements (e023, e031, e012) which have
    /// specific sign conventions relative to the standard normal vector (nx, ny, nz).
    /// The e031 basis element has opposite orientation from the standard e013 = e0∧e1∧e3,
    /// so the y-component of the reflection must be adjusted to match.
    ///
    /// This helper accounts for that convention difference to validate that our GA
    /// implementation produces geometrically correct reflections.
    fn reflect_point_through_origin_plane(
        point: &na::Point3<f64>,
        normal: &na::Vector3<f64>,
    ) -> na::Point3<f64> {
        // The GA stores plane normal with y-component sign flipped in e031
        // So from_plane_through_origin(nx, ny, nz) reflects as if normal is (nx, -ny, nz)
        let adjusted_normal = na::Vector3::new(normal.x, -normal.y, normal.z);
        let n_unit = adjusted_normal.normalize();
        let proj = point.coords.dot(&n_unit);
        na::Point3::from(point.coords - n_unit * (2.0 * proj))
    }

    proptest! {
        /// Verify flector reflection through origin planes matches manual formula.
        ///
        /// We test planes through the origin where the reflection formula is unambiguous:
        /// P' = P - 2 * (P · n_unit) * n_unit
        #[test]
        fn flector_reflection_through_origin(
            nx in -1.0f64..1.0,
            ny in -1.0f64..1.0,
            nz in -1.0f64..1.0,
            p in any::<UnitizedPoint<f64>>(),
        ) {
            // Skip degenerate normals
            let norm_len = (nx*nx + ny*ny + nz*nz).sqrt();
            prop_assume!(norm_len > 0.1);

            // Create a flector directly from plane through origin
            // This uses from_plane_through_origin which constructs the flector directly
            let flector = Flector::from_plane_through_origin(nx, ny, nz);

            // Reflect with flector
            let result_ga = flector.transform_point(&*p);

            // Reflect manually
            let na_p: na::Point3<f64> = (*p).try_into().unwrap();
            let na_normal = na::Vector3::new(nx, ny, nz);
            let result_na = reflect_point_through_origin_plane(&na_p, &na_normal);

            // Compare results
            prop_assert!(
                relative_eq!(result_ga.x(), result_na.x, epsilon = EPS, max_relative = EPS),
                "x mismatch: GA={} vs manual={}", result_ga.x(), result_na.x
            );
            prop_assert!(
                relative_eq!(result_ga.y(), result_na.y, epsilon = EPS, max_relative = EPS),
                "y mismatch: GA={} vs manual={}", result_ga.y(), result_na.y
            );
            prop_assert!(
                relative_eq!(result_ga.z(), result_na.z, epsilon = EPS, max_relative = EPS),
                "z mismatch: GA={} vs manual={}", result_ga.z(), result_na.z
            );
        }

        /// Verify that reflecting twice returns to the original point.
        #[test]
        fn flector_double_reflection_identity(
            nx in -1.0f64..1.0,
            ny in -1.0f64..1.0,
            nz in -1.0f64..1.0,
            d in -10.0f64..10.0,
            p in any::<UnitizedPoint<f64>>(),
        ) {
            // Skip degenerate normals
            let norm_len = (nx*nx + ny*ny + nz*nz).sqrt();
            prop_assume!(norm_len > 0.1);

            let plane = Plane::from_normal_and_distance(nx, ny, nz, d);
            let flector = Flector::from_plane(&plane);

            // Reflect twice
            let once = flector.transform_point(&*p);
            let twice = flector.transform_point(&once);

            // Should be back at original point
            prop_assert!(
                relative_eq!(p.x(), twice.x(), epsilon = EPS, max_relative = EPS),
                "x mismatch after double reflection: {} vs {}", p.x(), twice.x()
            );
            prop_assert!(
                relative_eq!(p.y(), twice.y(), epsilon = EPS, max_relative = EPS),
                "y mismatch after double reflection: {} vs {}", p.y(), twice.y()
            );
            prop_assert!(
                relative_eq!(p.z(), twice.z(), epsilon = EPS, max_relative = EPS),
                "z mismatch after double reflection: {} vs {}", p.z(), twice.z()
            );
        }
    }

    #[test]
    fn isometry_identity_to_motor() {
        let iso = na::Isometry3::<f64>::identity();
        let motor: Motor<f64> = iso.into();

        println!(
            "Converted motor: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            motor.s(),
            motor.e23(),
            motor.e31(),
            motor.e12(),
            motor.e01(),
            motor.e02(),
            motor.e03(),
            motor.e0123()
        );

        let p = Point::from_cartesian(1.0, 2.0, 3.0);
        let result = motor.transform_point(&p);

        println!(
            "Point ({}, {}, {}) -> ({}, {}, {})",
            p.x(),
            p.y(),
            p.z(),
            result.x(),
            result.y(),
            result.z()
        );

        assert!(relative_eq!(
            result.x(),
            p.x(),
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            result.y(),
            p.y(),
            epsilon = EPS,
            max_relative = EPS
        ));
        assert!(relative_eq!(
            result.z(),
            p.z(),
            epsilon = EPS,
            max_relative = EPS
        ));
    }
}
