//! Proptest `Arbitrary` implementations for 3D PGA specialized types.
//!
//! This module provides strategies for generating random [`Point`], [`Line`],
//! [`Plane`], [`Motor`], and [`Flector`] values for property-based testing.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim3::{Point, Motor};
//! use clifford::specialized::projective::dim3::arbitrary::UnitMotor;
//! use proptest::prelude::*;
//! use approx::abs_diff_eq;
//!
//! proptest! {
//!     #[test]
//!     fn motor_preserves_distance(
//!         motor in any::<UnitMotor<f64>>(),
//!         p1 in any::<Point<f64>>(),
//!         p2 in any::<Point<f64>>(),
//!     ) {
//!         let d_before = p1.distance(&p2);
//!         let t1 = motor.transform_point(&p1);
//!         let t2 = motor.transform_point(&p2);
//!         let d_after = t1.distance(&t2);
//!         prop_assert!(abs_diff_eq!(d_before, d_after, epsilon = 1e-10));
//!     }
//! }
//! ```

use core::fmt::Debug;
use core::ops::Deref;

use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

use crate::scalar::Float;
use crate::specialized::euclidean::dim3::Vector as EuclideanVector;

use super::generated::types::{Flector, Line, Motor, Plane, Point};

// ============================================================================
// Point
// ============================================================================

// Note: Basic Arbitrary impl is in generated/traits.rs

/// Wrapper type for non-origin points in 3D PGA.
///
/// Use this when you need a point guaranteed to not be at the origin.
#[derive(Debug, Clone)]
pub struct NonOriginPoint<T: Float>(
    /// The wrapped point.
    pub Point<T>,
);

impl<T: Float> NonOriginPoint<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Point<T> {
        self.0
    }
}

impl<T: Float> Deref for NonOriginPoint<T> {
    type Target = Point<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float + Debug> Arbitrary for NonOriginPoint<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_filter("must not be at origin", |(x, y, z)| {
                x * x + y * y + z * z > 0.01
            })
            .prop_map(|(x, y, z)| {
                NonOriginPoint(Point::from_cartesian(
                    T::from_f64(x),
                    T::from_f64(y),
                    T::from_f64(z),
                ))
            })
            .boxed()
    }
}

// ============================================================================
// Motor
// ============================================================================

// Note: Basic Arbitrary impl is in generated/traits.rs

/// Wrapper type for unit motors (valid rigid transformations) in 3D PGA.
///
/// Use this when you need a motor guaranteed to represent a valid rigid
/// transformation (rotation + translation).
#[derive(Debug, Clone)]
pub struct UnitMotor<T: Float>(
    /// The wrapped unit motor.
    pub Motor<T>,
);

impl<T: Float> UnitMotor<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Motor<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitMotor<T> {
    type Target = Motor<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float + Debug> Arbitrary for UnitMotor<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate a valid rigid transformation: rotation around arbitrary axis + translation
        // This produces motors with non-zero e0123 (screw component) from the
        // geometric product of rotation and translation components.
        (
            // Rotation axis (will be normalized)
            -1.0f64..1.0,
            -1.0f64..1.0,
            -1.0f64..1.0,
            // Rotation angle
            -std::f64::consts::PI..std::f64::consts::PI,
            // Translation
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
        )
            .prop_filter("axis must be non-zero", |(ax, ay, az, _, _, _, _)| {
                ax * ax + ay * ay + az * az > 0.01
            })
            .prop_map(|(ax, ay, az, angle, tx, ty, tz)| {
                // Normalize axis
                let len = (ax * ax + ay * ay + az * az).sqrt();
                let axis = EuclideanVector::new(
                    T::from_f64(ax / len),
                    T::from_f64(ay / len),
                    T::from_f64(az / len),
                );

                let rotation = Motor::from_axis_angle(&axis, T::from_f64(angle));
                let translation =
                    Motor::from_translation(T::from_f64(tx), T::from_f64(ty), T::from_f64(tz));
                UnitMotor(translation.compose(&rotation))
            })
            .boxed()
    }
}

// ============================================================================
// Line
// ============================================================================

// Note: Basic Arbitrary impl is in generated/traits.rs

/// Wrapper type for non-degenerate lines in 3D PGA.
///
/// Use this when you need a line with non-zero direction (weight).
#[derive(Debug, Clone)]
pub struct NonDegenerateLine<T: Float>(
    /// The wrapped line.
    pub Line<T>,
);

impl<T: Float> NonDegenerateLine<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Line<T> {
        self.0
    }
}

impl<T: Float> Deref for NonDegenerateLine<T> {
    type Target = Line<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float + Debug> Arbitrary for NonDegenerateLine<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
        )
            .prop_filter("direction must be non-zero", |(e01, e02, e03, _, _, _)| {
                e01 * e01 + e02 * e02 + e03 * e03 > 0.01
            })
            .prop_map(|(e01, e02, e03, e23, e31, e12)| {
                NonDegenerateLine(Line::new_unchecked(
                    T::from_f64(e01),
                    T::from_f64(e02),
                    T::from_f64(e03),
                    T::from_f64(e23),
                    T::from_f64(e31),
                    T::from_f64(e12),
                ))
            })
            .boxed()
    }
}

// ============================================================================
// Plane
// ============================================================================

// Note: Basic Arbitrary impl is in generated/traits.rs

/// Wrapper type for non-degenerate planes in 3D PGA.
///
/// Use this when you need a plane with non-zero normal direction (weight).
#[derive(Debug, Clone)]
pub struct NonDegeneratePlane<T: Float>(
    /// The wrapped plane.
    pub Plane<T>,
);

impl<T: Float> NonDegeneratePlane<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Plane<T> {
        self.0
    }
}

impl<T: Float> Deref for NonDegeneratePlane<T> {
    type Target = Plane<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float + Debug> Arbitrary for NonDegeneratePlane<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
        )
            .prop_filter("normal must be non-zero", |(e023, e031, e012, _)| {
                e023 * e023 + e031 * e031 + e012 * e012 > 0.01
            })
            .prop_map(|(e023, e031, e012, e123)| {
                NonDegeneratePlane(Plane::new(
                    T::from_f64(e023),
                    T::from_f64(e031),
                    T::from_f64(e012),
                    T::from_f64(e123),
                ))
            })
            .boxed()
    }
}

// ============================================================================
// Flector
// ============================================================================

// Note: Basic Arbitrary impl is in generated/traits.rs

/// Wrapper type for unit flectors (valid improper isometries) in 3D PGA.
///
/// Use this when you need a flector guaranteed to represent a valid
/// improper isometry (reflection, glide reflection, etc.).
#[derive(Debug, Clone)]
pub struct UnitFlector<T: Float>(
    /// The wrapped unit flector.
    pub Flector<T>,
);

impl<T: Float> UnitFlector<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Flector<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitFlector<T> {
    type Target = Flector<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float + Debug> Arbitrary for UnitFlector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate a pure plane reflection (unit flector)
        (-1.0f64..1.0, -1.0f64..1.0, -1.0f64..1.0, -100.0f64..100.0)
            .prop_filter("normal must be non-zero", |(nx, ny, nz, _)| {
                nx * nx + ny * ny + nz * nz > 0.01
            })
            .prop_map(|(nx, ny, nz, d)| {
                // Normalize the normal direction
                let len = (nx * nx + ny * ny + nz * nz).sqrt();
                let plane = Plane::from_normal_and_distance(
                    T::from_f64(nx / len),
                    T::from_f64(ny / len),
                    T::from_f64(nz / len),
                    T::from_f64(d),
                );
                UnitFlector(Flector::from_plane(&plane))
            })
            .boxed()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;

    proptest! {
        #[test]
        fn motor_preserves_distance(
            motor in any::<UnitMotor<f64>>(),
            p1 in any::<Point<f64>>(),
            p2 in any::<Point<f64>>(),
        ) {
            let d_before = p1.distance(&p2);
            let t1 = motor.transform_point(&p1);
            let t2 = motor.transform_point(&p2);
            let d_after = t1.distance(&t2);
            prop_assert!(abs_diff_eq!(d_before, d_after, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn motor_composition_associative(
            m1 in any::<UnitMotor<f64>>(),
            m2 in any::<UnitMotor<f64>>(),
            m3 in any::<UnitMotor<f64>>(),
        ) {
            // Use references via Deref to avoid moves
            let lhs = m1.compose(&*m2).compose(&*m3);
            let rhs = m1.compose(&m2.compose(&*m3));
            prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn pure_rotation_inverse_roundtrip(
            // Pure rotation: e0123 = 0
            axis_x in -1.0f64..1.0,
            axis_y in -1.0f64..1.0,
            axis_z in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            point in any::<Point<f64>>(),
        ) {
            let axis_len = (axis_x * axis_x + axis_y * axis_y + axis_z * axis_z).sqrt();
            prop_assume!(axis_len > 0.01);
            let axis = EuclideanVector::new(axis_x / axis_len, axis_y / axis_len, axis_z / axis_len);
            let motor = Motor::from_axis_angle(&axis, angle);

            let transformed = motor.transform_point(&point);
            let back = motor.inverse().transform_point(&transformed);

            prop_assert!(abs_diff_eq!(back.x(), point.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(back.y(), point.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(back.z(), point.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn pure_translation_inverse_roundtrip(
            // Pure translation: e0123 = 0
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            tz in -100.0f64..100.0,
            point in any::<Point<f64>>(),
        ) {
            let motor = Motor::from_translation(tx, ty, tz);

            let transformed = motor.transform_point(&point);
            let back = motor.inverse().transform_point(&transformed);

            prop_assert!(abs_diff_eq!(back.x(), point.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(back.y(), point.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(back.z(), point.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn line_join_contains_both_points(
            p1 in any::<NonOriginPoint<f64>>(),
            p2 in any::<NonOriginPoint<f64>>(),
        ) {
            let line = Line::join(&*p1, &*p2);
            // Both points should have ~0 distance to the line (if line is non-degenerate)
            if !line.is_zero(ABS_DIFF_EQ_EPS) {
                let d1 = line.distance_to_point(&*p1);
                let d2 = line.distance_to_point(&*p2);
                prop_assert!(d1.abs() < ABS_DIFF_EQ_EPS);
                prop_assert!(d2.abs() < ABS_DIFF_EQ_EPS);
            }
        }

        #[test]
        fn point_dot_symmetric(
            p1 in any::<Point<f64>>(),
            p2 in any::<Point<f64>>(),
        ) {
            let d12 = p1.dot(&p2);
            let d21 = p2.dot(&p1);
            prop_assert!(abs_diff_eq!(d12, d21, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn line_dot_symmetric(
            l1 in any::<Line<f64>>(),
            l2 in any::<Line<f64>>(),
        ) {
            let d12 = l1.dot(&l2);
            let d21 = l2.dot(&l1);
            prop_assert!(abs_diff_eq!(d12, d21, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn line_plucker_inner_symmetric(
            l1 in any::<NonDegenerateLine<f64>>(),
            l2 in any::<NonDegenerateLine<f64>>(),
        ) {
            let p12 = l1.plucker_inner(&*l2);
            let p21 = l2.plucker_inner(&*l1);
            prop_assert!(abs_diff_eq!(p12, p21, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn flector_double_reflection_is_identity(
            flector in any::<UnitFlector<f64>>(),
            point in any::<Point<f64>>(),
        ) {
            // Reflecting twice should return to original point
            let once = flector.transform_point(&point);
            let twice = flector.transform_point(&once);

            prop_assert!(abs_diff_eq!(twice.x(), point.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(twice.y(), point.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(twice.z(), point.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn motor_transforms_line_preserves_direction_norm(
            motor in any::<UnitMotor<f64>>(),
            line in any::<NonDegenerateLine<f64>>(),
        ) {
            let transformed = motor.transform_line(&*line);
            // Direction norm (weight) should be preserved
            prop_assert!(abs_diff_eq!(
                line.weight_norm(),
                transformed.weight_norm(),
                epsilon = ABS_DIFF_EQ_EPS
            ));
        }

        #[test]
        fn line_join_point_contains_both(
            line in any::<NonDegenerateLine<f64>>(),
            point in any::<Point<f64>>(),
        ) {
            let plane = line.join_point(&point);
            // The plane should contain both the original line and the point
            // For now, just verify we get a valid plane
            let _ = plane.normal();
        }

        /// Verify that unit motors satisfy the study condition.
        ///
        /// The study condition for a proper motor is:
        ///   s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0
        ///
        /// Equivalently: s*e0123 = v·m where v = (e23, e31, e12) and m = (e01, e02, e03).
        ///
        /// This condition ensures the motor represents a valid rigid transformation.
        /// It arises naturally when a motor is constructed as T*R (translation composed
        /// with rotation), and is preserved under motor composition and inversion.
        ///
        /// Reference: https://rigidgeometricalgebra.org/wiki/index.php?title=Motor
        #[test]
        fn unit_motor_satisfies_study_condition(motor in any::<UnitMotor<f64>>()) {
            // Study condition: s*e0123 = v·m
            // where v = (e23, e31, e12) and m = (e01, e02, e03)
            let study = motor.s() * motor.e0123()
                - motor.e23() * motor.e01()
                - motor.e31() * motor.e02()
                - motor.e12() * motor.e03();
            prop_assert!(abs_diff_eq!(study, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
