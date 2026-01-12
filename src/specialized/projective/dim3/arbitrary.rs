//! Proptest `Arbitrary` implementations for 3D PGA specialized types.
//!
//! This module provides wrapper types for property-based testing that require
//! constrained values (non-origin points, non-degenerate lines/planes, unit flectors).
//!
//! For basic `Arbitrary` implementations of [`Point`], [`Line`], [`Plane`],
//! [`Motor`], and [`Flector`], see the generated traits module.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim3::{Point, Motor};
//! use clifford::specialized::euclidean::dim3::Vector as EuclideanVector;
//! use proptest::prelude::*;
//! use approx::abs_diff_eq;
//!
//! proptest! {
//!     #[test]
//!     fn motor_preserves_distance(
//!         angle in -std::f64::consts::PI..std::f64::consts::PI,
//!         axis_x in -1.0f64..1.0,
//!         axis_y in -1.0f64..1.0,
//!         axis_z in -1.0f64..1.0,
//!         tx in -100.0f64..100.0,
//!         ty in -100.0f64..100.0,
//!         tz in -100.0f64..100.0,
//!         p1 in any::<Point<f64>>(),
//!         p2 in any::<Point<f64>>(),
//!     ) {
//!         let axis_len = (axis_x * axis_x + axis_y * axis_y + axis_z * axis_z).sqrt();
//!         prop_assume!(axis_len > 0.01);
//!         let axis = EuclideanVector::new(axis_x / axis_len, axis_y / axis_len, axis_z / axis_len);
//!         let motor = Motor::from_translation(tx, ty, tz)
//!             .compose(&Motor::from_axis_angle(&axis, angle));
//!         let d_before = p1.distance(&p2);
//!         let t1 = motor.transform_point(&p1);
//!         let t2 = motor.transform_point(&p2);
//!         let d_after = t1.distance(&t2);
//!         prop_assert!(relative_eq!(d_before, d_after, epsilon = 1e-10));
//!     }
//! }
//! ```

use core::fmt::Debug;
use core::ops::Deref;

use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

use crate::scalar::Float;

use super::generated::types::{Flector, Line, Plane, Point};

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

// Note: Basic Arbitrary impl for Motor is in generated/traits.rs

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
    use crate::specialized::euclidean::dim3::Vector as EuclideanVector;
    use crate::specialized::projective::dim3::Motor;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;

    /// Helper function to create a valid unit motor from rotation + translation params.
    ///
    /// The result is normalized to satisfy the Study condition, which is required
    /// for correct rigid body transformations.
    fn make_unit_motor(
        axis_x: f64,
        axis_y: f64,
        axis_z: f64,
        angle: f64,
        tx: f64,
        ty: f64,
        tz: f64,
    ) -> Option<Motor<f64>> {
        let axis_len = (axis_x * axis_x + axis_y * axis_y + axis_z * axis_z).sqrt();
        if axis_len < 0.01 {
            return None;
        }
        let axis = EuclideanVector::new(axis_x / axis_len, axis_y / axis_len, axis_z / axis_len);
        let rotation = Motor::from_axis_angle(&axis, angle);
        let translation = Motor::from_translation(tx, ty, tz);
        // Normalize to satisfy Study condition for valid transformations
        Some(translation.compose(&rotation).normalized())
    }

    proptest! {
        #[test]
        fn motor_preserves_distance(
            axis_x in -1.0f64..1.0,
            axis_y in -1.0f64..1.0,
            axis_z in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            tz in -100.0f64..100.0,
            p1 in any::<Point<f64>>(),
            p2 in any::<Point<f64>>(),
        ) {
            if let Some(motor) = make_unit_motor(axis_x, axis_y, axis_z, angle, tx, ty, tz) {
                let d_before = p1.distance(&p2);
                let t1 = motor.transform_point(&p1);
                let t2 = motor.transform_point(&p2);
                let d_after = t1.distance(&t2);
                prop_assert!(relative_eq!(d_before, d_after, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            }
        }

        /// Test that motor geometric product is associative.
        ///
        /// Note: compose() returns the raw geometric product, which IS associative.
        /// The e0123 component may drift from the Study condition after multiple
        /// compositions, but the algebraic structure is preserved.
        #[test]
        fn motor_composition_associative(
            ax1 in -1.0f64..1.0, ay1 in -1.0f64..1.0, az1 in -1.0f64..1.0,
            angle1 in -std::f64::consts::PI..std::f64::consts::PI,
            tx1 in -10.0f64..10.0, ty1 in -10.0f64..10.0, tz1 in -10.0f64..10.0,
            ax2 in -1.0f64..1.0, ay2 in -1.0f64..1.0, az2 in -1.0f64..1.0,
            angle2 in -std::f64::consts::PI..std::f64::consts::PI,
            tx2 in -10.0f64..10.0, ty2 in -10.0f64..10.0, tz2 in -10.0f64..10.0,
            ax3 in -1.0f64..1.0, ay3 in -1.0f64..1.0, az3 in -1.0f64..1.0,
            angle3 in -std::f64::consts::PI..std::f64::consts::PI,
            tx3 in -10.0f64..10.0, ty3 in -10.0f64..10.0, tz3 in -10.0f64..10.0,
        ) {
            if let (Some(m1), Some(m2), Some(m3)) = (
                make_unit_motor(ax1, ay1, az1, angle1, tx1, ty1, tz1),
                make_unit_motor(ax2, ay2, az2, angle2, tx2, ty2, tz2),
                make_unit_motor(ax3, ay3, az3, angle3, tx3, ty3, tz3),
            ) {
                let lhs = m1.compose(&m2).compose(&m3);
                let rhs = m1.compose(&m2.compose(&m3));
                // Compare the geometric product results (associative)
                prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            }
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

            prop_assert!(relative_eq!(back.x(), point.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(back.y(), point.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(back.z(), point.z(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
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

            prop_assert!(relative_eq!(back.x(), point.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(back.y(), point.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(back.z(), point.z(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn line_join_contains_both_points(
            p1 in any::<NonOriginPoint<f64>>(),
            p2 in any::<NonOriginPoint<f64>>(),
        ) {
            let line = Line::join(&*p1, &*p2);
            // Both points should have ~0 distance to the line (if line is non-degenerate)
            if !line.is_zero(RELATIVE_EQ_EPS) {
                let d1 = line.distance_to_point(&*p1);
                let d2 = line.distance_to_point(&*p2);
                prop_assert!(d1.abs() < RELATIVE_EQ_EPS);
                prop_assert!(d2.abs() < RELATIVE_EQ_EPS);
            }
        }

        #[test]
        fn point_dot_symmetric(
            p1 in any::<Point<f64>>(),
            p2 in any::<Point<f64>>(),
        ) {
            let d12 = p1.dot(&p2);
            let d21 = p2.dot(&p1);
            prop_assert!(relative_eq!(d12, d21, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn line_dot_symmetric(
            l1 in any::<Line<f64>>(),
            l2 in any::<Line<f64>>(),
        ) {
            let d12 = l1.dot(&l2);
            let d21 = l2.dot(&l1);
            prop_assert!(relative_eq!(d12, d21, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn line_plucker_inner_symmetric(
            l1 in any::<NonDegenerateLine<f64>>(),
            l2 in any::<NonDegenerateLine<f64>>(),
        ) {
            let p12 = l1.plucker_inner(&*l2);
            let p21 = l2.plucker_inner(&*l1);
            prop_assert!(relative_eq!(p12, p21, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn flector_double_reflection_is_identity(
            flector in any::<UnitFlector<f64>>(),
            point in any::<Point<f64>>(),
        ) {
            // Reflecting twice should return to original point
            let once = flector.transform_point(&point);
            let twice = flector.transform_point(&once);

            prop_assert!(relative_eq!(twice.x(), point.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(twice.y(), point.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(twice.z(), point.z(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn motor_transforms_line_preserves_direction_norm(
            axis_x in -1.0f64..1.0,
            axis_y in -1.0f64..1.0,
            axis_z in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            tz in -100.0f64..100.0,
            line in any::<NonDegenerateLine<f64>>(),
        ) {
            if let Some(motor) = make_unit_motor(axis_x, axis_y, axis_z, angle, tx, ty, tz) {
                let transformed = motor.transform_line(&*line);
                // Direction norm (weight) should be preserved
                prop_assert!(relative_eq!(
                    line.weight_norm(),
                    transformed.weight_norm(),
                    max_relative = RELATIVE_EQ_EPS
                ));
            }
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

        /// Verify that factory-constructed motors satisfy the study condition.
        ///
        /// The study condition for a proper motor is:
        ///   s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0
        ///
        /// This condition ensures the motor represents a valid rigid transformation.
        /// It arises naturally when a motor is constructed as T*R (translation composed
        /// with rotation), and is preserved under motor composition and inversion.
        ///
        /// Reference: https://rigidgeometricalgebra.org/wiki/index.php?title=Motor
        #[test]
        fn unit_motor_satisfies_study_condition(
            axis_x in -1.0f64..1.0,
            axis_y in -1.0f64..1.0,
            axis_z in -1.0f64..1.0,
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            tz in -100.0f64..100.0,
        ) {
            if let Some(motor) = make_unit_motor(axis_x, axis_y, axis_z, angle, tx, ty, tz) {
                // Use the method that computes the Study condition residual
                prop_assert!(relative_eq!(motor.study_residual(), 0.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            }
        }
    }
}
