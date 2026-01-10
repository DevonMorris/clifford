//! Proptest `Arbitrary` implementations for 3D PGA specialized types.
//!
//! This module provides strategies for generating random [`Point`] and
//! [`Motor`] values for property-based testing.
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

use super::types::{Motor, Point};

// ============================================================================
// Point
// ============================================================================

impl<T: Float + Debug> Arbitrary for Point<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate finite points with normalized w=1
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| Point::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}

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
                NonOriginPoint(Point::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            })
            .boxed()
    }
}

// ============================================================================
// Motor
// ============================================================================

impl<T: Float + Debug> Arbitrary for Motor<T> {
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
            -10.0f64..10.0,
            -10.0f64..10.0,
        )
            .prop_map(|(s, e23, e31, e12, e01, e02, e03, e0123)| {
                Motor::new(
                    T::from_f64(s),
                    T::from_f64(e23),
                    T::from_f64(e31),
                    T::from_f64(e12),
                    T::from_f64(e01),
                    T::from_f64(e02),
                    T::from_f64(e03),
                    T::from_f64(e0123),
                )
            })
            .boxed()
    }
}

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
    }
}
