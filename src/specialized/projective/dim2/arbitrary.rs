//! Proptest `Arbitrary` implementations for 2D PGA specialized types.
//!
//! This module provides wrapper types for property-based testing that require
//! constrained values (non-origin points, non-degenerate lines).
//!
//! For basic `Arbitrary` implementations of [`Point`], [`Line`], and [`Motor`],
//! see the types module.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim2::{Point, Motor};
//! use proptest::prelude::*;
//! use approx::abs_diff_eq;
//!
//! proptest! {
//!     #[test]
//!     fn motor_preserves_distance(
//!         angle in -std::f64::consts::PI..std::f64::consts::PI,
//!         tx in -100.0f64..100.0,
//!         ty in -100.0f64..100.0,
//!         p1 in any::<Point<f64>>(),
//!         p2 in any::<Point<f64>>(),
//!     ) {
//!         let motor = Motor::from_translation(tx, ty)
//!             .compose(&Motor::from_rotation(angle));
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

use super::types::{Line, Motor, Point};

// ============================================================================
// Point
// ============================================================================

impl<T: Float + Debug> Arbitrary for Point<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate finite points with normalized w=1
        (-100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y)| Point::new(T::from_f64(x), T::from_f64(y)))
            .boxed()
    }
}

/// Wrapper type for non-origin points in 2D PGA.
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
        (-100.0f64..100.0, -100.0f64..100.0)
            .prop_filter("must not be at origin", |(x, y)| x * x + y * y > 0.01)
            .prop_map(|(x, y)| NonOriginPoint(Point::new(T::from_f64(x), T::from_f64(y))))
            .boxed()
    }
}

// ============================================================================
// Line
// ============================================================================

impl<T: Float + Debug> Arbitrary for Line<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0)
            .prop_map(|(e12, e20, e01)| {
                Line::new(T::from_f64(e12), T::from_f64(e20), T::from_f64(e01))
            })
            .boxed()
    }
}

/// Wrapper type for non-degenerate lines in 2D PGA.
///
/// Use this when you need a line with non-zero normal direction.
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
        (-10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0)
            .prop_filter("normal must be non-zero", |(_, e20, e01)| {
                e20 * e20 + e01 * e01 > 0.01
            })
            .prop_map(|(e12, e20, e01)| {
                NonDegenerateLine(Line::new(
                    T::from_f64(e12),
                    T::from_f64(e20),
                    T::from_f64(e01),
                ))
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
        )
            .prop_map(|(s, e12, e20, e01)| {
                Motor::new_unchecked(
                    T::from_f64(s),
                    T::from_f64(e12),
                    T::from_f64(e20),
                    T::from_f64(e01),
                )
            })
            .boxed()
    }
}

// Note: UnitMotor wrapper type was removed. Use inline motor construction
// with Motor::from_rotation() and Motor::from_translation() in tests.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;

    /// Helper function to create a valid unit motor from rotation + translation.
    fn make_unit_motor(angle: f64, tx: f64, ty: f64) -> Motor<f64> {
        let rotation = Motor::from_rotation(angle);
        let translation = Motor::from_translation(tx, ty);
        translation.compose(&rotation)
    }

    proptest! {
        #[test]
        fn unit_motor_has_unit_norm(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
        ) {
            let motor = make_unit_motor(angle, tx, ty);
            // Unit motors should have norm ~1
            prop_assert!(relative_eq!(motor.norm(), 1.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn motor_preserves_distance(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            p1 in any::<Point<f64>>(),
            p2 in any::<Point<f64>>(),
        ) {
            let motor = make_unit_motor(angle, tx, ty);
            let d_before = p1.distance(&p2);
            let t1 = motor.transform_point(&p1);
            let t2 = motor.transform_point(&p2);
            let d_after = t1.distance(&t2);
            prop_assert!(relative_eq!(d_before, d_after, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn motor_inverse_roundtrip(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -100.0f64..100.0,
            ty in -100.0f64..100.0,
            point in any::<Point<f64>>(),
        ) {
            let motor = make_unit_motor(angle, tx, ty);
            let transformed = motor.transform_point(&point);
            let back = motor.inverse().transform_point(&transformed);

            prop_assert!(relative_eq!(back.x(), point.x(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(back.y(), point.y(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn motor_composition_associative(
            angle1 in -std::f64::consts::PI..std::f64::consts::PI,
            tx1 in -10.0f64..10.0, ty1 in -10.0f64..10.0,
            angle2 in -std::f64::consts::PI..std::f64::consts::PI,
            tx2 in -10.0f64..10.0, ty2 in -10.0f64..10.0,
            angle3 in -std::f64::consts::PI..std::f64::consts::PI,
            tx3 in -10.0f64..10.0, ty3 in -10.0f64..10.0,
        ) {
            let m1 = make_unit_motor(angle1, tx1, ty1);
            let m2 = make_unit_motor(angle2, tx2, ty2);
            let m3 = make_unit_motor(angle3, tx3, ty3);
            let lhs = m1.compose(&m2).compose(&m3);
            let rhs = m1.compose(&m2.compose(&m3));
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn join_contains_both_points(
            p1 in any::<Point<f64>>(),
            p2 in any::<Point<f64>>(),
        ) {
            let line = p1.join(&p2);
            // Both points should have ~0 distance to the line
            let d1 = line.distance_to_point(&p1).abs();
            let d2 = line.distance_to_point(&p2).abs();
            prop_assert!(d1 < RELATIVE_EQ_EPS || line.is_zero(RELATIVE_EQ_EPS));
            prop_assert!(d2 < RELATIVE_EQ_EPS || line.is_zero(RELATIVE_EQ_EPS));
        }
    }
}
