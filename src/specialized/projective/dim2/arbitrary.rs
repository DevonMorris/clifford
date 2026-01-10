//! Proptest `Arbitrary` implementations for 2D PGA specialized types.
//!
//! This module provides strategies for generating random [`Point`], [`Line`],
//! and [`Motor`] values for property-based testing.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::projective::dim2::{Point, Motor};
//! use clifford::specialized::projective::dim2::arbitrary::UnitMotor;
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
                Motor::new(
                    T::from_f64(s),
                    T::from_f64(e12),
                    T::from_f64(e20),
                    T::from_f64(e01),
                )
            })
            .boxed()
    }
}

/// Wrapper type for unit motors (valid rigid transformations) in 2D PGA.
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
        // Generate a valid rigid transformation: rotation + translation
        (
            -std::f64::consts::PI..std::f64::consts::PI, // angle
            -100.0f64..100.0,                            // tx
            -100.0f64..100.0,                            // ty
        )
            .prop_map(|(angle, tx, ty)| {
                let rotation = Motor::from_rotation(T::from_f64(angle));
                let translation = Motor::from_translation(T::from_f64(tx), T::from_f64(ty));
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
        fn unit_motor_has_unit_norm(motor in any::<UnitMotor<f64>>()) {
            // Unit motors should have norm ~1
            prop_assert!(abs_diff_eq!(motor.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        }

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
        fn motor_inverse_roundtrip(
            motor in any::<UnitMotor<f64>>(),
            point in any::<Point<f64>>(),
        ) {
            let transformed = motor.transform_point(&point);
            let back = motor.inverse().transform_point(&transformed);

            prop_assert!(abs_diff_eq!(back.x(), point.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(back.y(), point.y(), epsilon = ABS_DIFF_EQ_EPS));
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
        fn join_contains_both_points(
            p1 in any::<Point<f64>>(),
            p2 in any::<Point<f64>>(),
        ) {
            let line = p1.join(&p2);
            // Both points should have ~0 distance to the line
            let d1 = line.distance_to_point(&p1).abs();
            let d2 = line.distance_to_point(&p2).abs();
            prop_assert!(d1 < ABS_DIFF_EQ_EPS || line.is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(d2 < ABS_DIFF_EQ_EPS || line.is_zero(ABS_DIFF_EQ_EPS));
        }
    }
}
