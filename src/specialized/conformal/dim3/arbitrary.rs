//! Proptest `Arbitrary` implementations for 3D CGA types.
//!
//! This module provides [`Arbitrary`] implementations for property-based testing
//! with [`mod@proptest`]. It is available when either running tests or when the
//! `proptest-support` feature is enabled.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::conformal::dim3::Point;
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn point_is_always_null(p in any::<Point<f64>>()) {
//!         prop_assert!(p.is_null(1e-9));
//!     }
//! }
//! ```

use super::{Point, Sphere};
use crate::scalar::Float;
use core::fmt::Debug;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

// ============================================================================
// Point Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Point<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| Point::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}

// ============================================================================
// Sphere Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Sphere<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            0.01f64..100.0,
        )
            .prop_map(|(cx, cy, cz, r)| {
                Sphere::from_center_radius(
                    T::from_f64(cx),
                    T::from_f64(cy),
                    T::from_f64(cz),
                    T::from_f64(r),
                )
            })
            .boxed()
    }
}
