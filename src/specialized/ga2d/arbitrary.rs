//! Proptest `Arbitrary` implementations for 2D GA types.
//!
//! This module provides [`Arbitrary`] implementations for property-based testing
//! with [`proptest`]. It is available when either running tests or when the
//! `proptest-support` feature is enabled.
//!
//! # Wrapper Types
//!
//! For constrained values (non-zero, unit length), wrapper types are provided:
//!
//! | Type | Description |
//! |------|-------------|
//! | [`NonZeroVec2`] | Vec2 with non-zero magnitude |
//! | [`UnitVec2`] | Vec2 with unit magnitude |
//! | [`UnitRotor2`] | Rotor2 with unit magnitude |
//!
//! # Example
//!
//! ```
//! use clifford::specialized::ga2d::{Vec2, arbitrary::UnitVec2};
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn unit_vec_has_unit_length(v in any::<UnitVec2>()) {
//!         prop_assert!((v.0.norm() - 1.0).abs() < 1e-10);
//!     }
//! }
//! ```

use super::{Rotor2, Vec2};
use core::f64::consts::PI;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

impl Arbitrary for Vec2<f64> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0..100.0, -100.0..100.0)
            .prop_map(|(x, y)| Vec2::new(x, y))
            .boxed()
    }
}

/// Wrapper type for non-zero [`Vec2`].
///
/// Use this when you need a vector guaranteed to have non-zero magnitude,
/// for example when testing normalization or division.
#[derive(Debug, Clone, Copy)]
pub struct NonZeroVec2(
    /// The wrapped non-zero vector.
    pub Vec2<f64>,
);

impl Arbitrary for NonZeroVec2 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Vec2<f64>>()
            .prop_filter("non-zero vector", |v| v.norm_squared() > 1e-10)
            .prop_map(NonZeroVec2)
            .boxed()
    }
}

/// Wrapper type for unit [`Vec2`].
///
/// Use this when you need a vector guaranteed to have unit magnitude,
/// for example when testing rotations or reflections.
#[derive(Debug, Clone, Copy)]
pub struct UnitVec2(
    /// The wrapped unit vector.
    pub Vec2<f64>,
);

impl Arbitrary for UnitVec2 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<NonZeroVec2>()
            .prop_map(|v| UnitVec2(v.0.normalized()))
            .boxed()
    }
}

/// Wrapper type for unit [`Rotor2`].
///
/// Use this when you need a rotor guaranteed to have unit magnitude,
/// which is required for proper rotation operations.
#[derive(Debug, Clone, Copy)]
pub struct UnitRotor2(
    /// The wrapped unit rotor.
    pub Rotor2<f64>,
);

impl Arbitrary for UnitRotor2 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (0.0..2.0 * PI)
            .prop_map(|angle| UnitRotor2(Rotor2::from_angle(angle)))
            .boxed()
    }
}
