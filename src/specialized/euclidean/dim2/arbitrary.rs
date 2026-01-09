//! Proptest `Arbitrary` implementations for 2D GA types.
//!
//! This module provides [`Arbitrary`] implementations for property-based testing
//! with [`proptest`]. It is available when either running tests or when the
//! `proptest-support` feature is enabled.
//!
//! # Generic Over Float Type
//!
//! The `Arbitrary` trait is implemented for both `f32` and `f64` variants of all types.
//! Wrapper types are generic over `T: Float`, with blanket `Arbitrary` implementations
//! that work for any `T` where the underlying type implements `Arbitrary`.
//!
//! # Wrapper Types
//!
//! For constrained values (non-zero, unit length), wrapper types are provided:
//!
//! | Type | Description |
//! |------|-------------|
//! | [`NonZeroVec2<T>`] | Vec2 with non-zero magnitude |
//! | [`UnitVec2<T>`] | Vec2 with unit magnitude |
//! | [`UnitRotor2<T>`] | Rotor2 with unit magnitude |
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim2::{Vec2, arbitrary::UnitVec2};
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn unit_vec_has_unit_length(v in any::<UnitVec2<f64>>()) {
//!         prop_assert!((v.norm() - 1.0).abs() < 1e-10);
//!     }
//! }
//! ```

use super::{Bivec2, Rotor2, Vec2};
use crate::scalar::Float;
use core::fmt::Debug;
use core::ops::Deref;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

/// Minimum norm squared threshold for non-zero checks.
/// Using 1e-6 works for both f32 and f64.
const MIN_NORM_SQUARED: f64 = 1e-6;

// ============================================================================
// Vec2 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Vec2<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y)| Vec2::new(T::from_f64(x), T::from_f64(y)))
            .boxed()
    }
}

// ============================================================================
// Bivec2 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Bivec2<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0)
            .prop_map(|val| Bivec2::new(T::from_f64(val)))
            .boxed()
    }
}

// ============================================================================
// Rotor2 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Rotor2<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(s, xy)| Rotor2::new(T::from_f64(s), T::from_f64(xy)))
            .boxed()
    }
}

// ============================================================================
// NonZeroVec2
// ============================================================================

/// Wrapper type for non-zero [`Vec2`].
///
/// Use this when you need a vector guaranteed to have non-zero magnitude,
/// for example when testing normalization or division.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `Vec2<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct NonZeroVec2<T: Float>(
    /// The wrapped non-zero vector.
    pub Vec2<T>,
);

impl<T: Float> NonZeroVec2<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Vec2<T> {
        self.0
    }
}

impl<T: Float> Deref for NonZeroVec2<T> {
    type Target = Vec2<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Vec2<T>> for NonZeroVec2<T> {
    #[inline]
    fn as_ref(&self) -> &Vec2<T> {
        &self.0
    }
}

impl<T: Float> From<NonZeroVec2<T>> for Vec2<T> {
    #[inline]
    fn from(v: NonZeroVec2<T>) -> Self {
        v.0
    }
}

impl<T> Arbitrary for NonZeroVec2<T>
where
    T: Float + Debug,
    Vec2<T>: Arbitrary + Debug,
    <Vec2<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        let threshold = T::from_f64(MIN_NORM_SQUARED);
        any::<Vec2<T>>()
            .prop_filter("non-zero vector", move |v| v.norm_squared() > threshold)
            .prop_map(NonZeroVec2)
            .boxed()
    }
}

// ============================================================================
// UnitVec2
// ============================================================================

/// Wrapper type for unit [`Vec2`].
///
/// Use this when you need a vector guaranteed to have unit magnitude,
/// for example when testing rotations or reflections.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `NonZeroVec2<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct UnitVec2<T: Float>(
    /// The wrapped unit vector.
    pub Vec2<T>,
);

impl<T: Float> UnitVec2<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Vec2<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitVec2<T> {
    type Target = Vec2<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Vec2<T>> for UnitVec2<T> {
    #[inline]
    fn as_ref(&self) -> &Vec2<T> {
        &self.0
    }
}

impl<T: Float> From<UnitVec2<T>> for Vec2<T> {
    #[inline]
    fn from(v: UnitVec2<T>) -> Self {
        v.0
    }
}

impl<T> Arbitrary for UnitVec2<T>
where
    T: Float + Debug,
    NonZeroVec2<T>: Arbitrary,
    <NonZeroVec2<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<NonZeroVec2<T>>()
            .prop_map(|v| UnitVec2(v.0.normalized()))
            .boxed()
    }
}

// ============================================================================
// UnitRotor2
// ============================================================================

/// Wrapper type for unit [`Rotor2`].
///
/// Use this when you need a rotor guaranteed to have unit magnitude,
/// which is required for proper rotation operations.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float`.
#[derive(Debug, Clone, Copy)]
pub struct UnitRotor2<T: Float>(
    /// The wrapped unit rotor.
    pub Rotor2<T>,
);

impl<T: Float> UnitRotor2<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Rotor2<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitRotor2<T> {
    type Target = Rotor2<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Rotor2<T>> for UnitRotor2<T> {
    #[inline]
    fn as_ref(&self) -> &Rotor2<T> {
        &self.0
    }
}

impl<T: Float> From<UnitRotor2<T>> for Rotor2<T> {
    #[inline]
    fn from(r: UnitRotor2<T>) -> Self {
        r.0
    }
}

impl<T> Arbitrary for UnitRotor2<T>
where
    T: Float + Debug,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate angle in [0, 2π) by using f64 range and converting
        let two_pi = T::TWO * T::PI;
        (0.0f64..1.0)
            .prop_map(move |t| {
                let angle = T::from_f64(t) * two_pi;
                UnitRotor2(Rotor2::from_angle(angle))
            })
            .boxed()
    }
}
