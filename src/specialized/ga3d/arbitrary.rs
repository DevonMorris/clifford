//! Proptest `Arbitrary` implementations for 3D GA types.
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
//! | [`NonZeroVec3<T>`] | Vec3 with non-zero magnitude |
//! | [`UnitVec3<T>`] | Vec3 with unit magnitude |
//! | [`UnitBivec3<T>`] | Bivec3 with unit magnitude |
//! | [`UnitRotor3<T>`] | Rotor3 with unit magnitude |
//!
//! # Example
//!
//! ```
//! use clifford::specialized::ga3d::{Vec3, arbitrary::UnitRotor3};
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn rotation_preserves_norm(r in any::<UnitRotor3<f64>>(), v in any::<Vec3<f64>>()) {
//!         let rotated = r.rotate(v);
//!         prop_assert!((v.norm() - rotated.norm()).abs() < 1e-9);
//!     }
//! }
//! ```

use super::{Bivec3, Even3, Rotor3, Trivec3, Vec3};
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
// Vec3 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Vec3<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| Vec3::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}

// ============================================================================
// NonZeroVec3
// ============================================================================

/// Wrapper type for non-zero [`Vec3`].
///
/// Use this when you need a vector guaranteed to have non-zero magnitude,
/// for example when testing normalization or division.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `Vec3<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct NonZeroVec3<T: Float>(
    /// The wrapped non-zero vector.
    pub Vec3<T>,
);

impl<T: Float> NonZeroVec3<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Vec3<T> {
        self.0
    }
}

impl<T: Float> Deref for NonZeroVec3<T> {
    type Target = Vec3<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Vec3<T>> for NonZeroVec3<T> {
    #[inline]
    fn as_ref(&self) -> &Vec3<T> {
        &self.0
    }
}

impl<T: Float> From<NonZeroVec3<T>> for Vec3<T> {
    #[inline]
    fn from(v: NonZeroVec3<T>) -> Self {
        v.0
    }
}

impl<T> Arbitrary for NonZeroVec3<T>
where
    T: Float + Debug,
    Vec3<T>: Arbitrary + Debug,
    <Vec3<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        let threshold = T::from_f64(MIN_NORM_SQUARED);
        any::<Vec3<T>>()
            .prop_filter("non-zero vector", move |v| v.norm_squared() > threshold)
            .prop_map(NonZeroVec3)
            .boxed()
    }
}

// ============================================================================
// UnitVec3
// ============================================================================

/// Wrapper type for unit [`Vec3`].
///
/// Use this when you need a vector guaranteed to have unit magnitude,
/// for example when testing rotations or reflections.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `NonZeroVec3<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct UnitVec3<T: Float>(
    /// The wrapped unit vector.
    pub Vec3<T>,
);

impl<T: Float> UnitVec3<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Vec3<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitVec3<T> {
    type Target = Vec3<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Vec3<T>> for UnitVec3<T> {
    #[inline]
    fn as_ref(&self) -> &Vec3<T> {
        &self.0
    }
}

impl<T: Float> From<UnitVec3<T>> for Vec3<T> {
    #[inline]
    fn from(v: UnitVec3<T>) -> Self {
        v.0
    }
}

impl<T> Arbitrary for UnitVec3<T>
where
    T: Float + Debug,
    NonZeroVec3<T>: Arbitrary,
    <NonZeroVec3<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<NonZeroVec3<T>>()
            .prop_map(|v| UnitVec3(v.0.normalized()))
            .boxed()
    }
}

// ============================================================================
// Bivec3 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Bivec3<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(xy, xz, yz)| Bivec3::new(T::from_f64(xy), T::from_f64(xz), T::from_f64(yz)))
            .boxed()
    }
}

// ============================================================================
// UnitBivec3
// ============================================================================

/// Wrapper type for unit [`Bivec3`].
///
/// Use this when you need a bivector guaranteed to have unit magnitude,
/// for example when constructing rotors from angle-plane pairs.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `Bivec3<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct UnitBivec3<T: Float>(
    /// The wrapped unit bivector.
    pub Bivec3<T>,
);

impl<T: Float> UnitBivec3<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Bivec3<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitBivec3<T> {
    type Target = Bivec3<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Bivec3<T>> for UnitBivec3<T> {
    #[inline]
    fn as_ref(&self) -> &Bivec3<T> {
        &self.0
    }
}

impl<T: Float> From<UnitBivec3<T>> for Bivec3<T> {
    #[inline]
    fn from(b: UnitBivec3<T>) -> Self {
        b.0
    }
}

impl<T> Arbitrary for UnitBivec3<T>
where
    T: Float + Debug,
    Bivec3<T>: Arbitrary + Debug,
    <Bivec3<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        let threshold = T::from_f64(MIN_NORM_SQUARED);
        any::<Bivec3<T>>()
            .prop_filter("non-zero bivector", move |b| b.norm_squared() > threshold)
            .prop_map(|b| UnitBivec3(b.normalized()))
            .boxed()
    }
}

// ============================================================================
// Trivec3 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Trivec3<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0)
            .prop_map(|val| Trivec3::new(T::from_f64(val)))
            .boxed()
    }
}

// ============================================================================
// Rotor3 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Rotor3<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
        )
            .prop_map(|(s, xy, xz, yz)| {
                Rotor3::new(
                    T::from_f64(s),
                    Bivec3::new(T::from_f64(xy), T::from_f64(xz), T::from_f64(yz)),
                )
            })
            .boxed()
    }
}

// ============================================================================
// Even3 Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Even3<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
        )
            .prop_map(|(s, xy, xz, yz)| {
                Even3::new(
                    T::from_f64(s),
                    Bivec3::new(T::from_f64(xy), T::from_f64(xz), T::from_f64(yz)),
                )
            })
            .boxed()
    }
}

// ============================================================================
// UnitRotor3
// ============================================================================

/// Wrapper type for unit [`Rotor3`].
///
/// Use this when you need a rotor guaranteed to have unit magnitude,
/// which is required for proper rotation operations.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float + FloatConst` where `UnitBivec3<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct UnitRotor3<T: Float>(
    /// The wrapped unit rotor.
    pub Rotor3<T>,
);

impl<T: Float> UnitRotor3<T> {
    /// Creates a new unit rotor wrapper.
    #[inline]
    fn new(rotor: Rotor3<T>) -> Self {
        Self(rotor)
    }

    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Rotor3<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitRotor3<T> {
    type Target = Rotor3<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Rotor3<T>> for UnitRotor3<T> {
    #[inline]
    fn as_ref(&self) -> &Rotor3<T> {
        &self.0
    }
}

impl<T: Float> From<UnitRotor3<T>> for Rotor3<T> {
    #[inline]
    fn from(r: UnitRotor3<T>) -> Self {
        r.0
    }
}

impl<T> Arbitrary for UnitRotor3<T>
where
    T: Float + Debug,
    UnitBivec3<T>: Arbitrary,
    <UnitBivec3<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate angle in [0, 2π) and a unit bivector for the rotation plane
        let two_pi = T::TWO * T::PI;
        any::<UnitBivec3<T>>()
            .prop_flat_map(move |plane| {
                // We need to generate the angle separately since T may not support ranges directly
                // Use f64 range and convert
                (Just(plane), 0.0f64..1.0)
            })
            .prop_map(move |(plane, t)| {
                let angle = T::from_f64(t) * two_pi;
                UnitRotor3::new(Rotor3::from_angle_plane(angle, plane.0))
            })
            .boxed()
    }
}
