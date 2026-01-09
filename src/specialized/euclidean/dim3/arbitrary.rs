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
//! | [`NonZeroVector<T>`] | Vector with non-zero magnitude |
//! | [`UnitVector<T>`] | Vector with unit magnitude |
//! | [`UnitBivector<T>`] | Bivector with unit magnitude |
//! | [`UnitRotor<T>`] | Rotor with unit magnitude |
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim3::{Vector, arbitrary::UnitRotor};
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn rotation_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
//!         let rotated = r.rotate(v);
//!         prop_assert!((v.norm() - rotated.norm()).abs() < 1e-9);
//!     }
//! }
//! ```

use super::{Bivector, Even, Rotor, Trivector, Vector};
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
// Vector Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Vector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| Vector::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}

// ============================================================================
// NonZeroVector
// ============================================================================

/// Wrapper type for non-zero [`Vector`].
///
/// Use this when you need a vector guaranteed to have non-zero magnitude,
/// for example when testing normalization or division.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `Vector<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct NonZeroVector<T: Float>(
    /// The wrapped non-zero vector.
    pub Vector<T>,
);

impl<T: Float> NonZeroVector<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Vector<T> {
        self.0
    }
}

impl<T: Float> Deref for NonZeroVector<T> {
    type Target = Vector<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Vector<T>> for NonZeroVector<T> {
    #[inline]
    fn as_ref(&self) -> &Vector<T> {
        &self.0
    }
}

impl<T: Float> From<NonZeroVector<T>> for Vector<T> {
    #[inline]
    fn from(v: NonZeroVector<T>) -> Self {
        v.0
    }
}

impl<T> Arbitrary for NonZeroVector<T>
where
    T: Float + Debug,
    Vector<T>: Arbitrary + Debug,
    <Vector<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        let threshold = T::from_f64(MIN_NORM_SQUARED);
        any::<Vector<T>>()
            .prop_filter("non-zero vector", move |v| v.norm_squared() > threshold)
            .prop_map(NonZeroVector)
            .boxed()
    }
}

// ============================================================================
// UnitVector
// ============================================================================

/// Wrapper type for unit [`Vector`].
///
/// Use this when you need a vector guaranteed to have unit magnitude,
/// for example when testing rotations or reflections.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `NonZeroVector<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct UnitVector<T: Float>(
    /// The wrapped unit vector.
    pub Vector<T>,
);

impl<T: Float> UnitVector<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Vector<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitVector<T> {
    type Target = Vector<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Vector<T>> for UnitVector<T> {
    #[inline]
    fn as_ref(&self) -> &Vector<T> {
        &self.0
    }
}

impl<T: Float> From<UnitVector<T>> for Vector<T> {
    #[inline]
    fn from(v: UnitVector<T>) -> Self {
        v.0
    }
}

impl<T> Arbitrary for UnitVector<T>
where
    T: Float + Debug,
    NonZeroVector<T>: Arbitrary,
    <NonZeroVector<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<NonZeroVector<T>>()
            .prop_map(|v| UnitVector(v.0.normalized()))
            .boxed()
    }
}

// ============================================================================
// Bivector Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Bivector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(xy, xz, yz)| {
                Bivector::new(T::from_f64(xy), T::from_f64(xz), T::from_f64(yz))
            })
            .boxed()
    }
}

// ============================================================================
// UnitBivector
// ============================================================================

/// Wrapper type for unit [`Bivector`].
///
/// Use this when you need a bivector guaranteed to have unit magnitude,
/// for example when constructing rotors from angle-plane pairs.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float` where `Bivector<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct UnitBivector<T: Float>(
    /// The wrapped unit bivector.
    pub Bivector<T>,
);

impl<T: Float> UnitBivector<T> {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Bivector<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitBivector<T> {
    type Target = Bivector<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Bivector<T>> for UnitBivector<T> {
    #[inline]
    fn as_ref(&self) -> &Bivector<T> {
        &self.0
    }
}

impl<T: Float> From<UnitBivector<T>> for Bivector<T> {
    #[inline]
    fn from(b: UnitBivector<T>) -> Self {
        b.0
    }
}

impl<T> Arbitrary for UnitBivector<T>
where
    T: Float + Debug,
    Bivector<T>: Arbitrary + Debug,
    <Bivector<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        let threshold = T::from_f64(MIN_NORM_SQUARED);
        any::<Bivector<T>>()
            .prop_filter("non-zero bivector", move |b| b.norm_squared() > threshold)
            .prop_map(|b| UnitBivector(b.normalized()))
            .boxed()
    }
}

// ============================================================================
// Trivector Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Trivector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0)
            .prop_map(|val| Trivector::new(T::from_f64(val)))
            .boxed()
    }
}

// ============================================================================
// Rotor Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Rotor<T> {
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
                Rotor::new(
                    T::from_f64(s),
                    Bivector::new(T::from_f64(xy), T::from_f64(xz), T::from_f64(yz)),
                )
            })
            .boxed()
    }
}

// ============================================================================
// Even Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Even<T> {
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
                Even::new(
                    T::from_f64(s),
                    Bivector::new(T::from_f64(xy), T::from_f64(xz), T::from_f64(yz)),
                )
            })
            .boxed()
    }
}

// ============================================================================
// UnitRotor
// ============================================================================

/// Wrapper type for unit [`Rotor`].
///
/// Use this when you need a rotor guaranteed to have unit magnitude,
/// which is required for proper rotation operations.
///
/// Generic over the scalar type `T`. The `Arbitrary` implementation is available
/// for any `T: Float + FloatConst` where `UnitBivector<T>: Arbitrary`.
#[derive(Debug, Clone, Copy)]
pub struct UnitRotor<T: Float>(
    /// The wrapped unit rotor.
    pub Rotor<T>,
);

impl<T: Float> UnitRotor<T> {
    /// Creates a new unit rotor wrapper.
    #[inline]
    fn new(rotor: Rotor<T>) -> Self {
        Self(rotor)
    }

    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Rotor<T> {
        self.0
    }
}

impl<T: Float> Deref for UnitRotor<T> {
    type Target = Rotor<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Rotor<T>> for UnitRotor<T> {
    #[inline]
    fn as_ref(&self) -> &Rotor<T> {
        &self.0
    }
}

impl<T: Float> From<UnitRotor<T>> for Rotor<T> {
    #[inline]
    fn from(r: UnitRotor<T>) -> Self {
        r.0
    }
}

impl<T> Arbitrary for UnitRotor<T>
where
    T: Float + Debug,
    UnitBivector<T>: Arbitrary,
    <UnitBivector<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate angle in [0, 2π) and a unit bivector for the rotation plane
        let two_pi = T::TWO * T::PI;
        any::<UnitBivector<T>>()
            .prop_flat_map(move |plane| {
                // We need to generate the angle separately since T may not support ranges directly
                // Use f64 range and convert
                (Just(plane), 0.0f64..1.0)
            })
            .prop_map(move |(plane, t)| {
                let angle = T::from_f64(t) * two_pi;
                UnitRotor::new(Rotor::from_angle_plane(angle, plane.0))
            })
            .boxed()
    }
}
