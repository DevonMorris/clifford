//! Proptest `Arbitrary` wrapper types for 3D GA.
//!
//! This module provides wrapper types for constrained values in property-based testing.
//! The base types (Vector, Bivector, Rotor, etc.) have `Arbitrary` implementations
//! in the generated code.
//!
//! # Wrapper Types
//!
//! | Type | Description |
//! |------|-------------|
//! | [`UnitRotor<T>`] | Rotor with unit magnitude |
//! | [`NonZeroVector<T>`] | Vector with non-zero magnitude |
//! | [`UnitVector<T>`] | Vector with unit magnitude |
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim3::{Vector, arbitrary::UnitRotor};
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn unit_rotor_has_unit_norm(r in any::<UnitRotor<f64>>()) {
//!         prop_assert!((r.norm() - 1.0).abs() < 1e-10);
//!     }
//! }
//! ```

use super::{Rotor, Vector};
use crate::scalar::Float;
use core::fmt::Debug;
use core::ops::Deref;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

/// Minimum norm squared threshold for non-zero checks.
const MIN_NORM_SQUARED: f64 = 1e-6;

// ============================================================================
// UnitRotor
// ============================================================================

/// Wrapper type for unit [`Rotor`].
///
/// Use this when you need a rotor guaranteed to have unit magnitude,
/// which is required for valid rotations.
#[derive(Debug, Clone, Copy)]
pub struct UnitRotor<T: Float>(Rotor<T>);

impl<T: Float> UnitRotor<T> {
    /// Unwraps the inner rotor.
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
    fn from(wrapper: UnitRotor<T>) -> Self {
        wrapper.0
    }
}

impl<T> Arbitrary for UnitRotor<T>
where
    T: Float + Debug,
    Rotor<T>: Arbitrary + Debug,
    <Rotor<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Rotor<T>>()
            .prop_filter_map("rotor must be normalizable", |r| {
                let norm_sq = r.norm_squared();
                if norm_sq > T::from_f64(MIN_NORM_SQUARED) {
                    Some(UnitRotor(r.normalized()))
                } else {
                    None
                }
            })
            .boxed()
    }
}

// ============================================================================
// NonZeroVector
// ============================================================================

/// Wrapper type for non-zero [`Vector`].
///
/// Use this when you need a vector guaranteed to have non-zero magnitude.
#[derive(Debug, Clone, Copy)]
pub struct NonZeroVector<T: Float>(Vector<T>);

impl<T: Float> NonZeroVector<T> {
    /// Unwraps the inner vector.
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
    fn from(wrapper: NonZeroVector<T>) -> Self {
        wrapper.0
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
        any::<Vector<T>>()
            .prop_filter_map("vector must be non-zero", |v| {
                if v.norm_squared() > T::from_f64(MIN_NORM_SQUARED) {
                    Some(NonZeroVector(v))
                } else {
                    None
                }
            })
            .boxed()
    }
}

// ============================================================================
// UnitVector
// ============================================================================

/// Wrapper type for unit [`Vector`].
///
/// Use this when you need a vector guaranteed to have unit magnitude.
#[derive(Debug, Clone, Copy)]
pub struct UnitVector<T: Float>(Vector<T>);

impl<T: Float> UnitVector<T> {
    /// Unwraps the inner vector.
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
    fn from(wrapper: UnitVector<T>) -> Self {
        wrapper.0
    }
}

impl<T> Arbitrary for UnitVector<T>
where
    T: Float + Debug,
    Vector<T>: Arbitrary + Debug,
    <Vector<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Vector<T>>()
            .prop_filter_map("vector must be normalizable", |v| {
                let norm_sq = v.norm_squared();
                if norm_sq > T::from_f64(MIN_NORM_SQUARED) {
                    Some(UnitVector(v.normalized()))
                } else {
                    None
                }
            })
            .boxed()
    }
}
