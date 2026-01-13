//! Proptest `Arbitrary` implementations for generic [`Multivector`] types.
//!
//! This module provides [`Arbitrary`] implementations for property-based testing
//! with [`mod@proptest`]. It is available when either running tests or when the
//! `proptest-support` feature is enabled.
//!
//! # Generic Implementation
//!
//! The `Arbitrary` trait is implemented generically for `Multivector<T, S>` where
//! `T: Float` and `S: Signature`. This works for any float type and any signature,
//! generating random coefficients for all basis blades.
//!
//! # Wrapper Types
//!
//! For constrained values (vectors, non-zero vectors, unit vectors), wrapper types
//! are provided for 3D Euclidean space:
//!
//! | Type | Description |
//! |------|-------------|
//! | [`VectorE3`] | Grade-1 multivector (vector) |
//! | [`NonZeroVectorE3`] | Non-zero magnitude vector |
//! | [`UnitVectorE3`] | Unit magnitude vector |
//!
//! # Example
//!
//! ```
//! use clifford::algebra::{Multivector, arbitrary::VectorE3};
//! use clifford::signature::Euclidean3;
//! use approx::abs_diff_eq;
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn vector_wedge_is_antisymmetric(
//!         a in any::<VectorE3>(),
//!         b in any::<VectorE3>()
//!     ) {
//!         let ab = a.exterior(&*b);
//!         let ba = b.exterior(&*a);
//!         prop_assert!(abs_diff_eq!(ab, -&ba, epsilon = 1e-10));
//!     }
//! }
//! ```

use super::Multivector;
use crate::scalar::Float;
use crate::signature::{Euclidean3, Signature};
use core::fmt::Debug;
use core::ops::Deref;
use generic_array::ArrayLength;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

impl<T, S> Arbitrary for Multivector<T, S>
where
    T: Float + Debug,
    S: Signature,
    S::NumBlades: ArrayLength,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        let num_blades = S::num_blades();
        prop::collection::vec(-10.0f64..10.0, num_blades)
            .prop_map(|coeffs| {
                let t_coeffs: Vec<T> = coeffs.into_iter().map(|c| T::from_f64(c)).collect();
                Multivector::from_coeffs(&t_coeffs)
            })
            .boxed()
    }
}

/// Wrapper type for vectors (grade-1 multivectors) in 3D Euclidean space.
///
/// Use this when you need an arbitrary vector for property testing.
#[derive(Debug, Clone)]
pub struct VectorE3(
    /// The wrapped vector multivector.
    pub Multivector<f64, Euclidean3>,
);

impl VectorE3 {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Multivector<f64, Euclidean3> {
        self.0
    }
}

impl Deref for VectorE3 {
    type Target = Multivector<f64, Euclidean3>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AsRef<Multivector<f64, Euclidean3>> for VectorE3 {
    #[inline]
    fn as_ref(&self) -> &Multivector<f64, Euclidean3> {
        &self.0
    }
}

impl From<VectorE3> for Multivector<f64, Euclidean3> {
    #[inline]
    fn from(v: VectorE3) -> Self {
        v.0
    }
}

impl Arbitrary for VectorE3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        prop::array::uniform3(-10.0f64..10.0)
            .prop_map(|v| VectorE3(Multivector::vector(&v)))
            .boxed()
    }
}

/// Wrapper type for non-zero vectors in 3D Euclidean space.
///
/// Use this when you need a vector guaranteed to have non-zero magnitude,
/// for example when testing normalization or inverse operations.
#[derive(Debug, Clone)]
pub struct NonZeroVectorE3(
    /// The wrapped non-zero vector multivector.
    pub Multivector<f64, Euclidean3>,
);

impl NonZeroVectorE3 {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Multivector<f64, Euclidean3> {
        self.0
    }
}

impl Deref for NonZeroVectorE3 {
    type Target = Multivector<f64, Euclidean3>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AsRef<Multivector<f64, Euclidean3>> for NonZeroVectorE3 {
    #[inline]
    fn as_ref(&self) -> &Multivector<f64, Euclidean3> {
        &self.0
    }
}

impl From<NonZeroVectorE3> for Multivector<f64, Euclidean3> {
    #[inline]
    fn from(v: NonZeroVectorE3) -> Self {
        v.0
    }
}

impl Arbitrary for NonZeroVectorE3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        prop::array::uniform3(-10.0f64..10.0)
            .prop_filter("must be non-zero", |v| {
                v[0] * v[0] + v[1] * v[1] + v[2] * v[2] > 0.1
            })
            .prop_map(|v| NonZeroVectorE3(Multivector::vector(&v)))
            .boxed()
    }
}

/// Wrapper type for unit vectors in 3D Euclidean space.
///
/// Use this when you need a vector guaranteed to have unit magnitude,
/// for example when testing reflections or rotation axes.
#[derive(Debug, Clone)]
pub struct UnitVectorE3(
    /// The wrapped unit vector multivector.
    pub Multivector<f64, Euclidean3>,
);

impl UnitVectorE3 {
    /// Unwraps and returns the inner value.
    #[inline]
    pub fn into_inner(self) -> Multivector<f64, Euclidean3> {
        self.0
    }
}

impl Deref for UnitVectorE3 {
    type Target = Multivector<f64, Euclidean3>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AsRef<Multivector<f64, Euclidean3>> for UnitVectorE3 {
    #[inline]
    fn as_ref(&self) -> &Multivector<f64, Euclidean3> {
        &self.0
    }
}

impl From<UnitVectorE3> for Multivector<f64, Euclidean3> {
    #[inline]
    fn from(v: UnitVectorE3) -> Self {
        v.0
    }
}

impl Arbitrary for UnitVectorE3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<NonZeroVectorE3>()
            .prop_map(|v| UnitVectorE3(v.0.normalize().unwrap()))
            .boxed()
    }
}
