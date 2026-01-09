//! Proptest `Arbitrary` implementations for generic [`Multivector`] types.
//!
//! This module provides [`Arbitrary`] implementations for property-based testing
//! with [`proptest`]. It is available when either running tests or when the
//! `proptest-support` feature is enabled.
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
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn vector_wedge_is_antisymmetric(
//!         a in any::<VectorE3>(),
//!         b in any::<VectorE3>()
//!     ) {
//!         let ab = a.0.outer(&b.0);
//!         let ba = b.0.outer(&a.0);
//!         prop_assert!(ab.approx_eq(&(-&ba), 1e-10));
//!     }
//! }
//! ```

use super::Multivector;
use crate::signature::Euclidean3;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

impl Arbitrary for Multivector<f64, Euclidean3> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        prop::array::uniform8(-10.0f64..10.0)
            .prop_map(|coeffs| Multivector::from_coeffs(&coeffs))
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

impl Arbitrary for UnitVectorE3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<NonZeroVectorE3>()
            .prop_map(|v| UnitVectorE3(v.0.normalize().unwrap()))
            .boxed()
    }
}
