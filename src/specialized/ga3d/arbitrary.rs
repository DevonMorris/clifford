//! Proptest `Arbitrary` implementations for 3D GA types.
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
//! | [`NonZeroVec3`] | Vec3 with non-zero magnitude |
//! | [`UnitVec3`] | Vec3 with unit magnitude |
//! | [`UnitBivec3`] | Bivec3 with unit magnitude |
//! | [`UnitRotor3`] | Rotor3 with unit magnitude |
//!
//! # Example
//!
//! ```
//! use clifford::specialized::ga3d::{Vec3, arbitrary::UnitRotor3};
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn rotation_preserves_norm(r in any::<UnitRotor3>(), v in any::<Vec3<f64>>()) {
//!         let rotated = r.0.rotate(v);
//!         prop_assert!((v.norm() - rotated.norm()).abs() < 1e-9);
//!     }
//! }
//! ```

use super::{Bivec3, Rotor3, Vec3};
use core::f64::consts::PI;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

impl Arbitrary for Vec3<f64> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0..100.0, -100.0..100.0, -100.0..100.0)
            .prop_map(|(x, y, z)| Vec3::new(x, y, z))
            .boxed()
    }
}

/// Wrapper type for non-zero [`Vec3`].
///
/// Use this when you need a vector guaranteed to have non-zero magnitude,
/// for example when testing normalization or division.
#[derive(Debug, Clone, Copy)]
pub struct NonZeroVec3(
    /// The wrapped non-zero vector.
    pub Vec3<f64>,
);

impl Arbitrary for NonZeroVec3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Vec3<f64>>()
            .prop_filter("non-zero vector", |v| v.norm_squared() > 1e-10)
            .prop_map(NonZeroVec3)
            .boxed()
    }
}

/// Wrapper type for unit [`Vec3`].
///
/// Use this when you need a vector guaranteed to have unit magnitude,
/// for example when testing rotations or reflections.
#[derive(Debug, Clone, Copy)]
pub struct UnitVec3(
    /// The wrapped unit vector.
    pub Vec3<f64>,
);

impl Arbitrary for UnitVec3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<NonZeroVec3>()
            .prop_map(|v| UnitVec3(v.0.normalized()))
            .boxed()
    }
}

impl Arbitrary for Bivec3<f64> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0..100.0, -100.0..100.0, -100.0..100.0)
            .prop_map(|(xy, xz, yz)| Bivec3::new(xy, xz, yz))
            .boxed()
    }
}

/// Wrapper type for unit [`Bivec3`].
///
/// Use this when you need a bivector guaranteed to have unit magnitude,
/// for example when constructing rotors from angle-plane pairs.
#[derive(Debug, Clone, Copy)]
pub struct UnitBivec3(
    /// The wrapped unit bivector.
    pub Bivec3<f64>,
);

impl Arbitrary for UnitBivec3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Bivec3<f64>>()
            .prop_filter("non-zero bivector", |b| b.norm_squared() > 1e-10)
            .prop_map(|b| UnitBivec3(b.normalized()))
            .boxed()
    }
}

/// Wrapper type for unit [`Rotor3`].
///
/// Use this when you need a rotor guaranteed to have unit magnitude,
/// which is required for proper rotation operations.
#[derive(Debug, Clone, Copy)]
pub struct UnitRotor3(
    /// The wrapped unit rotor.
    pub Rotor3<f64>,
);

impl Arbitrary for UnitRotor3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (0.0..2.0 * PI, any::<UnitBivec3>())
            .prop_map(|(angle, plane)| UnitRotor3(Rotor3::from_angle_plane(angle, plane.0)))
            .boxed()
    }
}
