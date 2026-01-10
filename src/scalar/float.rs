//! Floating-point scalar type abstraction.
//!
//! This module provides the [`Float`] trait which abstracts over floating-point
//! types like `f32` and `f64`, allowing the library to be generic over precision.

use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use core::fmt::{Debug, Display};
use core::iter::{Product, Sum};
use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

/// Trait for floating-point scalar types used in geometric algebra computations.
///
/// This trait abstracts over `f32` and `f64` (and potentially other float types),
/// providing the necessary operations for geometric algebra while allowing users
/// to choose their desired precision.
///
/// # Mathematical Context
///
/// In geometric algebra, scalars form the grade-0 elements of the algebra.
/// All higher-grade elements (vectors, bivectors, etc.) have coefficients
/// that are scalars. This trait ensures those coefficients support the
/// arithmetic operations needed for geometric products and other operations.
///
/// # Design
///
/// This trait extends [`num_traits::Float`] to inherit standard floating-point
/// operations (`abs`, `sqrt`, `sin`, `cos`, etc.) and adds a few additional
/// requirements specific to geometric algebra:
/// - [`approx`] traits for floating-point comparisons
/// - [`Sum`] and [`Product`] for iterator operations
/// - Constants like `TWO` and `PI` commonly used in GA
///
/// # Example
///
/// ```
/// use clifford::scalar::Float;
///
/// fn compute_norm_squared<T: Float>(x: T, y: T, z: T) -> T {
///     x * x + y * y + z * z
/// }
///
/// let norm_sq: f64 = compute_norm_squared(1.0, 2.0, 3.0);
/// assert!((norm_sq - 14.0).abs() < f64::EPSILON);
/// ```
pub trait Float:
    num_traits::Float
    + Debug
    + Display
    + Default
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + Sum
    + Product
    + AbsDiffEq<Epsilon = Self>
    + RelativeEq
    + UlpsEq
    + 'static
{
    /// Two - useful for half-angle formulas in rotors.
    const TWO: Self;

    /// Pi constant (π ≈ 3.14159...).
    const PI: Self;

    /// Converts an `i8` to this float type.
    ///
    /// Used primarily for sign values (-1, 0, 1) in basis blade multiplication.
    fn from_i8(value: i8) -> Self;

    /// Converts a `usize` to this float type.
    ///
    /// Used for various index-based computations.
    fn from_usize(value: usize) -> Self;

    /// Converts an `f64` to this float type.
    ///
    /// Used for converting constants and test thresholds to the appropriate precision.
    /// Note: Converting from f64 to f32 may lose precision.
    fn from_f64(value: f64) -> Self;
}

impl Float for f32 {
    const TWO: Self = 2.0;
    const PI: Self = std::f32::consts::PI;

    #[inline]
    fn from_i8(value: i8) -> Self {
        value as Self
    }

    #[inline]
    fn from_usize(value: usize) -> Self {
        value as Self
    }

    #[inline]
    fn from_f64(value: f64) -> Self {
        value as Self
    }
}

impl Float for f64 {
    const TWO: Self = 2.0;
    const PI: Self = std::f64::consts::PI;

    #[inline]
    fn from_i8(value: i8) -> Self {
        value as Self
    }

    #[inline]
    fn from_usize(value: usize) -> Self {
        value as Self
    }

    #[inline]
    fn from_f64(value: f64) -> Self {
        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn f64_from_i8_roundtrip(x in -128i8..127) {
            let f = f64::from_i8(x);
            prop_assert_eq!(f as i8, x);
        }

        #[test]
        fn f32_from_i8_roundtrip(x in -128i8..127) {
            let f = f32::from_i8(x);
            prop_assert_eq!(f as i8, x);
        }
    }
}
