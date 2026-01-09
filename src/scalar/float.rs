//! Floating-point scalar type abstraction.
//!
//! This module provides the [`Float`] trait which abstracts over floating-point
//! types like `f32` and `f64`, allowing the library to be generic over precision.

use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use core::fmt::{Debug, Display};
use core::iter::{Product, Sum};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

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
    Copy
    + Clone
    + Default
    + PartialEq
    + PartialOrd
    + Debug
    + Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
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
    /// The additive identity (zero).
    const ZERO: Self;

    /// The multiplicative identity (one).
    const ONE: Self;

    /// Machine epsilon - the smallest value such that `1.0 + EPSILON != 1.0`.
    ///
    /// Used for approximate equality comparisons.
    const EPSILON: Self;

    /// Two - useful for half-angle formulas in rotors.
    const TWO: Self;

    /// Pi constant (π ≈ 3.14159...).
    const PI: Self;

    /// Returns the absolute value of `self`.
    fn abs(self) -> Self;

    /// Returns the square root of `self`.
    ///
    /// # Panics
    ///
    /// May panic or return NaN if `self` is negative.
    fn sqrt(self) -> Self;

    /// Returns the sine of `self` (in radians).
    fn sin(self) -> Self;

    /// Returns the cosine of `self` (in radians).
    fn cos(self) -> Self;

    /// Returns the arc tangent of `y / x`, using the signs to determine the quadrant.
    ///
    /// This is useful for converting between Cartesian and polar coordinates.
    fn atan2(self, x: Self) -> Self;

    /// Returns the arc cosine of `self` (in radians).
    ///
    /// The result is in the range `[0, π]`.
    fn acos(self) -> Self;

    /// Checks if `self` is approximately equal to `other` within `epsilon`.
    ///
    /// # Arguments
    ///
    /// * `other` - The value to compare against
    /// * `epsilon` - The maximum allowed difference
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::scalar::Float;
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 1.0 + 1e-12;
    /// assert!(a.approx_eq(b, 1e-10));
    /// ```
    fn approx_eq(self, other: Self, epsilon: Self) -> bool;

    /// Converts an `i8` to this float type.
    ///
    /// Used primarily for sign values (-1, 0, 1) in basis blade multiplication.
    fn from_i8(value: i8) -> Self;

    /// Converts a `usize` to this float type.
    ///
    /// Used for various index-based computations.
    fn from_usize(value: usize) -> Self;
}

impl Float for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const EPSILON: Self = f32::EPSILON;
    const TWO: Self = 2.0;
    const PI: Self = std::f32::consts::PI;

    #[inline]
    fn abs(self) -> Self {
        f32::abs(self)
    }

    #[inline]
    fn sqrt(self) -> Self {
        f32::sqrt(self)
    }

    #[inline]
    fn sin(self) -> Self {
        f32::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        f32::cos(self)
    }

    #[inline]
    fn atan2(self, x: Self) -> Self {
        f32::atan2(self, x)
    }

    #[inline]
    fn acos(self) -> Self {
        f32::acos(self)
    }

    #[inline]
    fn approx_eq(self, other: Self, epsilon: Self) -> bool {
        (self - other).abs() <= epsilon
    }

    #[inline]
    fn from_i8(value: i8) -> Self {
        value as Self
    }

    #[inline]
    fn from_usize(value: usize) -> Self {
        value as Self
    }
}

impl Float for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const EPSILON: Self = f64::EPSILON;
    const TWO: Self = 2.0;
    const PI: Self = std::f64::consts::PI;

    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }

    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }

    #[inline]
    fn sin(self) -> Self {
        f64::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        f64::cos(self)
    }

    #[inline]
    fn atan2(self, x: Self) -> Self {
        f64::atan2(self, x)
    }

    #[inline]
    fn acos(self) -> Self {
        f64::acos(self)
    }

    #[inline]
    fn approx_eq(self, other: Self, epsilon: Self) -> bool {
        (self - other).abs() <= epsilon
    }

    #[inline]
    fn from_i8(value: i8) -> Self {
        value as Self
    }

    #[inline]
    fn from_usize(value: usize) -> Self {
        value as Self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn f64_approx_eq_reflexive(x in -1000.0..1000.0f64) {
            prop_assert!(x.approx_eq(x, f64::EPSILON));
        }

        #[test]
        fn f64_from_i8_roundtrip(x in -128i8..127) {
            let f = f64::from_i8(x);
            prop_assert_eq!(f as i8, x);
        }

        #[test]
        fn f32_approx_eq_reflexive(x in -1000.0..1000.0f32) {
            prop_assert!(x.approx_eq(x, f32::EPSILON));
        }

        #[test]
        fn f32_from_i8_roundtrip(x in -128i8..127) {
            let f = f32::from_i8(x);
            prop_assert_eq!(f as i8, x);
        }
    }
}
