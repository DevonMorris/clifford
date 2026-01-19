//! Domain-specific extensions for Dual number types.
//!
//! This module adds transcendental functions and convenience methods
//! to the generated Dual type that enable automatic differentiation.
//!
//! # Automatic Differentiation
//!
//! Dual numbers have the form `a + bε` where `ε² = 0`. This nilpotent
//! property means that for any smooth function f:
//!
//! ```text
//! f(x + ε) = f(x) + f'(x)·ε
//! ```
//!
//! By evaluating `f(x + 1ε)` and extracting the dual part, we get
//! the exact derivative `f'(x)` with no truncation error!

use super::generated::types::Dual;
use crate::scalar::Float;
use std::ops::Div;

// ============================================================================
// Dual extensions
// ============================================================================

impl<T: Float> Dual<T> {
    /// Creates the multiplicative identity (1 + 0ε).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    ///
    /// let one = Dual::<f64>::one();
    /// assert_eq!(one.real(), 1.0);
    /// assert_eq!(one.dual(), 0.0);
    /// ```
    #[inline]
    pub fn one() -> Self {
        Self::new(T::one(), T::zero())
    }

    /// Creates a dual number from a real value with zero dual part.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    ///
    /// let d = Dual::from_real(3.0);
    /// assert_eq!(d.real(), 3.0);
    /// assert_eq!(d.dual(), 0.0);
    /// ```
    #[inline]
    pub fn from_real(real: T) -> Self {
        Self::new(real, T::zero())
    }

    /// Creates the infinitesimal unit ε.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    ///
    /// let eps = Dual::<f64>::epsilon();
    /// assert_eq!(eps.real(), 0.0);
    /// assert_eq!(eps.dual(), 1.0);
    /// ```
    #[inline]
    pub fn epsilon() -> Self {
        Self::new(T::zero(), T::one())
    }

    /// Creates a dual number for differentiation at point x.
    ///
    /// This creates `x + ε`, which when passed through a function f
    /// yields `f(x) + f'(x)ε`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    ///
    /// let d = Dual::variable(2.0);
    /// assert_eq!(d.real(), 2.0);
    /// assert_eq!(d.dual(), 1.0);
    /// ```
    #[inline]
    pub fn variable(x: T) -> Self {
        Self::new(x, T::one())
    }

    /// Evaluates f(x) and f'(x) simultaneously using dual number arithmetic.
    ///
    /// Given a function that operates on dual numbers, returns the value
    /// and derivative at point x.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    ///
    /// // Compute d/dx(x²) at x = 3
    /// let (value, derivative) = Dual::differentiate(3.0, |d| d * d);
    /// assert_eq!(value, 9.0);  // 3² = 9
    /// assert_eq!(derivative, 6.0);  // d/dx(x²) = 2x = 6
    /// ```
    #[inline]
    pub fn differentiate<F>(x: T, f: F) -> (T, T)
    where
        F: FnOnce(Dual<T>) -> Dual<T>,
    {
        let result = f(Self::variable(x));
        (result.real(), result.dual())
    }

    /// Returns the dual exponential e^d.
    ///
    /// For d = a + bε: e^d = e^a + b·e^a·ε
    ///
    /// This follows from f(x + ε) = f(x) + f'(x)ε where f(x) = e^x.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(e^x) at x = 0 is e^0 = 1
    /// let d = Dual::variable(0.0);
    /// let result = d.exp();
    /// assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 1.0, epsilon = 1e-10); // derivative
    /// ```
    #[inline]
    pub fn exp(&self) -> Self {
        let exp_real = self.real().exp();
        Self::new(exp_real, self.dual() * exp_real)
    }

    /// Returns the dual natural logarithm ln(d).
    ///
    /// For d = a + bε: ln(d) = ln(a) + (b/a)ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(ln(x)) at x = 2 is 1/2 = 0.5
    /// let d = Dual::variable(2.0);
    /// let result = d.ln();
    /// assert_relative_eq!(result.dual(), 0.5, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn ln(&self) -> Self {
        Self::new(self.real().ln(), self.dual() / self.real())
    }

    /// Returns the dual sine sin(d).
    ///
    /// For d = a + bε: sin(d) = sin(a) + b·cos(a)·ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // d/dx(sin(x)) at x = 0 is cos(0) = 1
    /// let d = Dual::variable(0.0);
    /// let result = d.sin();
    /// assert_relative_eq!(result.real(), 0.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 1.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn sin(&self) -> Self {
        Self::new(self.real().sin(), self.dual() * self.real().cos())
    }

    /// Returns the dual cosine cos(d).
    ///
    /// For d = a + bε: cos(d) = cos(a) - b·sin(a)·ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(cos(x)) at x = 0 is -sin(0) = 0
    /// let d = Dual::variable(0.0);
    /// let result = d.cos();
    /// assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 0.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn cos(&self) -> Self {
        Self::new(self.real().cos(), -self.dual() * self.real().sin())
    }

    /// Returns the dual tangent tan(d).
    ///
    /// For d = a + bε: tan(d) = tan(a) + b·sec²(a)·ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(tan(x)) at x = 0 is sec²(0) = 1
    /// let d = Dual::variable(0.0);
    /// let result = d.tan();
    /// assert_relative_eq!(result.real(), 0.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 1.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn tan(&self) -> Self {
        let cos_a = self.real().cos();
        let sec_sq = T::one() / (cos_a * cos_a);
        Self::new(self.real().tan(), self.dual() * sec_sq)
    }

    /// Returns the dual square root sqrt(d).
    ///
    /// For d = a + bε: sqrt(d) = sqrt(a) + (b / (2·sqrt(a)))·ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(sqrt(x)) at x = 4 is 1/(2·sqrt(4)) = 0.25
    /// let d = Dual::variable(4.0);
    /// let result = d.sqrt();
    /// assert_relative_eq!(result.real(), 2.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 0.25, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn sqrt(&self) -> Self {
        let sqrt_real = self.real().sqrt();
        Self::new(sqrt_real, self.dual() / (T::TWO * sqrt_real))
    }

    /// Returns the dual power d^n for integer n.
    ///
    /// For d = a + bε: d^n = a^n + n·a^(n-1)·b·ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(x³) at x = 2 is 3·2² = 12
    /// let d = Dual::variable(2.0);
    /// let result = d.powi(3);
    /// assert_relative_eq!(result.real(), 8.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 12.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn powi(&self, n: i32) -> Self {
        if n == 0 {
            return Self::one();
        }
        let a_pow_n = self.real().powi(n);
        let a_pow_nm1 = self.real().powi(n - 1);
        Self::new(a_pow_n, T::from_f64(f64::from(n)) * a_pow_nm1 * self.dual())
    }

    /// Returns the dual power d^p for real p.
    ///
    /// For d = a + bε: d^p = a^p + p·a^(p-1)·b·ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(x^1.5) at x = 4 is 1.5·4^0.5 = 3
    /// let d = Dual::variable(4.0);
    /// let result = d.powf(1.5);
    /// assert_relative_eq!(result.real(), 8.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 3.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn powf(&self, p: T) -> Self {
        let a_pow_p = self.real().powf(p);
        let a_pow_pm1 = self.real().powf(p - T::one());
        Self::new(a_pow_p, p * a_pow_pm1 * self.dual())
    }

    /// Returns the dual absolute value |d|.
    ///
    /// For d = a + bε: |d| = |a| + sign(a)·b·ε
    ///
    /// Note: Not differentiable at a = 0.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// let d = Dual::variable(-3.0);
    /// let result = d.abs();
    /// assert_relative_eq!(result.real(), 3.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), -1.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn abs(&self) -> Self {
        if self.real() >= T::zero() {
            *self
        } else {
            Self::new(-self.real(), -self.dual())
        }
    }
}

impl<T: Float> Div for Dual<T> {
    type Output = Self;

    /// Divides two dual numbers.
    ///
    /// For d1 = a + bε and d2 = c + dε:
    /// d1/d2 = (a/c) + ((b·c - a·d)/c²)ε
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::dual::Dual;
    /// use approx::assert_relative_eq;
    ///
    /// // d/dx(x/2) at x = 4 is 0.5
    /// let d1 = Dual::variable(4.0);
    /// let d2 = Dual::from_real(2.0);
    /// let result = d1 / d2;
    /// assert_relative_eq!(result.real(), 2.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.dual(), 0.5, epsilon = 1e-10);
    /// ```
    #[inline]
    fn div(self, other: Self) -> Self::Output {
        let c = other.real();
        let c_sq = c * c;
        Self::new(
            self.real() / c,
            (self.dual() * c - self.real() * other.dual()) / c_sq,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::{E, FRAC_PI_4};

    #[test]
    fn test_one() {
        let one = Dual::<f64>::one();
        assert_eq!(one.real(), 1.0);
        assert_eq!(one.dual(), 0.0);
    }

    #[test]
    fn test_variable() {
        let d = Dual::variable(3.0);
        assert_eq!(d.real(), 3.0);
        assert_eq!(d.dual(), 1.0);
    }

    #[test]
    fn test_differentiate_square() {
        let (value, derivative) = Dual::differentiate(3.0, |d| d * d);
        assert_eq!(value, 9.0);
        assert_eq!(derivative, 6.0);
    }

    #[test]
    fn test_differentiate_cubic() {
        let (value, derivative) = Dual::differentiate(2.0, |d| d * d * d);
        assert_eq!(value, 8.0);
        assert_eq!(derivative, 12.0);
    }

    #[test]
    fn test_exp() {
        let d = Dual::variable(0.0);
        let result = d.exp();
        assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 1.0, epsilon = 1e-10);

        let d2 = Dual::variable(1.0);
        let result2 = d2.exp();
        assert_relative_eq!(result2.real(), E, epsilon = 1e-10);
        assert_relative_eq!(result2.dual(), E, epsilon = 1e-10);
    }

    #[test]
    fn test_ln() {
        let d = Dual::variable(E);
        let result = d.ln();
        assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 1.0 / E, epsilon = 1e-10);
    }

    #[test]
    fn test_sin() {
        let d = Dual::variable(0.0);
        let result = d.sin();
        assert_relative_eq!(result.real(), 0.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cos() {
        let d = Dual::variable(0.0);
        let result = d.cos();
        assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_tan() {
        let d = Dual::variable(FRAC_PI_4);
        let result = d.tan();
        assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 2.0, epsilon = 1e-10); // sec²(π/4) = 2
    }

    #[test]
    fn test_sqrt() {
        let d = Dual::variable(4.0);
        let result = d.sqrt();
        assert_relative_eq!(result.real(), 2.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 0.25, epsilon = 1e-10);
    }

    #[test]
    fn test_powi() {
        let d = Dual::variable(2.0);
        let result = d.powi(3);
        assert_relative_eq!(result.real(), 8.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 12.0, epsilon = 1e-10);
    }

    #[test]
    fn test_powf() {
        let d = Dual::variable(4.0);
        let result = d.powf(1.5);
        assert_relative_eq!(result.real(), 8.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_division() {
        let d1 = Dual::variable(4.0);
        let d2 = Dual::from_real(2.0);
        let result = d1 / d2;
        assert_relative_eq!(result.real(), 2.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), 0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_chain_rule() {
        // d/dx(sin(x²)) at x = 1 should be 2x·cos(x²) = 2·cos(1)
        let (_, derivative) = Dual::differentiate(1.0, |d| (d * d).sin());
        assert_relative_eq!(derivative, 2.0 * 1.0_f64.cos(), epsilon = 1e-10);
    }

    #[test]
    fn test_abs() {
        let d = Dual::variable(-3.0);
        let result = d.abs();
        assert_relative_eq!(result.real(), 3.0, epsilon = 1e-10);
        assert_relative_eq!(result.dual(), -1.0, epsilon = 1e-10);
    }
}
