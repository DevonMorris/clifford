//! Domain-specific extensions for Complex number types.
//!
//! This module adds transcendental functions and convenience methods
//! to the generated Complex type that are essential for complex analysis.

use super::generated::types::Complex;
use crate::scalar::Float;
use std::ops::Div;

// ============================================================================
// Complex extensions
// ============================================================================

impl<T: Float> Complex<T> {
    /// Creates the multiplicative identity (1 + 0i).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    ///
    /// let one = Complex::<f64>::one();
    /// assert_eq!(one.real(), 1.0);
    /// assert_eq!(one.imag(), 0.0);
    /// ```
    #[inline]
    pub fn one() -> Self {
        Self::new(T::one(), T::zero())
    }

    /// Creates a purely real complex number.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    ///
    /// let z = Complex::from_real(3.0);
    /// assert_eq!(z.real(), 3.0);
    /// assert_eq!(z.imag(), 0.0);
    /// ```
    #[inline]
    pub fn from_real(real: T) -> Self {
        Self::new(real, T::zero())
    }

    /// Creates a purely imaginary complex number.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    ///
    /// let z = Complex::from_imag(3.0);
    /// assert_eq!(z.real(), 0.0);
    /// assert_eq!(z.imag(), 3.0);
    /// ```
    #[inline]
    pub fn from_imag(imag: T) -> Self {
        Self::new(T::zero(), imag)
    }

    /// Creates the imaginary unit i.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    ///
    /// let i = Complex::<f64>::i();
    /// assert_eq!(i.real(), 0.0);
    /// assert_eq!(i.imag(), 1.0);
    /// ```
    #[inline]
    pub fn i() -> Self {
        Self::new(T::zero(), T::one())
    }

    /// Creates a complex number from polar coordinates.
    ///
    /// `z = r * (cos(theta) + i*sin(theta)) = r * e^(i*theta)`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let z = Complex::from_polar(2.0, FRAC_PI_4);
    /// assert_relative_eq!(z.norm(), 2.0, epsilon = 1e-10);
    /// assert_relative_eq!(z.arg(), FRAC_PI_4, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn from_polar(r: T, theta: T) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }

    /// Returns the complex conjugate (a - bi for a + bi).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    ///
    /// let z = Complex::new(3.0, 4.0);
    /// let conj = z.conjugate();
    /// assert_eq!(conj.real(), 3.0);
    /// assert_eq!(conj.imag(), -4.0);
    /// ```
    #[inline]
    pub fn conjugate(&self) -> Self {
        Self::new(self.real(), -self.imag())
    }

    /// Returns the argument (phase angle) in radians.
    ///
    /// Returns the angle theta such that z = |z| * e^(i*theta).
    /// The result is in the range (-pi, pi].
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use std::f64::consts::FRAC_PI_4;
    /// use approx::assert_relative_eq;
    ///
    /// let z = Complex::new(1.0, 1.0);
    /// assert_relative_eq!(z.arg(), FRAC_PI_4, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn arg(&self) -> T {
        self.imag().atan2(self.real())
    }

    /// Returns the complex exponential e^z.
    ///
    /// For z = a + bi: e^z = e^a * (cos(b) + i*sin(b))
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    /// use std::f64::consts::PI;
    ///
    /// // Euler's identity: e^(i*pi) = -1
    /// let z = Complex::new(0.0, PI);
    /// let result = z.exp();
    /// assert_relative_eq!(result.real(), -1.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn exp(&self) -> Self {
        let exp_real = self.real().exp();
        Self::new(exp_real * self.imag().cos(), exp_real * self.imag().sin())
    }

    /// Returns the principal natural logarithm.
    ///
    /// For z = r * e^(i*theta): ln(z) = ln(r) + i*theta
    ///
    /// Note: This is the principal branch with imaginary part in (-pi, pi].
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    ///
    /// let z = Complex::new(1.0, 0.0);
    /// let result = z.ln();
    /// assert_relative_eq!(result.real(), 0.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn ln(&self) -> Self {
        Self::new(self.norm().ln(), self.arg())
    }

    /// Returns the complex sine.
    ///
    /// sin(z) = (e^(iz) - e^(-iz)) / (2i)
    ///
    /// For z = a + bi:
    /// sin(z) = sin(a)*cosh(b) + i*cos(a)*sinh(b)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let z = Complex::new(FRAC_PI_2, 0.0);
    /// let result = z.sin();
    /// assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn sin(&self) -> Self {
        Self::new(
            self.real().sin() * self.imag().cosh(),
            self.real().cos() * self.imag().sinh(),
        )
    }

    /// Returns the complex cosine.
    ///
    /// cos(z) = (e^(iz) + e^(-iz)) / 2
    ///
    /// For z = a + bi:
    /// cos(z) = cos(a)*cosh(b) - i*sin(a)*sinh(b)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    ///
    /// let z = Complex::new(0.0, 0.0);
    /// let result = z.cos();
    /// assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn cos(&self) -> Self {
        Self::new(
            self.real().cos() * self.imag().cosh(),
            -self.real().sin() * self.imag().sinh(),
        )
    }

    /// Returns the principal square root.
    ///
    /// For z = r * e^(i*theta): sqrt(z) = sqrt(r) * e^(i*theta/2)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    ///
    /// let z = Complex::new(0.0, 1.0); // i
    /// let result = z.sqrt();
    /// // sqrt(i) = (1 + i) / sqrt(2)
    /// let expected = 0.5_f64.sqrt();
    /// assert_relative_eq!(result.real(), expected + expected, epsilon = 1e-10);
    /// // Actually sqrt(i) = (1/sqrt(2)) + i*(1/sqrt(2))
    /// ```
    #[inline]
    pub fn sqrt(&self) -> Self {
        let r = self.norm();
        let theta = self.arg();
        Self::from_polar(r.sqrt(), theta / T::TWO)
    }

    /// Returns z raised to an integer power.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    ///
    /// let z = Complex::new(0.0, 1.0); // i
    /// let result = z.powi(2);
    /// // i^2 = -1
    /// assert_relative_eq!(result.real(), -1.0, epsilon = 1e-10);
    /// assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    /// ```
    #[inline]
    pub fn powi(&self, n: i32) -> Self {
        if n == 0 {
            return Self::one();
        }

        let mut result = Self::one();
        let mut base = *self;
        let mut exp = n.unsigned_abs();

        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            exp >>= 1;
        }

        if n < 0 { Self::one() / result } else { result }
    }
}

impl<T: Float> Div for Complex<T> {
    type Output = Self;

    /// Divides two complex numbers.
    ///
    /// z1 / z2 = z1 * conj(z2) / |z2|^2
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::complex::Complex;
    /// use approx::assert_relative_eq;
    ///
    /// let z1 = Complex::new(2.0, 3.0);
    /// let z2 = Complex::new(1.0, 1.0);
    /// let result = z1 / z2;
    /// // (2 + 3i) / (1 + i) = (2 + 3i)(1 - i) / 2 = (5 + i) / 2
    /// assert_relative_eq!(result.real(), 2.5, epsilon = 1e-10);
    /// assert_relative_eq!(result.imag(), 0.5, epsilon = 1e-10);
    /// ```
    #[inline]
    fn div(self, other: Self) -> Self::Output {
        let denom = other.norm_squared();
        let conj = other.conjugate();
        let numer = self * conj;
        Self::new(numer.real() / denom, numer.imag() / denom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::{E, FRAC_PI_2, FRAC_PI_4, PI};

    #[test]
    fn test_one() {
        let one = Complex::<f64>::one();
        assert_eq!(one.real(), 1.0);
        assert_eq!(one.imag(), 0.0);
    }

    #[test]
    fn test_from_polar() {
        let z = Complex::from_polar(2.0, FRAC_PI_4);
        assert_relative_eq!(z.norm(), 2.0, epsilon = 1e-10);
        assert_relative_eq!(z.arg(), FRAC_PI_4, epsilon = 1e-10);
    }

    #[test]
    fn test_conjugate() {
        let z = Complex::new(3.0, 4.0);
        let conj = z.conjugate();
        assert_eq!(conj.real(), 3.0);
        assert_eq!(conj.imag(), -4.0);
    }

    #[test]
    fn test_arg() {
        let z = Complex::new(1.0, 1.0);
        assert_relative_eq!(z.arg(), FRAC_PI_4, epsilon = 1e-10);

        let z2 = Complex::new(-1.0, 0.0);
        assert_relative_eq!(z2.arg(), PI, epsilon = 1e-10);
    }

    #[test]
    fn test_exp_euler() {
        // e^(i*pi) = -1
        let z = Complex::new(0.0, PI);
        let result = z.exp();
        assert_relative_eq!(result.real(), -1.0, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_exp_real() {
        // e^1 = e
        let z = Complex::new(1.0, 0.0);
        let result = z.exp();
        assert_relative_eq!(result.real(), E, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_ln() {
        let z = Complex::new(E, 0.0);
        let result = z.ln();
        assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_sin() {
        let z = Complex::new(FRAC_PI_2, 0.0);
        let result = z.sin();
        assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cos() {
        let z = Complex::new(0.0, 0.0);
        let result = z.cos();
        assert_relative_eq!(result.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_sqrt() {
        // sqrt(4) = 2
        let z = Complex::new(4.0, 0.0);
        let result = z.sqrt();
        assert_relative_eq!(result.real(), 2.0, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_sqrt_negative() {
        // sqrt(-1) = i
        let z = Complex::new(-1.0, 0.0);
        let result = z.sqrt();
        assert_relative_eq!(result.real(), 0.0, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_powi() {
        let z = Complex::new(0.0, 1.0); // i
        let result = z.powi(2);
        assert_relative_eq!(result.real(), -1.0, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.0, epsilon = 1e-10);

        let result4 = z.powi(4);
        assert_relative_eq!(result4.real(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(result4.imag(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_division() {
        let z1 = Complex::new(2.0, 3.0);
        let z2 = Complex::new(1.0, 1.0);
        let result = z1 / z2;
        assert_relative_eq!(result.real(), 2.5, epsilon = 1e-10);
        assert_relative_eq!(result.imag(), 0.5, epsilon = 1e-10);
    }
}
