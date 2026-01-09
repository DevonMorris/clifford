//! Type definitions for 2D geometric algebra.

use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::scalar::Float;

/// 2D vector (grade 1): `e₁`, `e₂`.
///
/// Represents a direction or position in 2D space.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Vec2<T: Float> {
    /// Coefficient of `e₁` (x-direction).
    pub x: T,
    /// Coefficient of `e₂` (y-direction).
    pub y: T,
}

impl<T: Float> Vec2<T> {
    /// Creates a new 2D vector.
    #[inline]
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }

    /// Creates the zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(T::ZERO, T::ZERO)
    }

    /// Creates the unit x-vector `e₁`.
    #[inline]
    pub fn unit_x() -> Self {
        Self::new(T::ONE, T::ZERO)
    }

    /// Creates the unit y-vector `e₂`.
    #[inline]
    pub fn unit_y() -> Self {
        Self::new(T::ZERO, T::ONE)
    }

    /// Returns the squared magnitude.
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    /// Returns the magnitude.
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized (unit length) vector.
    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.x / n, self.y / n)
    }

    /// Dot product: `a · b = a.x*b.x + a.y*b.y`.
    #[inline]
    pub fn dot(self, other: Self) -> T {
        self.x * other.x + self.y * other.y
    }

    /// Wedge product: `a ∧ b`.
    ///
    /// In 2D, this returns a bivector (the `e₁₂` coefficient).
    #[inline]
    pub fn wedge(self, other: Self) -> Bivec2<T> {
        Bivec2(self.x * other.y - self.y * other.x)
    }

    /// Perpendicular vector (90° counterclockwise rotation).
    #[inline]
    pub fn perp(&self) -> Self {
        Self::new(-self.y, self.x)
    }

    /// Geometric product of two vectors: `ab = a·b + a∧b`.
    #[inline]
    pub fn geometric(self, other: Self) -> Rotor2<T> {
        Rotor2 {
            s: self.dot(other),
            xy: self.x * other.y - self.y * other.x,
        }
    }
}

impl<T: Float> Default for Vec2<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// 2D bivector/pseudoscalar (grade 2): `e₁₂`.
///
/// In 2D, there's only one basis bivector, so this is a single coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Bivec2<T: Float>(
    /// Coefficient of `e₁₂`.
    pub T,
);

impl<T: Float> Bivec2<T> {
    /// Creates a new bivector.
    #[inline]
    pub fn new(value: T) -> Self {
        Self(value)
    }

    /// Creates the zero bivector.
    #[inline]
    pub fn zero() -> Self {
        Self(T::ZERO)
    }

    /// Creates the unit bivector `e₁₂`.
    #[inline]
    pub fn unit() -> Self {
        Self(T::ONE)
    }

    /// Returns the coefficient.
    #[inline]
    pub fn value(&self) -> T {
        self.0
    }

    /// Reverses the bivector: `(e₁₂)̃ = -e₁₂`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self(-self.0)
    }
}

impl<T: Float> Default for Bivec2<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// 2D rotor: scalar + bivector.
///
/// Represents a rotation in 2D. Equivalent to a complex number of unit magnitude.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Rotor2<T: Float> {
    /// Scalar part (grade 0).
    pub s: T,
    /// Bivector coefficient of `e₁₂` (grade 2).
    pub xy: T,
}

impl<T: Float> Rotor2<T> {
    /// Creates a new rotor.
    #[inline]
    pub fn new(s: T, xy: T) -> Self {
        Self { s, xy }
    }

    /// Creates the identity rotor (no rotation).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::ONE,
            xy: T::ZERO,
        }
    }

    /// Creates a rotor from a rotation angle.
    #[inline]
    pub fn from_angle(angle: T) -> Self {
        let half = angle / T::TWO;
        let cos_half = half.cos();
        let sin_half = half.sin();
        Self {
            s: cos_half,
            xy: sin_half,
        }
    }

    /// Creates a rotor that rotates vector `a` to vector `b`.
    #[inline]
    pub fn from_vectors(a: Vec2<T>, b: Vec2<T>) -> Self {
        let dot = a.dot(b);
        let wedge = a.x * b.y - a.y * b.x; // a ∧ b

        let norm = ((T::ONE + dot) * (T::ONE + dot) + wedge * wedge).sqrt();

        Self {
            s: (T::ONE + dot) / norm,
            xy: wedge / norm,
        }
    }

    /// Returns the squared magnitude.
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.s * self.s + self.xy * self.xy
    }

    /// Returns the magnitude.
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized (unit) rotor.
    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self {
            s: self.s / n,
            xy: self.xy / n,
        }
    }

    /// Returns the reverse: `R̃ = s - xy·e₁₂`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self {
            s: self.s,
            xy: -self.xy,
        }
    }

    /// Returns the inverse rotor.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        Self {
            s: self.s / norm_sq,
            xy: -self.xy / norm_sq,
        }
    }

    /// Applies this rotation to a vector: `v' = R̃ v R`.
    ///
    /// In 2D, both `R v R̃` and `R̃ v R` give the same result since the
    /// even subalgebra is commutative. We document `R̃ v R` for consistency
    /// with the 3D convention, which gives counterclockwise rotation.
    #[inline]
    pub fn rotate(&self, v: Vec2<T>) -> Vec2<T> {
        // Simplified for 2D: equivalent to rotation matrix
        let cos_theta = self.s * self.s - self.xy * self.xy;
        let sin_theta = T::TWO * self.s * self.xy;

        Vec2 {
            x: cos_theta * v.x - sin_theta * v.y,
            y: sin_theta * v.x + cos_theta * v.y,
        }
    }

    /// Composes two rotations: `R₂ ∘ R₁ = R₂ R₁`.
    #[inline]
    pub fn compose(&self, other: Self) -> Self {
        // Complex multiplication
        Self {
            s: self.s * other.s - self.xy * other.xy,
            xy: self.s * other.xy + self.xy * other.s,
        }
    }

    /// Linear interpolation (normalized).
    #[inline]
    pub fn lerp(&self, other: Self, t: T) -> Self {
        Self {
            s: self.s * (T::ONE - t) + other.s * t,
            xy: self.xy * (T::ONE - t) + other.xy * t,
        }
        .normalized()
    }

    /// Spherical linear interpolation.
    #[inline]
    pub fn slerp(&self, other: Self, t: T) -> Self {
        let dot = self.s * other.s + self.xy * other.xy;

        let dot = if dot > T::ONE {
            T::ONE
        } else if dot < -T::ONE {
            -T::ONE
        } else {
            dot
        };

        let theta = dot.acos();

        if theta.abs() < T::EPSILON {
            return self.lerp(other, t);
        }

        let sin_theta = theta.sin();
        let s1 = ((T::ONE - t) * theta).sin() / sin_theta;
        let s2 = (t * theta).sin() / sin_theta;

        Self {
            s: self.s * s1 + other.s * s2,
            xy: self.xy * s1 + other.xy * s2,
        }
    }

    /// Returns the rotation angle in radians.
    #[inline]
    pub fn angle(&self) -> T {
        self.xy.atan2(self.s) * T::TWO
    }
}

impl<T: Float> Default for Rotor2<T> {
    fn default() -> Self {
        Self::identity()
    }
}

// ============================================================================
// approx trait implementations for f32
// ============================================================================

impl AbsDiffEq for Vec2<f32> {
    type Epsilon = f32;

    fn default_epsilon() -> Self::Epsilon {
        f32::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f32::abs_diff_eq(&self.x, &other.x, epsilon) && f32::abs_diff_eq(&self.y, &other.y, epsilon)
    }
}

impl RelativeEq for Vec2<f32> {
    fn default_max_relative() -> Self::Epsilon {
        f32::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        f32::relative_eq(&self.x, &other.x, epsilon, max_relative)
            && f32::relative_eq(&self.y, &other.y, epsilon, max_relative)
    }
}

impl UlpsEq for Vec2<f32> {
    fn default_max_ulps() -> u32 {
        f32::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        f32::ulps_eq(&self.x, &other.x, epsilon, max_ulps)
            && f32::ulps_eq(&self.y, &other.y, epsilon, max_ulps)
    }
}

impl AbsDiffEq for Bivec2<f32> {
    type Epsilon = f32;

    fn default_epsilon() -> Self::Epsilon {
        f32::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f32::abs_diff_eq(&self.0, &other.0, epsilon)
    }
}

impl RelativeEq for Bivec2<f32> {
    fn default_max_relative() -> Self::Epsilon {
        f32::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        f32::relative_eq(&self.0, &other.0, epsilon, max_relative)
    }
}

impl UlpsEq for Bivec2<f32> {
    fn default_max_ulps() -> u32 {
        f32::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        f32::ulps_eq(&self.0, &other.0, epsilon, max_ulps)
    }
}

impl AbsDiffEq for Rotor2<f32> {
    type Epsilon = f32;

    fn default_epsilon() -> Self::Epsilon {
        f32::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f32::abs_diff_eq(&self.s, &other.s, epsilon)
            && f32::abs_diff_eq(&self.xy, &other.xy, epsilon)
    }
}

impl RelativeEq for Rotor2<f32> {
    fn default_max_relative() -> Self::Epsilon {
        f32::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        f32::relative_eq(&self.s, &other.s, epsilon, max_relative)
            && f32::relative_eq(&self.xy, &other.xy, epsilon, max_relative)
    }
}

impl UlpsEq for Rotor2<f32> {
    fn default_max_ulps() -> u32 {
        f32::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        f32::ulps_eq(&self.s, &other.s, epsilon, max_ulps)
            && f32::ulps_eq(&self.xy, &other.xy, epsilon, max_ulps)
    }
}

// ============================================================================
// approx trait implementations for f64
// ============================================================================

impl AbsDiffEq for Vec2<f64> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f64::abs_diff_eq(&self.x, &other.x, epsilon) && f64::abs_diff_eq(&self.y, &other.y, epsilon)
    }
}

impl RelativeEq for Vec2<f64> {
    fn default_max_relative() -> Self::Epsilon {
        f64::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        f64::relative_eq(&self.x, &other.x, epsilon, max_relative)
            && f64::relative_eq(&self.y, &other.y, epsilon, max_relative)
    }
}

impl UlpsEq for Vec2<f64> {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        f64::ulps_eq(&self.x, &other.x, epsilon, max_ulps)
            && f64::ulps_eq(&self.y, &other.y, epsilon, max_ulps)
    }
}

impl AbsDiffEq for Bivec2<f64> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f64::abs_diff_eq(&self.0, &other.0, epsilon)
    }
}

impl RelativeEq for Bivec2<f64> {
    fn default_max_relative() -> Self::Epsilon {
        f64::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        f64::relative_eq(&self.0, &other.0, epsilon, max_relative)
    }
}

impl UlpsEq for Bivec2<f64> {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        f64::ulps_eq(&self.0, &other.0, epsilon, max_ulps)
    }
}

impl AbsDiffEq for Rotor2<f64> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f64::abs_diff_eq(&self.s, &other.s, epsilon)
            && f64::abs_diff_eq(&self.xy, &other.xy, epsilon)
    }
}

impl RelativeEq for Rotor2<f64> {
    fn default_max_relative() -> Self::Epsilon {
        f64::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        f64::relative_eq(&self.s, &other.s, epsilon, max_relative)
            && f64::relative_eq(&self.xy, &other.xy, epsilon, max_relative)
    }
}

impl UlpsEq for Rotor2<f64> {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        f64::ulps_eq(&self.s, &other.s, epsilon, max_ulps)
            && f64::ulps_eq(&self.xy, &other.xy, epsilon, max_ulps)
    }
}
