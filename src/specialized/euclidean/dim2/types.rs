//! Type definitions for 2D geometric algebra.

use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::scalar::Float;

/// 2D vector (grade 1): `e₁`, `e₂`.
///
/// Represents a direction or position in 2D space.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Vector<T: Float> {
    /// Coefficient of `e₁` (x-direction).
    pub x: T,
    /// Coefficient of `e₂` (y-direction).
    pub y: T,
}

impl<T: Float> Vector<T> {
    /// Creates a new 2D vector.
    #[inline]
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }

    /// Creates the zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(T::zero(), T::zero())
    }

    /// Creates the unit x-vector `e₁`.
    #[inline]
    pub fn unit_x() -> Self {
        Self::new(T::one(), T::zero())
    }

    /// Creates the unit y-vector `e₂`.
    #[inline]
    pub fn unit_y() -> Self {
        Self::new(T::zero(), T::one())
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

    /// Scales the vector by a scalar: `s * v`.
    ///
    /// This is equivalent to `v * s` but allows scalar-first syntax,
    /// which is useful with generic scalar types like `Dual<f64>`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let v = Vector::new(1.0, 2.0);
    /// let scaled = v.scale(2.0);
    /// assert_eq!(scaled, Vector::new(2.0, 4.0));
    /// ```
    #[inline]
    pub fn scale(self, scalar: T) -> Self {
        Self::new(self.x * scalar, self.y * scalar)
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
    pub fn wedge(self, other: Self) -> Bivector<T> {
        Bivector(self.x * other.y - self.y * other.x)
    }

    /// Perpendicular vector (90° counterclockwise rotation).
    #[inline]
    pub fn perp(&self) -> Self {
        Self::new(-self.y, self.x)
    }

    /// Geometric product of two vectors: `ab = a·b + a∧b`.
    #[inline]
    pub fn geometric(self, other: Self) -> Rotor<T> {
        Rotor {
            s: self.dot(other),
            xy: self.x * other.y - self.y * other.x,
        }
    }
}

impl<T: Float> Default for Vector<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// 2D bivector/pseudoscalar (grade 2): `e₁₂`.
///
/// In 2D, there's only one basis bivector, so this is a single coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Bivector<T: Float>(
    /// Coefficient of `e₁₂`.
    pub T,
);

impl<T: Float> Bivector<T> {
    /// Creates a new bivector.
    #[inline]
    pub fn new(value: T) -> Self {
        Self(value)
    }

    /// Creates the zero bivector.
    #[inline]
    pub fn zero() -> Self {
        Self(T::zero())
    }

    /// Creates the unit bivector `e₁₂`.
    #[inline]
    pub fn unit() -> Self {
        Self(T::one())
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

impl<T: Float> Default for Bivector<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// 2D rotor: scalar + bivector.
///
/// Represents a rotation in 2D. Equivalent to a complex number of unit magnitude.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Rotor<T: Float> {
    /// Scalar part (grade 0).
    pub s: T,
    /// Bivector coefficient of `e₁₂` (grade 2).
    pub xy: T,
}

impl<T: Float> Rotor<T> {
    /// Creates a new rotor.
    #[inline]
    pub fn new(s: T, xy: T) -> Self {
        Self { s, xy }
    }

    /// Creates the identity rotor (no rotation).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::one(),
            xy: T::zero(),
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
    ///
    /// Uses `atan2` for numerical stability, especially when vectors are
    /// nearly parallel or anti-parallel.
    #[inline]
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        // Compute signed angle from a to b using atan2 for stability
        let angle_a = a.y.atan2(a.x);
        let angle_b = b.y.atan2(b.x);
        let angle = angle_b - angle_a;
        Self::from_angle(angle)
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
    pub fn rotate(&self, v: Vector<T>) -> Vector<T> {
        // Simplified for 2D: equivalent to rotation matrix
        let cos_theta = self.s * self.s - self.xy * self.xy;
        let sin_theta = T::TWO * self.s * self.xy;

        Vector {
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
            s: self.s * (T::one() - t) + other.s * t,
            xy: self.xy * (T::one() - t) + other.xy * t,
        }
        .normalized()
    }

    /// Spherical linear interpolation.
    #[inline]
    pub fn slerp(&self, other: Self, t: T) -> Self {
        let dot = self.s * other.s + self.xy * other.xy;

        let dot = if dot > T::one() {
            T::one()
        } else if dot < -T::one() {
            -T::one()
        } else {
            dot
        };

        let theta = dot.acos();

        if theta.abs() < T::epsilon() {
            return self.lerp(other, t);
        }

        let sin_theta = theta.sin();
        let s1 = ((T::one() - t) * theta).sin() / sin_theta;
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

impl<T: Float> Default for Rotor<T> {
    fn default() -> Self {
        Self::identity()
    }
}

// ============================================================================
// approx trait implementations (generic over Float)
// ============================================================================

impl<T: Float> AbsDiffEq for Vector<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.x, &other.x, epsilon) && T::abs_diff_eq(&self.y, &other.y, epsilon)
    }
}

impl<T: Float> RelativeEq for Vector<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.x, &other.x, epsilon, max_relative)
            && T::relative_eq(&self.y, &other.y, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Vector<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.x, &other.x, epsilon, max_ulps)
            && T::ulps_eq(&self.y, &other.y, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Bivector<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.0, &other.0, epsilon)
    }
}

impl<T: Float> RelativeEq for Bivector<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.0, &other.0, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Bivector<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.0, &other.0, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Rotor<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.s, &other.s, epsilon) && T::abs_diff_eq(&self.xy, &other.xy, epsilon)
    }
}

impl<T: Float> RelativeEq for Rotor<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.s, &other.s, epsilon, max_relative)
            && T::relative_eq(&self.xy, &other.xy, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Rotor<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.s, &other.s, epsilon, max_ulps)
            && T::ulps_eq(&self.xy, &other.xy, epsilon, max_ulps)
    }
}
