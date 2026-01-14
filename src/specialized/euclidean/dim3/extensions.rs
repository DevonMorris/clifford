//! Domain-specific extensions for 3D Euclidean GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to Euclidean 3D geometry.

use super::generated::types::{Bivector, Rotor, Trivector, Vector};
use crate::ops::{RightComplement, Transform, Wedge};
use crate::scalar::Float;

// ============================================================================
// Vector extensions
// ============================================================================

impl<T: Float> Vector<T> {
    /// Dot product (inner product): `a · b = a.x*b.x + a.y*b.y + a.z*b.z`.
    ///
    /// Returns the scalar part of the geometric product.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let a = Vector::new(1.0, 2.0, 3.0);
    /// let b = Vector::new(4.0, 5.0, 6.0);
    /// assert_eq!(a.dot(b), 32.0); // 1*4 + 2*5 + 3*6
    /// ```
    #[inline]
    pub fn dot(self, other: Self) -> T {
        self.x() * other.x() + self.y() * other.y() + self.z() * other.z()
    }

    /// Cross product: `a × b`.
    ///
    /// Computed as the Hodge dual of the wedge product.
    /// Returns a vector perpendicular to both inputs.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let x = Vector::<f64>::unit_x();
    /// let y = Vector::<f64>::unit_y();
    /// let z = x.cross(y);
    /// assert_eq!(z, Vector::unit_z());
    /// ```
    #[inline]
    pub fn cross(self, other: Self) -> Self {
        self.wedge(&other).dual()
    }

    /// Returns a normalized (unit length) version of this vector.
    ///
    /// This method divides by the norm without checking for zero.
    /// For a safe version, use `try_normalize()`.
    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.x() / n, self.y() / n, self.z() / n)
    }
}

// ============================================================================
// Bivector extensions
// ============================================================================

impl<T: Float> Bivector<T> {
    /// Computes the dual vector (Hodge star): `*B`.
    ///
    /// Maps bivector to vector: `*(e₁₂) = e₃`, `*(e₁₃) = -e₂`, `*(e₂₃) = e₁`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Vector, Bivector};
    ///
    /// let b = Bivector::<f64>::unit_bz();
    /// assert_eq!(b.dual(), Vector::unit_z());
    /// ```
    #[inline]
    pub fn dual(&self) -> Vector<T> {
        self.right_complement()
    }

    /// Returns a normalized (unit) bivector.
    ///
    /// This method divides by the norm without checking for zero.
    /// For a safe version, use `try_normalize()`.
    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.bz() / n, self.by() / n, self.bx() / n)
    }
}

// ============================================================================
// Trivector extensions
// ============================================================================

impl<T: Float> Trivector<T> {
    /// Creates the unit pseudoscalar `e₁₂₃`.
    #[inline]
    pub fn unit() -> Self {
        Self::unit_ps()
    }

    /// Returns the coefficient (alias for `ps()`).
    #[inline]
    pub fn value(&self) -> T {
        self.ps()
    }
}

// ============================================================================
// Rotor extensions (grades [0, 2]: scalar + bivector)
// ============================================================================

impl<T: Float> Rotor<T> {
    /// Returns the bivector part as a Bivector.
    #[inline]
    pub fn bivector(&self) -> Bivector<T> {
        Bivector::new(self.bz(), self.by(), self.bx())
    }

    /// Creates a rotor from a rotation angle and plane (bivector).
    ///
    /// The rotor `R = cos(θ/2) + sin(θ/2)B̂` rotates by angle `θ` in the plane `B`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// // 90° rotation in the xy-plane (around z-axis)
    /// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn from_angle_plane(angle: T, plane: Bivector<T>) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let b = plane.normalized();
        Self::new_unchecked(
            cos_half,
            sin_half * b.bz(),
            sin_half * b.by(),
            sin_half * b.bx(),
        )
    }

    /// Creates a rotor that rotates vector `a` to vector `b`.
    ///
    /// Both vectors should be non-zero. The rotation is the shortest
    /// arc from `a` to `b`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Vector};
    /// use approx::abs_diff_eq;
    ///
    /// let a = Vector::<f64>::unit_x();
    /// let b = Vector::<f64>::unit_y();
    /// let r = Rotor::from_vectors(a, b);
    /// let rotated = r.rotate(a);
    /// assert!(abs_diff_eq!(rotated.x(), b.x(), epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.y(), b.y(), epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.z(), b.z(), epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        // R = (1 + b*a) / |1 + b*a|
        // This gives the rotor for the shortest rotation from a to b
        let a_norm = a.normalized();
        let b_norm = b.normalized();

        // Compute 1 + b*a = 1 + b·a + b∧a
        let dot = b_norm.dot(a_norm);
        let wedge = b_norm.wedge(&a_norm);

        let r = Self::new_unchecked(T::one() + dot, wedge.bz(), wedge.by(), wedge.bx());
        r.normalize()
    }

    /// Returns the inverse rotor: `R⁻¹ = R̃ / |R|²`.
    ///
    /// For unit rotors, this is equivalent to the reverse.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        let rev = self.reverse();
        Self::new_unchecked(
            rev.s() / norm_sq,
            rev.bz() / norm_sq,
            rev.by() / norm_sq,
            rev.bx() / norm_sq,
        )
    }

    /// Applies this rotation to a vector via the sandwich product.
    ///
    /// For Euclidean GA (non-degenerate), the Transform trait uses sandwich.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn rotate(&self, v: Vector<T>) -> Vector<T> {
        self.transform(&v)
    }

    /// Composes two rotations: `R₂ ∘ R₁ = R₂ R₁`.
    ///
    /// The result applies `self` first, then `other`.
    #[inline]
    pub fn compose(&self, other: Self) -> Self {
        other * *self
    }

    /// Spherical linear interpolation between two rotors.
    ///
    /// Interpolates along the geodesic (great arc) between two rotors.
    /// Produces constant angular velocity.
    ///
    /// # Arguments
    ///
    /// * `other` - Target rotor
    /// * `t` - Interpolation parameter in [0, 1]
    #[inline]
    pub fn slerp(&self, other: Self, t: T) -> Self {
        let dot = self.s() * other.s()
            + self.bz() * other.bz()
            + self.by() * other.by()
            + self.bx() * other.bx();

        let dot = if dot > T::one() {
            T::one()
        } else if dot < -T::one() {
            -T::one()
        } else {
            dot
        };

        let theta = dot.acos();

        if theta.abs() < T::epsilon() {
            // Linear interpolation for small angles, then normalize
            return Self::new_unchecked(
                self.s() * (T::one() - t) + other.s() * t,
                self.bz() * (T::one() - t) + other.bz() * t,
                self.by() * (T::one() - t) + other.by() * t,
                self.bx() * (T::one() - t) + other.bx() * t,
            )
            .normalize();
        }

        let sin_theta = theta.sin();
        let s1 = ((T::one() - t) * theta).sin() / sin_theta;
        let s2 = (t * theta).sin() / sin_theta;

        // Spherical interpolation preserves unit norm
        Self::new_unchecked(
            self.s() * s1 + other.s() * s2,
            self.bz() * s1 + other.bz() * s2,
            self.by() * s1 + other.by() * s2,
            self.bx() * s1 + other.bx() * s2,
        )
    }
}
