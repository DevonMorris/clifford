//! Domain-specific extensions for 2D Euclidean GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to Euclidean 2D geometry.

use super::generated::products;
use super::generated::types::{Bivector, Rotor, Vector};
use crate::scalar::Float;

// ============================================================================
// Vector extensions
// ============================================================================

impl<T: Float> Vector<T> {
    /// Dot product (inner product): `a · b = a.x*b.x + a.y*b.y`.
    ///
    /// Returns the scalar part of the geometric product.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let a = Vector::new(1.0, 2.0);
    /// let b = Vector::new(3.0, 4.0);
    /// assert_eq!(a.dot(b), 11.0); // 1*3 + 2*4
    /// ```
    #[inline]
    pub fn dot(self, other: Self) -> T {
        self.x() * other.x() + self.y() * other.y()
    }

    /// Wedge product (outer product): `a ∧ b`.
    ///
    /// Returns the bivector (pseudoscalar in 2D) representing the
    /// signed area of the parallelogram spanned by the vectors.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::{Vector, Bivector};
    ///
    /// let a = Vector::<f64>::unit_x();
    /// let b = Vector::<f64>::unit_y();
    /// let ab = a.wedge(b);
    /// assert_eq!(ab.xy(), 1.0);
    /// ```
    #[inline]
    pub fn wedge(self, other: Self) -> Bivector<T> {
        products::exterior_vector_vector(&self, &other)
    }

    /// Geometric product of two vectors: `ab = a·b + a∧b`.
    ///
    /// Returns a rotor (scalar + bivector).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let a = Vector::new(1.0, 0.0);
    /// let b = Vector::new(0.0, 1.0);
    /// let ab = a.geometric(b);
    /// assert_eq!(ab.s(), 0.0);  // perpendicular, no dot product
    /// assert_eq!(ab.xy(), 1.0); // wedge product
    /// ```
    #[inline]
    pub fn geometric(self, other: Self) -> Rotor<T> {
        products::geometric_vector_vector(&self, &other)
    }

    /// Perpendicular vector (90° counterclockwise rotation).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let v = Vector::new(1.0, 0.0);
    /// let perp = v.perp();
    /// assert_eq!(perp.x(), 0.0);
    /// assert_eq!(perp.y(), 1.0);
    /// ```
    #[inline]
    pub fn perp(&self) -> Self {
        Self::new(-self.y(), self.x())
    }
}

// ============================================================================
// Bivector extensions
// ============================================================================

impl<T: Float> Bivector<T> {
    /// Creates the unit bivector `e₁₂`.
    ///
    /// Alias for `unit_xy()` for backward compatibility.
    #[inline]
    pub fn unit() -> Self {
        Self::unit_xy()
    }

    /// Returns the coefficient (alias for `xy()`).
    ///
    /// Provided for backward compatibility.
    #[inline]
    pub fn value(&self) -> T {
        self.xy()
    }
}

// ============================================================================
// Rotor extensions
// ============================================================================

impl<T: Float> Rotor<T> {
    /// Creates a rotor from a rotation angle.
    ///
    /// The rotor `R = cos(θ/2) + sin(θ/2)e₁₂` rotates by angle `θ`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::{Rotor, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// // 90° rotation
    /// let r = Rotor::from_angle(FRAC_PI_2);
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn from_angle(angle: T) -> Self {
        let half = angle / T::TWO;
        // Use new_unchecked since angle construction guarantees unit norm
        Self::new_unchecked(half.cos(), half.sin())
    }

    /// Creates a rotor that rotates vector `a` to vector `b`.
    ///
    /// Both vectors should be non-zero. The rotation is the shortest
    /// arc from `a` to `b`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim2::{Rotor, Vector};
    /// use approx::abs_diff_eq;
    ///
    /// let a = Vector::<f64>::unit_x();
    /// let b = Vector::<f64>::unit_y();
    /// let r = Rotor::from_vectors(a, b);
    /// let rotated = r.rotate(a);
    /// assert!(abs_diff_eq!(rotated.x(), b.x(), epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.y(), b.y(), epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        // Compute signed angle from a to b using atan2 for stability
        let angle_a = a.y().atan2(a.x());
        let angle_b = b.y().atan2(b.x());
        let angle = angle_b - angle_a;
        Self::from_angle(angle)
    }

    /// Returns the inverse rotor: `R⁻¹ = R̃ / |R|²`.
    ///
    /// For unit rotors, this is equivalent to the reverse.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        let rev = self.reverse();
        // Use new_unchecked since inverse preserves unit constraint
        Self::new_unchecked(rev.s() / norm_sq, rev.xy() / norm_sq)
    }

    /// Applies this rotation to a vector: `v' = R̃ v R`.
    ///
    /// In 2D, both `R v R̃` and `R̃ v R` give the same result since the
    /// even subalgebra is commutative. We document `R̃ v R` for consistency
    /// with the 3D convention, which gives counterclockwise rotation.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::{Rotor, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let r = Rotor::from_angle(FRAC_PI_2);
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn rotate(&self, _v: Vector<T>) -> Vector<T> {
        todo!("rotate needs generated sandwich product")
    }

    /// Composes two rotations: `R₂ ∘ R₁ = R₂ R₁`.
    ///
    /// The result applies `self` first, then `other`.
    #[inline]
    pub fn compose(&self, other: Self) -> Self {
        products::geometric_rotor_rotor(&other, self)
    }

    /// Linear interpolation (normalized).
    ///
    /// Interpolates linearly between two rotors and normalizes the result.
    ///
    /// # Arguments
    ///
    /// * `other` - Target rotor
    /// * `t` - Interpolation parameter in [0, 1]
    ///
    /// # Panics
    ///
    /// Panics if `self` and `other` are opposite rotors (s₁ ≈ -s₂, xy₁ ≈ -xy₂),
    /// as the linear interpolation passes through zero at t ≈ 0.5.
    /// Use [`slerp`](Self::slerp) for robust interpolation in all cases.
    #[inline]
    pub fn lerp(&self, other: Self, t: T) -> Self {
        Self::new_unchecked(
            self.s() * (T::one() - t) + other.s() * t,
            self.xy() * (T::one() - t) + other.xy() * t,
        )
        .normalize()
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
        let dot = self.s() * other.s() + self.xy() * other.xy();

        let dot = if dot > T::one() {
            T::one()
        } else if dot < -T::one() {
            -T::one()
        } else {
            dot
        };

        let theta = dot.acos();
        let sin_theta = theta.sin();

        // Fall back to lerp when sin_theta is near zero (rotors nearly identical or opposite)
        if sin_theta.abs() < T::epsilon() {
            return self.lerp(other, t);
        }
        let s1 = ((T::one() - t) * theta).sin() / sin_theta;
        let s2 = (t * theta).sin() / sin_theta;

        // Spherical interpolation preserves unit norm
        Self::new_unchecked(
            self.s() * s1 + other.s() * s2,
            self.xy() * s1 + other.xy() * s2,
        )
    }

    /// Returns the rotation angle in radians.
    #[inline]
    pub fn angle(&self) -> T {
        self.xy().atan2(self.s()) * T::TWO
    }
}

// ============================================================================
// Even type (alias for Rotor)
// ============================================================================

/// Even subalgebra element (scalar + bivector).
///
/// This is a type alias for [`Rotor`], since in 2D Euclidean GA
/// the even subalgebra elements and rotors have the same structure.
pub type Even<T> = Rotor<T>;
