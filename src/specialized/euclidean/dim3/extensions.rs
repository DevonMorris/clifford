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
    /// let b = Bivector::<f64>::unit_xy();
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
        Self::new(self.xy() / n, self.xz() / n, self.yz() / n)
    }
}

// ============================================================================
// Trivector extensions
// ============================================================================

impl<T: Float> Trivector<T> {
    /// Creates the unit pseudoscalar `e₁₂₃`.
    ///
    /// Alias for `unit_xyz()` for backward compatibility.
    #[inline]
    pub fn unit() -> Self {
        Self::unit_xyz()
    }

    /// Returns the coefficient (alias for `xyz()`).
    ///
    /// Provided for backward compatibility.
    #[inline]
    pub fn value(&self) -> T {
        self.xyz()
    }
}

// ============================================================================
// Rotor extensions (grades [1, 3]: vector + trivector)
// ============================================================================

impl<T: Float> Rotor<T> {
    /// Returns the vector part as a Vector.
    #[inline]
    pub fn v(&self) -> Vector<T> {
        Vector::new(self.x(), self.y(), self.z())
    }

    /// Returns the inverse rotor: `R⁻¹ = R̃ / |R|²`.
    ///
    /// For unit rotors, this is equivalent to the reverse.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        let rev = self.reverse();
        Self::new_unchecked(
            rev.x() / norm_sq,
            rev.y() / norm_sq,
            rev.z() / norm_sq,
            rev.xyz() / norm_sq,
        )
    }

    /// Applies this transformation to a vector via the antisandwich product.
    #[inline]
    pub fn rotate(&self, v: Vector<T>) -> Vector<T> {
        self.transform(&v)
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
        let dot = self.x() * other.x()
            + self.y() * other.y()
            + self.z() * other.z()
            + self.xyz() * other.xyz();

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
                self.x() * (T::one() - t) + other.x() * t,
                self.y() * (T::one() - t) + other.y() * t,
                self.z() * (T::one() - t) + other.z() * t,
                self.xyz() * (T::one() - t) + other.xyz() * t,
            )
            .normalize();
        }

        let sin_theta = theta.sin();
        let s1 = ((T::one() - t) * theta).sin() / sin_theta;
        let s2 = (t * theta).sin() / sin_theta;

        // Spherical interpolation preserves unit norm
        Self::new_unchecked(
            self.x() * s1 + other.x() * s2,
            self.y() * s1 + other.y() * s2,
            self.z() * s1 + other.z() * s2,
            self.xyz() * s1 + other.xyz() * s2,
        )
    }
}

