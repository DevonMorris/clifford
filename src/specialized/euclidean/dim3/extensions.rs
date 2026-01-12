//! Domain-specific extensions for 3D Euclidean GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to Euclidean 3D geometry.

use super::generated::products;
use super::generated::types::{Bivector, Rotor, Trivector, Vector};
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

    /// Wedge product (outer product): `a ∧ b`.
    ///
    /// Returns the bivector representing the oriented plane spanned by
    /// the two vectors.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Vector, Bivector};
    ///
    /// let a = Vector::<f64>::unit_x();
    /// let b = Vector::<f64>::unit_y();
    /// let ab = a.wedge(b);
    /// assert_eq!(ab, Bivector::unit_xy());
    /// ```
    #[inline]
    pub fn wedge(self, other: Self) -> Bivector<T> {
        products::exterior_vector_vector(&self, &other)
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
        self.wedge(other).dual()
    }

    /// Geometric product of two vectors: `ab = a·b + a∧b`.
    ///
    /// Returns a rotor (scalar + bivector).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let a = Vector::new(1.0, 0.0, 0.0);
    /// let b = Vector::new(0.0, 1.0, 0.0);
    /// let ab = a.geometric(b);
    /// assert_eq!(ab.s(), 0.0);  // perpendicular, no dot product
    /// ```
    #[inline]
    pub fn geometric(self, other: Self) -> Rotor<T> {
        products::geometric_vector_vector(&self, &other)
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
        Vector::new(self.yz(), -self.xz(), self.xy())
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
// Rotor extensions
// ============================================================================

impl<T: Float> Rotor<T> {
    /// Returns the bivector part as a Bivector.
    ///
    /// Provided for backward compatibility with code that used
    /// the previous struct layout.
    #[inline]
    pub fn b(&self) -> Bivector<T> {
        Bivector::new(self.xy(), self.xz(), self.yz())
    }

    /// Creates a rotor from an angle and rotation plane (unit bivector).
    ///
    /// The rotor `R = cos(θ/2) + sin(θ/2)B` rotates by angle `θ` in plane `B`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Bivector, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// // 90° rotation in xy-plane
    /// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn from_angle_plane(angle: T, plane: Bivector<T>) -> Self {
        let half = angle / T::TWO;
        let cos_half = half.cos();
        let sin_half = half.sin();
        // Use new_unchecked since angle/plane construction guarantees unit norm
        Self::new_unchecked(
            cos_half,
            plane.xy() * sin_half,
            plane.xz() * sin_half,
            plane.yz() * sin_half,
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
    /// ```
    #[inline]
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        let dot = a.dot(b);
        let wedge = a.wedge(b);

        let sum_sq = (T::one() + dot) * (T::one() + dot) + wedge.norm_squared();

        if sum_sq < T::epsilon() {
            // Vectors are anti-parallel, need to find perpendicular axis
            let perp = if a.x().abs() < a.y().abs() && a.x().abs() < a.z().abs() {
                Vector::unit_x()
            } else if a.y().abs() < a.z().abs() {
                Vector::unit_y()
            } else {
                Vector::unit_z()
            };
            let axis = a.cross(perp).normalized();
            let plane = a.wedge(axis).normalized();
            return Self::from_angle_plane(T::PI, plane);
        }

        let norm = sum_sq.sqrt();
        // Use new_unchecked since we're constructing from geometric product
        Self::new_unchecked(
            (T::one() + dot) / norm,
            wedge.xy() / norm,
            wedge.xz() / norm,
            wedge.yz() / norm,
        )
    }

    /// Returns the inverse rotor: `R⁻¹ = R̃ / |R|²`.
    ///
    /// For unit rotors, this is equivalent to the reverse.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        let rev = self.reverse();
        // Use new_unchecked since inverse preserves unit constraint
        Self::new_unchecked(
            rev.s() / norm_sq,
            rev.xy() / norm_sq,
            rev.xz() / norm_sq,
            rev.yz() / norm_sq,
        )
    }

    /// Applies this rotation to a vector: `v' = R̃ v R`.
    ///
    /// Note: Some sources use `R v R̃`. We use `R̃ v R` which gives counterclockwise
    /// rotation when looking along the rotation axis.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{Rotor, Bivector, Vector};
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
        let s = self.s();
        let bxy = self.xy();
        let bxz = self.xz();
        let byz = self.yz();

        // Intermediate: q = R̃ v (note: R̃ = s - b, so bivector signs flip)
        let qx = s * v.x() - bxy * v.y() - bxz * v.z();
        let qy = s * v.y() + bxy * v.x() - byz * v.z();
        let qz = s * v.z() + bxz * v.x() + byz * v.y();
        let qt = -bxy * v.z() + bxz * v.y() - byz * v.x();

        // Result: q R (only vector part survives)
        Vector::new(
            s * qx - bxy * qy - bxz * qz - byz * qt,
            s * qy + bxy * qx + bxz * qt - byz * qz,
            s * qz - bxy * qt + bxz * qx + byz * qy,
        )
    }

    /// Composes two rotations: `R₂ ∘ R₁ = R₂ R₁`.
    ///
    /// The result applies `self` first, then `other`.
    #[inline]
    pub fn compose(&self, other: Self) -> Self {
        products::geometric_rotor_rotor(&other, self)
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
            + self.xy() * other.xy()
            + self.xz() * other.xz()
            + self.yz() * other.yz();

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
                self.xy() * (T::one() - t) + other.xy() * t,
                self.xz() * (T::one() - t) + other.xz() * t,
                self.yz() * (T::one() - t) + other.yz() * t,
            )
            .normalize();
        }

        let sin_theta = theta.sin();
        let s1 = ((T::one() - t) * theta).sin() / sin_theta;
        let s2 = (t * theta).sin() / sin_theta;

        // Spherical interpolation preserves unit norm
        Self::new_unchecked(
            self.s() * s1 + other.s() * s2,
            self.xy() * s1 + other.xy() * s2,
            self.xz() * s1 + other.xz() * s2,
            self.yz() * s1 + other.yz() * s2,
        )
    }
}

// ============================================================================
// Even type (alias for Rotor)
// ============================================================================

/// Even subalgebra element (scalar + bivector).
///
/// This is a type alias for [`Rotor`], since in 3D Euclidean GA
/// the even subalgebra elements and rotors have the same structure.
pub type Even<T> = Rotor<T>;
