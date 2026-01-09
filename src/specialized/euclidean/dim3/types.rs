//! Type definitions for 3D geometric algebra.

use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::scalar::Float;

/// 3D vector (grade 1): `e₁`, `e₂`, `e₃`.
///
/// Represents a direction or position in 3D space.
///
/// # Components
///
/// - `x`: coefficient of `e₁` (x-direction)
/// - `y`: coefficient of `e₂` (y-direction)
/// - `z`: coefficient of `e₃` (z-direction)
///
/// # Example
///
/// ```
/// use clifford::specialized::euclidean::dim3::Vec3;
///
/// let v = Vec3::new(1.0, 2.0, 3.0);
/// assert_eq!(v.x, 1.0);
/// assert_eq!(v.y, 2.0);
/// assert_eq!(v.z, 3.0);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Vec3<T: Float> {
    /// Coefficient of `e₁` (x-direction).
    pub x: T,
    /// Coefficient of `e₂` (y-direction).
    pub y: T,
    /// Coefficient of `e₃` (z-direction).
    pub z: T,
}

impl<T: Float> Vec3<T> {
    /// Creates a new 3D vector.
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    /// Creates the zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(T::ZERO, T::ZERO, T::ZERO)
    }

    /// Creates the unit x-vector `e₁`.
    #[inline]
    pub fn unit_x() -> Self {
        Self::new(T::ONE, T::ZERO, T::ZERO)
    }

    /// Creates the unit y-vector `e₂`.
    #[inline]
    pub fn unit_y() -> Self {
        Self::new(T::ZERO, T::ONE, T::ZERO)
    }

    /// Creates the unit z-vector `e₃`.
    #[inline]
    pub fn unit_z() -> Self {
        Self::new(T::ZERO, T::ZERO, T::ONE)
    }

    /// Returns the squared magnitude: `x² + y² + z²`.
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns the magnitude (Euclidean norm): `√(x² + y² + z²)`.
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized (unit length) version of this vector.
    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.x / n, self.y / n, self.z / n)
    }

    /// Dot product (inner product): `a · b = a.x*b.x + a.y*b.y + a.z*b.z`.
    #[inline]
    pub fn dot(self, other: Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Wedge product (outer product): `a ∧ b`.
    #[inline]
    pub fn wedge(self, other: Self) -> Bivec3<T> {
        Bivec3 {
            xy: self.x * other.y - self.y * other.x,
            xz: self.x * other.z - self.z * other.x,
            yz: self.y * other.z - self.z * other.y,
        }
    }

    /// Cross product: `a × b`.
    #[inline]
    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Geometric product of two vectors: `ab = a·b + a∧b`.
    #[inline]
    pub fn geometric(self, other: Self) -> Even3<T> {
        Even3 {
            s: self.dot(other),
            b: self.wedge(other),
        }
    }
}

impl<T: Float> Default for Vec3<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// 3D bivector (grade 2): `e₁₂`, `e₁₃`, `e₂₃`.
///
/// Represents an oriented plane or rotation generator in 3D space.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Bivec3<T: Float> {
    /// Coefficient of `e₁₂` (xy-plane).
    pub xy: T,
    /// Coefficient of `e₁₃` (xz-plane).
    pub xz: T,
    /// Coefficient of `e₂₃` (yz-plane).
    pub yz: T,
}

impl<T: Float> Bivec3<T> {
    /// Creates a new bivector.
    #[inline]
    pub fn new(xy: T, xz: T, yz: T) -> Self {
        Self { xy, xz, yz }
    }

    /// Creates the zero bivector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(T::ZERO, T::ZERO, T::ZERO)
    }

    /// Creates the unit xy-plane bivector `e₁₂`.
    #[inline]
    pub fn unit_xy() -> Self {
        Self::new(T::ONE, T::ZERO, T::ZERO)
    }

    /// Creates the unit xz-plane bivector `e₁₃`.
    #[inline]
    pub fn unit_xz() -> Self {
        Self::new(T::ZERO, T::ONE, T::ZERO)
    }

    /// Creates the unit yz-plane bivector `e₂₃`.
    #[inline]
    pub fn unit_yz() -> Self {
        Self::new(T::ZERO, T::ZERO, T::ONE)
    }

    /// Returns the squared magnitude.
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.xy * self.xy + self.xz * self.xz + self.yz * self.yz
    }

    /// Returns the magnitude.
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized (unit) bivector.
    #[inline]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.xy / n, self.xz / n, self.yz / n)
    }

    /// Reverses the bivector: `B̃ = -B`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(-self.xy, -self.xz, -self.yz)
    }

    /// Computes the dual vector (Hodge star): `*B`.
    #[inline]
    pub fn dual(&self) -> Vec3<T> {
        Vec3::new(self.yz, -self.xz, self.xy)
    }
}

impl<T: Float> Default for Bivec3<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// 3D trivector/pseudoscalar (grade 3): `e₁₂₃`.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Trivec3<T: Float>(
    /// Coefficient of `e₁₂₃`.
    pub T,
);

impl<T: Float> Trivec3<T> {
    /// Creates a new trivector.
    #[inline]
    pub fn new(value: T) -> Self {
        Self(value)
    }

    /// Creates the zero trivector.
    #[inline]
    pub fn zero() -> Self {
        Self(T::ZERO)
    }

    /// Creates the unit pseudoscalar `e₁₂₃`.
    #[inline]
    pub fn unit() -> Self {
        Self(T::ONE)
    }

    /// Returns the coefficient.
    #[inline]
    pub fn value(&self) -> T {
        self.0
    }

    /// Reverses the trivector: `Ĩ = -I`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self(-self.0)
    }
}

impl<T: Float> Default for Trivec3<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// 3D rotor: scalar + bivector (even subalgebra).
///
/// Rotors represent rotations in geometric algebra. A unit rotor `R` with
/// `R R̃ = 1` rotates a vector `v` via the sandwich product: `v' = R̃ v R`.
/// This gives counterclockwise rotation when looking along the rotation axis.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Rotor3<T: Float> {
    /// Scalar part (grade 0).
    pub s: T,
    /// Bivector part (grade 2).
    pub b: Bivec3<T>,
}

impl<T: Float> Rotor3<T> {
    /// Creates a new rotor from scalar and bivector parts.
    #[inline]
    pub fn new(s: T, b: Bivec3<T>) -> Self {
        Self { s, b }
    }

    /// Creates the identity rotor (no rotation).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::ONE,
            b: Bivec3::zero(),
        }
    }

    /// Creates a rotor from an angle and rotation plane (unit bivector).
    ///
    /// The rotor `R = cos(θ/2) + sin(θ/2)B` rotates by angle `θ` in plane `B`.
    #[inline]
    pub fn from_angle_plane(angle: T, plane: Bivec3<T>) -> Self {
        let half = angle / T::TWO;
        let cos_half = half.cos();
        let sin_half = half.sin();
        Self {
            s: cos_half,
            b: Bivec3::new(
                plane.xy * sin_half,
                plane.xz * sin_half,
                plane.yz * sin_half,
            ),
        }
    }

    /// Creates a rotor that rotates vector `a` to vector `b`.
    #[inline]
    pub fn from_vectors(a: Vec3<T>, b: Vec3<T>) -> Self {
        let dot = a.dot(b);
        let wedge = a.wedge(b); // a ∧ b

        let sum_sq = (T::ONE + dot) * (T::ONE + dot) + wedge.norm_squared();

        if sum_sq < T::EPSILON {
            // Vectors are anti-parallel, need to find perpendicular axis
            let perp = if a.x.abs() < a.y.abs() && a.x.abs() < a.z.abs() {
                Vec3::unit_x()
            } else if a.y.abs() < a.z.abs() {
                Vec3::unit_y()
            } else {
                Vec3::unit_z()
            };
            let axis = a.cross(perp).normalized();
            let plane = a.wedge(axis).normalized();
            return Self::from_angle_plane(T::PI, plane);
        }

        let norm = sum_sq.sqrt();
        Self {
            s: (T::ONE + dot) / norm,
            b: Bivec3::new(wedge.xy / norm, wedge.xz / norm, wedge.yz / norm),
        }
    }

    /// Returns the squared magnitude: `R R̃ = s² + |b|²`.
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.s * self.s + self.b.norm_squared()
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
            b: Bivec3::new(self.b.xy / n, self.b.xz / n, self.b.yz / n),
        }
    }

    /// Returns the reverse: `R̃ = s - b`.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self {
            s: self.s,
            b: self.b.reverse(),
        }
    }

    /// Returns the inverse rotor: `R⁻¹ = R̃ / |R|²`.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        Self {
            s: self.s / norm_sq,
            b: Bivec3::new(
                -self.b.xy / norm_sq,
                -self.b.xz / norm_sq,
                -self.b.yz / norm_sq,
            ),
        }
    }

    /// Applies this rotation to a vector: `v' = R̃ v R`.
    ///
    /// Note: Some sources use `R v R̃`. We use `R̃ v R` which gives counterclockwise
    /// rotation when looking along the rotation axis.
    #[inline]
    pub fn rotate(&self, v: Vec3<T>) -> Vec3<T> {
        let s = self.s;
        let bxy = self.b.xy;
        let bxz = self.b.xz;
        let byz = self.b.yz;

        // Intermediate: q = R̃ v (note: R̃ = s - b, so bivector signs flip)
        let qx = s * v.x - bxy * v.y - bxz * v.z;
        let qy = s * v.y + bxy * v.x - byz * v.z;
        let qz = s * v.z + bxz * v.x + byz * v.y;
        let qt = -bxy * v.z + bxz * v.y - byz * v.x;

        // Result: q R (only vector part survives)
        Vec3 {
            x: s * qx - bxy * qy - bxz * qz - byz * qt,
            y: s * qy + bxy * qx + bxz * qt - byz * qz,
            z: s * qz - bxy * qt + bxz * qx + byz * qy,
        }
    }

    /// Composes two rotations: `R₂ ∘ R₁ = R₂ R₁`.
    #[inline]
    pub fn compose(&self, other: Self) -> Self {
        let s1 = self.s;
        let s2 = other.s;
        let b1 = &self.b;
        let b2 = &other.b;

        let s = s1 * s2 - (b1.xy * b2.xy + b1.xz * b2.xz + b1.yz * b2.yz);

        let bxy = s1 * b2.xy + s2 * b1.xy + (b1.xz * b2.yz - b1.yz * b2.xz);
        let bxz = s1 * b2.xz + s2 * b1.xz + (b1.yz * b2.xy - b1.xy * b2.yz);
        let byz = s1 * b2.yz + s2 * b1.yz + (b1.xy * b2.xz - b1.xz * b2.xy);

        Self {
            s,
            b: Bivec3::new(bxy, bxz, byz),
        }
    }

    /// Spherical linear interpolation between two rotors.
    #[inline]
    pub fn slerp(&self, other: Self, t: T) -> Self {
        let dot = self.s * other.s
            + self.b.xy * other.b.xy
            + self.b.xz * other.b.xz
            + self.b.yz * other.b.yz;

        let dot = if dot > T::ONE {
            T::ONE
        } else if dot < -T::ONE {
            -T::ONE
        } else {
            dot
        };

        let theta = dot.acos();

        if theta.abs() < T::EPSILON {
            return Self {
                s: self.s * (T::ONE - t) + other.s * t,
                b: Bivec3::new(
                    self.b.xy * (T::ONE - t) + other.b.xy * t,
                    self.b.xz * (T::ONE - t) + other.b.xz * t,
                    self.b.yz * (T::ONE - t) + other.b.yz * t,
                ),
            }
            .normalized();
        }

        let sin_theta = theta.sin();
        let s1 = ((T::ONE - t) * theta).sin() / sin_theta;
        let s2 = (t * theta).sin() / sin_theta;

        Self {
            s: self.s * s1 + other.s * s2,
            b: Bivec3::new(
                self.b.xy * s1 + other.b.xy * s2,
                self.b.xz * s1 + other.b.xz * s2,
                self.b.yz * s1 + other.b.yz * s2,
            ),
        }
    }
}

impl<T: Float> Default for Rotor3<T> {
    fn default() -> Self {
        Self::identity()
    }
}

/// Even subalgebra element (scalar + bivector).
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Even3<T: Float> {
    /// Scalar part (grade 0).
    pub s: T,
    /// Bivector part (grade 2).
    pub b: Bivec3<T>,
}

impl<T: Float> Even3<T> {
    /// Creates a new even multivector.
    #[inline]
    pub fn new(s: T, b: Bivec3<T>) -> Self {
        Self { s, b }
    }

    /// Creates zero.
    #[inline]
    pub fn zero() -> Self {
        Self {
            s: T::ZERO,
            b: Bivec3::zero(),
        }
    }

    /// Converts to a rotor (same representation).
    #[inline]
    pub fn to_rotor(self) -> Rotor3<T> {
        Rotor3 {
            s: self.s,
            b: self.b,
        }
    }
}

impl<T: Float> Default for Even3<T> {
    fn default() -> Self {
        Self::zero()
    }
}

// ============================================================================
// approx trait implementations (generic over Float)
// ============================================================================

impl<T: Float> AbsDiffEq for Vec3<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.x, &other.x, epsilon)
            && T::abs_diff_eq(&self.y, &other.y, epsilon)
            && T::abs_diff_eq(&self.z, &other.z, epsilon)
    }
}

impl<T: Float> RelativeEq for Vec3<T> {
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
            && T::relative_eq(&self.z, &other.z, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Vec3<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.x, &other.x, epsilon, max_ulps)
            && T::ulps_eq(&self.y, &other.y, epsilon, max_ulps)
            && T::ulps_eq(&self.z, &other.z, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Bivec3<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.xy, &other.xy, epsilon)
            && T::abs_diff_eq(&self.xz, &other.xz, epsilon)
            && T::abs_diff_eq(&self.yz, &other.yz, epsilon)
    }
}

impl<T: Float> RelativeEq for Bivec3<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.xy, &other.xy, epsilon, max_relative)
            && T::relative_eq(&self.xz, &other.xz, epsilon, max_relative)
            && T::relative_eq(&self.yz, &other.yz, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Bivec3<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.xy, &other.xy, epsilon, max_ulps)
            && T::ulps_eq(&self.xz, &other.xz, epsilon, max_ulps)
            && T::ulps_eq(&self.yz, &other.yz, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Trivec3<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.0, &other.0, epsilon)
    }
}

impl<T: Float> RelativeEq for Trivec3<T> {
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

impl<T: Float> UlpsEq for Trivec3<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.0, &other.0, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Rotor3<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.s, &other.s, epsilon)
            && Bivec3::abs_diff_eq(&self.b, &other.b, epsilon)
    }
}

impl<T: Float> RelativeEq for Rotor3<T> {
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
            && Bivec3::relative_eq(&self.b, &other.b, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Rotor3<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.s, &other.s, epsilon, max_ulps)
            && Bivec3::ulps_eq(&self.b, &other.b, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Even3<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.s, &other.s, epsilon)
            && Bivec3::abs_diff_eq(&self.b, &other.b, epsilon)
    }
}

impl<T: Float> RelativeEq for Even3<T> {
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
            && Bivec3::relative_eq(&self.b, &other.b, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Even3<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.s, &other.s, epsilon, max_ulps)
            && Bivec3::ulps_eq(&self.b, &other.b, epsilon, max_ulps)
    }
}
