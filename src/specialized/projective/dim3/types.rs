//! Type definitions for 3D Projective Geometric Algebra.
//!
//! # Convention Choices
//!
//! - **Basis ordering**: e₁, e₂, e₃, e₀ (Euclidean first, null last)
//! - **Rotation positive direction**: Counterclockwise when looking down the axis
//!   (right-hand rule)
//! - **Motor composition**: `a.compose(&b)` applies `a` first, then `b`
//! - **nalgebra correspondence**: `iso1 * iso2` applies `iso2` first, then `iso1`,
//!   so `motor2.compose(&motor1)` ≡ `iso1 * iso2`

use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::scalar::Float;

/// A point in 3D PGA (grade 1 vector).
///
/// In point-based PGA, a point is represented in homogeneous coordinates:
/// `P = x·e₁ + y·e₂ + z·e₃ + w·e₀`
///
/// where `(x/w, y/w, z/w)` are the Cartesian coordinates when `w ≠ 0`.
/// Points with `w = 0` represent ideal points (points at infinity).
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim3::Point;
///
/// // Point at (3, 4, 5)
/// let p = Point::new(3.0, 4.0, 5.0);
/// assert_eq!(p.x(), 3.0);
/// assert_eq!(p.y(), 4.0);
/// assert_eq!(p.z(), 5.0);
///
/// // Point at infinity in direction (1, 0, 0)
/// let ideal = Point::ideal(1.0, 0.0, 0.0);
/// assert!(ideal.is_ideal(1e-10));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Point<T: Float> {
    /// Coefficient of e₁ (x-coordinate × w).
    pub e1: T,
    /// Coefficient of e₂ (y-coordinate × w).
    pub e2: T,
    /// Coefficient of e₃ (z-coordinate × w).
    pub e3: T,
    /// Coefficient of e₀ (homogeneous weight).
    pub e0: T,
}

impl<T: Float> Point<T> {
    /// Creates a finite point at Cartesian coordinates (x, y, z).
    ///
    /// The homogeneous weight `w` is set to 1.
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        Self {
            e1: x,
            e2: y,
            e3: z,
            e0: T::one(),
        }
    }

    /// Creates a point from homogeneous coordinates.
    #[inline]
    pub fn from_homogeneous(e1: T, e2: T, e3: T, e0: T) -> Self {
        Self { e1, e2, e3, e0 }
    }

    /// Creates an ideal point (point at infinity) in the given direction.
    ///
    /// Ideal points have `w = 0` and represent directions rather than positions.
    #[inline]
    pub fn ideal(dx: T, dy: T, dz: T) -> Self {
        Self {
            e1: dx,
            e2: dy,
            e3: dz,
            e0: T::zero(),
        }
    }

    /// Origin point (0, 0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Returns the x-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn x(&self) -> T {
        self.e1 / self.e0
    }

    /// Returns the y-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn y(&self) -> T {
        self.e2 / self.e0
    }

    /// Returns the z-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn z(&self) -> T {
        self.e3 / self.e0
    }

    /// Returns the homogeneous weight.
    #[inline]
    pub fn w(&self) -> T {
        self.e0
    }

    /// Returns true if this is an ideal point (point at infinity).
    #[inline]
    pub fn is_ideal(&self, epsilon: T) -> bool {
        self.e0.abs() < epsilon
    }

    /// Returns true if this is a finite point (not at infinity).
    #[inline]
    pub fn is_finite(&self, epsilon: T) -> bool {
        self.e0.abs() >= epsilon
    }

    /// Normalizes the homogeneous coordinates so `w = 1` (if finite).
    ///
    /// Returns `None` if this is an ideal point.
    pub fn normalize(&self) -> Option<Self> {
        if self.e0.abs() < T::epsilon() {
            None
        } else {
            Some(Self {
                e1: self.e1 / self.e0,
                e2: self.e2 / self.e0,
                e3: self.e3 / self.e0,
                e0: T::one(),
            })
        }
    }

    /// Returns the Cartesian coordinates as a tuple, if finite.
    ///
    /// Returns `None` if this is an ideal point.
    #[inline]
    pub fn to_cartesian(&self) -> Option<(T, T, T)> {
        if self.e0.abs() < T::epsilon() {
            None
        } else {
            Some((self.e1 / self.e0, self.e2 / self.e0, self.e3 / self.e0))
        }
    }
}

impl<T: Float> Default for Point<T> {
    fn default() -> Self {
        Self::origin()
    }
}

/// A motor (rigid transformation) in 3D PGA.
///
/// Motors represent rigid body transformations (rotation + translation).
/// In 3D PGA, a motor is an even-grade element with 8 components:
/// - Scalar: `s`
/// - Bivector (6 components): `e₂₃, e₃₁, e₁₂` (rotation) and `e₀₁, e₀₂, e₀₃` (translation)
/// - Pseudoscalar: `e₀₁₂₃`
///
/// A motor transforms geometric objects via the sandwich product:
/// `X' = M X M̃`
///
/// # Construction
///
/// - [`Motor::identity`]: No transformation
/// - [`Motor::from_translation`]: Pure translation
/// - [`Motor::from_rotation_x`], [`Motor::from_rotation_y`], [`Motor::from_rotation_z`]:
///   Rotations around coordinate axes
/// - [`Motor::from_axis_angle`]: Rotation around arbitrary axis
/// - [`Motor::compose`]: Combine transformations
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim3::{Point, Motor};
/// use std::f64::consts::FRAC_PI_2;
/// use approx::abs_diff_eq;
///
/// // 90° rotation around Z axis
/// let rotation = Motor::from_rotation_z(FRAC_PI_2);
///
/// // Translation by (1, 2, 3)
/// let translation = Motor::from_translation(1.0, 2.0, 3.0);
///
/// // Transform a point
/// let p = Point::new(1.0, 0.0, 0.0);
/// let rotated = rotation.transform_point(&p);
/// assert!(abs_diff_eq!(rotated.x(), 0.0, epsilon = 1e-10));
/// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
/// assert!(abs_diff_eq!(rotated.z(), 0.0, epsilon = 1e-10));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Motor<T: Float> {
    /// Scalar part (cos(θ/2) for pure rotation).
    pub s: T,
    /// Coefficient of e₂₃ (rotation around x-axis).
    pub e23: T,
    /// Coefficient of e₃₁ (rotation around y-axis).
    pub e31: T,
    /// Coefficient of e₁₂ (rotation around z-axis).
    pub e12: T,
    /// Coefficient of e₀₁ (translation in x).
    pub e01: T,
    /// Coefficient of e₀₂ (translation in y).
    pub e02: T,
    /// Coefficient of e₀₃ (translation in z).
    pub e03: T,
    /// Coefficient of e₀₁₂₃ (pseudoscalar part).
    pub e0123: T,
}

impl<T: Float> Motor<T> {
    /// Creates a motor from all components.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub fn new(s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T, e0123: T) -> Self {
        Self {
            s,
            e23,
            e31,
            e12,
            e01,
            e02,
            e03,
            e0123,
        }
    }

    /// Creates the identity motor (no transformation).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::one(),
            e23: T::zero(),
            e31: T::zero(),
            e12: T::zero(),
            e01: T::zero(),
            e02: T::zero(),
            e03: T::zero(),
            e0123: T::zero(),
        }
    }

    /// Creates a pure translation motor.
    ///
    /// # Arguments
    ///
    /// * `dx` - Translation in x direction
    /// * `dy` - Translation in y direction
    /// * `dz` - Translation in z direction
    #[inline]
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self {
        // Translation motor: T = 1 + (d/2)·(dx·e₀₁ + dy·e₀₂ + dz·e₀₃)
        // The factor of 1/2 comes from the sandwich product doubling the translation
        Self {
            s: T::one(),
            e23: T::zero(),
            e31: T::zero(),
            e12: T::zero(),
            e01: dx / T::TWO,
            e02: dy / T::TWO,
            e03: dz / T::TWO,
            e0123: T::zero(),
        }
    }

    /// Creates a pure rotation motor around the Z axis.
    ///
    /// # Arguments
    ///
    /// * `angle` - Rotation angle in radians (counterclockwise when looking down Z)
    #[inline]
    pub fn from_rotation_z(angle: T) -> Self {
        let half = angle / T::TWO;
        Self {
            s: half.cos(),
            e23: T::zero(),
            e31: T::zero(),
            e12: half.sin(),
            e01: T::zero(),
            e02: T::zero(),
            e03: T::zero(),
            e0123: T::zero(),
        }
    }

    /// Creates a pure rotation motor around the X axis.
    ///
    /// # Arguments
    ///
    /// * `angle` - Rotation angle in radians (counterclockwise when looking down X)
    #[inline]
    pub fn from_rotation_x(angle: T) -> Self {
        let half = angle / T::TWO;
        Self {
            s: half.cos(),
            e23: half.sin(),
            e31: T::zero(),
            e12: T::zero(),
            e01: T::zero(),
            e02: T::zero(),
            e03: T::zero(),
            e0123: T::zero(),
        }
    }

    /// Creates a pure rotation motor around the Y axis.
    ///
    /// # Arguments
    ///
    /// * `angle` - Rotation angle in radians (counterclockwise when looking down Y)
    #[inline]
    pub fn from_rotation_y(angle: T) -> Self {
        let half = angle / T::TWO;
        Self {
            s: half.cos(),
            e23: T::zero(),
            e31: half.sin(),
            e12: T::zero(),
            e01: T::zero(),
            e02: T::zero(),
            e03: T::zero(),
            e0123: T::zero(),
        }
    }

    /// Creates a pure rotation motor around an arbitrary axis through the origin.
    ///
    /// # Arguments
    ///
    /// * `axis` - The rotation axis as (x, y, z). Must be normalized (unit length).
    /// * `angle` - Rotation angle in radians (counterclockwise when looking down the axis)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Point, Motor};
    /// use approx::abs_diff_eq;
    ///
    /// // Rotate around Z axis (same as from_rotation_z)
    /// let motor = Motor::from_axis_angle((0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
    /// let p = Point::new(1.0, 0.0, 0.0);
    /// let result = motor.transform_point(&p);
    ///
    /// assert!(abs_diff_eq!(result.x(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.y(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.z(), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn from_axis_angle(axis: (T, T, T), angle: T) -> Self {
        // Rotation motor: R = cos(θ/2) + sin(θ/2)·(ax·e₂₃ + ay·e₃₁ + az·e₁₂)
        // The axis components map to bivector components:
        // - x-axis rotation -> e₂₃ bivector
        // - y-axis rotation -> e₃₁ bivector
        // - z-axis rotation -> e₁₂ bivector
        let half = angle / T::TWO;
        let sin_half = half.sin();
        let (ax, ay, az) = axis;
        Self {
            s: half.cos(),
            e23: ax * sin_half,
            e31: ay * sin_half,
            e12: az * sin_half,
            e01: T::zero(),
            e02: T::zero(),
            e03: T::zero(),
            e0123: T::zero(),
        }
    }

    /// Composes two motors via geometric product: `self * other`.
    ///
    /// The composition order follows PGA convention where `self` is applied first,
    /// then `other`. So `a.compose(&b)` means "first apply transformation `a`, then `b`".
    ///
    /// **nalgebra correspondence**: nalgebra's `iso1 * iso2` applies `iso2` first, then `iso1`.
    /// So `motor2.compose(&motor1)` is equivalent to `iso1 * iso2`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Point, Motor};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// // First rotate 90° around Z, then translate by (1, 2, 3)
    /// let rotation = Motor::from_rotation_z(FRAC_PI_2);
    /// let translation = Motor::from_translation(1.0, 2.0, 3.0);
    /// let combined = rotation.compose(&translation);
    ///
    /// // (1,0,0) -> rotate to (0,1,0) -> translate to (1,3,3)
    /// let p = Point::new(1.0, 0.0, 0.0);
    /// let result = combined.transform_point(&p);
    ///
    /// assert!(abs_diff_eq!(result.x(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.y(), 3.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.z(), 3.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        // Motor multiplication in 3D PGA
        // M1 = s1 + B1 + I1 (scalar + bivector + pseudoscalar)
        // M2 = s2 + B2 + I2
        //
        // The even subalgebra multiplication follows specific rules.
        // We need to compute all grade-0, grade-2, and grade-4 parts.

        let s1 = self.s;
        let b23_1 = self.e23;
        let b31_1 = self.e31;
        let b12_1 = self.e12;
        let b01_1 = self.e01;
        let b02_1 = self.e02;
        let b03_1 = self.e03;
        let i1 = self.e0123;

        let s2 = other.s;
        let b23_2 = other.e23;
        let b31_2 = other.e31;
        let b12_2 = other.e12;
        let b01_2 = other.e01;
        let b02_2 = other.e02;
        let b03_2 = other.e03;
        let i2 = other.e0123;

        // Scalar part: s1*s2 - B1·B2
        // The Euclidean bivectors e23, e31, e12 square to -1, so their product contributes negatively
        // The null bivectors e01, e02, e03 square to 0 and don't contribute to scalar
        let s = s1 * s2 - b23_1 * b23_2 - b31_1 * b31_2 - b12_1 * b12_2;

        // Euclidean bivector parts: scalar*bivector + bivector*scalar + bivector×bivector
        //
        // The bivector multiplication rules in Cl(3,0,1) are:
        //   e23*e31 = -e12,  e31*e23 = +e12
        //   e31*e12 = -e23,  e12*e31 = +e23
        //   e12*e23 = -e31,  e23*e12 = +e31
        //
        // Therefore:
        //   e23 = s1*b23_2 + b23_1*s2 + b12_1*b31_2 - b31_1*b12_2
        //   e31 = s1*b31_2 + b31_1*s2 + b23_1*b12_2 - b12_1*b23_2
        //   e12 = s1*b12_2 + b12_1*s2 + b31_1*b23_2 - b23_1*b31_2
        let e23 = s1 * b23_2 + b23_1 * s2 + b12_1 * b31_2 - b31_1 * b12_2;
        let e31 = s1 * b31_2 + b31_1 * s2 + b23_1 * b12_2 - b12_1 * b23_2;
        let e12 = s1 * b12_2 + b12_1 * s2 + b31_1 * b23_2 - b23_1 * b31_2;

        // Combining all contributions for e01:
        // Note: e23*e0123 = e0123*e23 = -e01, so the pseudoscalar terms are negated
        let e01_new = s1 * b01_2 + b01_1 * s2 - b23_1 * i2 - i1 * b23_2 + b12_1 * b02_2
            - b02_1 * b12_2
            - b31_1 * b03_2
            + b03_1 * b31_2;

        // e02: similar pattern
        // Note: e31*e0123 = e0123*e31 = -e02, so the pseudoscalar terms are negated
        let e02_new = s1 * b02_2 + b02_1 * s2 - b31_1 * i2 - i1 * b31_2 + b23_1 * b03_2
            - b03_1 * b23_2
            - b12_1 * b01_2
            + b01_1 * b12_2;

        // e03: similar pattern
        // Note: e12*e0123 = e0123*e12 = -e03, so the pseudoscalar terms are negated
        let e03_new = s1 * b03_2 + b03_1 * s2 - b12_1 * i2 - i1 * b12_2 + b31_1 * b01_2
            - b01_1 * b31_2
            - b23_1 * b02_2
            + b02_1 * b23_2;

        // Pseudoscalar part: from bivector*bivector products that give grade 4
        // e23*e01 + e31*e02 + e12*e03 = e0123 (wedge contributions)
        // Plus s*I terms
        let e0123_new = s1 * i2
            + i1 * s2
            + b23_1 * b01_2
            + b01_1 * b23_2
            + b31_1 * b02_2
            + b02_1 * b31_2
            + b12_1 * b03_2
            + b03_1 * b12_2;

        Self {
            s,
            e23,
            e31,
            e12,
            e01: e01_new,
            e02: e02_new,
            e03: e03_new,
            e0123: e0123_new,
        }
    }

    /// Returns the reverse of the motor: `M̃`.
    ///
    /// For motors, the reverse negates the bivector and pseudoscalar parts.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self {
            s: self.s,
            e23: -self.e23,
            e31: -self.e31,
            e12: -self.e12,
            e01: -self.e01,
            e02: -self.e02,
            e03: -self.e03,
            e0123: -self.e0123,
        }
    }
}

impl<T: Float> Default for Motor<T> {
    fn default() -> Self {
        Self::identity()
    }
}

// ============================================================================
// approx trait implementations
// ============================================================================

impl<T: Float> AbsDiffEq for Point<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.e1, &other.e1, epsilon)
            && T::abs_diff_eq(&self.e2, &other.e2, epsilon)
            && T::abs_diff_eq(&self.e3, &other.e3, epsilon)
            && T::abs_diff_eq(&self.e0, &other.e0, epsilon)
    }
}

impl<T: Float> RelativeEq for Point<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.e1, &other.e1, epsilon, max_relative)
            && T::relative_eq(&self.e2, &other.e2, epsilon, max_relative)
            && T::relative_eq(&self.e3, &other.e3, epsilon, max_relative)
            && T::relative_eq(&self.e0, &other.e0, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Point<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.e1, &other.e1, epsilon, max_ulps)
            && T::ulps_eq(&self.e2, &other.e2, epsilon, max_ulps)
            && T::ulps_eq(&self.e3, &other.e3, epsilon, max_ulps)
            && T::ulps_eq(&self.e0, &other.e0, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Motor<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.s, &other.s, epsilon)
            && T::abs_diff_eq(&self.e23, &other.e23, epsilon)
            && T::abs_diff_eq(&self.e31, &other.e31, epsilon)
            && T::abs_diff_eq(&self.e12, &other.e12, epsilon)
            && T::abs_diff_eq(&self.e01, &other.e01, epsilon)
            && T::abs_diff_eq(&self.e02, &other.e02, epsilon)
            && T::abs_diff_eq(&self.e03, &other.e03, epsilon)
            && T::abs_diff_eq(&self.e0123, &other.e0123, epsilon)
    }
}

impl<T: Float> RelativeEq for Motor<T> {
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
            && T::relative_eq(&self.e23, &other.e23, epsilon, max_relative)
            && T::relative_eq(&self.e31, &other.e31, epsilon, max_relative)
            && T::relative_eq(&self.e12, &other.e12, epsilon, max_relative)
            && T::relative_eq(&self.e01, &other.e01, epsilon, max_relative)
            && T::relative_eq(&self.e02, &other.e02, epsilon, max_relative)
            && T::relative_eq(&self.e03, &other.e03, epsilon, max_relative)
            && T::relative_eq(&self.e0123, &other.e0123, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Motor<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.s, &other.s, epsilon, max_ulps)
            && T::ulps_eq(&self.e23, &other.e23, epsilon, max_ulps)
            && T::ulps_eq(&self.e31, &other.e31, epsilon, max_ulps)
            && T::ulps_eq(&self.e12, &other.e12, epsilon, max_ulps)
            && T::ulps_eq(&self.e01, &other.e01, epsilon, max_ulps)
            && T::ulps_eq(&self.e02, &other.e02, epsilon, max_ulps)
            && T::ulps_eq(&self.e03, &other.e03, epsilon, max_ulps)
            && T::ulps_eq(&self.e0123, &other.e0123, epsilon, max_ulps)
    }
}

/// A plane in 3D PGA (grade 3 trivector).
///
/// In point-based PGA, a plane is represented as a trivector:
/// `g = nx·e₀₂₃ + ny·e₀₃₁ + nz·e₀₁₂ + d·e₁₂₃`
///
/// This represents the implicit plane `nx·x + ny·y + nz·z + d = 0`
/// where `(nx, ny, nz)` is the normal direction.
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim3::Plane;
///
/// // XY plane (z = 0): normal (0, 0, 1), d = 0
/// let xy_plane = Plane::from_normal_and_distance(0.0, 0.0, 1.0, 0.0);
///
/// // Plane z = 5: normal (0, 0, 1), d = -5
/// let offset_plane = Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Plane<T: Float> {
    /// Coefficient of e₀₂₃ (x-component of normal).
    pub e023: T,
    /// Coefficient of e₀₃₁ (y-component of normal).
    pub e031: T,
    /// Coefficient of e₀₁₂ (z-component of normal).
    pub e012: T,
    /// Coefficient of e₁₂₃ (signed distance from origin, scaled by normal length).
    pub e123: T,
}

impl<T: Float> Plane<T> {
    /// Creates a plane from trivector coefficients.
    #[inline]
    pub fn new(e023: T, e031: T, e012: T, e123: T) -> Self {
        Self {
            e023,
            e031,
            e012,
            e123,
        }
    }

    /// Creates a plane from normal direction and signed distance.
    ///
    /// The plane equation is `nx·x + ny·y + nz·z + d = 0`.
    ///
    /// # Arguments
    ///
    /// * `nx, ny, nz` - Normal direction (should be unit length for geometric interpretation)
    /// * `d` - Signed distance parameter (negative of the distance from origin along normal)
    #[inline]
    pub fn from_normal_and_distance(nx: T, ny: T, nz: T, d: T) -> Self {
        Self {
            e023: nx,
            e031: ny,
            e012: nz,
            e123: d,
        }
    }

    /// Creates the XY plane (z = 0).
    #[inline]
    pub fn xy() -> Self {
        Self::from_normal_and_distance(T::zero(), T::zero(), T::one(), T::zero())
    }

    /// Creates the XZ plane (y = 0).
    #[inline]
    pub fn xz() -> Self {
        Self::from_normal_and_distance(T::zero(), T::one(), T::zero(), T::zero())
    }

    /// Creates the YZ plane (x = 0).
    #[inline]
    pub fn yz() -> Self {
        Self::from_normal_and_distance(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Returns the normal direction `(nx, ny, nz)`.
    #[inline]
    pub fn normal(&self) -> (T, T, T) {
        (self.e023, self.e031, self.e012)
    }

    /// Returns the signed distance parameter `d`.
    #[inline]
    pub fn distance(&self) -> T {
        self.e123
    }

    /// Returns the squared norm of the plane's weight (normal direction).
    #[inline]
    pub fn weight_norm_squared(&self) -> T {
        self.e023 * self.e023 + self.e031 * self.e031 + self.e012 * self.e012
    }

    /// Returns the weight norm (length of normal direction).
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    /// Returns a unitized plane (unit normal).
    pub fn unitized(&self) -> Self {
        let n = self.weight_norm();
        if n.abs() < T::epsilon() {
            *self
        } else {
            Self {
                e023: self.e023 / n,
                e031: self.e031 / n,
                e012: self.e012 / n,
                e123: self.e123 / n,
            }
        }
    }
}

impl<T: Float> Default for Plane<T> {
    fn default() -> Self {
        Self::xy()
    }
}

/// A flector (improper isometry) in 3D PGA.
///
/// Flectors represent orientation-reversing transformations: reflections,
/// glide reflections, rotoreflections, and inversions. A flector is an
/// odd-grade element combining a point (grade 1) and plane (grade 3):
///
/// `F = (px·e₁ + py·e₂ + pz·e₃ + pw·e₀) + (gx·e₀₂₃ + gy·e₀₃₁ + gz·e₀₁₂ + gw·e₁₂₃)`
///
/// The geometric constraint requires the point to lie in the plane.
///
/// A flector transforms geometric objects via the sandwich product with sign:
/// `X' = F X F̃` (note: some conventions use `X' = -F X F̃`)
///
/// # Construction
///
/// - [`Flector::from_plane`]: Pure reflection through a plane
/// - [`Flector::from_plane_through_origin`]: Reflection through plane containing origin
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim3::{Point, Plane, Flector};
/// use approx::abs_diff_eq;
///
/// // Reflection through the XY plane (z = 0)
/// let plane = Plane::xy();
/// let flector = Flector::from_plane(plane);
///
/// // Reflect a point above the plane
/// let p = Point::new(1.0, 2.0, 3.0);
/// let reflected = flector.transform_point(&p);
///
/// // z-coordinate is negated
/// assert!(abs_diff_eq!(reflected.x(), 1.0, epsilon = 1e-10));
/// assert!(abs_diff_eq!(reflected.y(), 2.0, epsilon = 1e-10));
/// assert!(abs_diff_eq!(reflected.z(), -3.0, epsilon = 1e-10));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Flector<T: Float> {
    // Grade 1 (point) components
    /// Coefficient of e₁.
    pub e1: T,
    /// Coefficient of e₂.
    pub e2: T,
    /// Coefficient of e₃.
    pub e3: T,
    /// Coefficient of e₀.
    pub e0: T,

    // Grade 3 (plane) components
    /// Coefficient of e₀₂₃.
    pub e023: T,
    /// Coefficient of e₀₃₁.
    pub e031: T,
    /// Coefficient of e₀₁₂.
    pub e012: T,
    /// Coefficient of e₁₂₃.
    pub e123: T,
}

impl<T: Float> Flector<T> {
    /// Creates a flector from all components.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub fn new(e1: T, e2: T, e3: T, e0: T, e023: T, e031: T, e012: T, e123: T) -> Self {
        Self {
            e1,
            e2,
            e3,
            e0,
            e023,
            e031,
            e012,
            e123,
        }
    }

    /// Creates a pure reflection flector from a plane.
    ///
    /// The resulting flector reflects points through the given plane.
    #[inline]
    pub fn from_plane(plane: Plane<T>) -> Self {
        Self {
            e1: T::zero(),
            e2: T::zero(),
            e3: T::zero(),
            e0: T::zero(),
            e023: plane.e023,
            e031: plane.e031,
            e012: plane.e012,
            e123: plane.e123,
        }
    }

    /// Creates a reflection through a plane passing through the origin.
    ///
    /// # Arguments
    ///
    /// * `nx, ny, nz` - Unit normal of the reflection plane
    #[inline]
    pub fn from_plane_through_origin(nx: T, ny: T, nz: T) -> Self {
        Self::from_plane(Plane::from_normal_and_distance(nx, ny, nz, T::zero()))
    }

    /// Creates a reflection through the XY plane (z = 0).
    #[inline]
    pub fn reflect_xy() -> Self {
        Self::from_plane(Plane::xy())
    }

    /// Creates a reflection through the XZ plane (y = 0).
    #[inline]
    pub fn reflect_xz() -> Self {
        Self::from_plane(Plane::xz())
    }

    /// Creates a reflection through the YZ plane (x = 0).
    #[inline]
    pub fn reflect_yz() -> Self {
        Self::from_plane(Plane::yz())
    }

    /// Returns the reverse of the flector: `F̃`.
    ///
    /// For odd-grade elements, the reverse negates grade-3 parts.
    /// Grade 1 is unchanged, grade 3 is negated.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self {
            e1: self.e1,
            e2: self.e2,
            e3: self.e3,
            e0: self.e0,
            e023: -self.e023,
            e031: -self.e031,
            e012: -self.e012,
            e123: -self.e123,
        }
    }

    /// Returns the point (grade 1) part of the flector.
    #[inline]
    pub fn point_part(&self) -> Point<T> {
        Point::from_homogeneous(self.e1, self.e2, self.e3, self.e0)
    }

    /// Returns the plane (grade 3) part of the flector.
    #[inline]
    pub fn plane_part(&self) -> Plane<T> {
        Plane::new(self.e023, self.e031, self.e012, self.e123)
    }

    /// Returns true if this is a pure plane reflection (no point component).
    #[inline]
    pub fn is_pure_reflection(&self, epsilon: T) -> bool {
        self.e1.abs() < epsilon
            && self.e2.abs() < epsilon
            && self.e3.abs() < epsilon
            && self.e0.abs() < epsilon
    }

    /// Returns the weight norm squared of the flector.
    ///
    /// For a unitized flector, this equals 1.
    #[inline]
    pub fn weight_norm_squared(&self) -> T {
        // Weight norm² = g_x² + g_y² + g_z² + p_w²
        self.e023 * self.e023 + self.e031 * self.e031 + self.e012 * self.e012 + self.e0 * self.e0
    }

    /// Returns the weight norm.
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    /// Returns a unitized flector (weight norm = 1).
    pub fn unitized(&self) -> Self {
        let n = self.weight_norm();
        if n.abs() < T::epsilon() {
            *self
        } else {
            Self {
                e1: self.e1 / n,
                e2: self.e2 / n,
                e3: self.e3 / n,
                e0: self.e0 / n,
                e023: self.e023 / n,
                e031: self.e031 / n,
                e012: self.e012 / n,
                e123: self.e123 / n,
            }
        }
    }
}

impl<T: Float> Default for Flector<T> {
    fn default() -> Self {
        Self::reflect_xy()
    }
}

// ============================================================================
// Plane approx trait implementations
// ============================================================================

impl<T: Float> AbsDiffEq for Plane<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.e023, &other.e023, epsilon)
            && T::abs_diff_eq(&self.e031, &other.e031, epsilon)
            && T::abs_diff_eq(&self.e012, &other.e012, epsilon)
            && T::abs_diff_eq(&self.e123, &other.e123, epsilon)
    }
}

impl<T: Float> RelativeEq for Plane<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.e023, &other.e023, epsilon, max_relative)
            && T::relative_eq(&self.e031, &other.e031, epsilon, max_relative)
            && T::relative_eq(&self.e012, &other.e012, epsilon, max_relative)
            && T::relative_eq(&self.e123, &other.e123, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Plane<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.e023, &other.e023, epsilon, max_ulps)
            && T::ulps_eq(&self.e031, &other.e031, epsilon, max_ulps)
            && T::ulps_eq(&self.e012, &other.e012, epsilon, max_ulps)
            && T::ulps_eq(&self.e123, &other.e123, epsilon, max_ulps)
    }
}

// ============================================================================
// Flector approx trait implementations
// ============================================================================

impl<T: Float> AbsDiffEq for Flector<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.e1, &other.e1, epsilon)
            && T::abs_diff_eq(&self.e2, &other.e2, epsilon)
            && T::abs_diff_eq(&self.e3, &other.e3, epsilon)
            && T::abs_diff_eq(&self.e0, &other.e0, epsilon)
            && T::abs_diff_eq(&self.e023, &other.e023, epsilon)
            && T::abs_diff_eq(&self.e031, &other.e031, epsilon)
            && T::abs_diff_eq(&self.e012, &other.e012, epsilon)
            && T::abs_diff_eq(&self.e123, &other.e123, epsilon)
    }
}

impl<T: Float> RelativeEq for Flector<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.e1, &other.e1, epsilon, max_relative)
            && T::relative_eq(&self.e2, &other.e2, epsilon, max_relative)
            && T::relative_eq(&self.e3, &other.e3, epsilon, max_relative)
            && T::relative_eq(&self.e0, &other.e0, epsilon, max_relative)
            && T::relative_eq(&self.e023, &other.e023, epsilon, max_relative)
            && T::relative_eq(&self.e031, &other.e031, epsilon, max_relative)
            && T::relative_eq(&self.e012, &other.e012, epsilon, max_relative)
            && T::relative_eq(&self.e123, &other.e123, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Flector<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.e1, &other.e1, epsilon, max_ulps)
            && T::ulps_eq(&self.e2, &other.e2, epsilon, max_ulps)
            && T::ulps_eq(&self.e3, &other.e3, epsilon, max_ulps)
            && T::ulps_eq(&self.e0, &other.e0, epsilon, max_ulps)
            && T::ulps_eq(&self.e023, &other.e023, epsilon, max_ulps)
            && T::ulps_eq(&self.e031, &other.e031, epsilon, max_ulps)
            && T::ulps_eq(&self.e012, &other.e012, epsilon, max_ulps)
            && T::ulps_eq(&self.e123, &other.e123, epsilon, max_ulps)
    }
}
