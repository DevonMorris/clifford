//! Type definitions for 2D Projective Geometric Algebra.
//!
//! # Convention Choices
//!
//! - **Basis ordering**: e₁, e₂, e₀ (Euclidean first, null last)
//! - **Motor composition**: `a.compose(&b)` applies `a` first, then `b`
//! - **nalgebra correspondence**: `iso1 * iso2` applies `iso2` first, then `iso1`,
//!   so `motor2.compose(&motor1)` ≡ `iso1 * iso2`

use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::scalar::Float;

/// A point in 2D PGA (grade 1 vector).
///
/// In point-based PGA, a point is represented in homogeneous coordinates:
/// `P = x·e₁ + y·e₂ + w·e₀`
///
/// where `(x/w, y/w)` are the Cartesian coordinates when `w ≠ 0`.
/// Points with `w = 0` represent ideal points (points at infinity).
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim2::Point;
///
/// // Point at (3, 4)
/// let p = Point::new(3.0, 4.0);
/// assert_eq!(p.x(), 3.0);
/// assert_eq!(p.y(), 4.0);
///
/// // Point at infinity in direction (1, 0)
/// let ideal = Point::ideal(1.0, 0.0);
/// assert!(ideal.is_ideal(1e-10));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Point<T: Float> {
    /// Coefficient of e₁ (x-coordinate × w).
    pub e1: T,
    /// Coefficient of e₂ (y-coordinate × w).
    pub e2: T,
    /// Coefficient of e₀ (homogeneous weight).
    pub e0: T,
}

impl<T: Float> Point<T> {
    /// Creates a finite point at Cartesian coordinates (x, y).
    ///
    /// The homogeneous weight `w` is set to 1.
    #[inline]
    pub fn new(x: T, y: T) -> Self {
        Self {
            e1: x,
            e2: y,
            e0: T::one(),
        }
    }

    /// Creates a point from homogeneous coordinates.
    #[inline]
    pub fn from_homogeneous(e1: T, e2: T, e0: T) -> Self {
        Self { e1, e2, e0 }
    }

    /// Creates an ideal point (point at infinity) in the given direction.
    ///
    /// Ideal points have `w = 0` and represent directions rather than positions.
    #[inline]
    pub fn ideal(dx: T, dy: T) -> Self {
        Self {
            e1: dx,
            e2: dy,
            e0: T::zero(),
        }
    }

    /// Origin point (0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero())
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
                e0: T::one(),
            })
        }
    }

    /// Returns the Cartesian coordinates as a tuple, if finite.
    ///
    /// Returns `None` if this is an ideal point.
    #[inline]
    pub fn to_cartesian(&self) -> Option<(T, T)> {
        if self.e0.abs() < T::epsilon() {
            None
        } else {
            Some((self.e1 / self.e0, self.e2 / self.e0))
        }
    }
}

impl<T: Float> Default for Point<T> {
    fn default() -> Self {
        Self::origin()
    }
}

/// A line in 2D PGA (grade 2 bivector).
///
/// In point-based PGA, a line is represented as a bivector:
/// `L = d·e₁₂ + a·e₂₀ + b·e₀₁`
///
/// This represents the implicit line `ax + by + d = 0` where `(a, b)` is
/// the normal direction.
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim2::{Point, Line};
///
/// // Line through two points
/// let p1 = Point::new(0.0, 0.0);
/// let p2 = Point::new(1.0, 1.0);
/// let line = p1.join(&p2);
///
/// // Line from implicit equation: x - y = 0
/// let line2 = Line::from_implicit(1.0, -1.0, 0.0);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Line<T: Float> {
    /// Coefficient of e₁₂ (related to signed distance from origin).
    pub e12: T,
    /// Coefficient of e₂₀ (x-component of normal).
    pub e20: T,
    /// Coefficient of e₀₁ (y-component of normal).
    pub e01: T,
}

impl<T: Float> Line<T> {
    /// Creates a line from bivector coefficients.
    #[inline]
    pub fn new(e12: T, e20: T, e01: T) -> Self {
        Self { e12, e20, e01 }
    }

    /// Creates the zero line (degenerate).
    #[inline]
    pub fn zero() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Creates a line from implicit equation `ax + by + d = 0`.
    ///
    /// The normal direction is `(a, b)`.
    #[inline]
    pub fn from_implicit(a: T, b: T, d: T) -> Self {
        Self {
            e12: d,
            e20: a,
            e01: b,
        }
    }

    /// Creates the x-axis (y = 0).
    #[inline]
    pub fn x_axis() -> Self {
        Self::from_implicit(T::zero(), T::one(), T::zero())
    }

    /// Creates the y-axis (x = 0).
    #[inline]
    pub fn y_axis() -> Self {
        Self::from_implicit(T::one(), T::zero(), T::zero())
    }

    /// Returns the normal direction `(a, b)` of the line `ax + by + d = 0`.
    #[inline]
    pub fn normal(&self) -> (T, T) {
        (self.e20, self.e01)
    }

    /// Returns the signed distance from the origin (scaled by normal length).
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        self.e12
    }

    /// Returns the squared norm of the line (sum of squared coefficients).
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.e12 * self.e12 + self.e20 * self.e20 + self.e01 * self.e01
    }

    /// Returns the norm of the line.
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized line (unit norm).
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        if n.abs() < T::epsilon() {
            *self
        } else {
            Self {
                e12: self.e12 / n,
                e20: self.e20 / n,
                e01: self.e01 / n,
            }
        }
    }

    /// Returns true if this line is degenerate (zero).
    #[inline]
    pub fn is_zero(&self, epsilon: T) -> bool {
        self.e12.abs() < epsilon && self.e20.abs() < epsilon && self.e01.abs() < epsilon
    }
}

impl<T: Float> Default for Line<T> {
    fn default() -> Self {
        Self::zero()
    }
}

/// A motor (rigid transformation) in 2D PGA.
///
/// Motors represent rigid body transformations (rotation + translation).
/// In 2D PGA, a motor is an even-grade element:
/// `M = s + d·e₁₂ + tx·e₂₀ + ty·e₀₁`
///
/// A motor transforms geometric objects via the sandwich product:
/// `X' = M X M̃`
///
/// # Construction
///
/// - [`Motor::identity`]: No transformation
/// - [`Motor::from_rotation`]: Pure rotation around origin
/// - [`Motor::from_translation`]: Pure translation
/// - [`Motor::compose`]: Combine transformations
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim2::{Point, Motor};
/// use std::f64::consts::FRAC_PI_2;
/// use approx::abs_diff_eq;
///
/// // 90° rotation around origin
/// let rotation = Motor::from_rotation(FRAC_PI_2);
///
/// // Translation by (1, 2)
/// let translation = Motor::from_translation(1.0, 2.0);
///
/// // Compose: first rotate, then translate (leftmost acts first)
/// let combined = rotation.compose(&translation);
///
/// // Transform a point
/// let p = Point::new(1.0, 0.0);
/// let result = combined.transform_point(&p);
///
/// // (1,0) rotated 90° becomes (0,1), then translated to (1,3)
/// assert!(abs_diff_eq!(result.x(), 1.0, epsilon = 1e-10));
/// assert!(abs_diff_eq!(result.y(), 3.0, epsilon = 1e-10));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Motor<T: Float> {
    /// Scalar part (cos(θ/2) for pure rotation).
    pub s: T,
    /// Coefficient of e₁₂ (sin(θ/2) for pure rotation).
    pub e12: T,
    /// Coefficient of e₂₀ (translation x / 2 for pure translation).
    pub e20: T,
    /// Coefficient of e₀₁ (translation y / 2 for pure translation).
    pub e01: T,
}

impl<T: Float> Motor<T> {
    /// Creates a motor from components.
    #[inline]
    pub fn new(s: T, e12: T, e20: T, e01: T) -> Self {
        Self { s, e12, e20, e01 }
    }

    /// Creates the identity motor (no transformation).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::one(),
            e12: T::zero(),
            e20: T::zero(),
            e01: T::zero(),
        }
    }

    /// Creates a pure rotation motor around the origin.
    ///
    /// # Arguments
    ///
    /// * `angle` - Rotation angle in radians (counterclockwise positive)
    #[inline]
    pub fn from_rotation(angle: T) -> Self {
        let half = angle / T::TWO;
        Self {
            s: half.cos(),
            e12: half.sin(),
            e20: T::zero(),
            e01: T::zero(),
        }
    }

    /// Creates a pure translation motor.
    ///
    /// # Arguments
    ///
    /// * `dx` - Translation in x direction
    /// * `dy` - Translation in y direction
    #[inline]
    pub fn from_translation(dx: T, dy: T) -> Self {
        // Translation: T = 1 + (d/2)·e₀i where e₀i are the ideal line bivectors
        // For translation by (dx, dy): T = 1 + (dx/2)·e₂₀ + (dy/2)·e₀₁
        Self {
            s: T::one(),
            e12: T::zero(),
            e20: dx / T::TWO,
            e01: dy / T::TWO,
        }
    }

    /// Creates a motor for rotation around an arbitrary point.
    ///
    /// # Arguments
    ///
    /// * `angle` - Rotation angle in radians
    /// * `center` - Center of rotation
    #[inline]
    pub fn from_rotation_around(angle: T, center: Point<T>) -> Self {
        // Translate center to origin, rotate, translate back
        let to_origin = Self::from_translation(-center.x(), -center.y());
        let rotation = Self::from_rotation(angle);
        let from_origin = Self::from_translation(center.x(), center.y());
        from_origin.compose(&rotation).compose(&to_origin)
    }

    /// Returns the reverse of the motor: `M̃`.
    ///
    /// For motors, the reverse negates the bivector parts.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self {
            s: self.s,
            e12: -self.e12,
            e20: -self.e20,
            e01: -self.e01,
        }
    }

    /// Returns the squared norm of the motor.
    ///
    /// For a unit motor, this equals 1.
    #[inline]
    pub fn norm_squared(&self) -> T {
        // In PGA, the norm is computed differently due to the null vector
        // For the rotation part: s² + e12²
        // The translation parts don't contribute to the squared norm in PGA
        self.s * self.s + self.e12 * self.e12
    }

    /// Returns the norm of the motor.
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Normalizes the motor to unit norm (rotation part only).
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        if n.abs() < T::epsilon() {
            *self
        } else {
            Self {
                s: self.s / n,
                e12: self.e12 / n,
                e20: self.e20 / n,
                e01: self.e01 / n,
            }
        }
    }

    /// Returns the inverse motor: `M⁻¹`.
    ///
    /// For a unit motor, this equals the reverse.
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        if norm_sq.abs() < T::epsilon() {
            return *self;
        }
        let rev = self.reverse();
        Self {
            s: rev.s / norm_sq,
            e12: rev.e12 / norm_sq,
            e20: rev.e20 / norm_sq,
            e01: rev.e01 / norm_sq,
        }
    }

    /// Composes two motors via geometric product: `self * other`.
    ///
    /// In PGA, the sandwich product `(self * other) P (self * other)⁻¹`
    /// applies `self` first, then `other`. This follows the PGA convention
    /// where the leftmost motor in the product acts first.
    ///
    /// To get "first A, then B" transformation, use `a.compose(&b)`.
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        // Motor multiplication in 2D PGA
        // (s₁ + b₁)(s₂ + b₂) where b represents the bivector parts
        Self {
            s: self.s * other.s - self.e12 * other.e12,
            e12: self.s * other.e12 + self.e12 * other.s,
            e20: self.s * other.e20 + self.e20 * other.s + self.e12 * other.e01
                - self.e01 * other.e12,
            e01: self.s * other.e01 + self.e01 * other.s - self.e12 * other.e20
                + self.e20 * other.e12,
        }
    }

    /// Linear interpolation between motors (normalized).
    #[inline]
    pub fn lerp(&self, other: &Self, t: T) -> Self {
        let one_minus_t = T::one() - t;
        Self {
            s: self.s * one_minus_t + other.s * t,
            e12: self.e12 * one_minus_t + other.e12 * t,
            e20: self.e20 * one_minus_t + other.e20 * t,
            e01: self.e01 * one_minus_t + other.e01 * t,
        }
        .normalized()
    }

    /// Returns the rotation angle in radians.
    #[inline]
    pub fn rotation_angle(&self) -> T {
        self.e12.atan2(self.s) * T::TWO
    }

    /// Returns the translation vector (dx, dy).
    #[inline]
    pub fn translation(&self) -> (T, T) {
        // For a pure translation: dx = 2*e20, dy = 2*e01
        // For a composed motor, this is more complex
        (self.e20 * T::TWO, self.e01 * T::TWO)
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
            && T::ulps_eq(&self.e0, &other.e0, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Line<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.e12, &other.e12, epsilon)
            && T::abs_diff_eq(&self.e20, &other.e20, epsilon)
            && T::abs_diff_eq(&self.e01, &other.e01, epsilon)
    }
}

impl<T: Float> RelativeEq for Line<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.e12, &other.e12, epsilon, max_relative)
            && T::relative_eq(&self.e20, &other.e20, epsilon, max_relative)
            && T::relative_eq(&self.e01, &other.e01, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Line<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.e12, &other.e12, epsilon, max_ulps)
            && T::ulps_eq(&self.e20, &other.e20, epsilon, max_ulps)
            && T::ulps_eq(&self.e01, &other.e01, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Motor<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.s, &other.s, epsilon)
            && T::abs_diff_eq(&self.e12, &other.e12, epsilon)
            && T::abs_diff_eq(&self.e20, &other.e20, epsilon)
            && T::abs_diff_eq(&self.e01, &other.e01, epsilon)
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
            && T::relative_eq(&self.e12, &other.e12, epsilon, max_relative)
            && T::relative_eq(&self.e20, &other.e20, epsilon, max_relative)
            && T::relative_eq(&self.e01, &other.e01, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Motor<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.s, &other.s, epsilon, max_ulps)
            && T::ulps_eq(&self.e12, &other.e12, epsilon, max_ulps)
            && T::ulps_eq(&self.e20, &other.e20, epsilon, max_ulps)
            && T::ulps_eq(&self.e01, &other.e01, epsilon, max_ulps)
    }
}
