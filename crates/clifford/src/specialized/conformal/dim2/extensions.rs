//! Domain-specific extensions for 2D Conformal GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to 2D conformal geometry.

use super::generated::types::{Circle, FlatPoint, Line, Motor, PointPair, RoundPoint};
use crate::scalar::Float;

// ============================================================================
// RoundPoint extensions
// ============================================================================

impl<T: Float> RoundPoint<T> {
    /// Creates a round point from Euclidean 2D coordinates.
    ///
    /// The point is embedded as a null vector:
    /// `P = x·e₁ + y·e₂ + e₃ + ½(x² + y²)·e₄`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::RoundPoint;
    ///
    /// let p = RoundPoint::from_euclidean(3.0, 4.0);
    /// assert_eq!(p.to_euclidean(), Some((3.0, 4.0)));
    /// ```
    #[inline]
    pub fn from_euclidean(x: T, y: T) -> Self {
        let half_sq = (x * x + y * y) / T::TWO;
        Self::new_unchecked(x, y, T::one(), half_sq)
    }

    /// Creates a round point from Euclidean coordinates with a weight.
    ///
    /// Useful for creating weighted point combinations.
    #[inline]
    pub fn from_euclidean_weighted(x: T, y: T, weight: T) -> Self {
        let half_sq = (x * x + y * y) / T::TWO;
        Self::new_unchecked(x * weight, y * weight, weight, half_sq * weight)
    }

    /// Origin point (0, 0).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::RoundPoint;
    ///
    /// let o = RoundPoint::<f64>::origin();
    /// assert_eq!(o.to_euclidean(), Some((0.0, 0.0)));
    /// ```
    #[inline]
    pub fn origin() -> Self {
        Self::new_unchecked(T::zero(), T::zero(), T::one(), T::zero())
    }

    /// Point at infinity.
    ///
    /// This is the e₄ basis vector, representing the point at infinity.
    #[inline]
    pub fn infinity() -> Self {
        Self::new_unchecked(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Extracts Euclidean coordinates from a finite round point.
    ///
    /// Returns `None` for points at infinity (o = 0).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::RoundPoint;
    ///
    /// let p = RoundPoint::from_euclidean(3.0_f64, 4.0);
    /// assert_eq!(p.to_euclidean(), Some((3.0, 4.0)));
    ///
    /// let inf = RoundPoint::<f64>::infinity();
    /// assert_eq!(inf.to_euclidean(), None);
    /// ```
    #[inline]
    pub fn to_euclidean(&self) -> Option<(T, T)> {
        if self.o().abs() < T::epsilon() {
            None
        } else {
            Some((self.x() / self.o(), self.y() / self.o()))
        }
    }

    /// Euclidean x-coordinate.
    ///
    /// # Panics
    ///
    /// May produce incorrect results for points at infinity (o = 0).
    #[inline]
    pub fn euclidean_x(&self) -> T {
        self.x() / self.o()
    }

    /// Euclidean y-coordinate.
    ///
    /// # Panics
    ///
    /// May produce incorrect results for points at infinity (o = 0).
    #[inline]
    pub fn euclidean_y(&self) -> T {
        self.y() / self.o()
    }

    /// Returns true if this is a point at infinity (o ≈ 0).
    #[inline]
    pub fn is_at_infinity(&self, epsilon: T) -> bool {
        self.o().abs() < epsilon
    }

    /// Squared Euclidean distance between two points.
    ///
    /// In CGA, the inner product of two null vectors gives:
    /// `P₁ · P₂ = -½|p₁ - p₂|²`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::RoundPoint;
    /// use approx::abs_diff_eq;
    ///
    /// let p1 = RoundPoint::<f64>::origin();
    /// let p2 = RoundPoint::from_euclidean(3.0, 4.0);
    /// assert!(abs_diff_eq!(p1.distance_squared(&p2), 25.0, epsilon = 1e-10));
    /// ```
    pub fn distance_squared(&self, other: &RoundPoint<T>) -> T {
        // For normalized points (o = 1):
        // P₁ · P₂ = x₁x₂ + y₁y₂ - o₁i₂ - i₁o₂
        // = x₁x₂ + y₁y₂ - ½(x₁² + y₁²) - ½(x₂² + y₂²)
        // = -½((x₁ - x₂)² + (y₁ - y₂)²)
        let inner = self.x() * other.x() + self.y() * other.y()
            - self.o() * other.i()
            - self.i() * other.o();
        -T::TWO * inner
    }

    /// Euclidean distance between two points.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::RoundPoint;
    /// use approx::abs_diff_eq;
    ///
    /// let p1 = RoundPoint::<f64>::origin();
    /// let p2 = RoundPoint::from_euclidean(3.0, 4.0);
    /// assert!(abs_diff_eq!(p1.distance(&p2), 5.0, epsilon = 1e-10));
    /// ```
    pub fn distance(&self, other: &RoundPoint<T>) -> T {
        self.distance_squared(other).abs().sqrt()
    }
}

// ============================================================================
// Circle extensions
// ============================================================================

impl<T: Float> Circle<T> {
    /// Creates a circle passing through three points.
    ///
    /// The circle is the outer product: `C = P₁ ∧ P₂ ∧ P₃`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::{Circle, RoundPoint};
    ///
    /// let p1 = RoundPoint::from_euclidean(1.0, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 1.0);
    /// let p3 = RoundPoint::from_euclidean(-1.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// // This circle is the unit circle centered at origin
    /// ```
    pub fn from_three_points(p1: &RoundPoint<T>, p2: &RoundPoint<T>, p3: &RoundPoint<T>) -> Self {
        use crate::ops::Wedge;
        p1.wedge(p2).wedge(p3)
    }

    /// Returns true if this represents a line (circle through infinity).
    ///
    /// A line is a circle with w = 0 (the e123 component).
    #[inline]
    pub fn is_line(&self, epsilon: T) -> bool {
        self.w().abs() < epsilon
    }
}

// ============================================================================
// Line extensions
// ============================================================================

impl<T: Float> Line<T> {
    /// Creates a line through two points.
    ///
    /// A line is a circle through infinity: `L = P₁ ∧ P₂ ∧ e∞`
    pub fn from_two_points(p1: &RoundPoint<T>, p2: &RoundPoint<T>) -> Self {
        use crate::ops::Wedge;
        let inf = RoundPoint::infinity();
        let circle: Circle<T> = p1.wedge(p2).wedge(&inf);
        // Extract line components from the circle (w should be 0)
        Self::new_unchecked(circle.cx(), circle.cy(), circle.r())
    }

    /// Creates a line from the equation `nx*x + ny*y = d`.
    ///
    /// The normal vector (nx, ny) points away from the origin if d > 0.
    #[inline]
    pub fn from_equation(nx: T, ny: T, d: T) -> Self {
        Self::new_unchecked(nx, ny, d)
    }
}

// ============================================================================
// FlatPoint extensions
// ============================================================================

impl<T: Float> FlatPoint<T> {
    /// Creates a flat point from Euclidean coordinates.
    ///
    /// A flat point represents a point using only the infinity-direction
    /// components of a point pair.
    #[inline]
    pub fn from_euclidean(x: T, y: T) -> Self {
        Self::new_unchecked(x, y, T::one())
    }

    /// Extracts Euclidean coordinates.
    #[inline]
    pub fn to_euclidean(&self) -> Option<(T, T)> {
        if self.oi().abs() < T::epsilon() {
            None
        } else {
            Some((self.ix() / self.oi(), self.iy() / self.oi()))
        }
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Identity motor (leaves all elements unchanged).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::Motor;
    ///
    /// let m = Motor::<f64>::identity();
    /// assert_eq!(m.s(), 1.0);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self::new_unchecked(
            T::one(),  // s
            T::zero(), // m
            T::zero(), // ox
            T::zero(), // ix
            T::zero(), // oy
            T::zero(), // iy
            T::zero(), // oi
            T::zero(), // ps
        )
    }

    /// Creates a translation motor.
    ///
    /// Translates by vector (dx, dy).
    ///
    /// In CGA, a translation is represented as:
    /// `T = 1 - ½(dx·e₁ + dy·e₂)·e∞ = 1 - ½dx·e₁₄ - ½dy·e₂₄`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::{Motor, RoundPoint};
    /// use clifford::ops::Transform;
    ///
    /// let m = Motor::from_translation(3.0, 4.0);
    /// let p = RoundPoint::<f64>::origin();
    /// let translated = m.transform(&p);
    /// // translated should be at (3, 4)
    /// ```
    pub fn from_translation(dx: T, dy: T) -> Self {
        let half = T::one() / T::TWO;
        Self::new_unchecked(
            T::one(),   // s = 1
            T::zero(),  // m = 0
            T::zero(),  // ox = 0
            -dx * half, // ix = -dx/2 (e14 component)
            T::zero(),  // oy = 0
            -dy * half, // iy = -dy/2 (e24 component)
            T::zero(),  // oi = 0
            T::zero(),  // ps = 0
        )
    }

    /// Creates a rotation motor around the origin.
    ///
    /// Rotates by angle (in radians) counterclockwise.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::{Motor, RoundPoint};
    /// use clifford::ops::Transform;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let m = Motor::from_rotation(FRAC_PI_2);
    /// let p = RoundPoint::from_euclidean(1.0, 0.0);
    /// let rotated = m.transform(&p);
    /// // rotated should be at approximately (0, 1)
    /// ```
    pub fn from_rotation(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            half.cos(),  // s
            -half.sin(), // m (rotation in e12 plane, negated for sandwich)
            T::zero(),   // ox
            T::zero(),   // ix
            T::zero(),   // oy
            T::zero(),   // iy
            T::zero(),   // oi
            T::zero(),   // ps
        )
    }

    /// Creates a uniform scaling (dilation) motor.
    ///
    /// Scales by factor `k` around the origin.
    ///
    /// In CGA, dilation is represented using hyperbolic functions
    /// in the origin-infinity plane (e₃₄).
    pub fn from_dilation(scale: T) -> Self {
        let half_log = scale.ln() / T::TWO;
        Self::new_unchecked(
            half_log.cosh(), // s
            T::zero(),       // m
            T::zero(),       // ox
            T::zero(),       // ix
            T::zero(),       // oy
            T::zero(),       // iy
            half_log.sinh(), // oi (dilation is in e34 plane)
            T::zero(),       // ps
        )
    }

    /// Composes two motors: applies `self` first, then `other`.
    ///
    /// The result is the geometric product: `other * self`
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        *other * *self
    }

    /// Returns the inverse motor.
    ///
    /// For a versor V, the inverse is V† / |V|².
    pub fn inverse(&self) -> Self {
        let rev = self.reverse();
        let norm_sq = self.norm_squared();
        if norm_sq.abs() < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            rev.s() / norm_sq,
            rev.m() / norm_sq,
            rev.ox() / norm_sq,
            rev.ix() / norm_sq,
            rev.oy() / norm_sq,
            rev.iy() / norm_sq,
            rev.oi() / norm_sq,
            rev.ps() / norm_sq,
        )
    }
}

// ============================================================================
// PointPair extensions
// ============================================================================

impl<T: Float> PointPair<T> {
    /// Creates a point pair from two points.
    ///
    /// The point pair is the outer product: `PP = P₁ ∧ P₂`
    pub fn from_points(p1: &RoundPoint<T>, p2: &RoundPoint<T>) -> Self {
        use crate::ops::Wedge;
        p1.wedge(p2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;

    #[test]
    fn round_point_euclidean_roundtrip() {
        let x = 3.0_f64;
        let y = 4.0_f64;
        let p = RoundPoint::from_euclidean(x, y);
        let (rx, ry) = p.to_euclidean().unwrap();

        assert!(relative_eq!(
            rx,
            x,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            ry,
            y,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn round_point_origin() {
        let o = RoundPoint::<f64>::origin();
        let (x, y) = o.to_euclidean().unwrap();
        assert!(relative_eq!(
            x,
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            y,
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn round_point_infinity() {
        let inf = RoundPoint::<f64>::infinity();
        assert!(inf.to_euclidean().is_none());
        assert!(inf.is_at_infinity(RELATIVE_EQ_EPS));
    }

    #[test]
    fn round_point_distance() {
        let p1 = RoundPoint::<f64>::origin();
        let p2 = RoundPoint::from_euclidean(3.0, 4.0);

        let dist = p1.distance(&p2);

        assert!(relative_eq!(
            dist,
            5.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_identity() {
        let m = Motor::<f64>::identity();
        assert!(relative_eq!(
            m.s(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            m.m(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
