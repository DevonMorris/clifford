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
    /// Creates a circle from center coordinates and radius.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::Circle;
    /// use approx::relative_eq;
    ///
    /// let circle = Circle::from_center_radius(1.0_f64, 2.0, 3.0);
    /// let (cx, cy) = circle.center().unwrap();
    /// let r = circle.radius().unwrap();
    ///
    /// assert!(relative_eq!(cx, 1.0, epsilon = 1e-10));
    /// assert!(relative_eq!(cy, 2.0, epsilon = 1e-10));
    /// assert!(relative_eq!(r, 3.0, epsilon = 1e-10));
    /// ```
    pub fn from_center_radius(cx: T, cy: T, radius: T) -> Self {
        // Construct by wedging three points on the circle
        let p1 = RoundPoint::from_euclidean(cx + radius, cy);
        let p2 = RoundPoint::from_euclidean(cx, cy + radius);
        let p3 = RoundPoint::from_euclidean(cx - radius, cy);
        Self::from_three_points(&p1, &p2, &p3)
    }

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
    /// A line is a circle with w = 0 (the e₁₂₃ component).
    /// When three points are collinear, their wedge product has w ≈ 0.
    #[inline]
    pub fn is_line(&self, epsilon: T) -> bool {
        self.w().abs() < epsilon
    }

    /// Extracts the Euclidean center coordinates of the circle.
    ///
    /// Returns `None` if this is a line (w ≈ 0), since lines have no finite center.
    ///
    /// # Mathematical Background
    ///
    /// In our Cl(3,1) orthonormal basis, the Circle trivector components encode
    /// the center as:
    /// - center_x = r / w (e₂₃₄ component divided by e₁₂₃ component)
    /// - center_y = -cy / w (negated e₁₃₄ component divided by e₁₂₃ component)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::{Circle, RoundPoint};
    /// use approx::relative_eq;
    ///
    /// let p1 = RoundPoint::from_euclidean(3.0_f64, 2.0);
    /// let p2 = RoundPoint::from_euclidean(2.0, 3.0);
    /// let p3 = RoundPoint::from_euclidean(1.0, 2.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// let (cx, cy) = circle.center().unwrap();
    ///
    /// assert!(relative_eq!(cx, 2.0, epsilon = 1e-10));
    /// assert!(relative_eq!(cy, 2.0, epsilon = 1e-10));
    /// ```
    pub fn center(&self) -> Option<(T, T)> {
        if self.w().abs() < T::epsilon() {
            return None; // This is a line, no finite center
        }
        let w = self.w();
        let center_x = self.r() / w;
        let center_y = -self.cy() / w;
        Some((center_x, center_y))
    }

    /// Extracts the radius of the circle.
    ///
    /// Returns `None` if this is a line (w ≈ 0), since lines have infinite radius.
    ///
    /// # Mathematical Background
    ///
    /// In our Cl(3,1) orthonormal basis, the radius is computed as:
    /// ```text
    /// radius² = 2 * cx / w + center_x² + center_y²
    /// ```
    /// where cx is the e₁₂₄ component, and center coordinates are extracted
    /// from other components.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::{Circle, RoundPoint};
    /// use approx::relative_eq;
    ///
    /// // Create unit circle at origin
    /// let p1 = RoundPoint::from_euclidean(1.0_f64, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 1.0);
    /// let p3 = RoundPoint::from_euclidean(-1.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// let radius = circle.radius().unwrap();
    ///
    /// assert!(relative_eq!(radius, 1.0, epsilon = 1e-10));
    /// ```
    pub fn radius(&self) -> Option<T> {
        if self.w().abs() < T::epsilon() {
            return None; // This is a line, infinite radius
        }
        let w = self.w();
        let center_x = self.r() / w;
        let center_y = -self.cy() / w;
        let radius_sq = T::TWO * self.cx() / w + center_x * center_x + center_y * center_y;
        // radius_sq should be non-negative for a valid circle
        if radius_sq < T::zero() {
            return None;
        }
        Some(radius_sq.sqrt())
    }

    /// Returns the curvature (1/radius) of the circle.
    ///
    /// Returns `None` if this is a line (curvature = 0) or if radius extraction fails.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::{Circle, RoundPoint};
    /// use approx::relative_eq;
    ///
    /// // Circle with radius 2
    /// let p1 = RoundPoint::from_euclidean(2.0_f64, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 2.0);
    /// let p3 = RoundPoint::from_euclidean(-2.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// let curvature = circle.curvature().unwrap();
    ///
    /// assert!(relative_eq!(curvature, 0.5, epsilon = 1e-10));
    /// ```
    pub fn curvature(&self) -> Option<T> {
        self.radius().map(|r| T::one() / r)
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
    pub fn compose(&self, other: Self) -> Self {
        other * *self
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

    #[test]
    fn test_circle_is_line_detection() {
        use crate::ops::Wedge;

        // Three collinear points should produce a line (w ≈ 0)
        let p1 = RoundPoint::from_euclidean(0.0_f64, 0.0);
        let p2 = RoundPoint::from_euclidean(1.0, 0.0);
        let p3 = RoundPoint::from_euclidean(2.0, 0.0);

        let result: Circle<f64> = p1.wedge(&p2).wedge(&p3);
        println!("\nCollinear points (0,0), (1,0), (2,0):");
        println!(
            "  w = {:.6}, is_line = {}",
            result.w(),
            result.is_line(1e-10)
        );
        assert!(
            result.is_line(1e-10),
            "Collinear points should produce a line (w=0)"
        );

        // Non-collinear points should produce a circle (w ≠ 0)
        let p1 = RoundPoint::from_euclidean(1.0_f64, 0.0);
        let p2 = RoundPoint::from_euclidean(0.0, 1.0);
        let p3 = RoundPoint::from_euclidean(-1.0, 0.0);

        let result: Circle<f64> = p1.wedge(&p2).wedge(&p3);
        println!("\nNon-collinear points (1,0), (0,1), (-1,0):");
        println!(
            "  w = {:.6}, is_line = {}",
            result.w(),
            result.is_line(1e-10)
        );
        assert!(
            !result.is_line(1e-10),
            "Non-collinear points should produce a circle (w≠0)"
        );
    }

    #[test]
    fn test_circle_center_radius_extraction() {
        // Test the Circle::center() and Circle::radius() methods
        let test_cases: Vec<(f64, f64, f64)> = vec![
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 2.0),
            (0.0, 0.0, 0.5),
            (1.0, 0.0, 1.0),
            (0.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (2.0, 3.0, 1.0),
            (1.0, 1.0, 2.0),
            (-1.0, -1.0, 1.5),
            (3.5, -2.7, 1.2),
        ];

        for (expected_cx, expected_cy, expected_r) in &test_cases {
            // Three points on circle
            let p1 = RoundPoint::from_euclidean(expected_cx + expected_r, *expected_cy);
            let p2 = RoundPoint::from_euclidean(*expected_cx, expected_cy + expected_r);
            let p3 = RoundPoint::from_euclidean(expected_cx - expected_r, *expected_cy);

            let circle = Circle::from_three_points(&p1, &p2, &p3);

            // Test center extraction
            let (calc_cx, calc_cy) = circle.center().expect("Circle should have finite center");
            assert!(
                relative_eq!(
                    calc_cx,
                    *expected_cx,
                    epsilon = RELATIVE_EQ_EPS,
                    max_relative = RELATIVE_EQ_EPS
                ),
                "Center x mismatch: expected {}, got {}",
                expected_cx,
                calc_cx
            );
            assert!(
                relative_eq!(
                    calc_cy,
                    *expected_cy,
                    epsilon = RELATIVE_EQ_EPS,
                    max_relative = RELATIVE_EQ_EPS
                ),
                "Center y mismatch: expected {}, got {}",
                expected_cy,
                calc_cy
            );

            // Test radius extraction
            let calc_r = circle.radius().expect("Circle should have finite radius");
            assert!(
                relative_eq!(
                    calc_r,
                    *expected_r,
                    epsilon = RELATIVE_EQ_EPS,
                    max_relative = RELATIVE_EQ_EPS
                ),
                "Radius mismatch: expected {}, got {}",
                expected_r,
                calc_r
            );

            // Test curvature
            let calc_curv = circle.curvature().expect("Circle should have curvature");
            let expected_curv = 1.0 / expected_r;
            assert!(
                relative_eq!(
                    calc_curv,
                    expected_curv,
                    epsilon = RELATIVE_EQ_EPS,
                    max_relative = RELATIVE_EQ_EPS
                ),
                "Curvature mismatch: expected {}, got {}",
                expected_curv,
                calc_curv
            );
        }
    }

    #[test]
    fn test_circle_line_has_no_center_or_radius() {
        // Collinear points form a line, which should return None for center/radius
        let p1 = RoundPoint::from_euclidean(0.0_f64, 0.0);
        let p2 = RoundPoint::from_euclidean(1.0, 0.0);
        let p3 = RoundPoint::from_euclidean(2.0, 0.0);

        let line = Circle::from_three_points(&p1, &p2, &p3);

        assert!(line.is_line(1e-10), "Collinear points should form a line");
        assert!(line.center().is_none(), "Line should have no finite center");
        assert!(line.radius().is_none(), "Line should have no finite radius");
        assert!(line.curvature().is_none(), "Line should have no curvature");
    }

    #[test]
    fn test_inverse_sandwich_computes_reflection_not_inversion() {
        use crate::ops::InverseSandwich;

        // The InverseSandwich operation C × P × C⁻¹ computes a REFLECTION
        // through the circle, NOT classic circle inversion.
        //
        // Classic circle inversion: P' = O + r²/|OP|² × (P - O)
        // For point (4,0) through circle at origin with r=2:
        // Expected classic inversion result: (1, 0) since 4 * 1 = 4 = 2²
        //
        // But InverseSandwich gives a different result because it's
        // a reflection operation, not inversive geometry.

        let inv_circle = Circle::from_center_radius(0.0_f64, 0.0, 2.0);
        let point = RoundPoint::from_euclidean(4.0, 0.0);

        let inverted = inv_circle.try_inverse_sandwich(&point);

        // The operation should succeed (not return None)
        assert!(
            inverted.is_some(),
            "InverseSandwich should return Some for valid circle and point"
        );

        // The result is NOT classic circle inversion
        let (ix, _iy) = inverted.unwrap().to_euclidean().unwrap();
        // Classic inversion would give 1.0, but reflection gives ~0.44
        assert!(
            (ix - 1.0).abs() > 0.1,
            "InverseSandwich is NOT classic circle inversion"
        );
    }
}
