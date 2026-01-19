//! Domain-specific extensions for 2D Conformal GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to 2D conformal geometry.
//!
//! # Basis Convention
//!
//! This module uses the **orthonormal conformal basis**:
//! - `e1`, `e2`: Euclidean basis vectors
//! - `ep` (e₊): Positive conformal basis (ep² = +1)
//! - `em` (e₋): Negative conformal basis (em² = -1)
//!
//! The null basis vectors relate to orthonormal as:
//! - `e₀` (origin) = (em - ep) / 2
//! - `e∞` (infinity) = em + ep
//!
//! Point embedding: For Euclidean (x, y), the orthonormal form is:
//! ```text
//! P = x·e₁ + y·e₂ + ½(x²+y²-1)·ep + ½(x²+y²+1)·em
//! ```

use super::generated::types::{Circle, FlatPoint, Line, Motor, PointPair, RoundPoint};
use crate::scalar::Float;

// ============================================================================
// RoundPoint extensions
// ============================================================================

impl<T: Float> RoundPoint<T> {
    /// Creates a round point from Euclidean 2D coordinates.
    ///
    /// The point is embedded using the orthonormal basis:
    /// ```text
    /// P = x·e₁ + y·e₂ + ½(x²+y²-1)·ep + ½(x²+y²+1)·em
    /// ```
    ///
    /// This produces a null vector (P·P = 0) representing the Euclidean point.
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
        let sq = x * x + y * y;
        let half = T::one() / T::TWO;
        // Orthonormal embedding:
        // ep = ½(x² + y² - 1)
        // em = ½(x² + y² + 1)
        let ep_coeff = (sq - T::one()) * half;
        let em_coeff = (sq + T::one()) * half;
        Self::new_unchecked(x, y, ep_coeff, em_coeff)
    }

    /// Creates a round point from Euclidean coordinates with a weight.
    ///
    /// Useful for creating weighted point combinations.
    #[inline]
    pub fn from_euclidean_weighted(x: T, y: T, weight: T) -> Self {
        let sq = x * x + y * y;
        let half = T::one() / T::TWO;
        let ep_coeff = (sq - T::one()) * half;
        let em_coeff = (sq + T::one()) * half;
        Self::new_unchecked(x * weight, y * weight, ep_coeff * weight, em_coeff * weight)
    }

    /// Origin point (0, 0).
    ///
    /// In orthonormal basis: P = -½·ep + ½·em = ½(em - ep) = e₀
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
        let half = T::one() / T::TWO;
        // Origin: x=0, y=0, so sq=0
        // ep = ½(0 - 1) = -½
        // em = ½(0 + 1) = ½
        Self::new_unchecked(T::zero(), T::zero(), -half, half)
    }

    /// Point at infinity.
    ///
    /// In orthonormal basis: e∞ = ep + em
    ///
    /// This represents the point at infinity in the conformal model.
    #[inline]
    pub fn infinity() -> Self {
        // e∞ = ep + em
        Self::new_unchecked(T::zero(), T::zero(), T::one(), T::one())
    }

    /// Extracts Euclidean coordinates from a finite round point.
    ///
    /// Returns `None` for points at infinity (where origin weight ≈ 0).
    ///
    /// The null basis origin weight is recovered as: o = em - ep
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
        // Recover null basis origin weight: o = em - ep
        let o = self.em() - self.ep();
        if o.abs() < T::epsilon() {
            None
        } else {
            Some((self.x() / o, self.y() / o))
        }
    }

    /// Euclidean x-coordinate.
    ///
    /// # Panics
    ///
    /// May produce incorrect results for points at infinity (o = 0).
    #[inline]
    pub fn euclidean_x(&self) -> T {
        let o = self.em() - self.ep();
        self.x() / o
    }

    /// Euclidean y-coordinate.
    ///
    /// # Panics
    ///
    /// May produce incorrect results for points at infinity (o = 0).
    #[inline]
    pub fn euclidean_y(&self) -> T {
        let o = self.em() - self.ep();
        self.y() / o
    }

    /// Returns the null basis origin weight (o = em - ep).
    ///
    /// For a normalized point embedding, this equals 1.
    /// For a point at infinity, this equals 0.
    #[inline]
    pub fn null_origin_weight(&self) -> T {
        self.em() - self.ep()
    }

    /// Returns the null basis infinity weight (i = (em + ep) / 2).
    ///
    /// For a normalized point at (x, y), this equals ½(x² + y²).
    #[inline]
    pub fn null_infinity_weight(&self) -> T {
        (self.em() + self.ep()) / T::TWO
    }

    /// Returns true if this is a point at infinity (o ≈ 0).
    #[inline]
    pub fn is_at_infinity(&self, epsilon: T) -> bool {
        self.null_origin_weight().abs() < epsilon
    }

    /// Squared Euclidean distance between two points.
    ///
    /// In CGA, the inner product of two null vectors gives:
    /// `P₁ · P₂ = -½|p₁ - p₂|²`
    ///
    /// Using the orthonormal metric (ep² = +1, em² = -1):
    /// `P · Q = x₁x₂ + y₁y₂ + ep₁·ep₂ - em₁·em₂`
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
        // Using the orthonormal metric:
        // P · Q = x₁x₂ + y₁y₂ + ep₁·ep₂·(+1) + em₁·em₂·(-1)
        //       = x₁x₂ + y₁y₂ + ep₁·ep₂ - em₁·em₂
        //
        // For null vectors representing Euclidean points:
        // P · Q = -½|p - q|²
        // So |p - q|² = -2(P · Q)
        let inner = self.x() * other.x() + self.y() * other.y() + self.ep() * other.ep()
            - self.em() * other.em();
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
    /// A line is detected when w = e12em, which makes the center/radius
    /// extraction formulas have a zero denominator.
    #[inline]
    pub fn is_line(&self, epsilon: T) -> bool {
        (self.w() - self.e12em()).abs() < epsilon
    }

    /// Extracts the Euclidean center coordinates of the circle.
    ///
    /// Returns `None` if this is a line (w = e12em), since lines have no finite center.
    ///
    /// # Mathematical Background
    ///
    /// For the orthonormal conformal basis, the center is extracted from the
    /// trivector components using empirically derived formulas:
    ///
    /// ```text
    /// cx = e2epem / (w - e12em)
    /// cy = -e1epem / (w - e12em)
    /// ```
    ///
    /// When w = e12em (denominator is zero), the trivector represents a line,
    /// which has no finite center.
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
        // Derived formula for orthonormal CGA basis:
        //   cx = e2epem / (w - e12em)
        //   cy = -e1epem / (w - e12em)
        //
        // When w = e12em, this is a line (no finite center).

        let denom = self.w() - self.e12em();

        if denom.abs() < T::epsilon() {
            // This is a line (w = e12em), no finite center
            return None;
        }

        let center_x = self.e2epem() / denom;
        let center_y = -self.e1epem() / denom;
        Some((center_x, center_y))
    }

    /// Extracts the radius of the circle.
    ///
    /// Returns `None` if this is a line (w = e12em), since lines have infinite radius.
    ///
    /// # Mathematical Background
    ///
    /// For the orthonormal conformal basis, the radius is extracted using:
    ///
    /// ```text
    /// r² = (e12em + w) / (e12em - w) + cx² + cy²
    /// ```
    ///
    /// This formula is derived from the structure of the wedge product P₁ ∧ P₂ ∧ P₃
    /// in the orthonormal basis.
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
        // Derived formula for orthonormal CGA basis:
        //   r² = (e12em + w) / (e12em - w) + cx² + cy²
        //
        // When w = e12em, the denominator (e12em - w) is zero → line (infinite radius).

        let w = self.w();
        let e12em = self.e12em();
        let denom = e12em - w; // Note: this is -(w - e12em)

        if denom.abs() < T::epsilon() {
            // This is a line (w = e12em), infinite radius
            return None;
        }

        // Extract center first
        let center_denom = w - e12em; // = -denom
        let center_x = self.e2epem() / center_denom;
        let center_y = -self.e1epem() / center_denom;

        // Compute radius squared
        let radius_sq = (e12em + w) / denom + center_x * center_x + center_y * center_y;

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

    /// Extracts line parameters (nx, ny, d) when this Circle represents a line.
    ///
    /// Returns `Some((nx, ny, d))` if this is a line (i.e., `is_line()` is true),
    /// where the line equation is `nx*x + ny*y = d`.
    ///
    /// Returns `None` if this is an actual circle, not a line.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim2::{Circle, RoundPoint};
    /// use clifford::ops::InverseSandwich;
    /// use approx::relative_eq;
    ///
    /// // Circle through the origin (inversion center)
    /// let inv_circle = Circle::from_center_radius(0.0_f64, 0.0, 2.0);
    /// let passing_circle = Circle::from_center_radius(2.0, 0.0, 2.0);
    ///
    /// // Invert the circle - should become a line
    /// if let Some(result) = inv_circle.try_inverse_sandwich(&passing_circle) {
    ///     if let Some((nx, ny, d)) = result.to_line_params(1e-10) {
    ///         // Result is a line
    ///         let len = (nx * nx + ny * ny).sqrt();
    ///         assert!(len > 0.0); // Normal is well-defined
    ///     }
    /// }
    /// ```
    pub fn to_line_params(&self, epsilon: T) -> Option<(T, T, T)> {
        if self.is_line(epsilon) {
            // Line components extracted from the circle representation
            // (same mapping as Line::from_two_points)
            Some((self.e12em(), self.e1epem(), self.e2epem()))
        } else {
            None
        }
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
        Self::new_unchecked(circle.e12em(), circle.e1epem(), circle.e2epem())
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
    /// A flat point represents a point using only the em-direction
    /// components of a point pair.
    #[inline]
    pub fn from_euclidean(x: T, y: T) -> Self {
        Self::new_unchecked(x, y, T::one())
    }

    /// Extracts Euclidean coordinates.
    #[inline]
    pub fn to_euclidean(&self) -> Option<(T, T)> {
        if self.epem().abs() < T::epsilon() {
            None
        } else {
            Some((self.e1em() / self.epem(), self.e2em() / self.epem()))
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
            T::zero(), // e1ep
            T::zero(), // e2ep
            T::zero(), // e1em
            T::zero(), // e2em
            T::zero(), // epem
            T::zero(), // ps
        )
    }

    /// Creates a translation motor.
    ///
    /// Translates by vector (dx, dy).
    ///
    /// In CGA, a translation is represented as:
    /// `T = 1 - ½(dx·e₁ + dy·e₂)·e∞ = 1 - ½dx·e₁em - ½dx·e₁ep - ½dy·e₂em - ½dy·e₂ep`
    ///
    /// Since e∞ = ep + em, e₁·e∞ = e₁ep + e₁em, etc.
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
        // Translation: T = 1 - ½·d·e∞ where d = dx·e₁ + dy·e₂
        // e₁·e∞ = e₁·(ep + em) = e₁ep + e₁em
        // So the bivector part is: -½dx·(e₁ep + e₁em) - ½dy·(e₂ep + e₂em)
        Self::new_unchecked(
            T::one(),   // s = 1
            T::zero(),  // m = 0
            -dx * half, // e1ep = -dx/2
            -dy * half, // e2ep = -dy/2
            -dx * half, // e1em = -dx/2
            -dy * half, // e2em = -dy/2
            T::zero(),  // epem = 0
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
            T::zero(),   // e1ep
            T::zero(),   // e2ep
            T::zero(),   // e1em
            T::zero(),   // e2em
            T::zero(),   // epem
            T::zero(),   // ps
        )
    }

    /// Creates a uniform scaling (dilation) motor.
    ///
    /// Scales by factor `k` around the origin.
    ///
    /// In CGA, dilation is represented using hyperbolic functions
    /// in the ep-em plane (e₃₄).
    pub fn from_dilation(scale: T) -> Self {
        let half_log = scale.ln() / T::TWO;
        Self::new_unchecked(
            half_log.cosh(), // s
            T::zero(),       // m
            T::zero(),       // e1ep
            T::zero(),       // e2ep
            T::zero(),       // e1em
            T::zero(),       // e2em
            half_log.sinh(), // epem (dilation is in ep-em plane)
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
            rev.e1ep() / norm_sq,
            rev.e2ep() / norm_sq,
            rev.e1em() / norm_sq,
            rev.e2em() / norm_sq,
            rev.epem() / norm_sq,
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
    fn round_point_is_null_vector() {
        // Verify the null vector property: P·P = 0
        let p = RoundPoint::from_euclidean(3.0_f64, 4.0);
        // Using orthonormal metric: P·P = x² + y² + ep² - em²
        let p_dot_p = p.x() * p.x() + p.y() * p.y() + p.ep() * p.ep() - p.em() * p.em();

        assert!(
            relative_eq!(p_dot_p, 0.0, epsilon = 1e-10),
            "Point should be a null vector: P·P = {} ≠ 0",
            p_dot_p
        );
    }

    #[test]
    fn round_point_origin_is_null_vector() {
        let o = RoundPoint::<f64>::origin();
        let o_dot_o = o.x() * o.x() + o.y() * o.y() + o.ep() * o.ep() - o.em() * o.em();

        assert!(
            relative_eq!(o_dot_o, 0.0, epsilon = 1e-10),
            "Origin should be a null vector: O·O = {} ≠ 0",
            o_dot_o
        );
    }

    #[test]
    fn round_point_infinity_is_null_vector() {
        let inf = RoundPoint::<f64>::infinity();
        let inf_dot_inf =
            inf.x() * inf.x() + inf.y() * inf.y() + inf.ep() * inf.ep() - inf.em() * inf.em();

        assert!(
            relative_eq!(inf_dot_inf, 0.0, epsilon = 1e-10),
            "Infinity should be a null vector: e∞·e∞ = {} ≠ 0",
            inf_dot_inf
        );
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
        assert!(
            result.is_line(1e-10),
            "Collinear points should produce a line (w=0), got w={}",
            result.w()
        );

        // Non-collinear points should produce a circle (w ≠ 0)
        let p1 = RoundPoint::from_euclidean(1.0_f64, 0.0);
        let p2 = RoundPoint::from_euclidean(0.0, 1.0);
        let p3 = RoundPoint::from_euclidean(-1.0, 0.0);

        let result: Circle<f64> = p1.wedge(&p2).wedge(&p3);
        assert!(
            !result.is_line(1e-10),
            "Non-collinear points should produce a circle (w≠0), got w={}",
            result.w()
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
                "Center x mismatch for ({}, {}, {}): expected {}, got {}",
                expected_cx,
                expected_cy,
                expected_r,
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
                "Center y mismatch for ({}, {}, {}): expected {}, got {}",
                expected_cx,
                expected_cy,
                expected_r,
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
                "Radius mismatch for ({}, {}, {}): expected {}, got {}",
                expected_cx,
                expected_cy,
                expected_r,
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
    fn test_inverse_sandwich_circle_inversion() {
        use crate::ops::InverseSandwich;

        // The InverseSandwich operation C × P × C⁻¹ computes circle inversion
        // in the orthonormal CGA basis.
        //
        // Classic circle inversion through a circle of radius r centered at origin:
        //   P' = (r² / |P|²) * P
        //
        // For a point at (4, 0) through a circle of radius 2:
        //   |P| = 4, r = 2, r² = 4
        //   P' = (4 / 16) * (4, 0) = (1, 0)

        let inv_circle = Circle::from_center_radius(0.0_f64, 0.0, 2.0);
        let point = RoundPoint::from_euclidean(4.0, 0.0);

        let inverted = inv_circle.try_inverse_sandwich(&point);

        // The operation should succeed (not return None)
        assert!(
            inverted.is_some(),
            "InverseSandwich should return Some for valid circle and point"
        );

        let (ix, iy) = inverted.unwrap().to_euclidean().unwrap();

        // With the correct orthonormal basis, InverseSandwich gives classic circle inversion
        assert!(
            relative_eq!(ix, 1.0, epsilon = 1e-10),
            "Inverted x should be 1.0, got {}",
            ix
        );
        assert!(
            relative_eq!(iy, 0.0, epsilon = 1e-10),
            "Inverted y should be 0.0, got {}",
            iy
        );
    }

    #[test]
    fn test_circle_from_three_points_various_positions() {
        use std::f64::consts::PI;

        // Test with three points at 0°, 120°, 240° around the circle
        let cx = 2.0_f64;
        let cy = 3.0;
        let r = 5.0;

        let angle1 = 0.0_f64;
        let angle2 = 2.0_f64 * PI / 3.0;
        let angle3 = 4.0_f64 * PI / 3.0;

        let p1 = RoundPoint::from_euclidean(cx + r * angle1.cos(), cy + r * angle1.sin());
        let p2 = RoundPoint::from_euclidean(cx + r * angle2.cos(), cy + r * angle2.sin());
        let p3 = RoundPoint::from_euclidean(cx + r * angle3.cos(), cy + r * angle3.sin());

        let circle = Circle::from_three_points(&p1, &p2, &p3);

        let (extracted_cx, extracted_cy) = circle.center().unwrap();
        let extracted_r = circle.radius().unwrap();

        assert!(
            relative_eq!(extracted_cx, cx, epsilon = 1e-10),
            "Center x: expected {}, got {}",
            cx,
            extracted_cx
        );
        assert!(
            relative_eq!(extracted_cy, cy, epsilon = 1e-10),
            "Center y: expected {}, got {}",
            cy,
            extracted_cy
        );
        assert!(
            relative_eq!(extracted_r, r, epsilon = 1e-10),
            "Radius: expected {}, got {}",
            r,
            extracted_r
        );
    }
}
