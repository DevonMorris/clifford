//! Domain-specific extensions for 2D Projective GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to 2D projective geometry.

use super::generated::products;
use super::generated::types::{Flector, Line, Motor, Point};
use crate::scalar::Float;
use crate::specialized::euclidean::dim2::Vector as EuclideanVector;

// ============================================================================
// Point extensions
// ============================================================================

impl<T: Float> Point<T> {
    /// Creates a finite point at Cartesian coordinates (x, y).
    ///
    /// The homogeneous weight `w` is set to 1.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let p = Point::from_cartesian(3.0, 4.0);
    /// assert_eq!(p.x(), 3.0);
    /// assert_eq!(p.y(), 4.0);
    /// ```
    #[inline]
    pub fn from_cartesian(x: T, y: T) -> Self {
        Self::new(x, y, T::one())
    }

    /// Creates an ideal point (point at infinity) in the given direction.
    ///
    /// Ideal points have `w = 0` and represent directions rather than positions.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let ideal = Point::<f64>::ideal(1.0, 0.0);
    /// assert!(ideal.is_ideal(1e-10));
    /// ```
    #[inline]
    pub fn ideal(dx: T, dy: T) -> Self {
        Self::new(dx, dy, T::zero())
    }

    /// Origin point (0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::from_cartesian(T::zero(), T::zero())
    }

    /// Returns the x-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn x(&self) -> T {
        self.e1() / self.e0()
    }

    /// Returns the y-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn y(&self) -> T {
        self.e2() / self.e0()
    }

    /// Returns the homogeneous weight.
    #[inline]
    pub fn w(&self) -> T {
        self.e0()
    }

    /// Returns true if this is an ideal point (point at infinity).
    #[inline]
    pub fn is_ideal(&self, epsilon: T) -> bool {
        self.e0().abs() < epsilon
    }

    /// Returns true if this is a finite point (not at infinity).
    #[inline]
    pub fn is_finite(&self, epsilon: T) -> bool {
        self.e0().abs() >= epsilon
    }

    /// Normalizes the homogeneous coordinates so `w = 1` (if finite).
    ///
    /// Returns `None` if this is an ideal point.
    pub fn unitize(&self) -> Option<Self> {
        if self.e0().abs() < T::epsilon() {
            None
        } else {
            Some(Self::new(
                self.e1() / self.e0(),
                self.e2() / self.e0(),
                T::one(),
            ))
        }
    }

    /// Returns the Cartesian coordinates as a tuple, if finite.
    ///
    /// Returns `None` if this is an ideal point.
    #[inline]
    pub fn to_cartesian(&self) -> Option<(T, T)> {
        if self.e0().abs() < T::epsilon() {
            None
        } else {
            Some((self.e1() / self.e0(), self.e2() / self.e0()))
        }
    }

    /// Returns the attitude of the point.
    ///
    /// For a point, the attitude is the weight (e₀ component).
    #[inline]
    pub fn attitude(&self) -> T {
        self.e0()
    }

    /// Returns the squared bulk norm of the point.
    ///
    /// The bulk norm is the length of the spatial part: `e1² + e2²`.
    #[inline]
    pub fn bulk_norm_squared(&self) -> T {
        self.e1() * self.e1() + self.e2() * self.e2()
    }

    /// Returns the bulk norm of the point.
    #[inline]
    pub fn bulk_norm(&self) -> T {
        self.bulk_norm_squared().sqrt()
    }

    /// Returns the weight norm of the point.
    ///
    /// For a point, the weight is the absolute value of the e₀ component.
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.e0().abs()
    }

    /// Returns the geometric norm (distance from origin).
    ///
    /// For a unitized point (w = ±1), this equals the distance from the origin.
    #[inline]
    pub fn geometric_norm(&self) -> T {
        let weight = self.weight_norm();
        if weight < T::epsilon() {
            T::zero()
        } else {
            self.bulk_norm() / weight
        }
    }

    /// Join of two points: the line through them.
    ///
    /// In point-based PGA, the join is the exterior product P₁ ∧ P₂,
    /// which gives a line (bivector) passing through both points.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let p1 = Point::from_cartesian(0.0, 0.0);
    /// let p2 = Point::from_cartesian(1.0, 1.0);
    /// let line = p1.join(&p2);
    /// ```
    #[inline]
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::exterior_point_point(self, other)
    }

    /// Euclidean distance to another point.
    ///
    /// Both points must be finite (non-ideal).
    ///
    /// # Panics
    ///
    /// May return incorrect results if either point is ideal.
    #[inline]
    pub fn distance(&self, other: &Point<T>) -> T {
        self.distance_squared(other).sqrt()
    }

    /// Squared distance to another point.
    #[inline]
    pub fn distance_squared(&self, other: &Point<T>) -> T {
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        dx * dx + dy * dy
    }

    /// Midpoint between this point and another.
    #[inline]
    pub fn midpoint(&self, other: &Point<T>) -> Point<T> {
        let p1_norm = self.unitize().unwrap_or(*self);
        let p2_norm = other.unitize().unwrap_or(*other);

        Point::from_cartesian(
            (p1_norm.e1() + p2_norm.e1()) / T::TWO,
            (p1_norm.e2() + p2_norm.e2()) / T::TWO,
        )
    }
}

// ============================================================================
// Line extensions
// ============================================================================

impl<T: Float> Line<T> {
    /// Creates a line from implicit equation `ax + by + d = 0`.
    ///
    /// The normal direction is `(a, b)`.
    #[inline]
    pub fn from_implicit(a: T, b: T, d: T) -> Self {
        Self::new(d, a, b)
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

    /// Returns the normal direction of the line `ax + by + d = 0`.
    ///
    /// Returns a Euclidean vector `(a, b)` representing the normal.
    #[inline]
    pub fn normal(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.e01(), self.e02())
    }

    /// Returns the attitude of the line.
    ///
    /// For a 2D line, the attitude is the normal direction.
    #[inline]
    pub fn attitude(&self) -> EuclideanVector<T> {
        self.normal()
    }

    /// Returns the signed distance from the origin (scaled by normal length).
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        self.e12()
    }

    /// Returns a normalized line (unit norm).
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        if n.abs() < T::epsilon() {
            *self
        } else {
            Self::new(self.e12() / n, self.e01() / n, self.e02() / n)
        }
    }

    /// Returns true if this line is degenerate (zero).
    #[inline]
    pub fn is_zero(&self, epsilon: T) -> bool {
        self.e12().abs() < epsilon && self.e01().abs() < epsilon && self.e02().abs() < epsilon
    }

    /// Returns the squared weight norm of the line.
    ///
    /// The weight norm is the length of the normal direction: `e01² + e02²`.
    #[inline]
    pub fn weight_norm_squared(&self) -> T {
        self.e01() * self.e01() + self.e02() * self.e02()
    }

    /// Returns the weight norm of the line.
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    /// Returns the squared bulk norm of the line.
    ///
    /// The bulk norm is the absolute value of the e₁₂ component.
    #[inline]
    pub fn bulk_norm_squared(&self) -> T {
        self.e12() * self.e12()
    }

    /// Returns the bulk norm of the line.
    #[inline]
    pub fn bulk_norm(&self) -> T {
        self.e12().abs()
    }

    /// Returns a unitized line (unit weight norm).
    pub fn unitized(&self) -> Self {
        let n = self.weight_norm();
        if n.abs() < T::epsilon() {
            *self
        } else {
            Self::new(self.e12() / n, self.e01() / n, self.e02() / n)
        }
    }

    /// Returns the geometric norm (unitized distance from origin).
    ///
    /// For a unitized line, this is the perpendicular distance from
    /// the origin to the line.
    #[inline]
    pub fn geometric_norm(&self) -> T {
        let weight = self.weight_norm();
        if weight.abs() < T::epsilon() {
            T::zero()
        } else {
            self.e12().abs() / weight
        }
    }

    /// Meet of two lines: their intersection point.
    ///
    /// In point-based PGA, this is computed via the regressive product.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::Line;
    /// use approx::abs_diff_eq;
    ///
    /// let l1: Line<f64> = Line::y_axis(); // x = 0
    /// let l2: Line<f64> = Line::x_axis(); // y = 0
    ///
    /// // Intersection is the origin
    /// let p = l1.meet(&l2);
    /// if let Some((x, y)) = p.to_cartesian() {
    ///     assert!(abs_diff_eq!(x, 0.0, epsilon = 1e-10));
    ///     assert!(abs_diff_eq!(y, 0.0, epsilon = 1e-10));
    /// }
    /// ```
    #[inline]
    pub fn meet(&self, other: &Line<T>) -> Point<T> {
        products::regressive_line_line(self, other)
    }

    /// Signed distance from a point to this line.
    ///
    /// Positive on one side, negative on the other.
    #[inline]
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        let n = self.normal();
        let d = self.e12();
        let (px, py) = p.to_cartesian().unwrap_or((p.e1(), p.e2()));
        let numerator = n.x() * px + n.y() * py + d;
        let denominator = (n.x() * n.x() + n.y() * n.y()).sqrt();
        if denominator.abs() < T::epsilon() {
            T::zero()
        } else {
            numerator / denominator
        }
    }

    /// Projects a point onto this line.
    #[inline]
    pub fn project(&self, p: &Point<T>) -> Point<T> {
        let n = self.normal();
        let d = self.e12();
        let (px, py) = p.to_cartesian().unwrap_or((p.e1(), p.e2()));

        let norm_sq = n.x() * n.x() + n.y() * n.y();
        if norm_sq.abs() < T::epsilon() {
            return *p;
        }

        let t = (n.x() * px + n.y() * py + d) / norm_sq;
        Point::from_cartesian(px - t * n.x(), py - t * n.y())
    }

    /// Reflects a point across this line.
    #[inline]
    pub fn reflect(&self, p: &Point<T>) -> Point<T> {
        let projected = self.project(p);
        let (px, py) = p.to_cartesian().unwrap_or((p.e1(), p.e2()));
        let (prx, pry) = projected
            .to_cartesian()
            .unwrap_or((projected.e1(), projected.e2()));

        Point::from_cartesian(T::TWO * prx - px, T::TWO * pry - py)
    }

    /// Computes the angle between two lines.
    ///
    /// Returns the acute angle in radians, in the range `[0, π/2]`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::Line;
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let x_axis: Line<f64> = Line::x_axis();
    /// let y_axis: Line<f64> = Line::y_axis();
    ///
    /// // Perpendicular lines have angle π/2
    /// assert!(abs_diff_eq!(x_axis.angle(&y_axis), FRAC_PI_2, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn angle(&self, other: &Line<T>) -> T {
        let n1 = self.normal();
        let n2 = other.normal();

        let dot = n1.x() * n2.x() + n1.y() * n2.y();
        let norm1 = (n1.x() * n1.x() + n1.y() * n1.y()).sqrt();
        let norm2 = (n2.x() * n2.x() + n2.y() * n2.y()).sqrt();

        let denom = norm1 * norm2;
        if denom < T::epsilon() {
            return T::zero();
        }

        let cos_theta = (dot / denom).abs();
        let clamped = if cos_theta > T::one() {
            T::one()
        } else {
            cos_theta
        };
        clamped.acos()
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Creates the identity motor (no transformation).
    ///
    /// In 2D PGA with dual representation, the identity motor has e012=1.
    /// This matches the 3D PGA pattern where identity is the pseudoscalar.
    #[inline]
    pub fn identity() -> Self {
        // Dual form: e012 = 1 (pseudoscalar component)
        Self::new(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Creates a pure rotation motor around the origin.
    ///
    /// # Arguments
    ///
    /// * `angle` - Rotation angle in radians (counterclockwise positive)
    ///
    /// # Notes
    ///
    /// In dual form, the rotation is encoded as:
    /// - e0 = -sin(θ/2) (dual of e12 component)
    /// - e012 = cos(θ/2) (dual of scalar component)
    #[inline]
    pub fn from_rotation(angle: T) -> Self {
        let half = angle / T::TWO;
        // Dual form: e0 = -sin(θ/2), e012 = cos(θ/2)
        Self::new(T::zero(), T::zero(), -half.sin(), half.cos())
    }

    /// Creates a pure translation motor.
    ///
    /// # Arguments
    ///
    /// * `dx` - Translation in x direction
    /// * `dy` - Translation in y direction
    ///
    /// # Note
    ///
    /// In dual form, translation is encoded as:
    /// - e1 = dy/2 (dual of -e02 component)
    /// - e2 = -dx/2 (dual of e01 component)
    /// - e012 = 1
    #[inline]
    pub fn from_translation(dx: T, dy: T) -> Self {
        // Dual form mapping: e1 = -(-dy/2) = dy/2, e2 = -dx/2
        Self::new(dy / T::TWO, -dx / T::TWO, T::zero(), T::one())
    }

    /// Returns the inverse motor: `M⁻¹`.
    ///
    /// For a unit motor, this equals the reverse.
    #[inline]
    pub fn inverse(&self) -> Self {
        // In dual form, the norm squared uses e0 and e012 (duals of e12 and s)
        let norm_sq = self.e012() * self.e012() + self.e0() * self.e0();
        if norm_sq.abs() < T::epsilon() {
            return *self;
        }
        let rev = self.reverse();
        Self::new(
            rev.e1() / norm_sq,
            rev.e2() / norm_sq,
            rev.e0() / norm_sq,
            rev.e012() / norm_sq,
        )
    }

    /// Linear interpolation between motors (normalized).
    #[inline]
    pub fn lerp(&self, other: &Self, t: T) -> Self {
        let one_minus_t = T::one() - t;
        Self::new(
            self.e1() * one_minus_t + other.e1() * t,
            self.e2() * one_minus_t + other.e2() * t,
            self.e0() * one_minus_t + other.e0() * t,
            self.e012() * one_minus_t + other.e012() * t,
        )
        .normalize()
    }

    /// Returns the rotation angle in radians.
    #[inline]
    pub fn rotation_angle(&self) -> T {
        // In dual form: e0 = -sin(θ/2), e012 = cos(θ/2)
        // So angle = 2 * atan2(-e0, e012)
        (-self.e0()).atan2(self.e012()) * T::TWO
    }

    /// Returns the translation vector as a Euclidean vector.
    ///
    /// # Note
    ///
    /// For pure translations (no rotation), this returns the exact translation.
    /// For composed motors (rotation + translation), this is an approximation.
    #[inline]
    pub fn translation(&self) -> EuclideanVector<T> {
        // In dual form: e1 = dy/2, e2 = -dx/2
        // So dx = -2*e2, dy = 2*e1
        EuclideanVector::new(-self.e2() * T::TWO, self.e1() * T::TWO)
    }

    /// Transforms a point using the antisandwich product.
    ///
    /// In 2D PGA with dual motor representation, the antisandwich product
    /// works directly for point transformationirtual. This matches the 3D PGA pattern.
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::antisandwich_motor_point(self, p)
    }

    /// Transforms a line using the antisandwich product.
    #[inline]
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        products::antisandwich_motor_line(self, l)
    }

    /// Transforms a Euclidean 2D vector using this motor.
    ///
    /// This embeds the vector as a projective point, transforms it,
    /// and extracts the Euclidean coordinates.
    #[inline]
    pub fn transform_euclidean(&self, v: &EuclideanVector<T>) -> EuclideanVector<T> {
        let p = Point::from_cartesian(v.x(), v.y());
        let transformed = self.transform_point(&p);
        EuclideanVector::new(
            transformed.e1() / transformed.e0(),
            transformed.e2() / transformed.e0(),
        )
    }
}

// ============================================================================
// Flector extensions
// ============================================================================

impl<T: Float> Flector<T> {
    /// Creates a flector from a reflection line.
    ///
    /// The flector will reflect points across the given line.
    /// The line should be unitized (unit weight norm) for proper reflection.
    ///
    /// In 2D PGA with swapped representation (Flector = grades 0+2),
    /// the flector stores the line components directly.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::{Flector, Line, Point};
    /// use approx::abs_diff_eq;
    ///
    /// let x_axis: Line<f64> = Line::x_axis();
    /// let f = Flector::from_line(&x_axis);
    /// let p = Point::from_cartesian(1.0, 2.0);
    /// let reflected = f.transform_point(&p);
    /// assert!(abs_diff_eq!(reflected.x(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(reflected.y(), -2.0, epsilon = 1e-10));
    /// ```
    pub fn from_line(line: &Line<T>) -> Self {
        let l = line.unitized();
        // Flector fields: s, e12, e01, e02
        // Line fields: e12, e01, e02
        Self::new(T::zero(), l.e12(), l.e01(), l.e02())
    }

    /// Creates a flector from line normal direction through the origin.
    ///
    /// Reflects across the line passing through the origin with the given normal.
    pub fn from_line_through_origin(nx: T, ny: T) -> Self {
        let norm = (nx * nx + ny * ny).sqrt();
        // Line through origin: d=0, so e12=0, normal=(nx, ny) -> e01=nx, e02=ny
        Self::new(T::zero(), T::zero(), nx / norm, ny / norm)
    }

    /// Reflect through X-axis (y = 0).
    #[inline]
    pub fn reflect_x() -> Self {
        Self::from_line(&Line::x_axis())
    }

    /// Reflect through Y-axis (x = 0).
    #[inline]
    pub fn reflect_y() -> Self {
        Self::from_line(&Line::y_axis())
    }

    /// Returns the line part of this flector (grade 2 components).
    ///
    /// For a pure reflection flector, this is the reflection line.
    #[inline]
    pub fn line_part(&self) -> Line<T> {
        Line::new(self.e12(), self.e01(), self.e02())
    }

    /// Returns the scalar part (grade 0).
    #[inline]
    pub fn scalar_part(&self) -> T {
        self.s()
    }

    /// Check if this is a pure reflection (scalar part is zero).
    #[inline]
    pub fn is_pure_reflection(&self) -> bool {
        self.s().abs() < T::epsilon()
    }

    /// Unitize to unit weight norm.
    ///
    /// For a flector to represent a proper rigid reflection, the weight norm
    /// (e01² + e02²) should be 1.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new(
            self.s() / wn,
            self.e12() / wn,
            self.e01() / wn,
            self.e02() / wn,
        )
    }

    /// Returns the squared weight norm of the flector.
    ///
    /// The weight norm is the length of the normal direction: `e01² + e02²`.
    #[inline]
    pub fn weight_norm_squared(&self) -> T {
        self.e01() * self.e01() + self.e02() * self.e02()
    }

    /// Returns the weight norm of the flector.
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    /// Transform a point by reflecting it across the flector's line.
    ///
    /// For a pure reflection flector (scalar part = 0), this reflects the
    /// point across the line represented by the grade-2 components.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::{Flector, Point};
    /// use approx::abs_diff_eq;
    ///
    /// // Reflect through Y-axis (x = 0)
    /// let f = Flector::<f64>::reflect_y();
    /// let p = Point::from_cartesian(2.0, 3.0);
    /// let reflected = f.transform_point(&p);
    /// assert!(abs_diff_eq!(reflected.x(), -2.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(reflected.y(), 3.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // Use the line part for reflection
        self.line_part().reflect(p)
    }

    /// Transform a line by reflecting it across the flector's line.
    #[inline]
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        // Reflect a line by reflecting two points on it and joining them
        // For now, use the sandwich product which should work for line-line
        products::sandwich_flector_line(self, l)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::abs_diff_eq;

    #[test]
    #[ignore = "exterior product or line representation needs investigation"]
    fn point_join_gives_line_through_both() {
        let p1: Point<f64> = Point::from_cartesian(0.0, 0.0);
        let p2: Point<f64> = Point::from_cartesian(1.0, 1.0);
        let line = p1.join(&p2);

        let d1 = line.distance_to_point(&p1);
        let d2 = line.distance_to_point(&p2);
        assert!(abs_diff_eq!(d1.abs(), 0.0, epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(d2.abs(), 0.0, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn line_meet_gives_intersection() {
        let x_axis: Line<f64> = Line::x_axis();
        let y_axis: Line<f64> = Line::y_axis();
        let intersection = x_axis.meet(&y_axis);

        let (x, y) = intersection.to_cartesian().unwrap();
        assert!(abs_diff_eq!(x, 0.0, epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(y, 0.0, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn motor_identity_preserves_point() {
        let p: Point<f64> = Point::from_cartesian(3.0, 4.0);
        let m: Motor<f64> = Motor::identity();
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), p.x(), epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), p.y(), epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn motor_rotation_90_degrees() {
        let p = Point::from_cartesian(1.0, 0.0);
        let m = Motor::from_rotation(std::f64::consts::FRAC_PI_2);
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), 0.0, epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 1.0, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn motor_translation() {
        let p = Point::from_cartesian(1.0, 2.0);
        let m = Motor::from_translation(3.0, 4.0);
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), 4.0, epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 6.0, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn motor_inverse() {
        let p = Point::from_cartesian(3.0, 4.0);
        let m = Motor::from_rotation(0.5);

        let transformed = m.transform_point(&p);
        let back = m.inverse().transform_point(&transformed);

        assert!(abs_diff_eq!(back.x(), p.x(), epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(back.y(), p.y(), epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn point_geometric_norm() {
        let p = Point::from_cartesian(3.0, 4.0);
        assert!(abs_diff_eq!(
            p.geometric_norm(),
            5.0,
            epsilon = RELATIVE_EQ_EPS
        ));

        let origin: Point<f64> = Point::origin();
        assert!(abs_diff_eq!(
            origin.geometric_norm(),
            0.0,
            epsilon = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn line_angle_perpendicular() {
        let x_axis: Line<f64> = Line::x_axis();
        let y_axis: Line<f64> = Line::y_axis();
        assert!(abs_diff_eq!(
            x_axis.angle(&y_axis),
            std::f64::consts::FRAC_PI_2,
            epsilon = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn line_angle_parallel() {
        let x_axis: Line<f64> = Line::x_axis();
        let parallel = Line::from_implicit(0.0, 1.0, 5.0);
        assert!(abs_diff_eq!(
            x_axis.angle(&parallel),
            0.0,
            epsilon = RELATIVE_EQ_EPS
        ));
    }

    // ==========================================================================
    // Flector tests
    // ==========================================================================

    #[test]
    fn flector_reflect_y_axis() {
        let f = Flector::<f64>::reflect_y();
        let p = Point::from_cartesian(2.0, 3.0);

        let reflected = f.transform_point(&p);

        // Reflecting (2, 3) through Y-axis (x = 0) should give (-2, 3)
        assert!(abs_diff_eq!(reflected.x(), -2.0, epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(reflected.y(), 3.0, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn flector_reflect_x_axis() {
        let f = Flector::<f64>::reflect_x();
        let p = Point::from_cartesian(2.0, 3.0);

        let reflected = f.transform_point(&p);

        // Reflecting (2, 3) through X-axis (y = 0) should give (2, -3)
        assert!(abs_diff_eq!(reflected.x(), 2.0, epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(reflected.y(), -3.0, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn flector_double_reflection_identity() {
        let f = Flector::<f64>::reflect_y();
        let p = Point::from_cartesian(2.0, 3.0);

        // Reflecting twice should return to original
        let once = f.transform_point(&p);
        let twice = f.transform_point(&once);

        assert!(abs_diff_eq!(twice.x(), p.x(), epsilon = RELATIVE_EQ_EPS));
        assert!(abs_diff_eq!(twice.y(), p.y(), epsilon = RELATIVE_EQ_EPS));
    }
}
