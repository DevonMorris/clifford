//! Operations for 2D Projective Geometric Algebra types.

use super::types::{Line, Motor, Point};
use crate::scalar::Float;

// ============================================================================
// Point operations
// ============================================================================

impl<T: Float> Point<T> {
    /// Join of two points: the line through them.
    ///
    /// In point-based PGA, the join is the regressive product P₁ ∨ P₂,
    /// which gives a line (bivector) passing through both points.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let p1 = Point::new(0.0, 0.0);
    /// let p2 = Point::new(1.0, 1.0);
    /// let line = p1.join(&p2);
    ///
    /// // Line x - y = 0
    /// let n = line.normal();
    /// // Normal should be proportional to (1, -1)
    /// ```
    #[inline]
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Regressive product of two vectors in 2D PGA
        // P₁ = x₁e₁ + y₁e₂ + w₁e₀
        // P₂ = x₂e₁ + y₂e₂ + w₂e₀
        //
        // In point-based PGA, join is: dual(dual(P₁) ∧ dual(P₂))
        // But we can compute directly:
        // e₁₂ component: x₁y₂ - y₁x₂
        // e₂₀ component: y₁w₂ - w₁y₂
        // e₀₁ component: w₁x₂ - x₁w₂
        Line {
            e12: self.e1 * other.e2 - self.e2 * other.e1,
            e20: self.e2 * other.e0 - self.e0 * other.e2,
            e01: self.e0 * other.e1 - self.e1 * other.e0,
        }
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
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        (dx * dx + dy * dy).sqrt()
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
        // Normalize both points first
        let p1_norm = self.normalize().unwrap_or(*self);
        let p2_norm = other.normalize().unwrap_or(*other);

        Point::new(
            (p1_norm.e1 + p2_norm.e1) / T::TWO,
            (p1_norm.e2 + p2_norm.e2) / T::TWO,
        )
    }
}

// ============================================================================
// Line operations
// ============================================================================

impl<T: Float> Line<T> {
    /// Meet of two lines: their intersection point.
    ///
    /// In point-based PGA, the meet is the outer product L₁ ∧ L₂,
    /// which gives the intersection point.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::Line;
    /// use approx::abs_diff_eq;
    ///
    /// // x = 0 (y-axis)
    /// let l1: Line<f64> = Line::y_axis();
    /// // y = 0 (x-axis)
    /// let l2: Line<f64> = Line::x_axis();
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
        // Outer product of two bivectors in 2D PGA gives a vector (point)
        // L₁ = d₁e₁₂ + a₁e₂₀ + b₁e₀₁
        // L₂ = d₂e₁₂ + a₂e₂₀ + b₂e₀₁
        //
        // The meet (outer product) for bivectors in grade-2:
        // e₁ coeff: e₂₀ ∧ e₀₁ gives e₂₀₀₁ = -e₀₀₁₂ = -e₁₂ (since e₀₀=0, this = 0)
        // Wait, need to reconsider...
        //
        // Actually in 2D PGA, meet of two lines (bivectors) is computed as:
        // The dual of the outer product of the duals
        // Or directly: P = L₁ × L₂ (cross-like operation)
        //
        // For lines: e₁ = a₁b₂ - b₁a₂, e₂ = d₁a₂ - a₁d₂, e₀ = b₁d₂ - d₁b₂
        // Wait, let me recalculate...
        //
        // Actually the correct formula for meet of two lines in 2D PGA:
        // e₁: e₂₀ ∧ e₀₁ component -> (a₁)(b₂) - (b₁)(a₂) for e₁
        // But bivector ∧ bivector in 2D doesn't directly give vectors...
        //
        // In point-based 2D PGA:
        // Lines are bivectors. To get their meet (intersection point),
        // we need the regressive product or work through the dual.
        //
        // Regressive product: L₁ ∨ L₂ = dual(dual(L₁) ∧ dual(L₂))
        //
        // Let's use explicit calculation:
        // e₁ component = e₂₀·e₀₁ - e₀₁·e₂₀ in the appropriate sense
        //
        // For two lines ax + by + d = 0 and a'x + b'y + d' = 0:
        // Intersection: x = (bd' - b'd)/(ab' - a'b), y = (a'd - ad')/(ab' - a'b)
        //
        // In homogeneous coords:
        // e₁ = bd' - b'd
        // e₂ = a'd - ad'
        // e₀ = ab' - a'b
        Point {
            e1: self.e01 * other.e12 - self.e12 * other.e01,
            e2: self.e12 * other.e20 - self.e20 * other.e12,
            e0: self.e20 * other.e01 - self.e01 * other.e20,
        }
    }

    /// Signed distance from a point to this line.
    ///
    /// Positive on one side, negative on the other.
    #[inline]
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        // For line ax + by + d = 0 and point (px, py):
        // Distance = (a*px + b*py + d) / sqrt(a² + b²)
        let n = self.normal();
        let d = self.e12;
        let (px, py) = p.to_cartesian().unwrap_or((p.e1, p.e2));
        let numerator = n.x * px + n.y * py + d;
        let denominator = (n.x * n.x + n.y * n.y).sqrt();
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
        let d = self.e12;
        let (px, py) = p.to_cartesian().unwrap_or((p.e1, p.e2));

        let norm_sq = n.x * n.x + n.y * n.y;
        if norm_sq.abs() < T::epsilon() {
            return *p;
        }

        // Projection formula: P' = P - ((a*px + b*py + d)/(a² + b²)) * (a, b)
        let t = (n.x * px + n.y * py + d) / norm_sq;
        Point::new(px - t * n.x, py - t * n.y)
    }

    /// Reflects a point across this line.
    #[inline]
    pub fn reflect(&self, p: &Point<T>) -> Point<T> {
        let projected = self.project(p);
        let (px, py) = p.to_cartesian().unwrap_or((p.e1, p.e2));
        let (prx, pry) = projected
            .to_cartesian()
            .unwrap_or((projected.e1, projected.e2));

        // Reflection: P' = 2*projection - P
        Point::new(T::TWO * prx - px, T::TWO * pry - py)
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
    ///
    /// // Parallel lines have angle 0
    /// let parallel = Line::from_implicit(0.0, 1.0, 5.0); // y = -5
    /// assert!(abs_diff_eq!(x_axis.angle(&parallel), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn angle(&self, other: &Line<T>) -> T {
        // Use normal vectors: cos(θ) = |n₁·n₂| / (|n₁||n₂|)
        let n1 = self.normal();
        let n2 = other.normal();

        let dot = n1.x * n2.x + n1.y * n2.y;
        let norm1 = (n1.x * n1.x + n1.y * n1.y).sqrt();
        let norm2 = (n2.x * n2.x + n2.y * n2.y).sqrt();

        let denom = norm1 * norm2;
        if denom < T::epsilon() {
            return T::zero();
        }

        // Use absolute value for acute angle
        let cos_theta = (dot / denom).abs();
        // Clamp to [-1, 1] to handle numerical errors
        let clamped = if cos_theta > T::one() {
            T::one()
        } else {
            cos_theta
        };
        clamped.acos()
    }
}

// ============================================================================
// Motor operations on geometric objects
// ============================================================================

impl<T: Float> Motor<T> {
    /// Transforms a point: `P' = M P M̃`.
    ///
    /// This applies the rigid transformation represented by the motor to the point.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim2::{Point, Motor};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let p = Point::new(1.0, 0.0);
    ///
    /// // 90° rotation
    /// let rotor = Motor::from_rotation(FRAC_PI_2);
    /// let rotated = rotor.transform_point(&p);
    /// assert!(abs_diff_eq!(rotated.x(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    ///
    /// // Translation
    /// let trans = Motor::from_translation(2.0, 3.0);
    /// let translated = trans.transform_point(&p);
    /// assert!(abs_diff_eq!(translated.x(), 3.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(translated.y(), 3.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // Sandwich product: M P M̃
        // For efficiency, we compute this directly rather than using generic multivector ops

        let s = self.s;
        let b12 = self.e12;
        let b20 = self.e20;
        let b01 = self.e01;

        let px = p.e1;
        let py = p.e2;
        let pw = p.e0;

        // Motor reverse: (s, -b12, -b20, -b01)
        // First compute M * P, then * M̃

        // M * P (motor times point)
        // Point is grade 1: e₁, e₂, e₀
        // Motor is grade 0+2: s, e₁₂, e₂₀, e₀₁

        // The full sandwich computation for 2D PGA:
        // This is optimized for the specific case

        // For pure rotation (b20 = b01 = 0):
        // x' = (s² - b12²)x - 2sb12·y
        // y' = 2sb12·x + (s² - b12²)y
        // w' = w

        // For general motor, we need the full computation
        let cos_full = s * s - b12 * b12;
        let sin_full = T::TWO * s * b12;

        // Rotation part
        let rx = cos_full * px - sin_full * py;
        let ry = sin_full * px + cos_full * py;

        // Translation part (acts on the weight)
        // The translation contributes: 2(s*b20 - b12*b01) to e1, 2(s*b01 + b12*b20) to e2
        let tx = T::TWO * (s * b20 - b12 * b01);
        let ty = T::TWO * (s * b01 + b12 * b20);

        Point {
            e1: rx + tx * pw,
            e2: ry + ty * pw,
            e0: pw,
        }
    }

    /// Transforms a line: `L' = M L M̃`.
    ///
    /// This applies the rigid transformation to the line.
    #[inline]
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        // For lines (bivectors), the sandwich product preserves the bivector grade
        // This is more complex to derive, so we use the general formula

        let s = self.s;
        let b12 = self.e12;
        let b20 = self.e20;
        let b01 = self.e01;

        let l12 = l.e12;
        let l20 = l.e20;
        let l01 = l.e01;

        // The e₁₂ component (rotation part) transforms as:
        // e₁₂' = (s² + b12²)e₁₂ - 2(s*b20 + b12*b01)*e₂₀ - 2(s*b01 - b12*b20)*e₀₁
        // Actually, the computation is simpler for e₁₂ which transforms trivially for
        // rotations but picks up translation terms

        // For rotation only (b20 = b01 = 0):
        // e₁₂' = e₁₂ (invariant!)
        // e₂₀' = cos(θ)*e₂₀ - sin(θ)*e₀₁
        // e₀₁' = sin(θ)*e₂₀ + cos(θ)*e₀₁

        let cos_full = s * s - b12 * b12;
        let sin_full = T::TWO * s * b12;

        // Rotation of the normal direction
        let new_e20 = cos_full * l20 - sin_full * l01;
        let new_e01 = sin_full * l20 + cos_full * l01;

        // Translation affects e₁₂ (the distance from origin)
        // d' = d + translation · normal
        let tx = T::TWO * (s * b20 - b12 * b01);
        let ty = T::TWO * (s * b01 + b12 * b20);
        let new_e12 = l12 + tx * l20 + ty * l01;

        Line {
            e12: new_e12,
            e20: new_e20,
            e01: new_e01,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;

    #[test]
    fn point_join_gives_line_through_both() {
        let p1: Point<f64> = Point::new(0.0, 0.0);
        let p2: Point<f64> = Point::new(1.0, 1.0);
        let line = p1.join(&p2);

        // Both points should be on the line (distance ~ 0)
        let d1 = line.distance_to_point(&p1);
        let d2 = line.distance_to_point(&p2);
        assert!(abs_diff_eq!(d1.abs(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d2.abs(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_meet_gives_intersection() {
        // x-axis and y-axis intersect at origin
        let x_axis: Line<f64> = Line::x_axis();
        let y_axis: Line<f64> = Line::y_axis();
        let intersection = x_axis.meet(&y_axis);

        let (x, y) = intersection.to_cartesian().unwrap();
        assert!(abs_diff_eq!(x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_identity_preserves_point() {
        let p: Point<f64> = Point::new(3.0, 4.0);
        let m: Motor<f64> = Motor::identity();
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_rotation_90_degrees() {
        let p = Point::new(1.0, 0.0);
        let m = Motor::from_rotation(std::f64::consts::FRAC_PI_2);
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_translation() {
        let p = Point::new(1.0, 2.0);
        let m = Motor::from_translation(3.0, 4.0);
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), 4.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 6.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_composition() {
        let p = Point::new(1.0, 0.0);

        // First rotate 90°, then translate by (1, 2)
        // In compose, leftmost motor acts first: rotation.compose(&translation) = R * T
        let rotation = Motor::from_rotation(std::f64::consts::FRAC_PI_2);
        let translation = Motor::from_translation(1.0, 2.0);
        let composed = rotation.compose(&translation);

        let result = composed.transform_point(&p);

        // (1,0) -> rotated to (0,1) -> translated to (1,3)
        assert!(abs_diff_eq!(result.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_inverse() {
        let p = Point::new(3.0, 4.0);
        let m = Motor::from_rotation(0.5).compose(&Motor::from_translation(1.0, 2.0));

        let transformed = m.transform_point(&p);
        let back = m.inverse().transform_point(&transformed);

        assert!(abs_diff_eq!(back.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(back.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Unary operation tests
    // ========================================================================

    #[test]
    fn point_geometric_norm() {
        // Point at (3, 4) has distance 5 from origin
        let p = Point::new(3.0, 4.0);
        assert!(abs_diff_eq!(
            p.geometric_norm(),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Origin has geometric norm 0
        let origin: Point<f64> = Point::origin();
        assert!(abs_diff_eq!(
            origin.geometric_norm(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn point_bulk_weight_norm() {
        let p = Point::new(3.0, 4.0);
        assert!(abs_diff_eq!(p.bulk_norm(), 5.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(
            p.weight_norm(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(p.attitude(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_geometric_norm() {
        // X axis at y=0 has distance 0 from origin
        let x_axis: Line<f64> = Line::x_axis();
        assert!(abs_diff_eq!(
            x_axis.geometric_norm(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Line y = 3 has distance 3 from origin
        let offset_line = Line::from_implicit(0.0, 1.0, -3.0);
        let unitized = offset_line.unitized();
        assert!(abs_diff_eq!(
            unitized.geometric_norm(),
            3.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_bulk_weight_norm() {
        // Y axis: normal (1, 0), distance 0
        let y_axis: Line<f64> = Line::y_axis();
        assert!(abs_diff_eq!(
            y_axis.weight_norm(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            y_axis.bulk_norm(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Attitude should be the normal direction
        let att = y_axis.attitude();
        assert!(abs_diff_eq!(att.x, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(att.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_reverse() {
        let line = Line::new(1.0, 2.0, 3.0);
        let rev = line.reverse();

        // Reverse negates all components
        assert!(abs_diff_eq!(rev.e12, -1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e20, -2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e01, -3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_unitized() {
        // Non-unit line
        let line: Line<f64> = Line::from_implicit(3.0, 4.0, 5.0);
        let unitized = line.unitized();

        // Weight norm should be 1
        assert!(abs_diff_eq!(
            unitized.weight_norm(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Normal should be unit length
        let n = unitized.normal();
        let norm = (n.x * n.x + n.y * n.y).sqrt();
        assert!(abs_diff_eq!(norm, 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Angle tests
    // ========================================================================

    #[test]
    fn line_angle_perpendicular() {
        let x_axis: Line<f64> = Line::x_axis();
        let y_axis: Line<f64> = Line::y_axis();
        assert!(abs_diff_eq!(
            x_axis.angle(&y_axis),
            std::f64::consts::FRAC_PI_2,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_angle_parallel() {
        let x_axis: Line<f64> = Line::x_axis();
        let parallel = Line::from_implicit(0.0, 1.0, 5.0); // y = -5
        assert!(abs_diff_eq!(
            x_axis.angle(&parallel),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_angle_symmetric() {
        let l1: Line<f64> = Line::from_implicit(1.0, 2.0, 3.0);
        let l2: Line<f64> = Line::from_implicit(2.0, -1.0, 1.0);
        assert!(abs_diff_eq!(
            l1.angle(&l2),
            l2.angle(&l1),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_angle_range() {
        let l1: Line<f64> = Line::from_implicit(1.0, 1.0, 0.0);
        let l2: Line<f64> = Line::from_implicit(1.0, -1.0, 0.0);
        let angle = l1.angle(&l2);
        // Angle should be in [0, π/2]
        assert!(angle >= 0.0);
        assert!(angle <= std::f64::consts::FRAC_PI_2 + ABS_DIFF_EQ_EPS);
    }
}
