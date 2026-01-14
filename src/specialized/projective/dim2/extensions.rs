//! Domain-specific extensions for 2D Projective GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to 2D projective geometry.

use super::generated::types::{Flector, Line, Motor, Point, Trivector};
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
    /// assert_eq!(p.cartesian_x(), 3.0);
    /// assert_eq!(p.cartesian_y(), 4.0);
    /// ```
    #[inline]
    pub fn from_cartesian(x: T, y: T) -> Self {
        Self::new(x, y, T::one())
    }

    /// Origin point (0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::from_cartesian(T::zero(), T::zero())
    }

    /// Returns the Cartesian x-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn cartesian_x(&self) -> T {
        self.x() / self.w()
    }

    /// Returns the Cartesian y-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn cartesian_y(&self) -> T {
        self.y() / self.w()
    }

    /// Returns the Cartesian coordinates as a tuple, if finite.
    ///
    /// Returns `None` if this is an ideal point.
    #[inline]
    pub fn to_cartesian(&self) -> Option<(T, T)> {
        if self.w().abs() < T::epsilon() {
            None
        } else {
            Some((self.x() / self.w(), self.y() / self.w()))
        }
    }

    /// Euclidean distance to another finite point.
    pub fn distance(&self, other: &Point<T>) -> T {
        self.distance_squared(other).sqrt()
    }

    /// Squared Euclidean distance.
    pub fn distance_squared(&self, other: &Point<T>) -> T {
        let dx = self.cartesian_x() - other.cartesian_x();
        let dy = self.cartesian_y() - other.cartesian_y();
        dx * dx + dy * dy
    }

    /// Midpoint between two finite points.
    pub fn midpoint(&self, other: &Point<T>) -> Point<T> {
        let two = T::TWO;
        Point::new(
            (self.x() * other.w() + other.x() * self.w()) / two,
            (self.y() * other.w() + other.y() * self.w()) / two,
            self.w() * other.w(),
        )
    }
}

// ============================================================================
// Line extensions
// ============================================================================

impl<T: Float> Line<T> {
    /// Creates a line through two points (their join/wedge product).
    #[inline]
    pub fn through_points(p: &Point<T>, q: &Point<T>) -> Self {
        crate::ops::Wedge::wedge(p, q)
    }

    /// Creates a line through a point in the given direction.
    pub fn from_point_and_direction(point: &Point<T>, direction: &EuclideanVector<T>) -> Self {
        // Ideal point (point at infinity) in the given direction
        let ideal = Point::new(direction.x(), direction.y(), T::zero());
        crate::ops::Wedge::wedge(point, &ideal)
    }

    /// Creates a line from normal (a, b) and distance c: ax + by + c = 0.
    ///
    /// The line equation is: mx*a + my*b + d*c = 0
    #[inline]
    pub fn from_equation(a: T, b: T, c: T) -> Self {
        Self::new(a, b, c)
    }

    /// X-axis (line y = 0).
    ///
    /// The x-axis (y = 0) has normal vector (0, 1), so it uses the my (e01) component.
    #[inline]
    pub fn x_axis() -> Self {
        Self::new(T::zero(), T::one(), T::zero())
    }

    /// Y-axis (line x = 0).
    ///
    /// The y-axis (x = 0) has normal vector (1, 0), so it uses the d (e02) component.
    #[inline]
    pub fn y_axis() -> Self {
        Self::new(T::zero(), T::zero(), T::one())
    }

    /// Normal vector (normal_x, normal_y).
    #[inline]
    pub fn normal(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.normal_x(), self.normal_y())
    }

    /// Distance from origin (dist component).
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        self.dist()
    }

    /// Check if line passes through origin.
    pub fn through_origin(&self, epsilon: T) -> bool {
        self.dist().abs() < epsilon
    }

    /// Check if two lines are parallel.
    pub fn is_parallel(&self, other: &Line<T>, epsilon: T) -> bool {
        // Lines are parallel if their normals are parallel
        let cross = self.normal_x() * other.normal_y() - self.normal_y() * other.normal_x();
        cross.abs() < epsilon
    }

    /// Angle between two lines (in radians).
    pub fn angle(&self, other: &Line<T>) -> T {
        use crate::norm::DegenerateNormed;
        let n1 = self.try_unitize().unwrap_or(*self);
        let n2 = other.try_unitize().unwrap_or(*other);
        let cos_angle = (n1.normal_x() * n2.normal_x() + n1.normal_y() * n2.normal_y())
            .abs()
            .min(T::one());
        cos_angle.acos()
    }

    /// Signed distance from a point to this line.
    ///
    /// The line should be unitized for accurate results.
    pub fn signed_distance(&self, point: &Point<T>) -> T {
        use crate::norm::DegenerateNormed;
        let l = self.try_unitize().unwrap_or(*self);
        (l.normal_x() * point.x() + l.normal_y() * point.y() + l.dist() * point.w()) / point.w()
    }

    /// Distance from a point to this line.
    pub fn distance_to_point(&self, point: &Point<T>) -> T {
        self.signed_distance(point).abs()
    }

    /// Returns true if the line is degenerate (zero normal).
    #[inline]
    pub fn is_zero(&self, epsilon: T) -> bool {
        use crate::norm::DegenerateNormed;
        self.weight_norm() < epsilon
    }

    /// Project a point onto this line.
    pub fn project_point(&self, point: &Point<T>) -> Point<T> {
        let dist = self.signed_distance(point);
        let n = self.normal();
        let n_norm_sq = n.x() * n.x() + n.y() * n.y();
        if n_norm_sq < T::epsilon() {
            return *point;
        }
        Point::new(
            point.x() - dist * n.x() * point.w() / n_norm_sq.sqrt(),
            point.y() - dist * n.y() * point.w() / n_norm_sq.sqrt(),
            point.w(),
        )
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Identity motor (leaves all elements unchanged).
    ///
    /// In 2D PGA with the odd subalgebra (grades 1 and 3), the identity motor is the
    /// pseudoscalar e012 = 1, which is the identity for the antiproduct.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim2::{Motor, Point};
    /// use clifford::ops::Transform;
    ///
    /// let m = Motor::<f64>::identity();
    /// let p = Point::from_cartesian(1.0, 2.0);
    /// let p2 = m.transform(&p);
    /// assert!((p2.cartesian_x() - p.cartesian_x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        // ps = 1 is the identity for the antiproduct
        Self::new(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Pure translation motor.
    ///
    /// Creates a motor that translates by the vector (dx, dy).
    ///
    /// In 2D PGA with Motor grades [1, 3], the antisandwich formula
    /// requires these sign conventions to achieve correct translation.
    pub fn from_translation(dx: T, dy: T) -> Self {
        let half = T::one() / T::TWO;
        Self::new(
            dy * half,  // x
            -dx * half, // y
            T::zero(),  // w
            T::one(),   // ps (identity part)
        )
    }

    /// Pure rotation motor around the origin.
    ///
    /// Creates a motor that rotates by the given angle (in radians).
    /// Positive angles rotate counter-clockwise (right-hand rule).
    pub fn from_rotation(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new(
            T::zero(),   // x
            T::zero(),   // y
            -half.sin(), // w
            half.cos(),  // ps
        )
    }

    /// Rotation around a point.
    ///
    /// Creates a motor that rotates by the given angle around the specified center point.
    pub fn from_rotation_around(center: &Point<T>, angle: T) -> Self {
        use crate::ops::Versor;
        let to_origin = Motor::from_translation(-center.cartesian_x(), -center.cartesian_y());
        let rotation = Motor::from_rotation(angle);
        let from_origin = Motor::from_translation(center.cartesian_x(), center.cartesian_y());

        // Compose: first translate to origin, then rotate, then translate back
        from_origin.compose(&rotation.compose(&to_origin))
    }

    /// Inverse motor.
    ///
    /// For a unit motor, the inverse is the antireverse.
    pub fn inverse(&self) -> Self {
        use crate::norm::DegenerateNormed;
        let norm_sq = self.weight_norm_squared();
        if norm_sq < T::epsilon() {
            return *self;
        }
        let rev = self.antireverse();
        Self::new(
            rev.x() / norm_sq,
            rev.y() / norm_sq,
            rev.w() / norm_sq,
            rev.ps() / norm_sq,
        )
    }

    /// Extract rotation angle.
    pub fn rotation_angle(&self) -> T {
        let cos_half = self.ps().min(T::one()).max(-T::one());
        let sin_half = self.w();
        sin_half.atan2(cos_half) * T::TWO
    }

    /// Extract translation vector.
    ///
    /// For a pure translation motor, this returns the translation vector.
    pub fn translation(&self) -> EuclideanVector<T> {
        // Inverse of from_translation encoding:
        // from_translation sets: x = dy/2, y = -dx/2
        // So: dx = -2*y, dy = 2*x
        EuclideanVector::new(-T::TWO * self.y(), T::TWO * self.x())
    }
}

// ============================================================================
// Flector extensions
// ============================================================================

impl<T: Float> Flector<T> {
    /// Create flector from reflection line.
    ///
    /// In 2D PGA with the even subalgebra (grades 0 and 2), flectors represent
    /// reflections through lines. The flector uses the same bivector components
    /// as the line.
    pub fn from_line(line: &Line<T>) -> Self {
        use crate::norm::DegenerateNormed;
        let l = line.try_unitize().unwrap_or(*line);
        Self::new(T::zero(), l.normal_x(), l.normal_y(), l.dist())
    }

    /// Reflect through y-axis (x = 0).
    #[inline]
    pub fn reflect_y_axis() -> Self {
        Self::from_line(&Line::y_axis())
    }

    /// Reflect through x-axis (y = 0).
    #[inline]
    pub fn reflect_x_axis() -> Self {
        Self::from_line(&Line::x_axis())
    }

    /// Reflect through line at angle θ from x-axis through origin.
    pub fn reflect_through_angle(angle: T) -> Self {
        let line = Line::new(angle.sin(), -angle.cos(), T::zero());
        Self::from_line(&line)
    }

    /// Line part (the grade-2 reflection line).
    #[inline]
    pub fn line_part(&self) -> Line<T> {
        Line::new(self.normal_x(), self.normal_y(), self.dist())
    }

    /// Check if this is a pure reflection (no scalar component).
    #[inline]
    pub fn is_pure_reflection(&self, epsilon: T) -> bool {
        self.s().abs() < epsilon
    }
}

// ============================================================================
// Trivector extensions
// ============================================================================

impl<T: Float> Trivector<T> {
    /// The unit pseudoscalar e₀₁₂.
    #[inline]
    pub fn unit() -> Self {
        Self::new(T::one())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::{Join, Meet, Transform};
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;

    #[test]
    fn motor_identity_preserves_point() {
        let p = Point::<f64>::from_cartesian(3.0, 4.0);
        let m = Motor::<f64>::identity();

        let result = m.transform(&p);

        assert!(relative_eq!(
            result.x(),
            p.x(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.y(),
            p.y(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.w(),
            p.w(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_translation_x() {
        let origin = Point::<f64>::origin();
        let t = Motor::<f64>::from_translation(2.0, 0.0);

        let result = t.transform(&origin);

        assert!(relative_eq!(
            result.cartesian_x(),
            2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.cartesian_y(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_translation_y() {
        let origin = Point::<f64>::origin();
        let t = Motor::<f64>::from_translation(0.0, 3.0);

        let result = t.transform(&origin);

        assert!(relative_eq!(
            result.cartesian_x(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.cartesian_y(),
            3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_rotation_90_degrees() {
        use std::f64::consts::FRAC_PI_2;

        let p = Point::<f64>::from_cartesian(1.0, 0.0);
        let r = Motor::<f64>::from_rotation(FRAC_PI_2);

        let result = r.transform(&p);

        // 90 degree rotation of (1, 0) should give (0, 1)
        assert!(relative_eq!(
            result.cartesian_x(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.cartesian_y(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn point_join_creates_line() {
        let p1 = Point::<f64>::origin();
        let p2 = Point::<f64>::from_cartesian(1.0, 0.0);

        let line = p1.join(&p2);

        // Line through origin and (1, 0) is the x-axis (y = 0)
        // The wedge product gives normal_y component (e01) for this line
        assert!(line.normal_y().abs() > 0.1);
        // And the line should pass through origin (normal_x = 0, dist = 0)
        assert!(line.normal_x().abs() < 0.1);
    }

    #[test]
    fn line_meet_finds_intersection() {
        let l1 = Line::<f64>::x_axis(); // y = 0
        let l2 = Line::<f64>::y_axis(); // x = 0

        let intersection = l1.meet(&l2);

        // Should be the origin (finite point, not ideal)
        let (x, y) = intersection.to_cartesian().expect("should be finite");
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
    fn flector_reflect_x_axis() {
        let f = Flector::<f64>::reflect_x_axis();
        let p = Point::<f64>::from_cartesian(1.0, 2.0);

        eprintln!(
            "Flector: s={}, normal_x={}, normal_y={}, dist={}",
            f.s(),
            f.normal_x(),
            f.normal_y(),
            f.dist()
        );
        eprintln!("Point: x={}, y={}, w={}", p.x(), p.y(), p.w());

        let result = f.transform(&p);

        eprintln!(
            "Result: x={}, y={}, w={}",
            result.x(),
            result.y(),
            result.w()
        );
        eprintln!(
            "Cartesian: x={}, y={}",
            result.cartesian_x(),
            result.cartesian_y()
        );

        // Reflecting (1, 2) through x-axis should give (1, -2)
        assert!(relative_eq!(
            result.cartesian_x(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.cartesian_y(),
            -2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn flector_reflect_y_axis() {
        let f = Flector::<f64>::reflect_y_axis();
        let p = Point::<f64>::from_cartesian(1.0, 2.0);

        let result = f.transform(&p);

        // Reflecting (1, 2) through y-axis should give (-1, 2)
        assert!(relative_eq!(
            result.cartesian_x(),
            -1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.cartesian_y(),
            2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn point_distance() {
        let p1 = Point::<f64>::origin();
        let p2 = Point::<f64>::from_cartesian(3.0, 4.0);

        let dist = p1.distance(&p2);

        assert!(relative_eq!(
            dist,
            5.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn point_midpoint() {
        let p1 = Point::<f64>::from_cartesian(0.0, 0.0);
        let p2 = Point::<f64>::from_cartesian(4.0, 6.0);

        let mid = p1.midpoint(&p2);

        assert!(relative_eq!(
            mid.cartesian_x(),
            2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            mid.cartesian_y(),
            3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
