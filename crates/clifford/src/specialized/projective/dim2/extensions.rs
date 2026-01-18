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
    /// Creates a line through a point in the given direction.
    pub fn from_point_and_direction(point: &Point<T>, direction: &EuclideanVector<T>) -> Self {
        // Ideal point (point at infinity) in the given direction
        let ideal = Point::new(direction.x(), direction.y(), T::zero());
        crate::ops::Wedge::wedge(point, &ideal)
    }

    /// Creates a line from the standard equation ax + by + c = 0.
    ///
    /// The coefficients map directly to internal representation:
    /// - nx = a (x-coefficient)
    /// - ny = b (y-coefficient)
    /// - d = c (constant term)
    #[inline]
    pub fn from_equation(a: T, b: T, c: T) -> Self {
        // new() takes (d, nx, ny) where the equation is: nx*x + ny*y + d = 0
        Self::new(c, a, b)
    }

    /// X-axis (line y = 0).
    ///
    /// The x-axis has equation 0*x + 1*y + 0 = 0.
    #[inline]
    pub fn x_axis() -> Self {
        Self::new(T::zero(), T::zero(), T::one())
    }

    /// Y-axis (line x = 0).
    ///
    /// The y-axis has equation 1*x + 0*y + 0 = 0.
    #[inline]
    pub fn y_axis() -> Self {
        Self::new(T::zero(), T::one(), T::zero())
    }

    /// Normal vector (a, b) where the line equation is ax + by + c = 0.
    ///
    /// Returns (nx, ny) directly.
    #[inline]
    pub fn normal(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.nx(), self.ny())
    }

    /// Distance from origin (d component, the constant term).
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        self.d()
    }

    /// Check if line passes through origin.
    pub fn through_origin(&self, epsilon: T) -> bool {
        self.d().abs() < epsilon
    }

    /// Check if two lines are parallel.
    pub fn is_parallel(&self, other: &Line<T>, epsilon: T) -> bool {
        // Lines are parallel if their normals are parallel.
        // Cross product of normals: nx1*ny2 - ny1*nx2
        let cross = self.nx() * other.ny() - self.ny() * other.nx();
        cross.abs() < epsilon
    }

    /// Angle between two lines (in radians).
    pub fn angle(&self, other: &Line<T>) -> T {
        use crate::norm::DegenerateNormed;
        let n1 = self.try_unitize().unwrap_or(*self);
        let n2 = other.try_unitize().unwrap_or(*other);
        // Dot product of normals (nx, ny)
        let cos_angle = (n1.nx() * n2.nx() + n1.ny() * n2.ny()).abs().min(T::one());
        cos_angle.acos()
    }

    /// Signed distance from a point to this line.
    ///
    /// The line should be unitized for accurate results.
    /// Uses the standard formula: (a*x + b*y + c) / sqrt(a² + b²)
    pub fn signed_distance(&self, point: &Point<T>) -> T {
        use crate::norm::DegenerateNormed;
        let l = self.try_unitize().unwrap_or(*self);
        // Line equation: nx*x + ny*y + d = 0
        // Distance = (nx*px + ny*py + d*w) / w
        (l.nx() * point.x() + l.ny() * point.y() + l.d() * point.w()) / point.w()
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
            dy * half,  // ty (controls y-translation)
            -dx * half, // tx (controls x-translation)
            T::zero(),  // r (rotation)
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
            T::zero(),   // ty
            T::zero(),   // tx
            -half.sin(), // r (rotation = -sin(θ/2))
            half.cos(),  // ps (= cos(θ/2))
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
            rev.ty() / norm_sq,
            rev.tx() / norm_sq,
            rev.r() / norm_sq,
            rev.ps() / norm_sq,
        )
    }

    /// Extract rotation angle.
    pub fn rotation_angle(&self) -> T {
        let cos_half = self.ps().min(T::one()).max(-T::one());
        let sin_half = self.r();
        sin_half.atan2(cos_half) * T::TWO
    }

    /// Extract translation vector.
    ///
    /// For a pure translation motor, this returns the translation vector.
    pub fn translation(&self) -> EuclideanVector<T> {
        // Inverse of from_translation encoding:
        // from_translation sets: ty = dy/2, tx = -dx/2
        // So: dx = -2*tx, dy = 2*ty
        EuclideanVector::new(-T::TWO * self.tx(), T::TWO * self.ty())
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
        // new() takes (s, d, nx, ny)
        Self::new(T::zero(), l.d(), l.nx(), l.ny())
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
        // Line through origin at angle θ from x-axis.
        // Normal is perpendicular to line direction, so normal = (sin(θ), -cos(θ))
        // For line equation a*x + b*y = 0: a = sin(θ), b = -cos(θ)
        let line = Line::new(T::zero(), angle.sin(), -angle.cos());
        Self::from_line(&line)
    }

    /// Line part (the grade-2 reflection line).
    #[inline]
    pub fn line_part(&self) -> Line<T> {
        // Line::new() takes (d, nx, ny)
        Line::new(self.d(), self.nx(), self.ny())
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
        // Standard equation: 0*x + 1*y + 0 = 0
        // With our convention: nx=0, ny=1, d=0
        assert!(relative_eq!(
            line.nx(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            line.ny(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        // And the line should pass through origin (d = 0)
        assert!(relative_eq!(
            line.d(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
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

        let result = f.transform(&p);

        // Reflecting (1, 2) through x-axis should give (1, -2)
        // Use Cartesian coordinates since antisandwich may return different
        // projective representatives (e.g., (-1, 2, -1) instead of (1, -2, 1))
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
        // Use Cartesian coordinates since antisandwich may return different
        // projective representatives (e.g., (1, -2, -1) instead of (-1, 2, 1))
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
