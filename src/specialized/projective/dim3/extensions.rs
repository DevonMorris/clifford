//! Domain-specific extensions for 3D Projective GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to 3D projective geometry.

use super::generated::products;
use super::generated::types::{Flector, Line, Motor, Plane, Point};
use crate::scalar::Float;
use crate::specialized::euclidean::dim3::Vector as EuclideanVector;

// ============================================================================
// Point extensions
// ============================================================================

impl<T: Float> Point<T> {
    /// Creates a finite point at Cartesian coordinates (x, y, z).
    ///
    /// The homogeneous weight `w` is set to 1.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Point;
    ///
    /// let p = Point::from_cartesian(3.0, 4.0, 5.0);
    /// assert_eq!(p.x(), 3.0);
    /// assert_eq!(p.y(), 4.0);
    /// assert_eq!(p.z(), 5.0);
    /// ```
    #[inline]
    pub fn from_cartesian(x: T, y: T, z: T) -> Self {
        Self::new(x, y, z, T::one())
    }

    /// Creates an ideal point (point at infinity) in the given direction.
    ///
    /// Ideal points have `w = 0` and represent directions rather than positions.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Point;
    ///
    /// let ideal = Point::<f64>::ideal(1.0, 0.0, 0.0);
    /// assert!(ideal.is_ideal(1e-10));
    /// ```
    #[inline]
    pub fn ideal(dx: T, dy: T, dz: T) -> Self {
        Self::new(dx, dy, dz, T::zero())
    }

    /// Origin point (0, 0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::from_cartesian(T::zero(), T::zero(), T::zero())
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

    /// Returns the z-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn z(&self) -> T {
        self.e3() / self.e0()
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
                self.e3() / self.e0(),
                T::one(),
            ))
        }
    }

    /// Returns the Cartesian coordinates as a tuple, if finite.
    ///
    /// Returns `None` if this is an ideal point.
    #[inline]
    pub fn to_cartesian(&self) -> Option<(T, T, T)> {
        if self.e0().abs() < T::epsilon() {
            None
        } else {
            Some((
                self.e1() / self.e0(),
                self.e2() / self.e0(),
                self.e3() / self.e0(),
            ))
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
    /// The bulk norm is the length of the spatial part: `e1² + e2² + e3²`.
    #[inline]
    pub fn bulk_norm_squared(&self) -> T {
        self.e1() * self.e1() + self.e2() * self.e2() + self.e3() * self.e3()
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

    /// Join with another point to create a line (regressive product).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Point;
    ///
    /// let p1 = Point::origin();
    /// let p2 = Point::from_cartesian(1.0, 0.0, 0.0);
    /// let line = p1.join(&p2);
    /// ```
    #[inline]
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::exterior_point_point(self, other)
    }

    /// Euclidean distance to another finite point.
    pub fn distance(&self, other: &Point<T>) -> T {
        self.distance_squared(other).sqrt()
    }

    /// Squared Euclidean distance.
    pub fn distance_squared(&self, other: &Point<T>) -> T {
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        let dz = self.z() - other.z();
        dx * dx + dy * dy + dz * dz
    }

    /// Midpoint between two finite points.
    pub fn midpoint(&self, other: &Point<T>) -> Point<T> {
        let two = T::TWO;
        Point::new(
            (self.e1() * other.e0() + other.e1() * self.e0()) / two,
            (self.e2() * other.e0() + other.e2() * self.e0()) / two,
            (self.e3() * other.e0() + other.e3() * self.e0()) / two,
            self.e0() * other.e0(),
        )
    }

    /// Inner product (dot product) with another point.
    #[inline]
    pub fn dot(&self, other: &Point<T>) -> T {
        self.e1() * other.e1() + self.e2() * other.e2() + self.e3() * other.e3()
    }

    /// Left contraction onto a line.
    #[inline]
    pub fn left_contract_line(&self, line: &Line<T>) -> Plane<T> {
        // P ⌋ L = P · L + P ∧ L (for grade-1 contracted with grade-2)
        // Result is grade 3 (plane)
        products::exterior_point_line(self, line)
    }
}

// ============================================================================
// Line extensions
// ============================================================================

impl<T: Float> Line<T> {
    /// Creates a line through two points (their join/wedge product).
    #[inline]
    pub fn join(p: &Point<T>, q: &Point<T>) -> Self {
        p.join(q)
    }

    /// Creates a line through a point in the given direction.
    pub fn from_point_and_direction(point: &Point<T>, direction: &EuclideanVector<T>) -> Self {
        let ideal = Point::ideal(direction.x(), direction.y(), direction.z());
        point.join(&ideal)
    }

    /// Creates a line from Plücker coordinates without validation.
    #[inline]
    pub fn from_plucker(direction: &EuclideanVector<T>, moment: &EuclideanVector<T>) -> Self {
        Self::new_unchecked(
            direction.x(),
            direction.y(),
            direction.z(),
            moment.x(),
            moment.y(),
            moment.z(),
        )
    }

    /// X-axis (line through origin along x).
    #[inline]
    pub fn x_axis() -> Self {
        Self::new_unchecked(
            T::one(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        )
    }

    /// Y-axis.
    #[inline]
    pub fn y_axis() -> Self {
        Self::new_unchecked(
            T::zero(),
            T::one(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        )
    }

    /// Z-axis.
    #[inline]
    pub fn z_axis() -> Self {
        Self::new_unchecked(
            T::zero(),
            T::zero(),
            T::one(),
            T::zero(),
            T::zero(),
            T::zero(),
        )
    }

    /// Direction vector (e01, e02, e03).
    #[inline]
    pub fn direction(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.e01(), self.e02(), self.e03())
    }

    /// Moment vector (e23, e31, e12).
    #[inline]
    pub fn moment(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.e23(), self.e31(), self.e12())
    }

    /// Weight norm (direction magnitude).
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.direction().norm()
    }

    /// Bulk norm (moment magnitude).
    #[inline]
    pub fn bulk_norm(&self) -> T {
        self.moment().norm()
    }

    /// Unitize to unit direction.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            self.e01() / wn,
            self.e02() / wn,
            self.e03() / wn,
            self.e23() / wn,
            self.e31() / wn,
            self.e12() / wn,
        )
    }

    /// Plücker condition residual: d · m.
    #[inline]
    pub fn plucker_residual(&self) -> T {
        self.e01() * self.e23() + self.e02() * self.e31() + self.e03() * self.e12()
    }

    /// Check if line satisfies Plücker condition.
    #[inline]
    pub fn satisfies_plucker_condition(&self, tolerance: T) -> bool {
        self.plucker_residual().abs() < tolerance
    }

    /// Check if line passes through origin.
    pub fn through_origin(&self) -> bool {
        self.moment().norm() < T::epsilon()
    }

    /// Check if two lines are parallel.
    pub fn is_parallel(&self, other: &Line<T>) -> bool {
        self.direction().cross(other.direction()).norm() < T::epsilon()
    }

    /// Inner product (dot product) with another line.
    #[inline]
    pub fn dot(&self, other: &Line<T>) -> T {
        self.e01() * other.e01()
            + self.e02() * other.e02()
            + self.e03() * other.e03()
            + self.e23() * other.e23()
            + self.e31() * other.e31()
            + self.e12() * other.e12()
    }

    /// Geometric norm.
    #[inline]
    pub fn geometric_norm(&self) -> T {
        let weight = self.weight_norm();
        if weight < T::epsilon() {
            T::zero()
        } else {
            self.bulk_norm() / weight
        }
    }

    /// Join with a point to create a plane.
    #[inline]
    pub fn join_point(&self, point: &Point<T>) -> Plane<T> {
        products::exterior_line_point(self, point)
    }

    /// Meet with a plane to find intersection point (regressive product).
    ///
    /// Returns the point where this line intersects the plane. If the line
    /// is parallel to the plane (no intersection), returns an ideal point.
    #[inline]
    pub fn meet_plane(&self, plane: &Plane<T>) -> Point<T> {
        products::regressive_line_plane(self, plane)
    }

    /// Angle between two lines (in radians).
    pub fn angle(&self, other: &Line<T>) -> T {
        let d1 = self.direction().normalized();
        let d2 = other.direction().normalized();
        let cos_angle = d1.dot(d2).abs().min(T::one());
        cos_angle.acos()
    }

    /// Shortest distance between two lines.
    pub fn distance(&self, other: &Line<T>) -> T {
        let d1 = self.direction();
        let d2 = other.direction();
        let m1 = self.moment();
        let m2 = other.moment();

        // Reciprocal product gives distance for normalized lines
        let cross = d1.cross(d2);
        let cross_norm = cross.norm();

        if cross_norm < T::epsilon() {
            // Parallel lines
            let l1_unit = self.unitized();
            let l2_unit = other.unitized();
            (l1_unit.moment() - l2_unit.moment()).norm()
        } else {
            // Skew lines
            let reciprocal = d1.dot(m2) + d2.dot(m1);
            reciprocal.abs() / cross_norm
        }
    }

    /// Returns true if the line is degenerate (zero direction).
    #[inline]
    pub fn is_zero(&self, epsilon: T) -> bool {
        self.weight_norm_squared() < epsilon * epsilon
    }

    /// Squared weight norm (direction magnitude squared).
    #[inline]
    pub fn weight_norm_squared(&self) -> T {
        self.e01() * self.e01() + self.e02() * self.e02() + self.e03() * self.e03()
    }

    /// Computes the Plücker inner product (used for testing intersection/parallelism).
    ///
    /// Two lines are:
    /// - Parallel or identical if the result is zero and they have parallel directions
    /// - Intersecting if the result is zero and they have non-parallel directions
    /// - Skew if the result is non-zero
    #[inline]
    pub fn plucker_inner(&self, other: &Line<T>) -> T {
        // direction1 · moment2 + direction2 · moment1
        self.e01() * other.e23()
            + self.e02() * other.e31()
            + self.e03() * other.e12()
            + other.e01() * self.e23()
            + other.e02() * self.e31()
            + other.e03() * self.e12()
    }

    /// Distance from a point to this line.
    ///
    /// The line should be unitized for accurate results.
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        let d1 = self.e01();
        let d2 = self.e02();
        let d3 = self.e03();
        let m1 = self.e23();
        let m2 = self.e31();
        let m3 = self.e12();

        let px = p.e1();
        let py = p.e2();
        let pz = p.e3();
        let pw = p.e0();

        // Cross product of direction and point position
        let cx = d2 * pz - d3 * py;
        let cy = d3 * px - d1 * pz;
        let cz = d1 * py - d2 * px;

        // Also include moment contribution
        let nx = cx + m1 * pw;
        let ny = cy + m2 * pw;
        let nz = cz + m3 * pw;

        // Direction magnitude (for normalization)
        let dir_norm = self.weight_norm();
        if dir_norm < T::epsilon() {
            return T::zero();
        }

        (nx * nx + ny * ny + nz * nz).sqrt() / (dir_norm * pw.abs())
    }

    /// Meet with a plane to find intersection point (alias for meet_plane).
    #[inline]
    pub fn meet(&self, plane: &Plane<T>) -> Point<T> {
        self.meet_plane(plane)
    }

    /// Closest point on this line to a given point.
    ///
    /// Projects the point onto the line and returns the closest point on the line.
    pub fn closest_point(&self, p: &Point<T>) -> Point<T> {
        let d = self.direction();
        let m = self.moment();
        let d_sq = d.x() * d.x() + d.y() * d.y() + d.z() * d.z();

        if d_sq < T::epsilon() {
            return Point::origin();
        }

        // A point on the line: P_line = d × m / |d|²
        let line_pt_x = (d.y() * m.z() - d.z() * m.y()) / d_sq;
        let line_pt_y = (d.z() * m.x() - d.x() * m.z()) / d_sq;
        let line_pt_z = (d.x() * m.y() - d.y() * m.x()) / d_sq;

        // Vector from line point to given point
        let px = p.x() - line_pt_x;
        let py = p.y() - line_pt_y;
        let pz = p.z() - line_pt_z;

        // Project onto direction
        let t = (px * d.x() + py * d.y() + pz * d.z()) / d_sq;

        Point::from_cartesian(
            line_pt_x + t * d.x(),
            line_pt_y + t * d.y(),
            line_pt_z + t * d.z(),
        )
    }
}

// ============================================================================
// Plane extensions
// ============================================================================

impl<T: Float> Plane<T> {
    /// Create plane from normal and distance.
    ///
    /// The plane equation is `n·x + d = 0`.
    pub fn from_normal_and_distance(nx: T, ny: T, nz: T, d: T) -> Self {
        Self::new(nx, ny, nz, d)
    }

    /// XY plane (z = 0).
    #[inline]
    pub fn xy() -> Self {
        Self::from_normal_and_distance(T::zero(), T::zero(), T::one(), T::zero())
    }

    /// XZ plane (y = 0).
    #[inline]
    pub fn xz() -> Self {
        Self::from_normal_and_distance(T::zero(), T::one(), T::zero(), T::zero())
    }

    /// YZ plane (x = 0).
    #[inline]
    pub fn yz() -> Self {
        Self::from_normal_and_distance(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Normal vector.
    #[inline]
    pub fn normal(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.e023(), self.e031(), self.e012())
    }

    /// Distance from origin (signed).
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        self.e123()
    }

    /// Weight norm (normal magnitude).
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.normal().norm()
    }

    /// Bulk norm.
    #[inline]
    pub fn bulk_norm(&self) -> T {
        self.e123().abs()
    }

    /// Unitize to unit normal.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new(
            self.e023() / wn,
            self.e031() / wn,
            self.e012() / wn,
            self.e123() / wn,
        )
    }

    /// Inner product with another plane.
    #[inline]
    pub fn dot(&self, other: &Plane<T>) -> T {
        self.e023() * other.e023() + self.e031() * other.e031() + self.e012() * other.e012()
    }

    /// Attitude (ideal line at infinity).
    #[inline]
    pub fn attitude(&self) -> T {
        self.e123()
    }

    /// Geometric norm.
    #[inline]
    pub fn geometric_norm(&self) -> T {
        let weight = self.weight_norm();
        if weight < T::epsilon() {
            T::zero()
        } else {
            self.bulk_norm() / weight
        }
    }

    /// Meet with another plane to find intersection line (regressive product).
    ///
    /// Returns the line where the two planes intersect. If the planes are
    /// parallel (no intersection), returns a line at infinity.
    #[inline]
    pub fn meet(&self, other: &Plane<T>) -> Line<T> {
        products::regressive_plane_plane(self, other)
    }

    /// Angle between two planes (in radians).
    pub fn angle(&self, other: &Plane<T>) -> T {
        let n1 = self.normal().normalized();
        let n2 = other.normal().normalized();
        let cos_angle = n1.dot(n2).abs().min(T::one());
        cos_angle.acos()
    }

    /// Angle to a line (in radians).
    pub fn angle_to_line(&self, line: &Line<T>) -> T {
        let n = self.normal().normalized();
        let d = line.direction().normalized();
        let sin_angle = n.dot(d).abs().min(T::one());
        sin_angle.asin()
    }

    /// Signed distance from a point to this plane.
    pub fn signed_distance(&self, point: &Point<T>) -> T {
        let p = self.unitized();
        (p.e023() * point.e1()
            + p.e031() * point.e2()
            + p.e012() * point.e3()
            + p.e123() * point.e0())
            / point.e0()
    }

    /// Project a point onto this plane.
    pub fn project_point(&self, point: &Point<T>) -> Point<T> {
        let dist = self.signed_distance(point);
        let n = self.normal().normalized();
        Point::new(
            point.e1() - dist * n.x() * point.e0(),
            point.e2() - dist * n.y() * point.e0(),
            point.e3() - dist * n.z() * point.e0(),
            point.e0(),
        )
    }

    /// Project a line onto this plane.
    pub fn project_line(&self, line: &Line<T>) -> Line<T> {
        let n = self.normal();
        let d = line.direction();

        // Project direction onto plane
        let dot = n.dot(d);
        let n_norm_sq = n.dot(n);
        let proj_d = if n_norm_sq < T::epsilon() {
            d
        } else {
            EuclideanVector::new(
                d.x() - dot * n.x() / n_norm_sq,
                d.y() - dot * n.y() / n_norm_sq,
                d.z() - dot * n.z() / n_norm_sq,
            )
        };

        // The projected line passes through the projection of any point on the original line
        // For simplicity, use the closest point to origin
        let m = line.moment();
        let d_norm_sq = d.dot(d);
        let closest_point = if d_norm_sq < T::epsilon() {
            Point::origin()
        } else {
            let cross = d.cross(m);
            Point::new(
                cross.x() / d_norm_sq,
                cross.y() / d_norm_sq,
                cross.z() / d_norm_sq,
                T::one(),
            )
        };

        let proj_point = self.project_point(&closest_point);
        Line::from_point_and_direction(&proj_point, &proj_d)
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Identity motor (leaves all elements unchanged).
    ///
    /// In PGA with the antisandwich product, the identity is the pseudoscalar 𝟙 = e₀₁₂₃.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Motor, Point};
    ///
    /// let m = Motor::<f64>::identity();
    /// let p = Point::from_cartesian(1.0, 2.0, 3.0);
    /// let p2 = m.transform_point(&p);
    /// assert!((p2.x() - p.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self::new_unchecked(
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::one(), // e0123 = 1
        )
    }

    /// Pure translation motor.
    ///
    /// Creates a motor that translates by the vector (dx, dy, dz).
    ///
    /// In PGA with our basis ordering, translation is encoded with specific sign conventions
    /// determined by the antisandwich formula.
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self {
        let half = T::one() / T::TWO;
        Self::new_unchecked(
            T::zero(),  // s
            dx * half,  // e23
            -dy * half, // e31 (negated due to basis ordering)
            dz * half,  // e12
            T::zero(),  // e01
            T::zero(),  // e02
            T::zero(),  // e03
            T::one(),   // e0123 (identity part)
        )
    }

    /// Pure rotation around x-axis through origin.
    ///
    /// In PGA, the rotation formula is R = l·sin(φ/2) + 𝟙·cos(φ/2)
    /// where l is the line (axis) and the antisandwich applies the rotation twice.
    pub fn from_rotation_x(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            T::zero(),  // s
            T::zero(),  // e23
            T::zero(),  // e31
            T::zero(),  // e12
            half.sin(), // e01 (x direction)
            T::zero(),  // e02
            T::zero(),  // e03
            half.cos(), // e0123
        )
    }

    /// Pure rotation around y-axis through origin.
    pub fn from_rotation_y(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            T::zero(),  // s
            T::zero(),  // e23
            T::zero(),  // e31
            T::zero(),  // e12
            T::zero(),  // e01
            half.sin(), // e02 (y direction)
            T::zero(),  // e03
            half.cos(), // e0123
        )
    }

    /// Pure rotation around z-axis through origin.
    pub fn from_rotation_z(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            T::zero(),  // s
            T::zero(),  // e23
            T::zero(),  // e31
            T::zero(),  // e12
            T::zero(),  // e01
            T::zero(),  // e02
            half.sin(), // e03 (z direction)
            half.cos(), // e0123
        )
    }

    /// Rotation around arbitrary axis through origin.
    ///
    /// The axis vector determines the rotation axis (will be normalized).
    /// The rotation follows the right-hand rule.
    pub fn from_axis_angle(axis: &EuclideanVector<T>, angle: T) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let axis_norm = axis.normalized();
        Self::new_unchecked(
            T::zero(),                // s
            T::zero(),                // e23
            T::zero(),                // e31
            T::zero(),                // e12
            sin_half * axis_norm.x(), // e01
            sin_half * axis_norm.y(), // e02
            sin_half * axis_norm.z(), // e03
            cos_half,                 // e0123
        )
    }

    /// Screw motion along a line.
    pub fn from_line(line: &Line<T>, angle: T, distance: T) -> Self {
        let line_unit = line.unitized();
        let half_angle = angle / T::TWO;
        let half_dist = distance / T::TWO;

        let (sin_a, cos_a) = (half_angle.sin(), half_angle.cos());

        let d = line_unit.direction();
        let m = line_unit.moment();

        Self::new_unchecked(
            cos_a,
            sin_a * d.x(),
            sin_a * d.y(),
            sin_a * d.z(),
            sin_a * m.x() + half_dist * cos_a * d.x(),
            sin_a * m.y() + half_dist * cos_a * d.y(),
            sin_a * m.z() + half_dist * cos_a * d.z(),
            -half_dist * sin_a,
        )
    }

    /// Compose motors: self then other.
    ///
    /// The result applies `self` first, then `other`.
    /// In PGA with the antisandwich transformation, motor composition uses
    /// the geometric antiproduct (∨) to properly combine transformations.
    pub fn compose(&self, other: &Motor<T>) -> Motor<T> {
        products::antigeometric_motor_motor(self, other)
    }

    /// Inverse motor.
    ///
    /// Uses the weight norm squared (rotation part: e01² + e02² + e03² + e0123²)
    /// for proper inversion. For a unit motor, this equals 1.
    pub fn inverse(&self) -> Self {
        use crate::norm::DegenerateNormed;
        let wn_sq = self.weight_norm_squared();
        let rev = self.reverse();
        Self::new_unchecked(
            rev.s() / wn_sq,
            rev.e23() / wn_sq,
            rev.e31() / wn_sq,
            rev.e12() / wn_sq,
            rev.e01() / wn_sq,
            rev.e02() / wn_sq,
            rev.e03() / wn_sq,
            rev.e0123() / wn_sq,
        )
    }

    /// Unitize to unit bulk norm (makes the rotor part have unit magnitude).
    ///
    /// For a motor to represent a proper rigid transformation, the bulk norm
    /// (rotor part: s² + e23² + e31² + e12²) should be 1.
    pub fn unitized(&self) -> Self {
        use crate::norm::DegenerateNormed;
        let bn = self.bulk_norm();
        if bn < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            self.s() / bn,
            self.e23() / bn,
            self.e31() / bn,
            self.e12() / bn,
            self.e01() / bn,
            self.e02() / bn,
            self.e03() / bn,
            self.e0123() / bn,
        )
    }

    /// Check if motor is unitized (bulk norm ≈ 1).
    ///
    /// A unitized motor has bulk_norm_squared ≈ 1, meaning the rotor part
    /// has unit magnitude.
    #[inline]
    pub fn is_unitized(&self, tolerance: T) -> bool {
        use crate::norm::DegenerateNormed;
        (self.bulk_norm_squared() - T::one()).abs() < tolerance
    }

    /// Geometric constraint residual.
    ///
    /// The geometric constraint requires: `s·e₀₁₂₃ + e₂₃·e₀₁ + e₃₁·e₀₂ + e₁₂·e₀₃ = 0`
    ///
    /// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint>
    #[inline]
    pub fn geometric_constraint_residual(&self) -> T {
        self.s() * self.e0123()
            + self.e23() * self.e01()
            + self.e31() * self.e02()
            + self.e12() * self.e03()
    }

    /// Check if motor satisfies the geometric constraint.
    ///
    /// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint>
    #[inline]
    pub fn satisfies_geometric_constraint(&self, tolerance: T) -> bool {
        self.geometric_constraint_residual().abs() < tolerance
    }

    /// Extract rotation angle.
    pub fn rotation_angle(&self) -> T {
        let cos_half = self.s().min(T::one()).max(-T::one());
        cos_half.acos() * T::TWO
    }

    /// Extract translation vector.
    ///
    /// For a pure translation motor, this returns the translation vector.
    /// For combined rotation-translation motors, this extracts the translational
    /// component based on the dual motor representation.
    ///
    /// **Note**: This is only exact for pure translation motors. For general
    /// motors (rotation + translation), use decomposition methods for accurate
    /// extraction.
    pub fn translation(&self) -> EuclideanVector<T> {
        // Inverse of from_translation encoding:
        // from_translation sets: e23 = dx/2, e31 = -dy/2, e12 = dz/2
        // So: dx = 2*e23, dy = -2*e31, dz = 2*e12
        EuclideanVector::new(
            T::TWO * self.e23(),
            -T::TWO * self.e31(),
            T::TWO * self.e12(),
        )
    }

    /// Commutator: [A, B] = AB - BA.
    #[inline]
    pub fn commutator(&self, other: &Motor<T>) -> Motor<T> {
        let ab = products::geometric_motor_motor(self, other);
        let ba = products::geometric_motor_motor(other, self);
        Motor::new_unchecked(
            ab.s() - ba.s(),
            ab.e23() - ba.e23(),
            ab.e31() - ba.e31(),
            ab.e12() - ba.e12(),
            ab.e01() - ba.e01(),
            ab.e02() - ba.e02(),
            ab.e03() - ba.e03(),
            ab.e0123() - ba.e0123(),
        )
    }

    /// Anticommutator: {A, B} = AB + BA.
    #[inline]
    pub fn anticommutator(&self, other: &Motor<T>) -> Motor<T> {
        let ab = products::geometric_motor_motor(self, other);
        let ba = products::geometric_motor_motor(other, self);
        Motor::new_unchecked(
            ab.s() + ba.s(),
            ab.e23() + ba.e23(),
            ab.e31() + ba.e31(),
            ab.e12() + ba.e12(),
            ab.e01() + ba.e01(),
            ab.e02() + ba.e02(),
            ab.e03() + ba.e03(),
            ab.e0123() + ba.e0123(),
        )
    }

    /// Transform a point using the antisandwich product.
    ///
    /// In PGA, transformations use the geometric antiproduct: P' = M ⊛ P ⊛ M̃
    /// where ⊛ is the geometric antiproduct and M̃ is the antireverse.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Motor, Point};
    ///
    /// // Translate a point
    /// let t = Motor::<f64>::from_translation(1.0, 0.0, 0.0);
    /// let p = Point::origin();
    /// let p2 = t.transform_point(&p);
    /// assert!((p2.x() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::antisandwich_motor_point(self, p)
    }

    /// Transform a line using the antisandwich product.
    ///
    /// In PGA, transformations use the geometric antiproduct: L' = M ⊛ L ⊛ M̃
    /// where ⊛ is the geometric antiproduct and M̃ is the antireverse.
    #[inline]
    pub fn transform_line(&self, line: &Line<T>) -> Line<T> {
        products::antisandwich_motor_line(self, line)
    }

    /// Transform a plane using the antisandwich product.
    ///
    /// In PGA, transformations use the geometric antiproduct: Π' = M ⊛ Π ⊛ M̃
    /// where ⊛ is the geometric antiproduct and M̃ is the antireverse.
    #[inline]
    pub fn transform_plane(&self, plane: &Plane<T>) -> Plane<T> {
        products::antisandwich_motor_plane(self, plane)
    }
}

// ============================================================================
// Flector extensions
// ============================================================================

impl<T: Float> Flector<T> {
    /// Create flector from reflection plane.
    pub fn from_plane(plane: &Plane<T>) -> Self {
        let p = plane.unitized();
        Self::new_unchecked(
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            p.e023(),
            p.e031(),
            p.e012(),
            p.e123(),
        )
    }

    /// Create reflection through plane at origin.
    pub fn from_plane_through_origin(nx: T, ny: T, nz: T) -> Self {
        let norm = (nx * nx + ny * ny + nz * nz).sqrt();
        Self::new_unchecked(
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            nx / norm,
            ny / norm,
            nz / norm,
            T::zero(),
        )
    }

    /// Reflect through XY plane.
    #[inline]
    pub fn reflect_xy() -> Self {
        Self::from_plane(&Plane::xy())
    }

    /// Reflect through XZ plane.
    #[inline]
    pub fn reflect_xz() -> Self {
        Self::from_plane(&Plane::xz())
    }

    /// Reflect through YZ plane.
    #[inline]
    pub fn reflect_yz() -> Self {
        Self::from_plane(&Plane::yz())
    }

    /// Point part (grade 1).
    #[inline]
    pub fn point_part(&self) -> Point<T> {
        Point::new(self.e1(), self.e2(), self.e3(), self.e0())
    }

    /// Plane part (grade 3).
    #[inline]
    pub fn plane_part(&self) -> Plane<T> {
        Plane::new(self.e023(), self.e031(), self.e012(), self.e123())
    }

    /// Check if this is a pure reflection (no point part).
    #[inline]
    pub fn is_pure_reflection(&self) -> bool {
        let pt = self.point_part();
        pt.bulk_norm() < T::epsilon() && pt.weight_norm() < T::epsilon()
    }

    /// Unitize to unit bulk norm.
    ///
    /// For a flector to represent a proper rigid reflection, the bulk norm
    /// (e1² + e2² + e3² + e123²) should be 1.
    pub fn unitized(&self) -> Self {
        use crate::norm::DegenerateNormed;
        let bn = self.bulk_norm();
        if bn < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            self.e1() / bn,
            self.e2() / bn,
            self.e3() / bn,
            self.e0() / bn,
            self.e023() / bn,
            self.e031() / bn,
            self.e012() / bn,
            self.e123() / bn,
        )
    }

    /// Compose two flectors (result is a motor).
    ///
    /// Since flector transformations use the antisandwich product, composition
    /// uses the geometric antiproduct to properly combine transformations.
    /// Two reflections compose to give a rotation around their intersection axis.
    #[inline]
    pub fn compose(&self, other: &Flector<T>) -> Motor<T> {
        products::antigeometric_flector_flector(self, other)
    }

    /// Transform a point using the antisandwich product.
    ///
    /// In PGA with our convention, flector transformations use the geometric antiproduct:
    /// `P' = F ⊛ P ⊛ F̃` where ⊛ is the geometric antiproduct and F̃ is the antireverse.
    ///
    /// For a pure reflection through a plane, this reflects the point across the plane.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Flector, Point};
    ///
    /// // Reflect through XY plane (z = 0)
    /// let f = Flector::<f64>::reflect_xy();
    /// let p = Point::from_cartesian(1.0, 2.0, 3.0);
    /// let p2 = f.transform_point(&p);
    /// // Point should be reflected in z
    /// assert!((p2.x() - 1.0).abs() < 1e-10);
    /// assert!((p2.y() - 2.0).abs() < 1e-10);
    /// assert!((p2.z() + 3.0).abs() < 1e-10); // z is negated
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::antisandwich_flector_point(self, p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;

    #[test]
    fn motor_identity_preserves_point() {
        let p = Point::<f64>::from_cartesian(3.0, 4.0, 5.0);
        let m = Motor::<f64>::identity();

        eprintln!(
            "Identity motor: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            m.s(),
            m.e23(),
            m.e31(),
            m.e12(),
            m.e01(),
            m.e02(),
            m.e03(),
            m.e0123()
        );
        eprintln!(
            "Input point: ({}, {}, {}, {})",
            p.e1(),
            p.e2(),
            p.e3(),
            p.e0()
        );

        let result = m.transform_point(&p);

        eprintln!(
            "Output point: ({}, {}, {}, {})",
            result.e1(),
            result.e2(),
            result.e3(),
            result.e0()
        );

        assert!(relative_eq!(
            result.e1(),
            p.e1(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e2(),
            p.e2(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e3(),
            p.e3(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e0(),
            p.e0(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_translation_x() {
        let origin = Point::<f64>::origin();
        let t = Motor::<f64>::from_translation(2.0, 0.0, 0.0);

        eprintln!(
            "Motor T: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            t.s(),
            t.e23(),
            t.e31(),
            t.e12(),
            t.e01(),
            t.e02(),
            t.e03(),
            t.e0123()
        );
        eprintln!(
            "Origin: ({}, {}, {}, {})",
            origin.e1(),
            origin.e2(),
            origin.e3(),
            origin.e0()
        );

        let result = t.transform_point(&origin);

        eprintln!(
            "Translated: ({}, {}, {}, {})",
            result.e1(),
            result.e2(),
            result.e3(),
            result.e0()
        );

        // Expected: (2, 0, 0) with w=1
        assert!(relative_eq!(
            result.e1(),
            2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e2(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e3(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e0(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_translation_y() {
        let origin = Point::<f64>::origin();
        let t = Motor::<f64>::from_translation(0.0, -19.26, 0.0);

        eprintln!(
            "Motor T: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            t.s(),
            t.e23(),
            t.e31(),
            t.e12(),
            t.e01(),
            t.e02(),
            t.e03(),
            t.e0123()
        );
        eprintln!(
            "Origin: ({}, {}, {}, {})",
            origin.e1(),
            origin.e2(),
            origin.e3(),
            origin.e0()
        );

        let result = t.transform_point(&origin);

        eprintln!(
            "Translated: ({}, {}, {}, {})",
            result.e1(),
            result.e2(),
            result.e3(),
            result.e0()
        );

        // Expected: (0, -19.26, 0) with w=1
        assert!(relative_eq!(
            result.e1(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e2(),
            -19.26,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e3(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_translation_z() {
        let origin = Point::<f64>::origin();
        let t = Motor::<f64>::from_translation(0.0, 0.0, 5.0);

        eprintln!(
            "Motor T: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            t.s(),
            t.e23(),
            t.e31(),
            t.e12(),
            t.e01(),
            t.e02(),
            t.e03(),
            t.e0123()
        );

        let result = t.transform_point(&origin);

        eprintln!(
            "Translated: ({}, {}, {}, {})",
            result.e1(),
            result.e2(),
            result.e3(),
            result.e0()
        );

        // Expected: (0, 0, 5) with w=1
        assert!(relative_eq!(
            result.e1(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e2(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.e3(),
            5.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn flector_reflect_xy_plane() {
        let f = Flector::<f64>::reflect_xy();
        let p = Point::<f64>::from_cartesian(1.0, 2.0, 3.0);

        eprintln!(
            "Flector: e1={}, e2={}, e3={}, e0={}, e023={}, e031={}, e012={}, e123={}",
            f.e1(),
            f.e2(),
            f.e3(),
            f.e0(),
            f.e023(),
            f.e031(),
            f.e012(),
            f.e123()
        );
        eprintln!(
            "Input point: ({}, {}, {}, {})",
            p.e1(),
            p.e2(),
            p.e3(),
            p.e0()
        );

        let result = f.transform_point(&p);

        eprintln!(
            "Output point: ({}, {}, {}, {})",
            result.e1(),
            result.e2(),
            result.e3(),
            result.e0()
        );

        // Reflecting (1, 2, 3) through XY plane should give (1, 2, -3)
        assert!(relative_eq!(
            result.x(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.y(),
            2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.z(),
            -3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn flector_reflect_yz_plane() {
        let f = Flector::<f64>::reflect_yz();
        let p = Point::<f64>::from_cartesian(1.0, 2.0, 3.0);

        let result = f.transform_point(&p);

        // Reflecting (1, 2, 3) through YZ plane should give (-1, 2, 3)
        assert!(relative_eq!(
            result.x(),
            -1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.y(),
            2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.z(),
            3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn flector_reflect_xz_plane() {
        let f = Flector::<f64>::reflect_xz();
        let p = Point::<f64>::from_cartesian(1.0, 2.0, 3.0);

        let result = f.transform_point(&p);

        // Reflecting (1, 2, 3) through XZ plane should give (1, -2, 3)
        assert!(relative_eq!(
            result.x(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.y(),
            -2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.z(),
            3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn flector_compose_two_reflections_gives_rotation() {
        // Composing two reflections through planes that meet at an angle gives a rotation
        // by twice that angle around their intersection axis

        let f1 = Flector::<f64>::reflect_xz(); // YZ normal, reflects y
        let f2 = Flector::<f64>::reflect_yz(); // XZ normal, reflects x

        // F1 * F2 should give a motor (rotation)
        let m = f1.compose(&f2);

        // Apply the motor to a point
        let p = Point::<f64>::from_cartesian(1.0, 0.0, 0.0);
        let result = m.transform_point(&p);

        eprintln!(
            "Composed motor: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            m.s(),
            m.e23(),
            m.e31(),
            m.e12(),
            m.e01(),
            m.e02(),
            m.e03(),
            m.e0123()
        );
        eprintln!(
            "Transformed point: ({}, {}, {})",
            result.x(),
            result.y(),
            result.z()
        );

        // This should be a 180-degree rotation around the z-axis
        // (1, 0, 0) -> (-1, 0, 0)
        assert!(relative_eq!(
            result.x(),
            -1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.y(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.z(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
