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
        // P ∧ Q for P = (px, py, pz, pw) and Q = (qx, qy, qz, qw)
        // Note: The generated exterior_point_point has incorrect sign/ordering mapping,
        // so we use the explicit formula from the PGA literature.
        Line::new_unchecked(
            self.e0() * other.e1() - self.e1() * other.e0(), // e01
            self.e0() * other.e2() - self.e2() * other.e0(), // e02
            self.e0() * other.e3() - self.e3() * other.e0(), // e03
            self.e2() * other.e3() - self.e3() * other.e2(), // e23
            self.e3() * other.e1() - self.e1() * other.e3(), // e31
            self.e1() * other.e2() - self.e2() * other.e1(), // e12
        )
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

    /// Left contraction onto a plane (returns scalar).
    #[inline]
    pub fn left_contract_plane(&self, plane: &Plane<T>) -> T {
        // P ⌋ Π for grade-1 and grade-3
        self.e1() * plane.e023()
            + self.e2() * plane.e031()
            + self.e3() * plane.e012()
            + self.e0() * plane.e123()
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
        // Regressive product L ∨ P computed via coordinate formula
        // Line: vx·e23 + vy·e31 + vz·e12 + mx·e01 + my·e02 + mz·e03
        // Plane: nx·e023 + ny·e031 + nz·e012 + d·e123
        let vx = self.e23();
        let vy = self.e31();
        let vz = self.e12();
        let mx = self.e01();
        let my = self.e02();
        let mz = self.e03();

        let nx = plane.e023();
        let ny = plane.e031();
        let nz = plane.e012();
        let d = plane.e123();

        // Point = (m × n + v·d, v · n)
        let e1 = my * nz - mz * ny + vx * d;
        let e2 = mz * nx - mx * nz + vy * d;
        let e3 = mx * ny - my * nx + vz * d;
        let e0 = vx * nx + vy * ny + vz * nz;

        Point::new(e1, e2, e3, e0)
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

    /// Left contraction onto a plane (returns point).
    ///
    /// Computes L ⌋ Π where L is a line (grade-2) and Π is a plane (grade-3).
    /// The result is a point (grade-1), which is the grade-1 part of the
    /// geometric product Line * Plane.
    #[inline]
    pub fn left_contract_plane(&self, plane: &Plane<T>) -> Point<T> {
        // Extract from geometric product Line * Plane -> Flector
        // The point (grade-1) part is the left contraction
        let e1 = -(plane.e023() * self.e03());
        let e2 = plane.e023() * self.e02();
        let e3 = -(plane.e023() * self.e01());
        let e0 = -(plane.e031() * self.e01())
            - (plane.e012() * self.e02())
            - (plane.e123() * self.e03());
        Point::new(e1, e2, e3, e0)
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
        // Regressive product P1 ∨ P2 computed via coordinate formula
        // Plane1: n1x·e023 + n1y·e031 + n1z·e012 + d1·e123
        // Plane2: n2x·e023 + n2y·e031 + n2z·e012 + d2·e123
        let n1x = self.e023();
        let n1y = self.e031();
        let n1z = self.e012();
        let d1 = self.e123();

        let n2x = other.e023();
        let n2y = other.e031();
        let n2z = other.e012();
        let d2 = other.e123();

        // Direction = n1 × n2 (stored in e01, e02, e03 for Line)
        let dir_x = n1y * n2z - n1z * n2y;
        let dir_y = n1z * n2x - n1x * n2z;
        let dir_z = n1x * n2y - n1y * n2x;

        // Moment from plane offsets (stored in e23, e31, e12 for Line)
        let mom_x = n1x * d2 - n2x * d1;
        let mom_y = n1y * d2 - n2y * d1;
        let mom_z = n1z * d2 - n2z * d1;

        // Line::new_unchecked(e01, e02, e03, e23, e31, e12)
        Line::new_unchecked(dir_x, dir_y, dir_z, mom_x, mom_y, mom_z)
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
    /// Pure translation motor.
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self {
        let half = T::one() / T::TWO;
        Self::new_unchecked(
            T::one(),
            T::zero(),
            T::zero(),
            T::zero(),
            dx * half,
            dy * half,
            dz * half,
            T::zero(),
        )
    }

    /// Pure rotation around x-axis through origin.
    pub fn from_rotation_x(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            half.cos(),
            half.sin(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        )
    }

    /// Pure rotation around y-axis through origin.
    pub fn from_rotation_y(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            half.cos(),
            T::zero(),
            half.sin(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        )
    }

    /// Pure rotation around z-axis through origin.
    pub fn from_rotation_z(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new_unchecked(
            half.cos(),
            T::zero(),
            T::zero(),
            half.sin(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        )
    }

    /// Rotation around arbitrary axis through origin.
    pub fn from_axis_angle(axis: &EuclideanVector<T>, angle: T) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let axis_norm = axis.normalized();
        Self::new_unchecked(
            cos_half,
            sin_half * axis_norm.x(),
            sin_half * axis_norm.y(),
            sin_half * axis_norm.z(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
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
    /// In PGA, the sandwich product `(A * B) * p * (A * B)̃` applies A first, then B.
    /// So for `compose(self, other)` to apply self first, we compute `self * other`.
    ///
    /// This uses an explicit formula derived from sympy to ensure correct composition.
    pub fn compose(&self, other: &Motor<T>) -> Motor<T> {
        // M1 = other, M2 = self → result = self * other
        let s1 = other.s();
        let b23_1 = other.e23();
        let b31_1 = other.e31();
        let b12_1 = other.e12();
        let d01_1 = other.e01();
        let d02_1 = other.e02();
        let d03_1 = other.e03();
        let i1 = other.e0123();

        let s2 = self.s();
        let b23_2 = self.e23();
        let b31_2 = self.e31();
        let b12_2 = self.e12();
        let d01_2 = self.e01();
        let d02_2 = self.e02();
        let d03_2 = self.e03();
        let i2 = self.e0123();

        // Geometric product: self * other = M2 * M1
        let s = -b12_1 * b12_2 - b23_1 * b23_2 - b31_1 * b31_2 + s1 * s2;
        let e23 = -b12_1 * b31_2 + b12_2 * b31_1 + b23_1 * s2 + b23_2 * s1;
        let e31 = b12_1 * b23_2 - b12_2 * b23_1 + b31_1 * s2 + b31_2 * s1;
        let e12 = b12_1 * s2 + b12_2 * s1 + b23_1 * b31_2 - b23_2 * b31_1;
        let e01 = -i1 * b23_2 - i2 * b23_1 + d01_1 * s2 + d01_2 * s1 + d02_1 * b12_2
            - d02_2 * b12_1
            - d03_1 * b31_2
            + d03_2 * b31_1;
        let e02 = -i1 * b31_2 - i2 * b31_1 - d01_1 * b12_2
            + d01_2 * b12_1
            + d02_1 * s2
            + d02_2 * s1
            + d03_1 * b23_2
            - d03_2 * b23_1;
        let e03 = -i1 * b12_2 - i2 * b12_1 + d01_1 * b31_2 - d01_2 * b31_1 - d02_1 * b23_2
            + d02_2 * b23_1
            + d03_1 * s2
            + d03_2 * s1;
        let e0123 = i1 * s2
            + i2 * s1
            + d01_1 * b23_2
            + d01_2 * b23_1
            + d02_1 * b31_2
            + d02_2 * b31_1
            + d03_1 * b12_2
            + d03_2 * b12_1;

        Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
    }

    /// Inverse motor.
    pub fn inverse(&self) -> Self {
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

    /// Weight norm squared (rotation part).
    #[inline]
    pub fn weight_norm_squared(&self) -> T {
        self.s() * self.s()
            + self.e23() * self.e23()
            + self.e31() * self.e31()
            + self.e12() * self.e12()
    }

    /// Weight norm.
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    /// Unitize to unit weight norm.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            self.s() / wn,
            self.e23() / wn,
            self.e31() / wn,
            self.e12() / wn,
            self.e01() / wn,
            self.e02() / wn,
            self.e03() / wn,
            self.e0123() / wn,
        )
    }

    /// Check if motor is unitized.
    #[inline]
    pub fn is_unitized(&self, tolerance: T) -> bool {
        (self.weight_norm_squared() - T::one()).abs() < tolerance
    }

    /// Study condition residual.
    #[inline]
    pub fn study_residual(&self) -> T {
        self.s() * self.e0123()
            + self.e23() * self.e01()
            + self.e31() * self.e02()
            + self.e12() * self.e03()
    }

    /// Check if motor satisfies Study condition.
    #[inline]
    pub fn satisfies_study_condition(&self, tolerance: T) -> bool {
        self.study_residual().abs() < tolerance
    }

    /// Extract rotation angle.
    pub fn rotation_angle(&self) -> T {
        let cos_half = self.s().min(T::one()).max(-T::one());
        cos_half.acos() * T::TWO
    }

    /// Extract translation vector.
    pub fn translation(&self) -> EuclideanVector<T> {
        let two = T::TWO;
        EuclideanVector::new(
            two * (self.s() * self.e01() + self.e31() * self.e03() - self.e12() * self.e02()
                + self.e23() * self.e0123()),
            two * (self.s() * self.e02() + self.e12() * self.e01() - self.e23() * self.e03()
                + self.e31() * self.e0123()),
            two * (self.s() * self.e03() + self.e23() * self.e02() - self.e31() * self.e01()
                + self.e12() * self.e0123()),
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

    /// Transform a point by this motor.
    ///
    /// Uses the PGA-specific transformation formula from Rigid Geometric Algebra.
    /// The standard antisandwich formula fails due to e0² = 0, so we use the
    /// explicit coordinate formula:
    ///
    /// ```text
    /// a = v × p_xyz + pw * m
    /// p'_xyz = p_xyz + 2(s * a + v × a - e0123 * pw * v)
    /// p'_w = pw
    /// ```
    ///
    /// where v = (e23, e31, e12), m = (e01, e02, e03).
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        let two = T::TWO;

        // Motor components
        let s = self.s();
        let vx = self.e23();
        let vy = self.e31();
        let vz = self.e12();
        let mx = self.e01();
        let my = self.e02();
        let mz = self.e03();
        let mw = self.e0123();

        // Point components
        let px = p.e1();
        let py = p.e2();
        let pz = p.e3();
        let pw = p.e0();

        // a = v × p + pw * m
        let ax = vy * pz - vz * py + pw * mx;
        let ay = vz * px - vx * pz + pw * my;
        let az = vx * py - vy * px + pw * mz;

        // v × a
        let vxa_x = vy * az - vz * ay;
        let vxa_y = vz * ax - vx * az;
        let vxa_z = vx * ay - vy * ax;

        // p'_xyz = p_xyz + 2(s * a + v × a - mw * pw * v)
        let px_new = px + two * (s * ax + vxa_x - mw * pw * vx);
        let py_new = py + two * (s * ay + vxa_y - mw * pw * vy);
        let pz_new = pz + two * (s * az + vxa_z - mw * pw * vz);

        // p'_w = pw (unchanged)
        Point::new(px_new, py_new, pz_new, pw)
    }

    /// Transform a line via antisandwich product.
    #[inline]
    pub fn transform_line(&self, line: &Line<T>) -> Line<T> {
        products::antisandwich_motor_line(self, line)
    }

    /// Transform a plane via antisandwich product.
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

    /// Weight norm squared.
    #[inline]
    pub fn weight_norm_squared(&self) -> T {
        self.e1() * self.e1()
            + self.e2() * self.e2()
            + self.e3() * self.e3()
            + self.e023() * self.e023()
            + self.e031() * self.e031()
            + self.e012() * self.e012()
    }

    /// Weight norm.
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.weight_norm_squared().sqrt()
    }

    /// Unitize to unit weight norm.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            self.e1() / wn,
            self.e2() / wn,
            self.e3() / wn,
            self.e0() / wn,
            self.e023() / wn,
            self.e031() / wn,
            self.e012() / wn,
            self.e123() / wn,
        )
    }

    /// Compose two flectors (result is a motor).
    #[inline]
    pub fn compose(&self, other: &Flector<T>) -> Motor<T> {
        products::geometric_flector_flector(self, other)
    }

    /// Transform a point via antisandwich product.
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::antisandwich_flector_point(self, p)
    }

    /// Transform a line via antisandwich product.
    #[inline]
    pub fn transform_line(&self, line: &Line<T>) -> Line<T> {
        products::antisandwich_flector_line(self, line)
    }

    /// Transform a plane via antisandwich product.
    #[inline]
    pub fn transform_plane(&self, plane: &Plane<T>) -> Plane<T> {
        products::antisandwich_flector_plane(self, plane)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::abs_diff_eq;

    const EPS: f64 = 1e-10;

    #[test]
    fn antisandwich_translation_works() {
        // Create a pure translation motor
        let motor = Motor::from_translation(10.0, 20.0, 30.0);
        let p = Point::from_cartesian(1.0, 2.0, 3.0);

        // Debug: print motor components
        eprintln!(
            "Motor: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            motor.s(),
            motor.e23(),
            motor.e31(),
            motor.e12(),
            motor.e01(),
            motor.e02(),
            motor.e03(),
            motor.e0123()
        );
        eprintln!(
            "Point: e1={}, e2={}, e3={}, e0={}",
            p.e1(),
            p.e2(),
            p.e3(),
            p.e0()
        );

        // Use generated antisandwich transformation
        let result = motor.transform_point(&p);

        // Debug: print result components
        eprintln!(
            "Result: e1={}, e2={}, e3={}, e0={}",
            result.e1(),
            result.e2(),
            result.e3(),
            result.e0()
        );

        // Expected: point translated by (10, 20, 30)
        assert!(
            abs_diff_eq!(result.x(), 11.0, epsilon = EPS),
            "x: expected 11.0, got {}",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 22.0, epsilon = EPS),
            "y: expected 22.0, got {}",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), 33.0, epsilon = EPS),
            "z: expected 33.0, got {}",
            result.z()
        );
    }

    #[test]
    fn antisandwich_rotation_works() {
        use std::f64::consts::FRAC_PI_2;

        // Create a 90° rotation around Z axis
        let motor = Motor::from_rotation_z(FRAC_PI_2);
        let p = Point::from_cartesian(1.0, 0.0, 0.0);

        // Debug: print motor components
        eprintln!(
            "Motor: s={}, e23={}, e31={}, e12={}, e01={}, e02={}, e03={}, e0123={}",
            motor.s(),
            motor.e23(),
            motor.e31(),
            motor.e12(),
            motor.e01(),
            motor.e02(),
            motor.e03(),
            motor.e0123()
        );
        eprintln!(
            "Point: e1={}, e2={}, e3={}, e0={}",
            p.e1(),
            p.e2(),
            p.e3(),
            p.e0()
        );

        // Use generated antisandwich transformation
        let result = motor.transform_point(&p);

        // Debug: print result components
        eprintln!(
            "Result: e1={}, e2={}, e3={}, e0={}",
            result.e1(),
            result.e2(),
            result.e3(),
            result.e0()
        );

        // Expected: (1, 0, 0) rotated 90° around Z -> (0, 1, 0)
        assert!(
            abs_diff_eq!(result.x(), 0.0, epsilon = EPS),
            "x: expected 0.0, got {}",
            result.x()
        );
        assert!(
            abs_diff_eq!(result.y(), 1.0, epsilon = EPS),
            "y: expected 1.0, got {}",
            result.y()
        );
        assert!(
            abs_diff_eq!(result.z(), 0.0, epsilon = EPS),
            "z: expected 0.0, got {}",
            result.z()
        );
    }
}
