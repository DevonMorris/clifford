//! Domain-specific extensions for 3D Projective GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to 3D projective geometry.

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

    /// Returns the Cartesian z-coordinate (requires `w ≠ 0`).
    ///
    /// # Panics
    ///
    /// Division by zero if `w = 0` (ideal point).
    #[inline]
    pub fn cartesian_z(&self) -> T {
        self.z() / self.w()
    }

    /// Returns true if this is an ideal point (point at infinity).
    #[inline]
    pub fn is_ideal(&self, epsilon: T) -> bool {
        self.w().abs() < epsilon
    }

    /// Returns true if this is a finite point (not at infinity).
    #[inline]
    pub fn is_finite(&self, epsilon: T) -> bool {
        self.w().abs() >= epsilon
    }

    /// Normalizes the homogeneous coordinates so `w = 1` (if finite).
    ///
    /// Returns `None` if this is an ideal point.
    pub fn unitize(&self) -> Option<Self> {
        if self.w().abs() < T::epsilon() {
            None
        } else {
            Some(Self::new(
                self.x() / self.w(),
                self.y() / self.w(),
                self.z() / self.w(),
                T::one(),
            ))
        }
    }

    /// Returns the Cartesian coordinates as a tuple, if finite.
    ///
    /// Returns `None` if this is an ideal point.
    #[inline]
    pub fn to_cartesian(&self) -> Option<(T, T, T)> {
        if self.w().abs() < T::epsilon() {
            None
        } else {
            Some((
                self.x() / self.w(),
                self.y() / self.w(),
                self.z() / self.w(),
            ))
        }
    }

    /// Returns the attitude of the point.
    ///
    /// For a point, the attitude is the weight (e₀ component).
    #[inline]
    pub fn attitude(&self) -> T {
        self.w()
    }

    /// Returns the squared bulk norm of the point.
    ///
    /// The bulk norm is the length of the spatial part: `e1² + e2² + e3²`.
    #[inline]
    pub fn bulk_norm_squared(&self) -> T {
        self.x() * self.x() + self.y() * self.y() + self.z() * self.z()
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
        self.w().abs()
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
        crate::ops::Wedge::wedge(self, other)
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
            (self.x() * other.w() + other.x() * self.w()) / two,
            (self.y() * other.w() + other.y() * self.w()) / two,
            (self.z() * other.w() + other.z() * self.w()) / two,
            self.w() * other.w(),
        )
    }

    /// Inner product (dot product) with another point.
    #[inline]
    pub fn dot(&self, other: &Point<T>) -> T {
        self.x() * other.x() + self.y() * other.y() + self.z() * other.z()
    }

    /// Left contraction onto a line.
    #[inline]
    pub fn left_contract_line(&self, line: &Line<T>) -> Plane<T> {
        // P ⌋ L = P · L + P ∧ L (for grade-1 contracted with grade-2)
        // Result is grade 3 (plane)
        crate::ops::Wedge::wedge(self, line)
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
    ///
    /// # Field mapping
    ///
    /// Plücker coordinates map to the lexicographic field order as:
    /// - direction (e01, e02, e03) → (dir_x, dir_y, dir_z) at positions 2, 4, 5
    /// - moment (e23, e31, e12) → (moment_x, -moment_y, moment_z) at positions 3, 1, 0
    ///
    /// Note: moment_y is negated due to the e31 sign convention.
    #[inline]
    pub fn from_plucker(direction: &EuclideanVector<T>, moment: &EuclideanVector<T>) -> Self {
        // Lexicographic field order: [moment_z, moment_y, dir_x, moment_x, dir_y, dir_z]
        // Position 0: e12 = moment_z, Position 1: e13=e31 = -moment_y
        // Position 2: e14=e01 = dir_x, Position 3: e23 = moment_x
        // Position 4: e24=e02 = dir_y, Position 5: e34=e03 = dir_z
        Self::new_unchecked(
            moment.z(),     // position 0: moment_z (e12)
            -moment.y(),    // position 1: moment_y (e13=e31, negated for sign)
            direction.x(),  // position 2: dir_x (e14=e01)
            moment.x(),     // position 3: moment_x (e23)
            direction.y(),  // position 4: dir_y (e24=e02)
            direction.z(),  // position 5: dir_z (e34=e03)
        )
    }

    /// X-axis (line through origin along x).
    #[inline]
    pub fn x_axis() -> Self {
        // Direction = (1, 0, 0), Moment = (0, 0, 0)
        // Field order: [moment_z, moment_y, dir_x, moment_x, dir_y, dir_z]
        Self::new_unchecked(
            T::zero(),  // moment_z
            T::zero(),  // moment_y
            T::one(),   // dir_x = 1
            T::zero(),  // moment_x
            T::zero(),  // dir_y
            T::zero(),  // dir_z
        )
    }

    /// Y-axis.
    #[inline]
    pub fn y_axis() -> Self {
        // Direction = (0, 1, 0), Moment = (0, 0, 0)
        Self::new_unchecked(
            T::zero(),  // moment_z
            T::zero(),  // moment_y
            T::zero(),  // dir_x
            T::zero(),  // moment_x
            T::one(),   // dir_y = 1
            T::zero(),  // dir_z
        )
    }

    /// Z-axis.
    #[inline]
    pub fn z_axis() -> Self {
        // Direction = (0, 0, 1), Moment = (0, 0, 0)
        Self::new_unchecked(
            T::zero(),  // moment_z
            T::zero(),  // moment_y
            T::zero(),  // dir_x
            T::zero(),  // moment_x
            T::zero(),  // dir_y
            T::one(),   // dir_z = 1
        )
    }

    /// Direction vector (e01, e02, e03).
    ///
    /// Returns the direction part of the Plücker coordinates.
    #[inline]
    pub fn direction(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.dir_x(), self.dir_y(), self.dir_z())
    }

    /// Moment vector (e23, e31, e12).
    ///
    /// Returns the moment part of the Plücker coordinates.
    /// Note: moment_y is negated to account for the e31 sign convention.
    #[inline]
    pub fn moment(&self) -> EuclideanVector<T> {
        EuclideanVector::new(self.moment_x(), -self.moment_y(), self.moment_z())
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
            self.moment_x() / wn,
            self.moment_y() / wn,
            self.moment_z() / wn,
            self.dir_x() / wn,
            self.dir_y() / wn,
            self.dir_z() / wn,
        )
    }

    /// Plücker condition residual: d · m.
    #[inline]
    pub fn plucker_residual(&self) -> T {
        self.moment_x() * self.dir_x()
            + self.moment_y() * self.dir_y()
            + self.moment_z() * self.dir_z()
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
        self.moment_x() * other.moment_x()
            + self.moment_y() * other.moment_y()
            + self.moment_z() * other.moment_z()
            + self.dir_x() * other.dir_x()
            + self.dir_y() * other.dir_y()
            + self.dir_z() * other.dir_z()
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
        crate::ops::Wedge::wedge(self, point)
    }

    /// Meet with a plane to find intersection point (regressive product).
    ///
    /// Returns the point where this line intersects the plane. If the line
    /// is parallel to the plane (no intersection), returns an ideal point.
    #[inline]
    pub fn meet_plane(&self, plane: &Plane<T>) -> Point<T> {
        crate::ops::Antiwedge::antiwedge(self, plane)
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
        self.moment_x() * self.moment_x()
            + self.moment_y() * self.moment_y()
            + self.moment_z() * self.moment_z()
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
        self.moment_x() * other.dir_x()
            + self.moment_y() * other.dir_y()
            + self.moment_z() * other.dir_z()
            + other.moment_x() * self.dir_x()
            + other.moment_y() * self.dir_y()
            + other.moment_z() * self.dir_z()
    }

    /// Distance from a point to this line.
    ///
    /// The line should be unitized for accurate results.
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        let d1 = self.moment_x();
        let d2 = self.moment_y();
        let d3 = self.moment_z();
        let m1 = self.dir_x();
        let m2 = self.dir_y();
        let m3 = self.dir_z();

        let px = p.x();
        let py = p.y();
        let pz = p.z();
        let pw = p.w();

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
    /// The plane equation is `n·x + dist = 0`.
    pub fn from_normal_and_distance(cart_nx: T, cart_ny: T, cart_nz: T, distance: T) -> Self {
        // PGA trivector basis mapping:
        // - nx (field) = e123 (ideal plane component)
        // - ny (field) = e012 = Cartesian nz (plane perpendicular to z)
        // - nz (field) = e031 = Cartesian ny (plane perpendicular to y), negated for e31 convention
        // - dist (field) = e023 = Cartesian nx (plane perpendicular to x)
        Self::new(distance, cart_nz, -cart_ny, cart_nx)
    }

    /// XY plane (z = 0).
    #[inline]
    pub fn xy() -> Self {
        // Normal (0, 0, 1) in Cartesian: d=0, nz=1, ny=0, nx=0
        Self::new(T::zero(), T::one(), T::zero(), T::zero())
    }

    /// XZ plane (y = 0).
    #[inline]
    pub fn xz() -> Self {
        // Normal (0, 1, 0) in Cartesian: d=0, nz=0, ny=-1 (e031 sign), nx=0
        Self::new(T::zero(), T::zero(), -T::one(), T::zero())
    }

    /// YZ plane (x = 0).
    #[inline]
    pub fn yz() -> Self {
        // Normal (1, 0, 0) in Cartesian: d=0, nz=0, ny=0, nx=1
        Self::new(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Normal vector (Cartesian).
    ///
    /// Inverse of from_normal_and_distance mapping.
    #[inline]
    pub fn normal(&self) -> EuclideanVector<T> {
        // PGA field → Cartesian normal mapping:
        // nx (e023) → cart_nx, -ny (e031) → cart_ny, nz (e012) → cart_nz
        EuclideanVector::new(self.nx(), -self.ny(), self.nz())
    }

    /// Distance from origin (signed).
    ///
    /// Stored in the dist field (e123 ideal component).
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        self.dist()
    }

    /// Weight norm (normal magnitude).
    #[inline]
    pub fn weight_norm(&self) -> T {
        self.normal().norm()
    }

    /// Bulk norm (distance component magnitude).
    #[inline]
    pub fn bulk_norm(&self) -> T {
        self.dist().abs()
    }

    /// Unitize to unit normal.
    pub fn unitized(&self) -> Self {
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }
        // Field order is [dist, nz, ny, nx]
        Self::new(
            self.dist() / wn,
            self.nz() / wn,
            self.ny() / wn,
            self.nx() / wn,
        )
    }

    /// Inner product with another plane.
    #[inline]
    pub fn dot(&self, other: &Plane<T>) -> T {
        self.nx() * other.nx() + self.ny() * other.ny() + self.nz() * other.nz()
    }

    /// Attitude (ideal line at infinity).
    #[inline]
    pub fn attitude(&self) -> T {
        self.dist()
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
        crate::ops::Antiwedge::antiwedge(self, other)
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
        (p.nx() * point.x() + p.ny() * point.y() + p.nz() * point.z() + p.dist() * point.w())
            / point.w()
    }

    /// Project a point onto this plane.
    pub fn project_point(&self, point: &Point<T>) -> Point<T> {
        let dist = self.signed_distance(point);
        let n = self.normal().normalized();
        Point::new(
            point.x() - dist * n.x() * point.w(),
            point.y() - dist * n.y() * point.w(),
            point.z() - dist * n.z() * point.w(),
            point.w(),
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
    /// In PGA with the sandwich product, the identity is the scalar s = 1.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Motor, Point};
    /// use clifford::ops::Transform;
    ///
    /// let m = Motor::<f64>::identity();
    /// let p = Point::from_cartesian(1.0, 2.0, 3.0);
    /// let p2 = m.transform(&p);
    /// assert!((p2.x() - p.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        // In PGA with antiproduct-based composition and antisandwich-based
        // transforms, the identity motor is the pseudoscalar e0123=1, not
        // the scalar s=1. This is because:
        // - compose() uses antiproduct: M1 ⊛ M2 = ∁(∁M1 × ∁M2)
        // - transform() uses antisandwich: M ⊛ X ⊛ antirev(M)
        // The pseudoscalar is the identity for the antiproduct.
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
    /// Translation motor in RGA convention: T = τₓe₂₃ - τᵧe₃₁ + τᵤe₁₂ + e₀₁₂₃
    /// where τ = (dx/2, dy/2, dz/2) is half the displacement.
    /// Note: e31 has opposite sign due to cyclic ordering (e31 = -e13).
    ///
    /// This uses e0123 (antiscalar) as the identity, matching rotation motors,
    /// so they compose correctly under antiproduct.
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self {
        let half = T::one() / T::TWO;
        // Motor fields: [s, tz, ty, tx, rx, ry, rz, ps]
        // Translation uses tz, ty, tx (positions 1-3)
        Self::new_unchecked(
            T::zero(),  // s
            dz * half,  // tz - z-translation
            -dy * half, // ty - y-translation (negated due to e13 sign)
            dx * half,  // tx - x-translation
            T::zero(),  // rx
            T::zero(),  // ry
            T::zero(),  // rz
            T::one(),   // ps (identity)
        )
    }

    /// Pure rotation around x-axis through origin.
    pub fn from_rotation_x(angle: T) -> Self {
        let half = angle / T::TWO;
        // Motor fields: [s, tz, ty, tx, rx, ry, rz, ps]
        // Rotation around x uses rx (position 4)
        Self::new_unchecked(
            T::zero(),   // s
            T::zero(),   // tz
            T::zero(),   // ty
            T::zero(),   // tx
            -half.sin(), // rx - rotation around x-axis
            T::zero(),   // ry
            T::zero(),   // rz
            half.cos(),  // ps = cos(θ/2)
        )
    }

    /// Pure rotation around y-axis through origin.
    pub fn from_rotation_y(angle: T) -> Self {
        let half = angle / T::TWO;
        // Motor fields: [s, tz, ty, tx, rx, ry, rz, ps]
        // Rotation around y uses ry (position 5)
        Self::new_unchecked(
            T::zero(),   // s
            T::zero(),   // tz
            T::zero(),   // ty
            T::zero(),   // tx
            T::zero(),   // rx
            -half.sin(), // ry - rotation around y-axis
            T::zero(),   // rz
            half.cos(),  // ps = cos(θ/2)
        )
    }

    /// Pure rotation around z-axis through origin.
    ///
    /// Uses positions 4-6 (bx, ty, tz) for rotation encoding, matching original.
    pub fn from_rotation_z(angle: T) -> Self {
        let half = angle / T::TWO;
        // Motor fields (lexicographic): [s, bz(e12), by(e13), tx(e14), bx(e23), ty(e24), tz(e34), ps(e0123)]
        // Rotation around z uses position 6 (tz = e34)
        Self::new_unchecked(
            T::zero(),   // s
            T::zero(),   // bz (e12)
            T::zero(),   // by (e13)
            T::zero(),   // tx (e14)
            T::zero(),   // bx (e23)
            T::zero(),   // ty (e24)
            -half.sin(), // tz (e34) - rotation around z-axis
            half.cos(),  // ps (e0123) = cos(θ/2)
        )
    }

    /// Rotation around arbitrary axis through origin.
    ///
    /// The axis vector determines the rotation axis (will be normalized).
    /// The rotation follows the right-hand rule.
    ///
    /// Uses positions 4-6 (bx, ty, tz) for rotation encoding.
    pub fn from_axis_angle(axis: &EuclideanVector<T>, angle: T) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let axis_norm = axis.normalized();
        // Motor fields (lexicographic): [s, bz(e12), by(e13), tx(e14), bx(e23), ty(e24), tz(e34), ps(e0123)]
        // Rotation uses positions 4, 5, 6 (bx, ty, tz)
        Self::new_unchecked(
            T::zero(),                 // s
            T::zero(),                 // bz (e12)
            T::zero(),                 // by (e13)
            T::zero(),                 // tx (e14)
            -sin_half * axis_norm.x(), // bx (e23) - rotation around x-axis
            -sin_half * axis_norm.y(), // ty (e24) - rotation around y-axis
            -sin_half * axis_norm.z(), // tz (e34) - rotation around z-axis
            cos_half,                  // ps (e0123) = cos(θ/2)
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
            sin_a * m.x() + half_dist * cos_a * d.x(), // e23 (translation)
            sin_a * m.y() + half_dist * cos_a * d.y(), // e31 (translation)
            sin_a * m.z() + half_dist * cos_a * d.z(), // e12 (translation)
            sin_a * d.x(),                             // e01 (rotation)
            sin_a * d.y(),                             // e02 (rotation)
            sin_a * d.z(),                             // e03 (rotation)
            -half_dist * sin_a,
        )
    }

    /// Inverse motor.
    ///
    /// Uses the weight norm squared (rotation part: e01² + e02² + e03² + e0123²)
    /// for proper inversion. For a unit motor, this equals 1.
    pub fn inverse(&self) -> Self {
        use crate::norm::DegenerateNormed;
        let wn_sq = self.weight_norm_squared();
        let rev = self.reverse();
        // new_unchecked signature: (s, bz, by, tx, bx, ty, tz, ps)
        Self::new_unchecked(
            rev.s() / wn_sq,
            rev.tz() / wn_sq,
            rev.ty() / wn_sq,
            rev.tx() / wn_sq,
            rev.rx() / wn_sq,
            rev.ry() / wn_sq,
            rev.rz() / wn_sq,
            rev.ps() / wn_sq,
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
        // new_unchecked signature: (s, bz, by, tx, bx, ty, tz, ps)
        Self::new_unchecked(
            self.s() / bn,
            self.tz() / bn,
            self.ty() / bn,
            self.tx() / bn,
            self.rx() / bn,
            self.ry() / bn,
            self.rz() / bn,
            self.ps() / bn,
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
        self.s() * self.ps() + self.rx() * self.tx() + self.ty() * self.ry() + self.tz() * self.rz()
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
        // from_translation sets: bx = dz/2, by = -dy/2, bz = dx/2
        // So: dx = 2*bz, dy = -2*by, dz = 2*bx
        EuclideanVector::new(T::TWO * self.tz(), -T::TWO * self.ty(), T::TWO * self.rx())
    }
}

// ============================================================================
// Flector extensions
// ============================================================================

impl<T: Float> Flector<T> {
    /// Create flector from reflection plane.
    pub fn from_plane(plane: &Plane<T>) -> Self {
        let p = plane.unitized();
        // Flector field order: [px, py, pz, pw, dist, nz, ny, nx]
        Self::new_unchecked(
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            p.dist(),
            p.nz(),
            p.ny(),
            p.nx(),
        )
    }

    /// Create reflection through plane at origin with Cartesian normal.
    pub fn from_plane_through_origin(cart_nx: T, cart_ny: T, cart_nz: T) -> Self {
        let norm = (cart_nx * cart_nx + cart_ny * cart_ny + cart_nz * cart_nz).sqrt();
        // Convert Cartesian to PGA: nx=cart_nx, ny=-cart_ny (e031 sign), nz=cart_nz, dist=0
        // Flector field order: [px, py, pz, pw, dist, nz, ny, nx]
        Self::new_unchecked(
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),       // dist
            cart_nz / norm,  // nz
            -cart_ny / norm, // ny (negated for e031)
            cart_nx / norm,  // nx
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
        Point::new(self.px(), self.py(), self.pz(), self.pw())
    }

    /// Plane part (grade 3).
    #[inline]
    pub fn plane_part(&self) -> Plane<T> {
        // Plane::new expects [dist, nz, ny, nx]
        Plane::new(self.dist(), self.nz(), self.ny(), self.nx())
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
        // Flector field order: [px, py, pz, pw, dist, nz, ny, nx]
        Self::new_unchecked(
            self.px() / bn,
            self.py() / bn,
            self.pz() / bn,
            self.pw() / bn,
            self.dist() / bn,
            self.nz() / bn,
            self.ny() / bn,
            self.nx() / bn,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::Transform;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;

    #[test]
    fn motor_identity_preserves_point() {
        let p = Point::<f64>::from_cartesian(3.0, 4.0, 5.0);
        let m = Motor::<f64>::identity();

        eprintln!(
            "Identity motor: s={}, bx={}, by={}, bz={}, tx={}, ty={}, tz={}, ps={}",
            m.s(),
            m.rx(),
            m.ty(),
            m.tz(),
            m.tx(),
            m.ry(),
            m.rz(),
            m.ps()
        );
        eprintln!("Input point: ({}, {}, {}, {})", p.x(), p.y(), p.z(), p.w());

        let result = m.transform(&p);

        eprintln!(
            "Output point: ({}, {}, {}, {})",
            result.x(),
            result.y(),
            result.z(),
            result.w()
        );

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
            result.z(),
            p.z(),
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
        let t = Motor::<f64>::from_translation(2.0, 0.0, 0.0);

        eprintln!(
            "Motor T: s={}, bx={}, by={}, bz={}, tx={}, ty={}, tz={}, ps={}",
            t.s(),
            t.rx(),
            t.ty(),
            t.tz(),
            t.tx(),
            t.ry(),
            t.rz(),
            t.ps()
        );
        eprintln!(
            "Origin: ({}, {}, {}, {})",
            origin.x(),
            origin.y(),
            origin.z(),
            origin.w()
        );

        let result = t.transform(&origin);

        eprintln!(
            "Translated: ({}, {}, {}, {})",
            result.x(),
            result.y(),
            result.z(),
            result.w()
        );

        // Expected: (2, 0, 0) with w=1
        assert!(relative_eq!(
            result.x(),
            2.0,
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
        assert!(relative_eq!(
            result.w(),
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
            "Motor T: s={}, bx={}, by={}, bz={}, tx={}, ty={}, tz={}, ps={}",
            t.s(),
            t.rx(),
            t.ty(),
            t.tz(),
            t.tx(),
            t.ry(),
            t.rz(),
            t.ps()
        );
        eprintln!(
            "Origin: ({}, {}, {}, {})",
            origin.x(),
            origin.y(),
            origin.z(),
            origin.w()
        );

        let result = t.transform(&origin);

        eprintln!(
            "Translated: ({}, {}, {}, {})",
            result.x(),
            result.y(),
            result.z(),
            result.w()
        );

        // Expected: (0, -19.26, 0) with w=1
        assert!(relative_eq!(
            result.x(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.y(),
            -19.26,
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

    #[test]
    fn motor_translation_z() {
        let origin = Point::<f64>::origin();
        let t = Motor::<f64>::from_translation(0.0, 0.0, 5.0);

        eprintln!(
            "Motor T: s={}, bx={}, by={}, bz={}, tx={}, ty={}, tz={}, ps={}",
            t.s(),
            t.rx(),
            t.ty(),
            t.tz(),
            t.tx(),
            t.ry(),
            t.rz(),
            t.ps()
        );

        let result = t.transform(&origin);

        eprintln!(
            "Translated: ({}, {}, {}, {})",
            result.x(),
            result.y(),
            result.z(),
            result.w()
        );

        // Expected: (0, 0, 5) with w=1
        assert!(relative_eq!(
            result.x(),
            0.0,
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
            "Flector: px={}, py={}, pz={}, pw={}, nx={}, ny={}, nz={}, dist={}",
            f.px(),
            f.py(),
            f.pz(),
            f.pw(),
            f.nx(),
            f.ny(),
            f.nz(),
            f.dist()
        );
        eprintln!("Input point: ({}, {}, {}, {})", p.x(), p.y(), p.z(), p.w());

        let result = f.transform(&p);

        eprintln!(
            "Output point: ({}, {}, {}, {})",
            result.x(),
            result.y(),
            result.z(),
            result.w()
        );

        // Reflecting (1, 2, 3) through XY plane should give (1, 2, -3)
        // Use Cartesian coordinates (normalized by w) since homogeneous coords may have sign flip
        assert!(relative_eq!(
            result.cartesian_x(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.cartesian_y(),
            2.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            result.cartesian_z(),
            -3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn flector_reflect_yz_plane() {
        let f = Flector::<f64>::reflect_yz();
        let p = Point::<f64>::from_cartesian(1.0, 2.0, 3.0);

        let result = f.transform(&p);

        // Reflecting (1, 2, 3) through YZ plane should give (-1, 2, 3)
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
        assert!(relative_eq!(
            result.cartesian_z(),
            3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn flector_reflect_xz_plane() {
        let f = Flector::<f64>::reflect_xz();
        let p = Point::<f64>::from_cartesian(1.0, 2.0, 3.0);

        let result = f.transform(&p);

        // Reflecting (1, 2, 3) through XZ plane should give (1, -2, 3)
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
        assert!(relative_eq!(
            result.cartesian_z(),
            3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
