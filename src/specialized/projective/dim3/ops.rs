//! Operations for 3D Projective Geometric Algebra types.

use super::types::{Flector, Line, Motor, Plane, Point};
use crate::scalar::Float;

// ============================================================================
// Point operations
// ============================================================================

impl<T: Float> Point<T> {
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
        let dz = self.z() - other.z();
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Squared distance to another point.
    #[inline]
    pub fn distance_squared(&self, other: &Point<T>) -> T {
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        let dz = self.z() - other.z();
        dx * dx + dy * dy + dz * dz
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
            (p1_norm.e3 + p2_norm.e3) / T::TWO,
        )
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
    /// use clifford::specialized::projective::dim3::{Point, Motor};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let p = Point::new(1.0, 0.0, 0.0);
    ///
    /// // 90° rotation around Z axis
    /// let rotor = Motor::from_rotation_z(FRAC_PI_2);
    /// let rotated = rotor.transform_point(&p);
    /// assert!(abs_diff_eq!(rotated.x(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.y(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(rotated.z(), 0.0, epsilon = 1e-10));
    ///
    /// // Translation
    /// let trans = Motor::from_translation(2.0, 3.0, 4.0);
    /// let translated = trans.transform_point(&p);
    /// assert!(abs_diff_eq!(translated.x(), 3.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(translated.y(), 3.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(translated.z(), 4.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // Sandwich product: M P M̃
        // Using the unified formula from Rigid Geometric Algebra:
        // a = v × p + pw * m
        // p' = p + 2(s * a + v × a + e0123 * pw * v)
        //
        // where v = (e23, e31, e12) is the rotation bivector,
        // m = (e01, e02, e03) is the translation bivector.

        let s = self.s;
        let b23 = self.e23;
        let b31 = self.e31;
        let b12 = self.e12;
        let b01 = self.e01;
        let b02 = self.e02;
        let b03 = self.e03;
        let i = self.e0123;

        let px = p.e1;
        let py = p.e2;
        let pz = p.e3;
        let pw = p.e0;

        let two = T::TWO;

        // Compute intermediate vector a = v × p + pw * m
        let ax = b31 * pz - b12 * py + pw * b01;
        let ay = b12 * px - b23 * pz + pw * b02;
        let az = b23 * py - b31 * px + pw * b03;

        // Compute v × a
        let vxa_x = b31 * az - b12 * ay;
        let vxa_y = b12 * ax - b23 * az;
        let vxa_z = b23 * ay - b31 * ax;

        // Final transformation: p' = p + 2(s * a + v × a + e0123 * pw * v)
        // Note: sign of e0123 term depends on convention; we use + here
        Point {
            e1: px + two * (s * ax + vxa_x + i * pw * b23),
            e2: py + two * (s * ay + vxa_y + i * pw * b31),
            e3: pz + two * (s * az + vxa_z + i * pw * b12),
            e0: pw,
        }
    }

    /// Transforms a line: `L' = M L M̃`.
    ///
    /// This applies the rigid transformation represented by the motor to the line.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Line, Motor};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// // X axis
    /// let line = Line::x_axis();
    ///
    /// // 90° rotation around Z axis
    /// let rotation = Motor::from_rotation_z(FRAC_PI_2);
    /// let rotated = rotation.transform_line(&line);
    ///
    /// // X axis rotated 90° around Z becomes Y axis
    /// let d = rotated.direction();
    /// assert!(abs_diff_eq!(d.x, 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(d.y, 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(d.z, 0.0, epsilon = 1e-10));
    /// ```
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        // For a line L = d + m (direction bivector + moment bivector)
        // The transformation is L' = M L M̃
        //
        // For a pure rotation R, both direction and moment transform as vectors:
        // d' = R d R̃ (rotate direction)
        // m' = R m R̃ (rotate moment)
        //
        // For translation, the direction is unchanged but the moment changes.

        let s = self.s;
        let b23 = self.e23;
        let b31 = self.e31;
        let b12 = self.e12;
        let b01 = self.e01;
        let b02 = self.e02;
        let b03 = self.e03;

        // Line components
        let d1 = l.e01; // direction x
        let d2 = l.e02; // direction y
        let d3 = l.e03; // direction z
        let m1 = l.e23; // moment x
        let m2 = l.e31; // moment y
        let m3 = l.e12; // moment z

        let two = T::TWO;

        // Rotate the direction (like a vector): d' = R d R̃
        // Using Rodrigues formula: d' = d + 2s(v × d) + 2(v × (v × d))
        // where v = (b23, b31, b12) is the rotation bivector

        // v × d
        let vxd_x = b31 * d3 - b12 * d2;
        let vxd_y = b12 * d1 - b23 * d3;
        let vxd_z = b23 * d2 - b31 * d1;

        // v × (v × d)
        let vxvxd_x = b31 * vxd_z - b12 * vxd_y;
        let vxvxd_y = b12 * vxd_x - b23 * vxd_z;
        let vxvxd_z = b23 * vxd_y - b31 * vxd_x;

        let d1_new = d1 + two * (s * vxd_x + vxvxd_x);
        let d2_new = d2 + two * (s * vxd_y + vxvxd_y);
        let d3_new = d3 + two * (s * vxd_z + vxvxd_z);

        // Rotate the moment: m' = R m R̃
        let vxm_x = b31 * m3 - b12 * m2;
        let vxm_y = b12 * m1 - b23 * m3;
        let vxm_z = b23 * m2 - b31 * m1;

        let vxvxm_x = b31 * vxm_z - b12 * vxm_y;
        let vxvxm_y = b12 * vxm_x - b23 * vxm_z;
        let vxvxm_z = b23 * vxm_y - b31 * vxm_x;

        let m1_rot = m1 + two * (s * vxm_x + vxvxm_x);
        let m2_rot = m2 + two * (s * vxm_y + vxvxm_y);
        let m3_rot = m3 + two * (s * vxm_z + vxvxm_z);

        // Translation contribution to moment: m' += t × d'
        // where t = (b01, b02, b03) is the translation
        let txd_x = b02 * d3_new - b03 * d2_new;
        let txd_y = b03 * d1_new - b01 * d3_new;
        let txd_z = b01 * d2_new - b02 * d1_new;

        Line {
            e01: d1_new,
            e02: d2_new,
            e03: d3_new,
            e23: m1_rot + two * txd_x,
            e31: m2_rot + two * txd_y,
            e12: m3_rot + two * txd_z,
        }
    }
}

// ============================================================================
// Line operations
// ============================================================================

impl<T: Float> Line<T> {
    /// Computes the meet of this line with a plane.
    ///
    /// Returns the point where the line intersects the plane.
    /// If the line is parallel to the plane, returns an ideal point (point at infinity).
    ///
    /// # Formula
    ///
    /// The meet is computed as the regressive product: `P = L ∨ G`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Line, Plane, Point};
    /// use approx::abs_diff_eq;
    ///
    /// // Line along Z axis through origin
    /// let line = Line::z_axis();
    ///
    /// // XY plane at z = 5
    /// let plane = Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);
    ///
    /// let intersection = line.meet(&plane);
    ///
    /// // Should intersect at (0, 0, 5)
    /// assert!(abs_diff_eq!(intersection.x(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(intersection.y(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(intersection.z(), 5.0, epsilon = 1e-10));
    /// ```
    pub fn meet(&self, plane: &Plane<T>) -> Point<T> {
        // The meet/regressive product L ∨ G gives a point where the line
        // intersects the plane.
        //
        // Using the Plücker coordinate formula for line-plane intersection:
        // P = (s·d + n×m, n·d)
        //
        // Where:
        // - d = line direction (e₀₁, e₀₂, e₀₃)
        // - m = line moment (e₂₃, e₃₁, e₁₂)
        // - n = plane normal (e₀₂₃, e₀₃₁, e₀₁₂)
        // - s = plane signed distance = -e₁₂₃ (our convention uses n·x + d = 0)
        //
        // For our basis, this gives:
        // e₁: -g₄·d₁ + g₂·m₃ - g₃·m₂
        // e₂: -g₄·d₂ + g₃·m₁ - g₁·m₃
        // e₃: -g₄·d₃ + g₁·m₂ - g₂·m₁
        // e₀: g₁·d₁ + g₂·d₂ + g₃·d₃

        let d1 = self.e01;
        let d2 = self.e02;
        let d3 = self.e03;
        let m1 = self.e23;
        let m2 = self.e31;
        let m3 = self.e12;

        let g1 = plane.e023; // nx
        let g2 = plane.e031; // ny
        let g3 = plane.e012; // nz
        let g4 = plane.e123; // -s (plane distance in our convention)

        Point {
            e1: -g4 * d1 + g2 * m3 - g3 * m2,
            e2: -g4 * d2 + g3 * m1 - g1 * m3,
            e3: -g4 * d3 + g1 * m2 - g2 * m1,
            e0: g1 * d1 + g2 * d2 + g3 * d3,
        }
    }

    /// Computes the signed distance from a point to this line.
    ///
    /// The line should be unitized for correct distance values.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Line, Point};
    /// use approx::abs_diff_eq;
    ///
    /// // Z axis (unitized)
    /// let line = Line::z_axis();
    ///
    /// // Point at (3, 4, 0) - distance 5 from Z axis
    /// let p = Point::new(3.0, 4.0, 0.0);
    /// assert!(abs_diff_eq!(line.distance_to_point(&p), 5.0, epsilon = 1e-10));
    /// ```
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        // For a unitized line L and point P:
        // distance = |L ∨ P| / |d|
        // where L ∨ P gives a plane and we take its weight norm

        let d1 = self.e01;
        let d2 = self.e02;
        let d3 = self.e03;
        let m1 = self.e23;
        let m2 = self.e31;
        let m3 = self.e12;

        let px = p.e1;
        let py = p.e2;
        let pz = p.e3;
        let pw = p.e0;

        // The regressive product L ∨ P gives a plane
        // We compute the weight norm of this plane divided by |d| and |pw|

        // Cross product of direction and point position
        let cx = d2 * pz - d3 * py;
        let cy = d3 * px - d1 * pz;
        let cz = d1 * py - d2 * px;

        // Add moment contribution
        let nx = cx + m1 * pw;
        let ny = cy + m2 * pw;
        let nz = cz + m3 * pw;

        let d_norm = self.weight_norm();
        (nx * nx + ny * ny + nz * nz).sqrt() / (d_norm * pw.abs())
    }

    /// Computes the angle between two lines.
    ///
    /// Returns the acute angle in radians, in the range `[0, π/2]`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Line;
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
    /// assert!(abs_diff_eq!(x_axis.angle(&x_axis), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn angle(&self, other: &Line<T>) -> T {
        // Use direction vectors: cos(θ) = |d₁·d₂| / (|d₁||d₂|)
        let d1 = self.direction();
        let d2 = other.direction();

        let dot = d1.x * d2.x + d1.y * d2.y + d1.z * d2.z;
        let norm1 = (d1.x * d1.x + d1.y * d1.y + d1.z * d1.z).sqrt();
        let norm2 = (d2.x * d2.x + d2.y * d2.y + d2.z * d2.z).sqrt();

        let denom = norm1 * norm2;
        if denom < T::epsilon() {
            return T::zero();
        }

        // Use absolute value for acute angle
        let cos_theta = (dot / denom).abs();
        // Clamp to [0, 1] to handle numerical errors
        let clamped = if cos_theta > T::one() {
            T::one()
        } else {
            cos_theta
        };
        clamped.acos()
    }

    /// Computes the closest point on this line to a given point.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Line, Point};
    /// use approx::abs_diff_eq;
    ///
    /// // Z axis
    /// let line = Line::z_axis();
    ///
    /// // Point at (3, 4, 5)
    /// let p = Point::new(3.0, 4.0, 5.0);
    /// let closest = line.closest_point(&p);
    ///
    /// // Closest point on Z axis should be (0, 0, 5)
    /// assert!(abs_diff_eq!(closest.x(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(closest.y(), 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(closest.z(), 5.0, epsilon = 1e-10));
    /// ```
    pub fn closest_point(&self, p: &Point<T>) -> Point<T> {
        // Project point onto line
        // P_closest = P_line + ((P - P_line) · d) * d / |d|²
        // where P_line is any point on the line

        let d = self.direction();
        let m = self.moment();
        let d_sq = d.x * d.x + d.y * d.y + d.z * d.z;

        if d_sq < T::epsilon() {
            return Point::origin();
        }

        // A point on the line: P_line = d × m / |d|² (when line doesn't pass through origin)
        // For line through origin (m = 0), use origin
        let line_pt_x = (d.y * m.z - d.z * m.y) / d_sq;
        let line_pt_y = (d.z * m.x - d.x * m.z) / d_sq;
        let line_pt_z = (d.x * m.y - d.y * m.x) / d_sq;

        // Vector from line point to given point
        let px = p.x() - line_pt_x;
        let py = p.y() - line_pt_y;
        let pz = p.z() - line_pt_z;

        // Project onto direction
        let t = (px * d.x + py * d.y + pz * d.z) / d_sq;

        Point::new(
            line_pt_x + t * d.x,
            line_pt_y + t * d.y,
            line_pt_z + t * d.z,
        )
    }
}

// ============================================================================
// Plane operations
// ============================================================================

impl<T: Float> Plane<T> {
    /// Computes the signed distance from a point to this plane.
    ///
    /// Positive distance means the point is on the side of the normal.
    /// The plane should be unitized for correct distance values.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Plane, Point};
    /// use approx::abs_diff_eq;
    ///
    /// let plane = Plane::xy(); // z = 0 plane
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// assert!(abs_diff_eq!(plane.signed_distance(&p), 3.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn signed_distance(&self, p: &Point<T>) -> T {
        // For plane nx*x + ny*y + nz*z + d = 0,
        // signed distance = (nx*px + ny*py + nz*pz + d*pw) / (pw * |n|)
        let numerator = self.e023 * p.e1 + self.e031 * p.e2 + self.e012 * p.e3 + self.e123 * p.e0;
        numerator / (p.e0 * self.weight_norm())
    }

    /// Returns true if the point lies on this plane (within epsilon).
    #[inline]
    pub fn contains_point(&self, p: &Point<T>, epsilon: T) -> bool {
        self.signed_distance(p).abs() < epsilon
    }

    /// Computes the dihedral angle between two planes.
    ///
    /// Returns the angle in radians, in the range `[0, π]`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Plane;
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let xy_plane: Plane<f64> = Plane::xy();
    /// let xz_plane: Plane<f64> = Plane::xz();
    ///
    /// // XY and XZ planes are perpendicular
    /// assert!(abs_diff_eq!(xy_plane.angle(&xz_plane), FRAC_PI_2, epsilon = 1e-10));
    ///
    /// // Same plane has angle 0
    /// assert!(abs_diff_eq!(xy_plane.angle(&xy_plane), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn angle(&self, other: &Plane<T>) -> T {
        // Use normal vectors: cos(θ) = (n₁·n₂) / (|n₁||n₂|)
        let n1 = self.normal();
        let n2 = other.normal();

        let dot = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
        let norm1 = (n1.x * n1.x + n1.y * n1.y + n1.z * n1.z).sqrt();
        let norm2 = (n2.x * n2.x + n2.y * n2.y + n2.z * n2.z).sqrt();

        let denom = norm1 * norm2;
        if denom < T::epsilon() {
            return T::zero();
        }

        let cos_theta = dot / denom;
        // Clamp to [-1, 1] to handle numerical errors
        let clamped = if cos_theta > T::one() {
            T::one()
        } else if cos_theta < -T::one() {
            -T::one()
        } else {
            cos_theta
        };
        clamped.acos()
    }

    /// Computes the angle between this plane and a line.
    ///
    /// Returns the angle in radians, in the range `[0, π/2]`.
    /// This is the angle between the line and the plane surface,
    /// not the angle between the line and the normal.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Plane, Line};
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::abs_diff_eq;
    ///
    /// let xy_plane: Plane<f64> = Plane::xy();
    /// let z_axis: Line<f64> = Line::z_axis();
    /// let x_axis: Line<f64> = Line::x_axis();
    ///
    /// // Z axis is perpendicular to XY plane (angle π/2)
    /// assert!(abs_diff_eq!(xy_plane.angle_to_line(&z_axis), FRAC_PI_2, epsilon = 1e-10));
    ///
    /// // X axis lies in XY plane (angle 0)
    /// assert!(abs_diff_eq!(xy_plane.angle_to_line(&x_axis), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn angle_to_line(&self, line: &Line<T>) -> T {
        // sin(θ) = |n·d| / (|n||d|)
        // where θ is the angle between line and plane surface
        let n = self.normal();
        let d = line.direction();

        let dot = (n.x * d.x + n.y * d.y + n.z * d.z).abs();
        let norm_n = (n.x * n.x + n.y * n.y + n.z * n.z).sqrt();
        let norm_d = (d.x * d.x + d.y * d.y + d.z * d.z).sqrt();

        let denom = norm_n * norm_d;
        if denom < T::epsilon() {
            return T::zero();
        }

        let sin_theta = dot / denom;
        // Clamp to [0, 1] to handle numerical errors
        let clamped = if sin_theta > T::one() {
            T::one()
        } else {
            sin_theta
        };
        clamped.asin()
    }

    /// Computes the meet (intersection) of this plane with another plane.
    ///
    /// Returns the line where the two planes intersect.
    /// If the planes are parallel, returns a degenerate line (zero direction).
    ///
    /// # Formula
    ///
    /// The meet is computed as the regressive product: `L = G₁ ∨ G₂`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Plane, Line};
    /// use approx::abs_diff_eq;
    ///
    /// // XY plane (z = 0) and XZ plane (y = 0) intersect along X axis
    /// let xy_plane: Plane<f64> = Plane::xy();
    /// let xz_plane: Plane<f64> = Plane::xz();
    ///
    /// let line = xy_plane.meet(&xz_plane);
    ///
    /// // Direction should be along X axis (or its negative)
    /// let d = line.direction();
    /// assert!(abs_diff_eq!(d.y, 0.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(d.z, 0.0, epsilon = 1e-10));
    /// assert!(d.x.abs() > 0.9); // dx = ±1
    /// ```
    pub fn meet(&self, other: &Plane<T>) -> Line<T> {
        // The regressive product of two planes gives a line
        // G₁ = n₁e₀₂₃ + n₂e₀₃₁ + n₃e₀₁₂ + d₁e₁₂₃
        // G₂ = m₁e₀₂₃ + m₂e₀₃₁ + m₃e₀₁₂ + d₂e₁₂₃
        //
        // The result is the cross product of normals (direction) and
        // a moment computed from the distance terms.

        let n1 = self.e023;
        let n2 = self.e031;
        let n3 = self.e012;
        let d1 = self.e123;

        let m1 = other.e023;
        let m2 = other.e031;
        let m3 = other.e012;
        let d2 = other.e123;

        // Direction is n × m (cross product of normals)
        let dir_x = n2 * m3 - n3 * m2;
        let dir_y = n3 * m1 - n1 * m3;
        let dir_z = n1 * m2 - n2 * m1;

        // Moment encodes the position of the line
        // m = d₁ * n₂ - d₂ * n₁ (for each component pair)
        let mom_x = d1 * m1 - d2 * n1;
        let mom_y = d1 * m2 - d2 * n2;
        let mom_z = d1 * m3 - d2 * n3;

        Line {
            e01: dir_x,
            e02: dir_y,
            e03: dir_z,
            e23: mom_x,
            e31: mom_y,
            e12: mom_z,
        }
    }
}

// ============================================================================
// Flector operations on geometric objects
// ============================================================================

impl<T: Float> Flector<T> {
    /// Transforms a point via the flector sandwich product: `P' = F P F̃`.
    ///
    /// For a pure plane reflection, this reflects the point through the plane.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Plane, Point, Flector};
    /// use approx::abs_diff_eq;
    ///
    /// // Reflection through XY plane (z = 0)
    /// let flector = Flector::reflect_xy();
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// let reflected = flector.transform_point(&p);
    ///
    /// assert!(abs_diff_eq!(reflected.x(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(reflected.y(), 2.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(reflected.z(), -3.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // For a pure plane reflection g, the transformation is:
        // P' = g P g (since g² = -|g|² for a plane through origin)
        //
        // For a general flector F = p + g (point + plane), the formula is more complex.
        // We use the explicit formula from geometric algebra.
        //
        // For pure reflection through plane with normal (nx, ny, nz) and distance d:
        // The reflected point is: P' = P - 2 * (P · n + d) * n
        //
        // In PGA terms with the sandwich product F P F̃:

        let px = p.e1;
        let py = p.e2;
        let pz = p.e3;
        let pw = p.e0;

        // Plane components (grade 3)
        let gx = self.e023; // normal x
        let gy = self.e031; // normal y
        let gz = self.e012; // normal z
        let gw = self.e123; // distance parameter

        // Point components of flector (grade 1)
        let fx = self.e1;
        let fy = self.e2;
        let fz = self.e3;
        let fw = self.e0;

        // For a pure plane reflection (fx = fy = fz = fw = 0):
        // The formula simplifies to standard reflection formula.
        //
        // For a general flector, we compute the full sandwich product.
        // This involves the geometric antiproduct in PGA.

        // Simplified formula for pure plane reflection:
        // P' = P - 2 * ((P · g) / |g|²) * g
        // where P · g = px*gx + py*gy + pz*gz + pw*gw

        let g_norm_sq = gx * gx + gy * gy + gz * gz;

        if g_norm_sq < T::epsilon() {
            // Degenerate plane, return original point
            return *p;
        }

        // For pure reflection (common case), use optimized formula
        if self.is_pure_reflection(T::epsilon()) {
            // Dot product of point with plane normal + distance term
            let dot = px * gx + py * gy + pz * gz + pw * gw;
            let factor = T::TWO * dot / g_norm_sq;

            return Point {
                e1: px - factor * gx,
                e2: py - factor * gy,
                e3: pz - factor * gz,
                e0: pw,
            };
        }

        // General flector case (point + plane)
        // This is a rotoreflection or glide reflection
        // Full sandwich product computation would go here
        // For now, we handle just the pure reflection case and
        // use an approximation for the general case

        // The general formula involves computing F P F̃ where F = point + plane
        // This is complex and involves the geometric antiproduct

        // Simplified: treat as pure plane reflection plus point contribution
        // This is an approximation; full formula would need proper antiproduct
        let dot = px * gx + py * gy + pz * gz + pw * gw;
        let factor = T::TWO * dot / g_norm_sq;

        // Add point contribution (this creates glide/rotoreflection effect)
        let scale = T::TWO / g_norm_sq;

        Point {
            e1: px - factor * gx + scale * (fy * gz - fz * gy + fw * gx),
            e2: py - factor * gy + scale * (fz * gx - fx * gz + fw * gy),
            e3: pz - factor * gz + scale * (fx * gy - fy * gx + fw * gz),
            e0: pw,
        }
    }

    /// Composes two flectors via geometric product.
    ///
    /// The composition of two reflections yields a rotation (motor).
    /// This returns a Motor representing the composed transformation.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Flector, Motor, Point};
    /// use approx::abs_diff_eq;
    ///
    /// // Two reflections through perpendicular planes = 180° rotation
    /// let f1 = Flector::reflect_yz(); // x = 0
    /// let f2 = Flector::reflect_xz(); // y = 0
    /// let motor = f1.compose(&f2);
    ///
    /// // Apply to a point
    /// let p = Point::new(1.0, 1.0, 0.0);
    /// let result = motor.transform_point(&p);
    ///
    /// // 180° rotation around Z axis: (1,1,0) -> (-1,-1,0)
    /// assert!(abs_diff_eq!(result.x(), -1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(result.y(), -1.0, epsilon = 1e-10));
    /// ```
    pub fn compose(&self, other: &Flector<T>) -> Motor<T> {
        // Flector * Flector = Motor (odd * odd = even)
        // For pure plane reflections g1 and g2:
        // g1 * g2 = (n1 · n2) + (n1 × n2) as a rotor

        let g1x = self.e023;
        let g1y = self.e031;
        let g1z = self.e012;
        let g1w = self.e123;

        let g2x = other.e023;
        let g2y = other.e031;
        let g2z = other.e012;
        let g2w = other.e123;

        // For pure plane reflections, the product gives a rotor
        // Scalar part: g1 · g2 = g1x*g2x + g1y*g2y + g1z*g2z
        // Bivector part: g1 × g2 (cross product of normals)

        let s = g1x * g2x + g1y * g2y + g1z * g2z;

        // Bivector components from cross product
        // e23: rotation around x (from y×z)
        // e31: rotation around y (from z×x)
        // e12: rotation around z (from x×y)
        let e23 = g1y * g2z - g1z * g2y;
        let e31 = g1z * g2x - g1x * g2z;
        let e12 = g1x * g2y - g1y * g2x;

        // Translation components from distance terms
        let e01 = g1w * g2x - g1x * g2w;
        let e02 = g1w * g2y - g1y * g2w;
        let e03 = g1w * g2z - g1z * g2w;

        // Pseudoscalar from distance terms
        let e0123 = g1w * g2w;

        Motor::new(s, e23, e31, e12, e01, e02, e03, e0123)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::specialized::euclidean::dim3::Vector as EuclideanVector;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;

    #[test]
    fn motor_identity_preserves_point() {
        let p: Point<f64> = Point::new(3.0, 4.0, 5.0);
        let m: Motor<f64> = Motor::identity();
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), p.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_rotation_90_degrees_z() {
        let p = Point::new(1.0, 0.0, 0.0);
        let m = Motor::from_rotation_z(std::f64::consts::FRAC_PI_2);
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_translation() {
        let p = Point::new(1.0, 2.0, 3.0);
        let m = Motor::from_translation(3.0, 4.0, 5.0);
        let result = m.transform_point(&p);

        assert!(abs_diff_eq!(result.x(), 4.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 6.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 8.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_composition() {
        let p = Point::new(1.0, 0.0, 0.0);

        // First rotate 90° around Z, then translate by (1, 2, 3)
        let rotation = Motor::from_rotation_z(std::f64::consts::FRAC_PI_2);
        let translation = Motor::from_translation(1.0, 2.0, 3.0);
        let composed = rotation.compose(&translation);

        let result = composed.transform_point(&p);

        // (1,0,0) -> rotated to (0,1,0) -> translated to (1,3,3)
        assert!(abs_diff_eq!(result.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn point_distance() {
        let p1: Point<f64> = Point::new(0.0, 0.0, 0.0);
        let p2: Point<f64> = Point::new(3.0, 4.0, 0.0);
        assert!(abs_diff_eq!(
            p1.distance(&p2),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn point_midpoint() {
        let p1: Point<f64> = Point::new(0.0, 0.0, 0.0);
        let p2: Point<f64> = Point::new(4.0, 6.0, 8.0);
        let mid = p1.midpoint(&p2);

        assert!(abs_diff_eq!(mid.x(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(mid.y(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(mid.z(), 4.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Flector tests
    // ========================================================================

    #[test]
    fn flector_reflect_xy_plane() {
        let f: Flector<f64> = Flector::reflect_xy();
        let p = Point::new(1.0, 2.0, 3.0);
        let result = f.transform_point(&p);

        // XY plane reflection negates z
        assert!(abs_diff_eq!(result.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), -3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flector_reflect_xz_plane() {
        let f: Flector<f64> = Flector::reflect_xz();
        let p = Point::new(1.0, 2.0, 3.0);
        let result = f.transform_point(&p);

        // XZ plane reflection negates y
        assert!(abs_diff_eq!(result.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), -2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flector_reflect_yz_plane() {
        let f: Flector<f64> = Flector::reflect_yz();
        let p = Point::new(1.0, 2.0, 3.0);
        let result = f.transform_point(&p);

        // YZ plane reflection negates x
        assert!(abs_diff_eq!(result.x(), -1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flector_double_reflection_is_identity() {
        let f: Flector<f64> = Flector::reflect_xy();
        let p = Point::new(1.0, 2.0, 3.0);

        // Reflecting twice should return to original
        let once = f.transform_point(&p);
        let twice = f.transform_point(&once);

        assert!(abs_diff_eq!(twice.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(twice.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(twice.z(), p.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flector_composition_two_perpendicular_reflections() {
        // Two reflections through perpendicular planes = 180° rotation
        let f1: Flector<f64> = Flector::reflect_yz(); // x = 0
        let f2: Flector<f64> = Flector::reflect_xz(); // y = 0
        let motor = f1.compose(&f2);

        let p = Point::new(1.0, 1.0, 0.0);
        let result = motor.transform_point(&p);

        // 180° rotation around Z axis: (1,1,0) -> (-1,-1,0)
        assert!(abs_diff_eq!(result.x(), -1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.y(), -1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(result.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_signed_distance() {
        let plane: Plane<f64> = Plane::xy(); // z = 0
        let p_above = Point::new(0.0, 0.0, 5.0);
        let p_below = Point::new(0.0, 0.0, -3.0);
        let p_on = Point::new(1.0, 2.0, 0.0);

        assert!(abs_diff_eq!(
            plane.signed_distance(&p_above),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            plane.signed_distance(&p_below),
            -3.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            plane.signed_distance(&p_on),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // ========================================================================
    // Line tests
    // ========================================================================

    #[test]
    fn line_join_two_points() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let line = Line::join(&p1, &p2);

        // Direction should be (1, 0, 0)
        let d = line.direction();
        assert!(abs_diff_eq!(d.x, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));

        // Moment should be zero (line through origin)
        let m = line.moment();
        assert!(abs_diff_eq!(m.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(m.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(m.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_join_offset_points() {
        // Line through (0, 1, 0) and (1, 1, 0) - parallel to X axis at y=1
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(1.0, 1.0, 0.0);
        let line = Line::join(&p1, &p2);

        // Direction should be (1, 0, 0)
        let d = line.direction();
        assert!(abs_diff_eq!(d.x, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));

        // Moment should be non-zero (line not through origin)
        // For this line, moment = direction × point_on_line = (1,0,0) × (0,1,0) = (0,0,1)
        // But we need to check the actual formula
        assert!(!line.through_origin(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_from_point_and_direction() {
        let p = Point::new(1.0, 2.0, 3.0);
        let line = Line::from_point_and_direction(&p, &EuclideanVector::new(0.0, 0.0, 1.0));

        // Direction should be (0, 0, 1)
        let d = line.direction();
        assert!(abs_diff_eq!(d.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_axis_through_origin() {
        let x_axis = Line::x_axis();
        let y_axis = Line::y_axis();
        let z_axis = Line::z_axis();

        assert!(x_axis.through_origin(ABS_DIFF_EQ_EPS));
        assert!(y_axis.through_origin(ABS_DIFF_EQ_EPS));
        assert!(z_axis.through_origin(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_meet_plane() {
        // Z axis meets XY plane at z=5 at point (0, 0, 5)
        let line = Line::z_axis();
        let plane = Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);

        let intersection = line.meet(&plane);
        assert!(abs_diff_eq!(
            intersection.x(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            intersection.y(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            intersection.z(),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_meet_plane_origin() {
        // Z axis meets XY plane at origin
        let line: Line<f64> = Line::z_axis();
        let plane: Plane<f64> = Plane::xy();

        let intersection = line.meet(&plane);
        assert!(abs_diff_eq!(
            intersection.x(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            intersection.y(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            intersection.z(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_meet_plane() {
        // XY plane meets XZ plane along X axis
        let xy_plane: Plane<f64> = Plane::xy();
        let xz_plane: Plane<f64> = Plane::xz();

        let line = xy_plane.meet(&xz_plane);

        // Direction should be along X axis
        let d = line.direction();
        assert!(abs_diff_eq!(d.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(d.x.abs() > 0.9); // dx = ±1

        // Line should pass through origin
        assert!(line.through_origin(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_distance_to_point() {
        // Z axis, point at (3, 4, 0) - distance 5
        let line: Line<f64> = Line::z_axis();
        let p = Point::new(3.0, 4.0, 0.0);

        assert!(abs_diff_eq!(
            line.distance_to_point(&p),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_closest_point() {
        // Z axis, closest point to (3, 4, 5) is (0, 0, 5)
        let line: Line<f64> = Line::z_axis();
        let p = Point::new(3.0, 4.0, 5.0);
        let closest = line.closest_point(&p);

        assert!(abs_diff_eq!(closest.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(closest.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(closest.z(), 5.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_transform_line_rotation() {
        // Rotate X axis 90° around Z to get Y axis
        let line: Line<f64> = Line::x_axis();
        let rotation = Motor::from_rotation_z(std::f64::consts::FRAC_PI_2);
        let rotated = rotation.transform_line(&line);

        let d = rotated.direction();
        assert!(abs_diff_eq!(d.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.y, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_transform_line_translation() {
        // Translate Z axis by (1, 0, 0)
        let line: Line<f64> = Line::z_axis();
        let translation = Motor::from_translation(1.0, 0.0, 0.0);
        let translated = translation.transform_line(&line);

        // Direction should still be (0, 0, 1)
        let d = translated.direction();
        assert!(abs_diff_eq!(d.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 1.0, epsilon = ABS_DIFF_EQ_EPS));

        // But line should no longer pass through origin
        assert!(!translated.through_origin(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_parallel_check() {
        let line1: Line<f64> = Line::z_axis();
        let line2 = Line::from_point_and_direction(
            &Point::new(1.0, 0.0, 0.0),
            &EuclideanVector::new(0.0, 0.0, 1.0),
        );

        // These lines are parallel (both in Z direction)
        assert!(line1.is_parallel(&line2, ABS_DIFF_EQ_EPS));

        // But they should intersect (Plücker inner = 0 for parallel lines)
        assert!(line1.intersects(&line2, ABS_DIFF_EQ_EPS));

        // X and Y axes are not parallel
        let x_axis = Line::x_axis();
        let y_axis = Line::y_axis();
        assert!(!x_axis.is_parallel(&y_axis, ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_skew_check() {
        // Create two skew lines
        let line1 = Line::z_axis(); // Z axis through origin
        let line2 = Line::from_point_and_direction(
            &Point::new(1.0, 0.0, 0.0),
            &EuclideanVector::new(0.0, 1.0, 0.0),
        ); // Y direction at x=1

        // These lines are skew (non-parallel, non-intersecting)
        assert!(!line1.is_parallel(&line2, ABS_DIFF_EQ_EPS));
        assert!(!line1.intersects(&line2, ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Unary operation tests
    // ========================================================================

    #[test]
    fn point_geometric_norm() {
        // Point at (3, 4, 0) has distance 5 from origin
        let p = Point::new(3.0, 4.0, 0.0);
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
        let p = Point::new(3.0, 4.0, 0.0);
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
        // Z axis through origin has geometric norm 0
        let z_axis: Line<f64> = Line::z_axis();
        assert!(abs_diff_eq!(
            z_axis.geometric_norm(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Line at x=3, parallel to z has geometric norm 3
        let offset_line = Line::from_point_and_direction(
            &Point::new(3.0, 0.0, 0.0),
            &EuclideanVector::new(0.0, 0.0, 1.0),
        );
        let unitized = offset_line.unitized();
        assert!(abs_diff_eq!(
            unitized.geometric_norm(),
            3.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_bulk_weight_norm() {
        // Z axis: direction (0,0,1), moment (0,0,0)
        let z_axis: Line<f64> = Line::z_axis();
        assert!(abs_diff_eq!(
            z_axis.weight_norm(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            z_axis.bulk_norm(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Attitude should be the direction
        let att = z_axis.attitude();
        assert!(abs_diff_eq!(att.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(att.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(att.z, 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_reverse() {
        let line = Line::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        let rev = line.reverse();

        // Reverse negates all components
        assert!(abs_diff_eq!(rev.e01, -1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e02, -2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e03, -3.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e23, -4.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e31, -5.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e12, -6.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_geometric_norm() {
        // XY plane at z=0 has distance 0 from origin
        let xy: Plane<f64> = Plane::xy();
        assert!(abs_diff_eq!(
            xy.geometric_norm(),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Plane z=5 has distance 5 from origin
        let offset = Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);
        assert!(abs_diff_eq!(
            offset.geometric_norm(),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_reverse() {
        let plane = Plane::new(1.0, 2.0, 3.0, 4.0);
        let rev = plane.reverse();

        // Reverse negates all components
        assert!(abs_diff_eq!(rev.e023, -1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e031, -2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e012, -3.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rev.e123, -4.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_attitude() {
        let plane = Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);
        let att = plane.attitude();

        assert!(abs_diff_eq!(att.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(att.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(att.z, 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_inverse() {
        // Rotation motor
        let rotation = Motor::from_rotation_z(std::f64::consts::FRAC_PI_4);
        let inv = rotation.inverse();
        let identity = rotation.compose(&inv);

        // Should give identity
        assert!(abs_diff_eq!(identity.s, 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(identity.e12, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_unitized() {
        // Create non-unit motor by scaling
        let rotation = Motor::from_rotation_z(std::f64::consts::FRAC_PI_4);
        let scaled = Motor::new(
            rotation.s * 2.0,
            rotation.e23 * 2.0,
            rotation.e31 * 2.0,
            rotation.e12 * 2.0,
            rotation.e01 * 2.0,
            rotation.e02 * 2.0,
            rotation.e03 * 2.0,
            rotation.e0123 * 2.0,
        );

        // Unitize and check norm
        let unitized = scaled.unitized();
        assert!(abs_diff_eq!(
            unitized.weight_norm(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn motor_is_unitized() {
        let rotation = Motor::from_rotation_z(std::f64::consts::FRAC_PI_4);
        assert!(rotation.is_unitized(ABS_DIFF_EQ_EPS));

        let identity = Motor::identity();
        assert!(identity.is_unitized(ABS_DIFF_EQ_EPS));
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
        let z_axis: Line<f64> = Line::z_axis();
        let parallel = Line::from_point_and_direction(
            &Point::new(1.0, 0.0, 0.0),
            &EuclideanVector::new(0.0, 0.0, 1.0),
        );
        assert!(abs_diff_eq!(
            z_axis.angle(&parallel),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_angle_symmetric() {
        let x_axis: Line<f64> = Line::x_axis();
        let y_axis: Line<f64> = Line::y_axis();
        assert!(abs_diff_eq!(
            x_axis.angle(&y_axis),
            y_axis.angle(&x_axis),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_angle_perpendicular() {
        let xy_plane: Plane<f64> = Plane::xy();
        let xz_plane: Plane<f64> = Plane::xz();
        assert!(abs_diff_eq!(
            xy_plane.angle(&xz_plane),
            std::f64::consts::FRAC_PI_2,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_angle_same() {
        let xy_plane: Plane<f64> = Plane::xy();
        assert!(abs_diff_eq!(
            xy_plane.angle(&xy_plane),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_angle_symmetric() {
        let xy_plane: Plane<f64> = Plane::xy();
        let xz_plane: Plane<f64> = Plane::xz();
        assert!(abs_diff_eq!(
            xy_plane.angle(&xz_plane),
            xz_plane.angle(&xy_plane),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_angle_to_line_perpendicular() {
        let xy_plane: Plane<f64> = Plane::xy();
        let z_axis: Line<f64> = Line::z_axis();
        assert!(abs_diff_eq!(
            xy_plane.angle_to_line(&z_axis),
            std::f64::consts::FRAC_PI_2,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_angle_to_line_parallel() {
        let xy_plane: Plane<f64> = Plane::xy();
        let x_axis: Line<f64> = Line::x_axis();
        assert!(abs_diff_eq!(
            xy_plane.angle_to_line(&x_axis),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }
}
