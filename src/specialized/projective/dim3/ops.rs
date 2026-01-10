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

    /// Computes the inner product of two points.
    ///
    /// In PGA, the inner product of two grade-1 elements (points) gives a scalar.
    /// For normalized points, this equals the dot product of their Euclidean parts
    /// plus the product of their weights.
    ///
    /// # Formula
    ///
    /// For points `P₁ = (e₁, e₂, e₃, e₀)` and `P₂ = (e₁', e₂', e₃', e₀')`:
    /// `P₁ · P₂ = e₁·e₁' + e₂·e₂' + e₃·e₃'`
    ///
    /// Note: The `e₀` components don't contribute because `e₀² = 0` in PGA.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Point;
    /// use approx::abs_diff_eq;
    ///
    /// let p1 = Point::new(1.0, 0.0, 0.0);
    /// let p2 = Point::new(1.0, 0.0, 0.0);
    /// // Same direction -> dot = 1
    /// assert!(abs_diff_eq!(p1.dot(&p2), 1.0, epsilon = 1e-10));
    ///
    /// let p3 = Point::new(0.0, 1.0, 0.0);
    /// // Perpendicular -> dot = 0
    /// assert!(abs_diff_eq!(p1.dot(&p3), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn dot(&self, other: &Point<T>) -> T {
        // In PGA with signature (+,+,+,0), e₀² = 0
        // So only the Euclidean (e₁, e₂, e₃) parts contribute
        self.e1 * other.e1 + self.e2 * other.e2 + self.e3 * other.e3
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

    /// Computes the commutator of two motors.
    ///
    /// `[M₁, M₂] = M₁M₂ - M₂M₁`
    ///
    /// The commutator measures the non-commutativity of two transformations.
    /// For rotations around the same axis, the commutator is zero.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Motor;
    /// use std::f64::consts::FRAC_PI_4;
    /// use approx::abs_diff_eq;
    ///
    /// // Rotations around the same axis commute
    /// let r1 = Motor::from_rotation_z(FRAC_PI_4);
    /// let r2 = Motor::from_rotation_z(FRAC_PI_4 / 2.0);
    /// let comm = r1.commutator(&r2);
    ///
    /// // Result should be close to zero
    /// assert!(abs_diff_eq!(comm.s, 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn commutator(&self, other: &Motor<T>) -> Motor<T> {
        // [M₁, M₂] = M₁M₂ - M₂M₁
        let m1m2 = self.compose(other);
        let m2m1 = other.compose(self);

        Motor {
            s: m1m2.s - m2m1.s,
            e23: m1m2.e23 - m2m1.e23,
            e31: m1m2.e31 - m2m1.e31,
            e12: m1m2.e12 - m2m1.e12,
            e01: m1m2.e01 - m2m1.e01,
            e02: m1m2.e02 - m2m1.e02,
            e03: m1m2.e03 - m2m1.e03,
            e0123: m1m2.e0123 - m2m1.e0123,
        }
    }

    /// Computes the anticommutator of two motors.
    ///
    /// `{M₁, M₂} = M₁M₂ + M₂M₁`
    ///
    /// The anticommutator extracts the symmetric part of the product.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Motor;
    /// use std::f64::consts::FRAC_PI_4;
    /// use approx::abs_diff_eq;
    ///
    /// // Same motor anticommuted with itself = 2 * M²
    /// let m = Motor::from_rotation_z(FRAC_PI_4);
    /// let anti = m.anticommutator(&m);
    /// let m_sq = m.compose(&m);
    ///
    /// assert!(abs_diff_eq!(anti.s, 2.0 * m_sq.s, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn anticommutator(&self, other: &Motor<T>) -> Motor<T> {
        // {M₁, M₂} = M₁M₂ + M₂M₁
        let m1m2 = self.compose(other);
        let m2m1 = other.compose(self);

        Motor {
            s: m1m2.s + m2m1.s,
            e23: m1m2.e23 + m2m1.e23,
            e31: m1m2.e31 + m2m1.e31,
            e12: m1m2.e12 + m2m1.e12,
            e01: m1m2.e01 + m2m1.e01,
            e02: m1m2.e02 + m2m1.e02,
            e03: m1m2.e03 + m2m1.e03,
            e0123: m1m2.e0123 + m2m1.e0123,
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

    /// Computes the perpendicular distance between two lines.
    ///
    /// For skew lines, returns the shortest distance between them.
    /// For intersecting lines, returns 0.
    /// For parallel lines, returns the constant perpendicular distance.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Line, Point};
    /// use clifford::specialized::euclidean::dim3::Vector;
    /// use approx::abs_diff_eq;
    ///
    /// // Z axis and a parallel line at (3, 4, 0)
    /// let z_axis: Line<f64> = Line::z_axis();
    /// let parallel = Line::from_point_and_direction(
    ///     &Point::new(3.0, 4.0, 0.0),
    ///     &Vector::new(0.0, 0.0, 1.0),
    /// );
    ///
    /// // Distance should be 5 (sqrt(3² + 4²))
    /// assert!(abs_diff_eq!(z_axis.distance(&parallel), 5.0, epsilon = 1e-10));
    ///
    /// // Intersecting lines have distance 0
    /// let x_axis: Line<f64> = Line::x_axis();
    /// let y_axis: Line<f64> = Line::y_axis();
    /// assert!(abs_diff_eq!(x_axis.distance(&y_axis), 0.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn distance(&self, other: &Line<T>) -> T {
        let d1 = self.direction();
        let d2 = other.direction();

        // Cross product of directions
        let cross_x = d1.y * d2.z - d1.z * d2.y;
        let cross_y = d1.z * d2.x - d1.x * d2.z;
        let cross_z = d1.x * d2.y - d1.y * d2.x;
        let cross_norm = (cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).sqrt();

        if cross_norm < T::epsilon() {
            // Lines are parallel - compute distance using closest point approach
            // Get a point on each line and compute perpendicular distance
            let m1 = self.moment();
            let m2 = other.moment();
            let d1_sq = d1.x * d1.x + d1.y * d1.y + d1.z * d1.z;

            if d1_sq < T::epsilon() {
                return T::zero();
            }

            // Point on line 1: d1 × m1 / |d1|²
            let p1_x = (d1.y * m1.z - d1.z * m1.y) / d1_sq;
            let p1_y = (d1.z * m1.x - d1.x * m1.z) / d1_sq;
            let p1_z = (d1.x * m1.y - d1.y * m1.x) / d1_sq;

            // Point on line 2: d2 × m2 / |d2|²
            let d2_sq = d2.x * d2.x + d2.y * d2.y + d2.z * d2.z;
            if d2_sq < T::epsilon() {
                return T::zero();
            }
            let p2_x = (d2.y * m2.z - d2.z * m2.y) / d2_sq;
            let p2_y = (d2.z * m2.x - d2.x * m2.z) / d2_sq;
            let p2_z = (d2.x * m2.y - d2.y * m2.x) / d2_sq;

            // Vector between points
            let v_x = p2_x - p1_x;
            let v_y = p2_y - p1_y;
            let v_z = p2_z - p1_z;

            // Project onto perpendicular direction (perpendicular to d1)
            // For parallel lines, we need the component perpendicular to d1
            let proj = (v_x * d1.x + v_y * d1.y + v_z * d1.z) / d1_sq;
            let perp_x = v_x - proj * d1.x;
            let perp_y = v_y - proj * d1.y;
            let perp_z = v_z - proj * d1.z;

            (perp_x * perp_x + perp_y * perp_y + perp_z * perp_z).sqrt()
        } else {
            // Skew or intersecting lines
            // Distance = |d1·m2 + d2·m1| / |d1×d2|
            let plucker = self.plucker_inner(other);
            plucker.abs() / cross_norm
        }
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

    /// Computes the inner product of two lines.
    ///
    /// In PGA, the inner product of two grade-2 elements (lines) gives a scalar.
    /// This is related to the Plücker inner product and measures the
    /// "reciprocity" of the two lines.
    ///
    /// # Formula
    ///
    /// For lines `L₁ = d₁ + m₁` and `L₂ = d₂ + m₂` (direction + moment):
    /// `L₁ · L₂ = d₁·m₂ + d₂·m₁`
    ///
    /// This equals zero when the lines intersect (including parallel lines).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Line, Point};
    /// use clifford::specialized::euclidean::dim3::Vector;
    /// use approx::abs_diff_eq;
    ///
    /// // Intersecting lines have dot = 0
    /// let x_axis: Line<f64> = Line::x_axis();
    /// let y_axis: Line<f64> = Line::y_axis();
    /// assert!(abs_diff_eq!(x_axis.dot(&y_axis), 0.0, epsilon = 1e-10));
    ///
    /// // Skew lines have non-zero dot
    /// let z_axis: Line<f64> = Line::z_axis();
    /// let skew = Line::from_point_and_direction(
    ///     &Point::new(1.0, 0.0, 0.0),
    ///     &Vector::new(0.0, 1.0, 0.0),
    /// );
    /// // Non-zero for skew lines
    /// assert!(z_axis.dot(&skew).abs() > 0.0);
    /// ```
    #[inline]
    pub fn dot(&self, other: &Line<T>) -> T {
        // The inner product of two bivectors in PGA
        // This is the Plücker inner product: d₁·m₂ + d₂·m₁
        self.plucker_inner(other)
    }

    /// Computes the plane containing this line and a point.
    ///
    /// This is the join (exterior product) `L ∧ P`, giving the unique plane
    /// that contains both the line and the point.
    ///
    /// If the point lies on the line, returns a degenerate plane.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Line, Plane, Point};
    /// use approx::abs_diff_eq;
    ///
    /// // Z axis
    /// let line: Line<f64> = Line::z_axis();
    ///
    /// // Point at (1, 0, 0) - not on the line
    /// let p = Point::new(1.0, 0.0, 0.0);
    /// let plane = line.join_point(&p);
    ///
    /// // Plane should be the XZ plane (normal in Y direction)
    /// let n = plane.normal();
    /// assert!(abs_diff_eq!(n.x, 0.0, epsilon = 1e-10));
    /// assert!(n.y.abs() > 0.9); // ny = ±1
    /// assert!(abs_diff_eq!(n.z, 0.0, epsilon = 1e-10));
    /// ```
    pub fn join_point(&self, p: &Point<T>) -> Plane<T> {
        // L ∧ P: line (bivector) ∧ point (vector) = trivector (plane)
        //
        // Line L = d₁e₀₁ + d₂e₀₂ + d₃e₀₃ + m₁e₂₃ + m₂e₃₁ + m₃e₁₂
        // Point P = p₁e₁ + p₂e₂ + p₃e₃ + p₀e₀
        //
        // Computing each exterior product term:
        // e₀₁ ∧ e₂ = e₀₁₂,  e₀₁ ∧ e₃ = -e₀₃₁
        // e₀₂ ∧ e₁ = -e₀₁₂, e₀₂ ∧ e₃ = e₀₂₃
        // e₀₃ ∧ e₁ = e₀₃₁,  e₀₃ ∧ e₂ = -e₀₂₃
        // e₂₃ ∧ e₀ = e₀₂₃,  e₂₃ ∧ e₁ = e₁₂₃
        // e₃₁ ∧ e₀ = -e₀₃₁, e₃₁ ∧ e₂ = -e₁₂₃
        // e₁₂ ∧ e₀ = e₀₁₂,  e₁₂ ∧ e₃ = e₁₂₃
        //
        // Collecting terms:
        // e₀₂₃: d₂p₃ - d₃p₂ + m₁p₀
        // e₀₃₁: d₃p₁ - d₁p₃ - m₂p₀
        // e₀₁₂: d₁p₂ - d₂p₁ + m₃p₀
        // e₁₂₃: m₁p₁ - m₂p₂ + m₃p₃

        let d1 = self.e01;
        let d2 = self.e02;
        let d3 = self.e03;
        let m1 = self.e23;
        let m2 = self.e31;
        let m3 = self.e12;

        let p0 = p.e0;
        let p1 = p.e1;
        let p2 = p.e2;
        let p3 = p.e3;

        Plane {
            e023: d2 * p3 - d3 * p2 + m1 * p0,
            e031: d3 * p1 - d1 * p3 - m2 * p0,
            e012: d1 * p2 - d2 * p1 + m3 * p0,
            e123: m1 * p1 - m2 * p2 + m3 * p3,
        }
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

    /// Computes the inner product of two planes.
    ///
    /// In PGA, the inner product of two grade-3 elements (planes) gives a scalar.
    /// This measures the alignment of the two planes' normals.
    ///
    /// # Formula
    ///
    /// For planes `G₁ = (n₁, d₁)` and `G₂ = (n₂, d₂)` (normal + distance):
    /// `G₁ · G₂ = n₁·n₂`
    ///
    /// The `e₁₂₃` (distance) components don't contribute because `e₁₂₃² = 0` in PGA.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::Plane;
    /// use approx::abs_diff_eq;
    ///
    /// let xy_plane: Plane<f64> = Plane::xy();
    /// let xz_plane: Plane<f64> = Plane::xz();
    ///
    /// // Perpendicular planes have dot = 0
    /// assert!(abs_diff_eq!(xy_plane.dot(&xz_plane), 0.0, epsilon = 1e-10));
    ///
    /// // Same plane has dot = 1 (for unit normals)
    /// assert!(abs_diff_eq!(xy_plane.dot(&xy_plane), 1.0, epsilon = 1e-10));
    /// ```
    #[inline]
    pub fn dot(&self, other: &Plane<T>) -> T {
        // Inner product of two trivectors in PGA
        // Only the normal parts (e023, e031, e012) contribute
        // e123² = 0 in the degenerate metric
        self.e023 * other.e023 + self.e031 * other.e031 + self.e012 * other.e012
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

    /// Projects a point onto this plane.
    ///
    /// Returns the closest point on the plane to the given point.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Plane, Point};
    /// use approx::abs_diff_eq;
    ///
    /// // XY plane (z = 0)
    /// let plane: Plane<f64> = Plane::xy();
    ///
    /// // Point at (1, 2, 5)
    /// let p = Point::new(1.0, 2.0, 5.0);
    /// let projected = plane.project_point(&p);
    ///
    /// // Projected point should be (1, 2, 0)
    /// assert!(abs_diff_eq!(projected.x(), 1.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(projected.y(), 2.0, epsilon = 1e-10));
    /// assert!(abs_diff_eq!(projected.z(), 0.0, epsilon = 1e-10));
    /// ```
    pub fn project_point(&self, p: &Point<T>) -> Point<T> {
        // P' = P - (signed_distance) * n
        // where n is the unit normal and signed_distance is the distance from P to plane
        let n = self.normal();
        let n_norm_sq = n.x * n.x + n.y * n.y + n.z * n.z;

        if n_norm_sq < T::epsilon() {
            return *p;
        }

        // Signed distance = (n·p + d) / |n| where p is the point position
        // For our representation: (e023*p1 + e031*p2 + e012*p3 + e123*p0) / (p0 * |n|)
        let numerator = self.e023 * p.e1 + self.e031 * p.e2 + self.e012 * p.e3 + self.e123 * p.e0;
        let dist = numerator / (p.e0 * n_norm_sq.sqrt());

        Point::new(
            p.x() - dist * n.x / n_norm_sq.sqrt(),
            p.y() - dist * n.y / n_norm_sq.sqrt(),
            p.z() - dist * n.z / n_norm_sq.sqrt(),
        )
    }

    /// Projects a line onto this plane.
    ///
    /// Returns the orthogonal projection of the line onto the plane.
    /// If the line is perpendicular to the plane, returns a degenerate line.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::projective::dim3::{Plane, Line, Point};
    /// use clifford::specialized::euclidean::dim3::Vector;
    /// use approx::abs_diff_eq;
    ///
    /// // XY plane (z = 0)
    /// let plane: Plane<f64> = Plane::xy();
    ///
    /// // Line through (0, 0, 1) in direction (1, 1, 1) - diagonal line
    /// let line = Line::from_point_and_direction(
    ///     &Point::new(0.0, 0.0, 1.0),
    ///     &Vector::new(1.0, 1.0, 1.0),
    /// );
    /// let projected = plane.project_line(&line);
    ///
    /// // Projected direction should be (1, 1, 0) (normalized)
    /// let d = projected.direction();
    /// // Z component should be 0
    /// assert!(abs_diff_eq!(d.z, 0.0, epsilon = 1e-10));
    /// ```
    pub fn project_line(&self, line: &Line<T>) -> Line<T> {
        // Project both the direction and a point on the line
        // Then construct the projected line from the projected point and direction

        let n = self.normal();
        let n_norm_sq = n.x * n.x + n.y * n.y + n.z * n.z;

        if n_norm_sq < T::epsilon() {
            return *line;
        }

        let n_norm = n_norm_sq.sqrt();

        // Get direction and project onto plane
        let d = line.direction();
        let d_dot_n = d.x * n.x + d.y * n.y + d.z * n.z;
        let projected_dx = d.x - d_dot_n * n.x / n_norm_sq;
        let projected_dy = d.y - d_dot_n * n.y / n_norm_sq;
        let projected_dz = d.z - d_dot_n * n.z / n_norm_sq;

        // Get a point on the line and project it
        let m = line.moment();
        let d_sq = d.x * d.x + d.y * d.y + d.z * d.z;
        if d_sq < T::epsilon() {
            return Line::zero();
        }

        // Point on line: d × m / |d|²
        let pt_x = (d.y * m.z - d.z * m.y) / d_sq;
        let pt_y = (d.z * m.x - d.x * m.z) / d_sq;
        let pt_z = (d.x * m.y - d.y * m.x) / d_sq;

        // Project point onto plane
        let pt_dot_n = pt_x * n.x + pt_y * n.y + pt_z * n.z;
        let plane_dist = self.e123 / n_norm;
        let signed_dist = (pt_dot_n + plane_dist) / n_norm;

        let proj_pt_x = pt_x - signed_dist * n.x / n_norm;
        let proj_pt_y = pt_y - signed_dist * n.y / n_norm;
        let proj_pt_z = pt_z - signed_dist * n.z / n_norm;

        // Construct projected line from point and direction
        let proj_point = Point::new(proj_pt_x, proj_pt_y, proj_pt_z);
        let proj_direction = crate::specialized::euclidean::dim3::Vector::new(
            projected_dx,
            projected_dy,
            projected_dz,
        );

        Line::from_point_and_direction(&proj_point, &proj_direction)
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

    // ========================================================================
    // Line-Line distance tests
    // ========================================================================

    #[test]
    fn line_distance_intersecting() {
        // X axis and Y axis intersect at origin
        let x_axis: Line<f64> = Line::x_axis();
        let y_axis: Line<f64> = Line::y_axis();
        assert!(abs_diff_eq!(
            x_axis.distance(&y_axis),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_distance_parallel() {
        // Z axis and a parallel line at (3, 4, 0)
        let z_axis: Line<f64> = Line::z_axis();
        let parallel = Line::from_point_and_direction(
            &Point::new(3.0, 4.0, 0.0),
            &EuclideanVector::new(0.0, 0.0, 1.0),
        );
        // Distance should be 5 (sqrt(3² + 4²))
        assert!(abs_diff_eq!(
            z_axis.distance(&parallel),
            5.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_distance_skew() {
        // Z axis and a line in Y direction at (1, 0, 0)
        let z_axis: Line<f64> = Line::z_axis();
        let skew = Line::from_point_and_direction(
            &Point::new(1.0, 0.0, 0.0),
            &EuclideanVector::new(0.0, 1.0, 0.0),
        );
        // Distance should be 1 (perpendicular distance between skew lines)
        assert!(abs_diff_eq!(
            z_axis.distance(&skew),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_distance_symmetric() {
        // Distance should be symmetric
        let z_axis: Line<f64> = Line::z_axis();
        let skew = Line::from_point_and_direction(
            &Point::new(2.0, 0.0, 0.0),
            &EuclideanVector::new(0.0, 1.0, 0.0),
        );
        assert!(abs_diff_eq!(
            z_axis.distance(&skew),
            skew.distance(&z_axis),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_distance_skew_diagonal() {
        // Two skew lines at 45 degrees
        // Line 1: through origin, direction (1, 0, 0)
        let line1: Line<f64> = Line::x_axis();
        // Line 2: through (0, 0, 1), direction (0, 1, 0)
        let line2 = Line::from_point_and_direction(
            &Point::new(0.0, 0.0, 1.0),
            &EuclideanVector::new(0.0, 1.0, 0.0),
        );
        // Distance between these skew lines is 1 (the z separation)
        assert!(abs_diff_eq!(
            line1.distance(&line2),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // ========================================================================
    // Line::join_point tests
    // ========================================================================

    #[test]
    fn line_join_point_z_axis() {
        // Z axis joined with point (1, 0, 0) gives XZ plane (y = 0)
        let z_axis: Line<f64> = Line::z_axis();
        let p = Point::new(1.0, 0.0, 0.0);
        let plane = z_axis.join_point(&p);

        // Normal should be in Y direction
        let n = plane.normal();
        assert!(abs_diff_eq!(n.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(n.y.abs() > 0.9);
        assert!(abs_diff_eq!(n.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_join_point_x_axis() {
        // X axis joined with point (0, 0, 1) gives XZ plane (y = 0)
        let x_axis: Line<f64> = Line::x_axis();
        let p = Point::new(0.0, 0.0, 1.0);
        let plane = x_axis.join_point(&p);

        // Normal should be in Y direction
        let n = plane.normal();
        assert!(abs_diff_eq!(n.x, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(n.y.abs() > 0.9);
        assert!(abs_diff_eq!(n.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_join_point_contains_line() {
        // The resulting plane should contain the original line
        let z_axis: Line<f64> = Line::z_axis();
        let p = Point::new(1.0, 0.0, 0.0);
        let plane = z_axis.join_point(&p);

        // Points on Z axis should lie in the plane
        let p_on_line = Point::new(0.0, 0.0, 5.0);
        assert!(plane.contains_point(&p_on_line, ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_join_point_contains_point() {
        // The resulting plane should contain the given point
        let z_axis: Line<f64> = Line::z_axis();
        let p = Point::new(1.0, 2.0, 0.0);
        let plane = z_axis.join_point(&p);

        // The point itself should lie in the plane
        assert!(plane.contains_point(&p, ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_join_point_on_line_degenerate() {
        // Joining a line with a point on the line gives degenerate plane
        let z_axis: Line<f64> = Line::z_axis();
        let p = Point::new(0.0, 0.0, 5.0); // On the line
        let plane = z_axis.join_point(&p);

        // Normal should be degenerate (zero or very small)
        let n = plane.normal();
        let norm = (n.x * n.x + n.y * n.y + n.z * n.z).sqrt();
        assert!(norm < ABS_DIFF_EQ_EPS);
    }

    // ========================================================================
    // Plane projection tests
    // ========================================================================

    #[test]
    fn plane_project_point_xy() {
        // XY plane (z = 0)
        let plane: Plane<f64> = Plane::xy();
        let p = Point::new(1.0, 2.0, 5.0);
        let projected = plane.project_point(&p);

        assert!(abs_diff_eq!(projected.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(projected.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(projected.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_project_point_xz() {
        // XZ plane (y = 0)
        let plane: Plane<f64> = Plane::xz();
        let p = Point::new(1.0, 5.0, 3.0);
        let projected = plane.project_point(&p);

        assert!(abs_diff_eq!(projected.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(projected.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(projected.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_project_point_offset() {
        // Plane z = 2
        let plane = Plane::from_normal_and_distance(0.0, 0.0, 1.0, -2.0);
        let p = Point::new(1.0, 2.0, 5.0);
        let projected = plane.project_point(&p);

        assert!(abs_diff_eq!(projected.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(projected.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(projected.z(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_project_point_on_plane() {
        // Point already on plane should not move
        let plane: Plane<f64> = Plane::xy();
        let p = Point::new(3.0, 4.0, 0.0);
        let projected = plane.project_point(&p);

        assert!(abs_diff_eq!(
            projected.x(),
            p.x(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            projected.y(),
            p.y(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            projected.z(),
            p.z(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_project_line_parallel() {
        // Line parallel to XY plane should project to itself (in XY)
        let plane: Plane<f64> = Plane::xy();
        let line = Line::from_point_and_direction(
            &Point::new(0.0, 0.0, 5.0),
            &EuclideanVector::new(1.0, 0.0, 0.0),
        );
        let projected = plane.project_line(&line);

        // Direction should stay along X axis
        let d = projected.direction();
        assert!(abs_diff_eq!(d.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(d.x.abs() > 0.9);
    }

    #[test]
    fn plane_project_line_diagonal() {
        // Diagonal line through origin
        let plane: Plane<f64> = Plane::xy();
        let line =
            Line::from_point_and_direction(&Point::origin(), &EuclideanVector::new(1.0, 1.0, 1.0));
        let projected = plane.project_line(&line);

        // Projected direction should be (1, 1, 0)
        let d = projected.direction();
        assert!(abs_diff_eq!(d.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_project_line_in_plane() {
        // Line already in plane should not change
        let plane: Plane<f64> = Plane::xy();
        let line: Line<f64> = Line::x_axis();
        let projected = plane.project_line(&line);

        let d = projected.direction();
        assert!(abs_diff_eq!(d.y, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(d.z, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(d.x.abs() > 0.9);
    }

    // ========================================================================
    // Dot product tests
    // ========================================================================

    #[test]
    fn point_dot_same() {
        let p: Point<f64> = Point::new(1.0, 0.0, 0.0);
        assert!(abs_diff_eq!(p.dot(&p), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn point_dot_perpendicular() {
        let p1: Point<f64> = Point::new(1.0, 0.0, 0.0);
        let p2: Point<f64> = Point::new(0.0, 1.0, 0.0);
        assert!(abs_diff_eq!(p1.dot(&p2), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn point_dot_symmetric() {
        let p1 = Point::new(1.0, 2.0, 3.0);
        let p2 = Point::new(4.0, 5.0, 6.0);
        assert!(abs_diff_eq!(
            p1.dot(&p2),
            p2.dot(&p1),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_dot_intersecting() {
        // Intersecting lines have dot = 0
        let x_axis: Line<f64> = Line::x_axis();
        let y_axis: Line<f64> = Line::y_axis();
        assert!(abs_diff_eq!(
            x_axis.dot(&y_axis),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn line_dot_skew() {
        // Skew lines have non-zero dot
        let z_axis: Line<f64> = Line::z_axis();
        let skew = Line::from_point_and_direction(
            &Point::new(1.0, 0.0, 0.0),
            &EuclideanVector::new(0.0, 1.0, 0.0),
        );
        assert!(z_axis.dot(&skew).abs() > ABS_DIFF_EQ_EPS);
    }

    #[test]
    fn line_dot_symmetric() {
        let l1: Line<f64> = Line::x_axis();
        let l2 = Line::from_point_and_direction(
            &Point::new(0.0, 1.0, 0.0),
            &EuclideanVector::new(0.0, 0.0, 1.0),
        );
        assert!(abs_diff_eq!(
            l1.dot(&l2),
            l2.dot(&l1),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_dot_perpendicular() {
        let xy_plane: Plane<f64> = Plane::xy();
        let xz_plane: Plane<f64> = Plane::xz();
        assert!(abs_diff_eq!(
            xy_plane.dot(&xz_plane),
            0.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_dot_same() {
        let xy_plane: Plane<f64> = Plane::xy();
        assert!(abs_diff_eq!(
            xy_plane.dot(&xy_plane),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_dot_symmetric() {
        let p1: Plane<f64> = Plane::xy();
        let p2 = Plane::from_normal_and_distance(1.0, 1.0, 0.0, 0.0);
        assert!(abs_diff_eq!(
            p1.dot(&p2),
            p2.dot(&p1),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // ========================================================================
    // Motor commutator/anticommutator tests
    // ========================================================================

    #[test]
    fn motor_commutator_same_axis() {
        // Rotations around the same axis commute
        let r1 = Motor::from_rotation_z(std::f64::consts::FRAC_PI_4);
        let r2 = Motor::from_rotation_z(std::f64::consts::FRAC_PI_4 / 2.0);
        let comm = r1.commutator(&r2);

        // All components should be approximately zero
        assert!(abs_diff_eq!(comm.s, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(comm.e23, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(comm.e31, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(comm.e12, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_commutator_different_axes() {
        // Rotations around different axes don't commute
        let rx = Motor::from_rotation_x(std::f64::consts::FRAC_PI_4);
        let ry = Motor::from_rotation_y(std::f64::consts::FRAC_PI_4);
        let comm = rx.commutator(&ry);

        // Commutator should be non-zero
        let norm =
            comm.s * comm.s + comm.e23 * comm.e23 + comm.e31 * comm.e31 + comm.e12 * comm.e12;
        assert!(norm > ABS_DIFF_EQ_EPS);
    }

    #[test]
    fn motor_commutator_antisymmetric() {
        // [A, B] = -[B, A]
        let m1 = Motor::from_rotation_x(0.3);
        let m2 = Motor::from_rotation_y(0.5);
        let comm_ab = m1.commutator(&m2);
        let comm_ba = m2.commutator(&m1);

        assert!(abs_diff_eq!(
            comm_ab.s,
            -comm_ba.s,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            comm_ab.e23,
            -comm_ba.e23,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            comm_ab.e31,
            -comm_ba.e31,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            comm_ab.e12,
            -comm_ba.e12,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn motor_anticommutator_self() {
        // {M, M} = 2 * M²
        let m = Motor::from_rotation_z(std::f64::consts::FRAC_PI_4);
        let anti = m.anticommutator(&m);
        let m_sq = m.compose(&m);

        assert!(abs_diff_eq!(
            anti.s,
            2.0 * m_sq.s,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            anti.e12,
            2.0 * m_sq.e12,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn motor_anticommutator_symmetric() {
        // {A, B} = {B, A}
        let m1 = Motor::from_rotation_x(0.3);
        let m2 = Motor::from_rotation_y(0.5);
        let anti_ab = m1.anticommutator(&m2);
        let anti_ba = m2.anticommutator(&m1);

        assert!(abs_diff_eq!(
            anti_ab.s,
            anti_ba.s,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            anti_ab.e23,
            anti_ba.e23,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            anti_ab.e31,
            anti_ba.e31,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            anti_ab.e12,
            anti_ba.e12,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn motor_translations_commute() {
        // Pure translations commute
        let t1 = Motor::from_translation(1.0, 0.0, 0.0);
        let t2 = Motor::from_translation(0.0, 1.0, 0.0);
        let comm = t1.commutator(&t2);

        assert!(abs_diff_eq!(comm.s, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(comm.e01, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(comm.e02, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(comm.e03, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}
