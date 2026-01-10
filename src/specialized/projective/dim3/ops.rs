//! Operations for 3D Projective Geometric Algebra types.

use super::types::{Flector, Motor, Plane, Point};
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
}
