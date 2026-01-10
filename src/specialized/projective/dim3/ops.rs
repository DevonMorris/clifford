//! Operations for 3D Projective Geometric Algebra types.

use super::types::{Motor, Point};
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
}
