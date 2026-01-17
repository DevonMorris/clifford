//! Rerun visualization support for 3D Projective GA types.
//!
//! This module provides conversions from clifford's 3D PGA types to
//! Rerun's visualization primitives.
//!
//! Enable with feature `rerun-0_28`.
//!
//! # Conversions
//!
//! | clifford | rerun | Notes |
//! |----------|-------|-------|
//! | [`Point<f32>`] | [`rerun::Position3D`] | Cartesian extraction |
//! | [`Motor<f32>`] | [`rerun::Transform3D`] | Full rigid transform |
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::projective::dim3::{Point, Motor};
//! use std::f32::consts::FRAC_PI_4;
//!
//! let rec = rerun::RecordingStreamBuilder::new("clifford_demo").connect_tcp()?;
//!
//! // Log a PGA point
//! let point = Point::from_cartesian(1.0_f32, 2.0, 3.0);
//! rec.log("point", &rerun::Points3D::new([point]))?;
//!
//! // Log a motor transform
//! let motor = Motor::from_rotation_z(FRAC_PI_4);
//! rec.log("transform", &rerun::Transform3D::from(motor))?;
//! ```

use rerun_0_28 as rerun;

use super::{Motor, Plane, Point};

// ============================================================================
// Point -> Position3D
// ============================================================================

impl From<Point<f32>> for rerun::Position3D {
    /// Converts a PGA point to a Rerun 3D position.
    ///
    /// The point's homogeneous coordinates are normalized to extract
    /// Cartesian coordinates `(x/w, y/w, z/w)`.
    ///
    /// # Panics
    ///
    /// May produce NaN if the point is ideal (w = 0).
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::Point;
    ///
    /// let p = Point::from_cartesian(1.0_f32, 2.0, 3.0);
    /// let pos: rerun::Position3D = p.into();
    /// ```
    #[inline]
    fn from(p: Point<f32>) -> Self {
        rerun::Position3D::new(p.x(), p.y(), p.z())
    }
}

// ============================================================================
// Motor -> Transform3D
// ============================================================================

impl From<Motor<f32>> for rerun::Transform3D {
    /// Converts a PGA motor to a Rerun 3D transform.
    ///
    /// The motor represents a rigid transformation (rotation + translation).
    ///
    /// # Convention
    ///
    /// The motor's rotation part maps to a quaternion:
    /// - `w = ps` (pseudoscalar)
    /// - `x = -rx` (rotation around x, negated for quaternion convention)
    /// - `y = -ry` (rotation around y, negated)
    /// - `z = -rz` (rotation around z, negated)
    ///
    /// The translation is extracted by transforming the origin point.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::Motor;
    /// use std::f32::consts::FRAC_PI_2;
    ///
    /// let motor = Motor::from_rotation_z(FRAC_PI_2);
    /// let transform: rerun::Transform3D = motor.into();
    /// ```
    #[inline]
    fn from(motor: Motor<f32>) -> Self {
        use super::Point;
        use crate::ops::Transform;

        let m = motor.unitized();

        // Extract rotation quaternion
        // Motor semantic fields: s, tz, ty, tx, rx, ry, rz, ps
        // Quaternion: (w, x, y, z) = (ps, -rx, -ry, -rz)
        // Negate bivector components to match quaternion rotation convention
        let quat = rerun::Quaternion::from_xyzw([-m.rx(), -m.ry(), -m.rz(), m.ps()]);
        let rotation = rerun::components::RotationQuat(quat);

        // Extract translation by transforming the origin
        // This correctly handles the interaction between rotation and translation
        let origin = Point::origin();
        let transformed = m.transform(&origin);
        let translation = rerun::Vec3D::new(
            transformed.cartesian_x(),
            transformed.cartesian_y(),
            transformed.cartesian_z(),
        );

        rerun::Transform3D::from_translation_rotation(translation, rotation)
    }
}

// ============================================================================
// Plane -> Vec3D (normal direction)
// ============================================================================

impl From<Plane<f32>> for rerun::Vec3D {
    /// Converts a PGA plane to a Rerun 3D vector (the normal direction).
    ///
    /// This extracts the plane's normal direction, useful for visualizing
    /// the plane's orientation as an arrow.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::Plane;
    ///
    /// let plane = Plane::xy();  // Normal is (0, 0, 1)
    /// let normal: rerun::Vec3D = plane.into();
    /// ```
    #[inline]
    fn from(plane: Plane<f32>) -> Self {
        let n = plane.normal();
        rerun::Vec3D::new(n.x(), n.y(), n.z())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    use crate::specialized::projective::dim3::BulkMotor;

    /// Epsilon for f32 comparisons.
    const EPS: f32 = 1e-4;

    proptest! {
        #[test]
        fn point_to_position_preserves_cartesian(
            x in -100.0f32..100.0,
            y in -100.0f32..100.0,
            z in -100.0f32..100.0,
        ) {
            let p = Point::from_cartesian(x, y, z);
            let pos: rerun::Position3D = p.into();
            prop_assert!(abs_diff_eq!(pos.x(), x, epsilon = EPS));
            prop_assert!(abs_diff_eq!(pos.y(), y, epsilon = EPS));
            prop_assert!(abs_diff_eq!(pos.z(), z, epsilon = EPS));
        }

        #[test]
        fn plane_to_vec3d_is_normal(
            nx in -1.0f32..1.0,
            ny in -1.0f32..1.0,
            nz in -1.0f32..1.0,
            d in -10.0f32..10.0,
        ) {
            // Skip near-zero normals
            if nx*nx + ny*ny + nz*nz < 0.01 {
                return Ok(());
            }
            let plane = Plane::from_normal_and_distance(nx, ny, nz, d);
            let v: rerun::Vec3D = plane.into();
            prop_assert!(abs_diff_eq!(v.x(), nx, epsilon = EPS));
            prop_assert!(abs_diff_eq!(v.y(), ny, epsilon = EPS));
            prop_assert!(abs_diff_eq!(v.z(), nz, epsilon = EPS));
        }

        #[test]
        fn motor_to_transform3d_does_not_panic(m in any::<BulkMotor<f32>>()) {
            let _t: rerun::Transform3D = (*m).into();
        }
    }

    #[test]
    fn motor_pure_translation_converts() {
        let motor = Motor::from_translation(1.0_f32, 2.0, 3.0);
        // Just verify conversion works - Transform3D internals are opaque in SDK mode
        let _t: rerun::Transform3D = motor.into();
    }

    #[test]
    fn motor_pure_rotation_converts() {
        let motor = Motor::from_rotation_z(std::f32::consts::FRAC_PI_2);
        // Just verify conversion works - Transform3D internals are opaque in SDK mode
        let _t: rerun::Transform3D = motor.into();
    }
}
