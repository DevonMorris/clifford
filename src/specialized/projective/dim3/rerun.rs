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
//! let point = Point::new(1.0_f32, 2.0, 3.0);
//! rec.log("point", &rerun::Points3D::new([point]))?;
//!
//! // Log a motor transform
//! let motor = Motor::from_rotation_z(FRAC_PI_4)
//!     .compose(&Motor::from_translation(1.0, 0.0, 0.0));
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
    /// let p = Point::new(1.0_f32, 2.0, 3.0);
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
    /// The motor's rotation part `(s, e23, e31, e12)` maps to a quaternion:
    /// - `w = s` (scalar)
    /// - `x = e23` (rotation around x)
    /// - `y = e31` (rotation around y)
    /// - `z = e12` (rotation around z)
    ///
    /// The translation is extracted from the motor's translation components.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim3::Motor;
    /// use std::f32::consts::FRAC_PI_2;
    ///
    /// let motor = Motor::from_rotation_z(FRAC_PI_2)
    ///     .compose(&Motor::from_translation(1.0, 2.0, 3.0));
    /// let transform: rerun::Transform3D = motor.into();
    /// ```
    #[inline]
    fn from(motor: Motor<f32>) -> Self {
        let m = motor.unitized();

        // Extract rotation quaternion from the rotor part
        // Motor components: s, e23, e31, e12
        // Quaternion: w = s, x = e23, y = e31, z = e12
        let quat = rerun::Quaternion::from_xyzw([m.e23(), m.e31(), m.e12(), m.s()]);
        let rotation = rerun::components::RotationQuat(quat);

        // Extract translation
        // For a motor M = R + d, the translation is: t = 2 * (e01, e02, e03)
        // But for a composed motor, we need to compute the actual translation
        // using: t = 2 * (s*d - b×d) where s is scalar, b is bivector, d is null bivector
        //
        // Simplified for PGA: t = 2 * (s*e01 + e12*e02 - e31*e03, ...)
        // For a unit motor where the translation was encoded with 1/2 factor:
        let tx = 2.0 * (m.s() * m.e01() + m.e12() * m.e02() - m.e31() * m.e03());
        let ty = 2.0 * (m.s() * m.e02() + m.e23() * m.e03() - m.e12() * m.e01());
        let tz = 2.0 * (m.s() * m.e03() + m.e31() * m.e01() - m.e23() * m.e02());

        let translation = rerun::Vec3D::new(tx, ty, tz);

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

    use crate::specialized::projective::dim3::arbitrary::UnitMotor;

    /// Epsilon for f32 comparisons.
    const EPS: f32 = 1e-4;

    proptest! {
        #[test]
        fn point_to_position_preserves_cartesian(
            x in -100.0f32..100.0,
            y in -100.0f32..100.0,
            z in -100.0f32..100.0,
        ) {
            let p = Point::new(x, y, z);
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
        fn motor_to_transform3d_does_not_panic(m in any::<UnitMotor<f32>>()) {
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
