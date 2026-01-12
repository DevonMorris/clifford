//! Rerun visualization support for 2D Projective GA types.
//!
//! This module provides conversions from clifford's 2D PGA types to
//! Rerun's visualization primitives.
//!
//! Enable with feature `rerun-0_28`.
//!
//! # Conversions
//!
//! | clifford | rerun | Notes |
//! |----------|-------|-------|
//! | [`Point<f32>`] | [`rerun::Position2D`] | Cartesian extraction |
//! | [`Line<f32>`] | [`rerun::Vec2D`] | Normal direction |
//! | [`Motor<f32>`] | [`rerun::Transform3D`] | 2D transform embedded in 3D |
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::projective::dim2::{Point, Motor};
//! use std::f32::consts::FRAC_PI_4;
//!
//! let rec = rerun::RecordingStreamBuilder::new("clifford_demo").connect_tcp()?;
//!
//! // Log a PGA point
//! let point = Point::new(1.0_f32, 2.0);
//! rec.log("point", &rerun::Points2D::new([point]))?;
//!
//! // Log a motor transform
//! let motor = Motor::from_rotation(FRAC_PI_4)
//!     .compose(&Motor::from_translation(1.0, 0.0));
//! rec.log("transform", &rerun::Transform3D::from(motor))?;
//! ```

use rerun_0_28 as rerun;

use super::{Line, Motor, Point};

// ============================================================================
// Point -> Position2D
// ============================================================================

impl From<Point<f32>> for rerun::Position2D {
    /// Converts a PGA point to a Rerun 2D position.
    ///
    /// The point's homogeneous coordinates are normalized to extract
    /// Cartesian coordinates `(x/w, y/w)`.
    ///
    /// # Panics
    ///
    /// May produce NaN if the point is ideal (w = 0).
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim2::Point;
    ///
    /// let p = Point::new(1.0_f32, 2.0);
    /// let pos: rerun::Position2D = p.into();
    /// ```
    #[inline]
    fn from(p: Point<f32>) -> Self {
        rerun::Position2D::new(p.x(), p.y())
    }
}

// ============================================================================
// Line -> Vec2D (normal direction)
// ============================================================================

impl From<Line<f32>> for rerun::Vec2D {
    /// Converts a PGA line to a Rerun 2D vector (the normal direction).
    ///
    /// This extracts the line's normal direction, useful for visualizing
    /// the line's orientation as an arrow.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim2::Line;
    ///
    /// let line = Line::x_axis();  // Normal is (0, 1)
    /// let normal: rerun::Vec2D = line.into();
    /// ```
    #[inline]
    fn from(line: Line<f32>) -> Self {
        let n = line.normal();
        rerun::Vec2D::new(n.x(), n.y())
    }
}

// ============================================================================
// Motor -> Transform3D (2D transform embedded in 3D)
// ============================================================================

impl From<Motor<f32>> for rerun::Transform3D {
    /// Converts a 2D PGA motor to a Rerun 3D transform.
    ///
    /// The 2D rigid transformation is embedded in 3D space:
    /// - Rotation around the z-axis
    /// - Translation in the xy-plane
    ///
    /// # Convention
    ///
    /// The motor's rotation `e₁₂` maps to z-axis rotation (in the xy-plane).
    /// The translation is extracted from the motor's bivector components.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::projective::dim2::Motor;
    /// use std::f32::consts::FRAC_PI_2;
    ///
    /// let motor = Motor::from_rotation(FRAC_PI_2)
    ///     .compose(&Motor::from_translation(1.0, 2.0));
    /// let transform: rerun::Transform3D = motor.into();
    /// ```
    #[inline]
    fn from(motor: Motor<f32>) -> Self {
        let m = motor.normalized();

        // Extract rotation as z-axis rotation (quaternion)
        // For 2D rotation: q = (0, 0, sin(θ/2), cos(θ/2))
        // Motor has s = cos(θ/2), e12 = sin(θ/2)
        let quat = rerun::Quaternion::from_xyzw([0.0, 0.0, m.e12(), m.s()]);
        let rotation = rerun::components::RotationQuat(quat);

        // Extract translation (similar to 3D motor)
        // For a 2D motor M = s + e12 + e20 + e01:
        // Translation t = 2 * (s*e20 - e12*e01, s*e01 + e12*e20)
        let tx = 2.0 * (m.s() * m.e20() - m.e12() * m.e01());
        let ty = 2.0 * (m.s() * m.e01() + m.e12() * m.e20());

        let translation = rerun::Vec3D::new(tx, ty, 0.0);

        rerun::Transform3D::from_translation_rotation(translation, rotation)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    /// Epsilon for f32 comparisons.
    const EPS: f32 = 1e-4;

    proptest! {
        #[test]
        fn point_to_position_preserves_cartesian(
            x in -100.0f32..100.0,
            y in -100.0f32..100.0,
        ) {
            let p = Point::new(x, y);
            let pos: rerun::Position2D = p.into();
            prop_assert!(abs_diff_eq!(pos.x(), x, epsilon = EPS));
            prop_assert!(abs_diff_eq!(pos.y(), y, epsilon = EPS));
        }

        #[test]
        fn line_to_vec2d_is_normal(
            nx in -1.0f32..1.0,
            ny in -1.0f32..1.0,
            d in -10.0f32..10.0,
        ) {
            // Skip near-zero normals
            if nx*nx + ny*ny < 0.01 {
                return Ok(());
            }
            let line = Line::from_implicit(nx, ny, d);
            let v: rerun::Vec2D = line.into();
            prop_assert!(abs_diff_eq!(v.x(), nx, epsilon = EPS));
            prop_assert!(abs_diff_eq!(v.y(), ny, epsilon = EPS));
        }

        #[test]
        fn motor_to_transform3d_does_not_panic(
            angle in -std::f32::consts::PI..std::f32::consts::PI,
            tx in -100.0f32..100.0,
            ty in -100.0f32..100.0,
        ) {
            let rotation = Motor::from_rotation(angle);
            let translation = Motor::from_translation(tx, ty);
            let motor = translation.compose(&rotation);
            let _t: rerun::Transform3D = motor.into();
        }
    }

    #[test]
    fn motor_pure_translation_converts() {
        let motor = Motor::from_translation(1.0_f32, 2.0);
        // Just verify conversion works - Transform3D internals are opaque in SDK mode
        let _t: rerun::Transform3D = motor.into();
    }

    #[test]
    fn motor_pure_rotation_converts() {
        let motor = Motor::from_rotation(std::f32::consts::FRAC_PI_2);
        // Just verify conversion works - Transform3D internals are opaque in SDK mode
        let _t: rerun::Transform3D = motor.into();
    }
}
