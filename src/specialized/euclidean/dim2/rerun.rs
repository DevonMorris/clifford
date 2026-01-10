//! Rerun visualization support for 2D Euclidean types.
//!
//! This module provides conversions from clifford's 2D Euclidean types to
//! Rerun's visualization primitives.
//!
//! Enable with feature `rerun-0_28`.
//!
//! # Conversions
//!
//! | clifford | rerun | Notes |
//! |----------|-------|-------|
//! | [`Vector<f32>`] | [`rerun::Vec2D`] | Direction/displacement |
//! | [`AsPosition<Vector<f32>>`] | [`rerun::Position2D`] | Point in space |
//! | [`Rotor<f32>`] | [`rerun::components::RotationAxisAngle`] | Rotation around z-axis |
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::euclidean::dim2::{Vector, Rotor};
//! use clifford::specialized::visualization::AsPosition;
//! use std::f32::consts::FRAC_PI_4;
//!
//! let rec = rerun::RecordingStreamBuilder::new("clifford_demo").connect_tcp()?;
//!
//! // Log a 2D point
//! let point = Vector::new(1.0_f32, 2.0);
//! rec.log("point", &rerun::Points2D::new([AsPosition(point)]))?;
//!
//! // Log a 2D direction arrow
//! let direction = Vector::unit_x();
//! rec.log("arrow", &rerun::Arrows2D::from_vectors([direction]))?;
//! ```

use rerun_0_28 as rerun;

use super::{Rotor, Vector};
use crate::specialized::visualization::AsPosition;

// ============================================================================
// Vector -> Vec2D (direction/displacement)
// ============================================================================

impl From<Vector<f32>> for rerun::Vec2D {
    /// Converts a 2D vector to a Rerun 2D vector.
    ///
    /// Use this for directions, displacements, or velocities.
    /// For points in space, use [`AsPosition<Vector<f32>>`] instead.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Vector;
    ///
    /// let v = Vector::new(1.0_f32, 2.0);
    /// let rerun_v: rerun::Vec2D = v.into();
    /// rec.log("arrow", &rerun::Arrows2D::from_vectors([rerun_v]))?;
    /// ```
    #[inline]
    fn from(v: Vector<f32>) -> Self {
        rerun::Vec2D::new(v.x(), v.y())
    }
}

// ============================================================================
// AsPosition<Vector> -> Position2D (point in space)
// ============================================================================

impl From<AsPosition<Vector<f32>>> for rerun::Position2D {
    /// Converts a vector wrapped as a position to a Rerun 2D position.
    ///
    /// Use this for points in space (positions).
    /// For directions/displacements, convert directly to [`rerun::Vec2D`].
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Vector;
    /// use clifford::specialized::visualization::AsPosition;
    ///
    /// let point = Vector::new(1.0_f32, 2.0);
    /// rec.log("point", &rerun::Points2D::new([AsPosition(point)]))?;
    /// ```
    #[inline]
    fn from(v: AsPosition<Vector<f32>>) -> Self {
        rerun::Position2D::new(v.0.x(), v.0.y())
    }
}

// ============================================================================
// Rotor -> RotationAxisAngle (rotation around z-axis)
// ============================================================================

impl From<Rotor<f32>> for rerun::components::RotationAxisAngle {
    /// Converts a 2D rotor to a Rerun rotation (around the z-axis).
    ///
    /// Since 2D rotations are around the z-axis (perpendicular to the xy-plane),
    /// we express this as an axis-angle rotation with axis (0, 0, 1).
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use std::f32::consts::FRAC_PI_2;
    ///
    /// let rotor = Rotor::from_angle(FRAC_PI_2);
    /// let rot: rerun::components::RotationAxisAngle = rotor.into();
    /// ```
    #[inline]
    fn from(rotor: Rotor<f32>) -> Self {
        let r = rotor.normalized();
        let angle = r.angle();
        // Z-axis rotation
        rerun::components::RotationAxisAngle::new(rerun::Vec3D::new(0.0, 0.0, 1.0), angle)
    }
}

// ============================================================================
// Rotor -> Transform3D (2D rotation embedded in 3D)
// ============================================================================

impl From<Rotor<f32>> for rerun::Transform3D {
    /// Converts a 2D rotor to a Rerun 3D transform (rotation around z-axis).
    ///
    /// This embeds the 2D rotation in 3D space as a rotation around the z-axis.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use std::f32::consts::FRAC_PI_4;
    ///
    /// let rotor = Rotor::from_angle(FRAC_PI_4);
    /// rec.log("rotation_2d", &rerun::Transform3D::from(rotor))?;
    /// ```
    #[inline]
    fn from(rotor: Rotor<f32>) -> Self {
        let rotation: rerun::components::RotationAxisAngle = rotor.into();
        rerun::Transform3D::from_rotation(rotation)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::abs_diff_eq;
    use proptest::prelude::*;
    use std::f32::consts::{FRAC_PI_2, PI};

    use crate::specialized::euclidean::dim2::arbitrary::UnitRotor;

    /// Epsilon for f32 comparisons.
    const EPS: f32 = 1e-5;

    proptest! {
        #[test]
        fn vector_to_vec2d_preserves_components(
            x in -100.0f32..100.0,
            y in -100.0f32..100.0,
        ) {
            let v = Vector::new(x, y);
            let rerun_v: rerun::Vec2D = v.into();
            prop_assert!(abs_diff_eq!(rerun_v.x(), x, epsilon = EPS));
            prop_assert!(abs_diff_eq!(rerun_v.y(), y, epsilon = EPS));
        }

        #[test]
        fn as_position_to_position2d_preserves_components(
            x in -100.0f32..100.0,
            y in -100.0f32..100.0,
        ) {
            let v = Vector::new(x, y);
            let pos: rerun::Position2D = AsPosition(v).into();
            prop_assert!(abs_diff_eq!(pos.x(), x, epsilon = EPS));
            prop_assert!(abs_diff_eq!(pos.y(), y, epsilon = EPS));
        }

        #[test]
        fn rotor_to_axis_angle_preserves_angle(r in any::<UnitRotor<f32>>()) {
            let original_angle = r.angle();
            let rot: rerun::components::RotationAxisAngle = (*r).into();

            // The angle should be preserved (possibly with 2π periodicity)
            let rerun_angle = rot.angle.radians();
            let diff = (original_angle - rerun_angle).abs();
            let diff_mod = diff.min((diff - 2.0 * PI).abs());
            prop_assert!(diff_mod < EPS, "Angle mismatch: {} vs {}", original_angle, rerun_angle);
        }

        #[test]
        fn rotor_to_transform3d_does_not_panic(r in any::<UnitRotor<f32>>()) {
            let _t: rerun::Transform3D = (*r).into();
        }
    }

    #[test]
    fn rotor_90_degrees() {
        let rotor = Rotor::from_angle(FRAC_PI_2);
        let rot: rerun::components::RotationAxisAngle = rotor.into();

        // Should be 90° around z-axis
        assert!(abs_diff_eq!(rot.angle.radians(), FRAC_PI_2, epsilon = EPS));
        assert!(abs_diff_eq!(rot.axis.x(), 0.0, epsilon = EPS));
        assert!(abs_diff_eq!(rot.axis.y(), 0.0, epsilon = EPS));
        assert!(abs_diff_eq!(rot.axis.z(), 1.0, epsilon = EPS));
    }
}
