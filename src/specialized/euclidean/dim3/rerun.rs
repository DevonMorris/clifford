//! Rerun visualization support for 3D Euclidean types.
//!
//! This module provides conversions from clifford's 3D Euclidean types to
//! Rerun's visualization primitives.
//!
//! Enable with feature `rerun-0_28`.
//!
//! # Conversions
//!
//! | clifford | rerun | Notes |
//! |----------|-------|-------|
//! | [`Vector<f32>`] | [`rerun::Vec3D`] | Direction/displacement |
//! | [`AsPosition<Vector<f32>>`] | [`rerun::Position3D`] | Point in space |
//! | [`Bivector<f32>`] | [`rerun::Vec3D`] | Via Hodge dual (rotation axis) |
//! | [`Rotor<f32>`] | [`rerun::components::RotationQuat`] | Quaternion representation |
//! | [`Rotor<f32>`] | [`rerun::Transform3D`] | Pure rotation transform |
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::euclidean::dim3::{Vector, Bivector, Rotor};
//! use clifford::specialized::visualization::AsPosition;
//! use std::f32::consts::FRAC_PI_4;
//!
//! let rec = rerun::RecordingStreamBuilder::new("clifford_demo").connect_tcp()?;
//!
//! // Log a point
//! let point = Vector::new(1.0_f32, 2.0, 3.0);
//! rec.log("point", &rerun::Points3D::new([AsPosition(point)]))?;
//!
//! // Log a direction arrow
//! let direction = Vector::unit_x();
//! rec.log("arrow", &rerun::Arrows3D::from_vectors([direction]))?;
//!
//! // Log a rotation
//! let rotor = Rotor::from_angle_plane(FRAC_PI_4, Bivector::unit_xy());
//! rec.log("rotation", &rerun::Transform3D::from(rotor))?;
//! ```

use rerun_0_28 as rerun;

use super::{Bivector, Rotor, Vector};
use crate::specialized::visualization::AsPosition;

// ============================================================================
// Vector -> Vec3D (direction/displacement)
// ============================================================================

impl From<Vector<f32>> for rerun::Vec3D {
    /// Converts a 3D vector to a Rerun 3D vector.
    ///
    /// Use this for directions, displacements, or velocities.
    /// For points in space, use [`AsPosition<Vector<f32>>`] instead.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let v = Vector::new(1.0_f32, 2.0, 3.0);
    /// let rerun_v: rerun::Vec3D = v.into();
    /// rec.log("arrow", &rerun::Arrows3D::from_vectors([rerun_v]))?;
    /// ```
    #[inline]
    fn from(v: Vector<f32>) -> Self {
        rerun::Vec3D::new(v.x(), v.y(), v.z())
    }
}

// ============================================================================
// AsPosition<Vector> -> Position3D (point in space)
// ============================================================================

impl From<AsPosition<Vector<f32>>> for rerun::Position3D {
    /// Converts a vector wrapped as a position to a Rerun 3D position.
    ///
    /// Use this for points in space (positions).
    /// For directions/displacements, convert directly to [`rerun::Vec3D`].
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Vector;
    /// use clifford::specialized::visualization::AsPosition;
    ///
    /// let point = Vector::new(1.0_f32, 2.0, 3.0);
    /// rec.log("point", &rerun::Points3D::new([AsPosition(point)]))?;
    /// ```
    #[inline]
    fn from(v: AsPosition<Vector<f32>>) -> Self {
        rerun::Position3D::new(v.0.x(), v.0.y(), v.0.z())
    }
}

// ============================================================================
// Bivector -> Vec3D (via Hodge dual)
// ============================================================================

impl From<Bivector<f32>> for rerun::Vec3D {
    /// Converts a bivector to a Rerun 3D vector via the Hodge dual.
    ///
    /// # Mathematical Correspondence
    ///
    /// In 3D, bivectors and vectors are dual to each other:
    /// - `e₂₃` (yz-plane) ↔ `e₁` (x-axis)
    /// - `e₁₃` (xz-plane) ↔ `-e₂` (negative y-axis)
    /// - `e₁₂` (xy-plane) ↔ `e₃` (z-axis)
    ///
    /// This is useful for visualizing rotation planes as their corresponding
    /// rotation axes.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Bivector;
    ///
    /// // The xy-plane bivector becomes the z-axis
    /// let b = Bivector::unit_xy();
    /// let axis: rerun::Vec3D = b.into();
    /// // axis ≈ (0, 0, 1)
    /// ```
    #[inline]
    fn from(b: Bivector<f32>) -> Self {
        // dual: Bivector(xy, xz, yz) -> Vector(yz, -xz, xy)
        rerun::Vec3D::new(b.yz(), -b.xz(), b.xy())
    }
}

// ============================================================================
// Rotor -> RotationQuat
// ============================================================================

impl From<Rotor<f32>> for rerun::components::RotationQuat {
    /// Converts a 3D rotor to a Rerun rotation quaternion.
    ///
    /// # Convention
    ///
    /// Clifford rotor `R = s + xy·e₁₂ + xz·e₁₃ + yz·e₂₃` maps to
    /// quaternion `(x, y, z, w) = (yz, -xz, xy, s)`.
    ///
    /// This follows the standard correspondence between bivector planes
    /// and quaternion imaginary units:
    /// - `e₂₃` ↔ `i` (rotation in yz-plane, around x-axis)
    /// - `e₁₃` ↔ `-j` (rotation in xz-plane, around y-axis, note sign)
    /// - `e₁₂` ↔ `k` (rotation in xy-plane, around z-axis)
    ///
    /// # Normalization
    ///
    /// The rotor is normalized before conversion to ensure a valid
    /// unit quaternion.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::{Rotor, Bivector};
    /// use std::f32::consts::FRAC_PI_2;
    ///
    /// let rotor = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let quat: rerun::components::RotationQuat = rotor.into();
    /// ```
    #[inline]
    fn from(rotor: Rotor<f32>) -> Self {
        let r = rotor.normalized();
        // Mapping: (x, y, z, w) = (yz, -xz, xy, s)
        let quat = rerun::Quaternion::from_xyzw([r.b().yz(), -r.b().xz(), r.b().xy(), r.s()]);
        rerun::components::RotationQuat(quat)
    }
}

// ============================================================================
// Rotor -> Transform3D (pure rotation)
// ============================================================================

impl From<Rotor<f32>> for rerun::Transform3D {
    /// Converts a 3D rotor to a Rerun 3D transform (pure rotation).
    ///
    /// The resulting transform contains only the rotation component;
    /// translation and scale are identity.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::{Rotor, Bivector};
    /// use std::f32::consts::FRAC_PI_4;
    ///
    /// let rotor = Rotor::from_angle_plane(FRAC_PI_4, Bivector::unit_xy());
    /// rec.log("rotation", &rerun::Transform3D::from(rotor))?;
    /// ```
    #[inline]
    fn from(rotor: Rotor<f32>) -> Self {
        let quat: rerun::components::RotationQuat = rotor.into();
        rerun::Transform3D::from_rotation(quat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    use crate::specialized::euclidean::dim3::arbitrary::UnitRotor;

    /// Epsilon for f32 comparisons.
    const EPS: f32 = 1e-5;

    proptest! {
        #[test]
        fn vector_to_vec3d_preserves_components(
            x in -100.0f32..100.0,
            y in -100.0f32..100.0,
            z in -100.0f32..100.0,
        ) {
            let v = Vector::new(x, y, z);
            let rerun_v: rerun::Vec3D = v.into();
            prop_assert!(abs_diff_eq!(rerun_v.x(), x, epsilon = EPS));
            prop_assert!(abs_diff_eq!(rerun_v.y(), y, epsilon = EPS));
            prop_assert!(abs_diff_eq!(rerun_v.z(), z, epsilon = EPS));
        }

        #[test]
        fn as_position_to_position3d_preserves_components(
            x in -100.0f32..100.0,
            y in -100.0f32..100.0,
            z in -100.0f32..100.0,
        ) {
            let v = Vector::new(x, y, z);
            let pos: rerun::Position3D = AsPosition(v).into();
            prop_assert!(abs_diff_eq!(pos.x(), x, epsilon = EPS));
            prop_assert!(abs_diff_eq!(pos.y(), y, epsilon = EPS));
            prop_assert!(abs_diff_eq!(pos.z(), z, epsilon = EPS));
        }

        #[test]
        fn bivector_to_vec3d_matches_dual(
            xy in -100.0f32..100.0,
            xz in -100.0f32..100.0,
            yz in -100.0f32..100.0,
        ) {
            let b = Bivector::new(xy, xz, yz);
            let rerun_v: rerun::Vec3D = b.into();
            let dual = b.dual();

            prop_assert!(abs_diff_eq!(rerun_v.x(), dual.x(), epsilon = EPS));
            prop_assert!(abs_diff_eq!(rerun_v.y(), dual.y(), epsilon = EPS));
            prop_assert!(abs_diff_eq!(rerun_v.z(), dual.z(), epsilon = EPS));
        }

        #[test]
        fn rotor_to_quaternion_does_not_panic(r in any::<UnitRotor<f32>>()) {
            let _quat: rerun::components::RotationQuat = (*r).into();
        }

        #[test]
        fn rotor_to_transform3d_does_not_panic(r in any::<UnitRotor<f32>>()) {
            let _t: rerun::Transform3D = (*r).into();
        }
    }

    #[test]
    fn rotor_z_rotation_converts() {
        // 90° rotation around z-axis (in xy-plane)
        let rotor = Rotor::from_angle_plane(std::f32::consts::FRAC_PI_2, Bivector::unit_xy());
        // Just verify conversion works - Quaternion internals are opaque in SDK mode
        let _quat: rerun::components::RotationQuat = rotor.into();
    }
}
