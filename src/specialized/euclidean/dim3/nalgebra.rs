//! nalgebra interoperability for 3D types.
//!
//! This module provides bidirectional conversions between clifford's 3D types
//! and nalgebra's equivalent types.
//!
//! Enable with feature `nalgebra-0_33` or `nalgebra-0_34`.
//!
//! # Conversions
//!
//! | clifford | nalgebra | Notes |
//! |----------|----------|-------|
//! | [`Vector<T>`] | [`na::Vector3<T>`] | Direct component mapping |
//! | [`Bivector<T>`] | [`na::Vector3<T>`] | Via Hodge dual |
//! | [`Bivector<T>`] | [`na::Matrix3<T>`] | Antisymmetric (skew-symmetric) matrix |
//! | [`Rotor<T>`] | [`na::UnitQuaternion<T>`] | Rotation representation |
//! | [`Rotor<T>`] | [`na::Rotation3<T>`] | Via quaternion conversion |

use core::fmt;

#[cfg(feature = "nalgebra-0_32")]
use nalgebra_0_32 as na;
#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use crate::scalar::Float;

use super::{Bivector, Rotor, Vector};

// ============================================================================
// Error type
// ============================================================================

/// Error when converting from nalgebra types to clifford types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NalgebraConversionError {
    /// Matrix is not antisymmetric (skew-symmetric).
    NotAntisymmetric,
}

impl fmt::Display for NalgebraConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotAntisymmetric => write!(f, "matrix is not antisymmetric"),
        }
    }
}

impl std::error::Error for NalgebraConversionError {}

// ============================================================================
// Vector <-> Vector3
// ============================================================================

impl<T: Float + na::Scalar> From<na::Vector3<T>> for Vector<T> {
    /// Converts a nalgebra 3D vector to a clifford 3D vector.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Vector;
    /// use nalgebra::Vector3;
    ///
    /// let na_v = Vector3::new(1.0, 2.0, 3.0);
    /// let v: Vector<f64> = na_v.into();
    /// assert_eq!(v.x, 1.0);
    /// assert_eq!(v.y, 2.0);
    /// assert_eq!(v.z, 3.0);
    /// ```
    #[inline]
    fn from(v: na::Vector3<T>) -> Self {
        Vector::new(v.x, v.y, v.z)
    }
}

impl<T: Float + na::Scalar> From<Vector<T>> for na::Vector3<T> {
    /// Converts a clifford 3D vector to a nalgebra 3D vector.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Vector;
    /// use nalgebra::Vector3;
    ///
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// let na_v: Vector3<f64> = v.into();
    /// assert_eq!(na_v.x, 1.0);
    /// assert_eq!(na_v.y, 2.0);
    /// assert_eq!(na_v.z, 3.0);
    /// ```
    #[inline]
    fn from(v: Vector<T>) -> Self {
        na::Vector3::new(v.x(), v.y(), v.z())
    }
}

// ============================================================================
// Bivector <-> Vector3 (via Hodge dual)
// ============================================================================

impl<T: Float + na::Scalar> From<Bivector<T>> for na::Vector3<T> {
    /// Converts a 3D bivector to its dual vector (Hodge star).
    ///
    /// # Mathematical Correspondence
    ///
    /// In 3D, bivectors and vectors are dual to each other via the Hodge star:
    /// - `e₂₃` (yz-plane) ↔ `e₁` (x-axis)
    /// - `e₁₃` (xz-plane) ↔ `-e₂` (negative y-axis)
    /// - `e₁₂` (xy-plane) ↔ `e₃` (z-axis)
    ///
    /// This is the same operation as [`Bivector::dual()`].
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use nalgebra::Vector3;
    ///
    /// // The xy-plane bivector is dual to the z-axis
    /// let b = Bivector::unit_xy();
    /// let v: Vector3<f64> = b.into();
    /// assert!((v.z - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(b: Bivector<T>) -> Self {
        // dual: Bivector(xy, xz, yz) -> Vector(yz, -xz, xy)
        na::Vector3::new(b.yz(), -b.xz(), b.xy())
    }
}

impl<T: Float + na::Scalar> From<na::Vector3<T>> for Bivector<T> {
    /// Converts a vector to its dual bivector (inverse Hodge star).
    ///
    /// # Mathematical Correspondence
    ///
    /// This is the inverse of the Hodge dual operation:
    /// - `e₁` (x-axis) ↔ `e₂₃` (yz-plane)
    /// - `e₂` (y-axis) ↔ `-e₁₃` (negative xz-plane)
    /// - `e₃` (z-axis) ↔ `e₁₂` (xy-plane)
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use nalgebra::Vector3;
    ///
    /// // The z-axis is dual to the xy-plane bivector
    /// let v = Vector3::new(0.0, 0.0, 1.0);
    /// let b: Bivector<f64> = v.into();
    /// assert!((b.xy - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(v: na::Vector3<T>) -> Self {
        // undual: Vector(x, y, z) -> Bivector(xy=z, xz=-y, yz=x)
        Bivector::new(v.z, -v.y, v.x)
    }
}

// ============================================================================
// Bivector <-> Matrix3 (antisymmetric matrix)
// ============================================================================

impl<T: Float + na::Scalar> From<Bivector<T>> for na::Matrix3<T> {
    /// Converts a bivector to its antisymmetric (skew-symmetric) matrix representation.
    ///
    /// # Mathematical Correspondence
    ///
    /// The bivector `B = xy·e₁₂ + xz·e₁₃ + yz·e₂₃` maps to the antisymmetric matrix:
    ///
    /// ```text
    /// [  0  -xy -xz ]
    /// [ xy   0  -yz ]
    /// [ xz  yz   0  ]
    /// ```
    ///
    /// This matrix `[B]` satisfies `[B] · v = dual(B) × v` for any vector `v`,
    /// where `dual(B)` is the Hodge dual of the bivector.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use nalgebra::Matrix3;
    ///
    /// let b = Bivector::new(1.0, 2.0, 3.0); // xy=1, xz=2, yz=3
    /// let m: Matrix3<f64> = b.into();
    ///
    /// // Verify antisymmetry: M = -M^T
    /// assert!((m[(0,1)] + m[(1,0)]).abs() < 1e-10);
    /// assert!((m[(0,2)] + m[(2,0)]).abs() < 1e-10);
    /// assert!((m[(1,2)] + m[(2,1)]).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(b: Bivector<T>) -> Self {
        na::Matrix3::new(
            T::zero(),
            -b.xy(),
            -b.xz(), // row 0
            b.xy(),
            T::zero(),
            -b.yz(), // row 1
            b.xz(),
            b.yz(),
            T::zero(), // row 2
        )
    }
}

impl<T: Float + na::Scalar> TryFrom<na::Matrix3<T>> for Bivector<T> {
    type Error = NalgebraConversionError;

    /// Extracts a bivector from an antisymmetric matrix.
    ///
    /// Returns an error if the matrix is not antisymmetric within tolerance.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use nalgebra::Matrix3;
    ///
    /// let m = Matrix3::new(
    ///     0.0, -1.0, -2.0,
    ///     1.0,  0.0, -3.0,
    ///     2.0,  3.0,  0.0,
    /// );
    /// let b: Bivector<f64> = m.try_into().unwrap();
    /// assert!((b.xy - 1.0).abs() < 1e-10);
    /// assert!((b.xz - 2.0).abs() < 1e-10);
    /// assert!((b.yz - 3.0).abs() < 1e-10);
    /// ```
    fn try_from(m: na::Matrix3<T>) -> Result<Self, Self::Error> {
        // Check antisymmetry: m[i,j] = -m[j,i]
        let tol = T::from_f64(1e-10);

        // Check diagonal is zero
        if m[(0, 0)].abs() > tol || m[(1, 1)].abs() > tol || m[(2, 2)].abs() > tol {
            return Err(NalgebraConversionError::NotAntisymmetric);
        }

        // Check off-diagonal antisymmetry
        if (m[(0, 1)] + m[(1, 0)]).abs() > tol
            || (m[(0, 2)] + m[(2, 0)]).abs() > tol
            || (m[(1, 2)] + m[(2, 1)]).abs() > tol
        {
            return Err(NalgebraConversionError::NotAntisymmetric);
        }

        // Extract bivector components from upper triangle (negated)
        // Matrix: [0, -xy, -xz; xy, 0, -yz; xz, yz, 0]
        Ok(Bivector::new(-m[(0, 1)], -m[(0, 2)], -m[(1, 2)]))
    }
}

// ============================================================================
// Rotor <-> UnitQuaternion
// ============================================================================

impl<T: Float + na::RealField> From<Rotor<T>> for na::UnitQuaternion<T> {
    /// Converts a 3D rotor to a nalgebra unit quaternion.
    ///
    /// # Mathematical Correspondence
    ///
    /// A rotor `R = s + xy·e₁₂ + xz·e₁₃ + yz·e₂₃` maps to quaternion
    /// `q = w + i·i + j·j + k·k` where:
    ///
    /// - `w = s` (scalar part)
    /// - `i = yz` (rotation in yz-plane ↔ around x-axis)
    /// - `j = -xz` (rotation in xz-plane ↔ around y-axis, note sign)
    /// - `k = xy` (rotation in xy-plane ↔ around z-axis)
    ///
    /// The sign difference for `j` arises from the handedness conventions:
    /// bivector `e₁₃` represents the xz-plane with orientation `e₁ ∧ e₃`,
    /// while quaternion `j` rotates around the positive y-axis.
    ///
    /// # Normalization
    ///
    /// The input rotor is normalized before conversion to ensure a valid
    /// unit quaternion.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::{Bivector, Rotor};
    /// use nalgebra::UnitQuaternion;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // 90° rotation around z-axis
    /// let rotor = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let q: UnitQuaternion<f64> = rotor.into();
    ///
    /// // Verify the quaternion represents the same rotation
    /// let axis = q.axis().unwrap();
    /// assert!((axis.z.abs() - 1.0).abs() < 1e-10); // rotation around z
    /// ```
    #[inline]
    fn from(rotor: Rotor<T>) -> Self {
        let r = rotor.normalize();
        // Mapping: (w, i, j, k) = (s, yz, -xz, xy)
        let q = na::Quaternion::new(r.s(), r.b().yz(), -r.b().xz(), r.b().xy());
        na::UnitQuaternion::new_normalize(q)
    }
}

impl<T: Float + na::RealField> From<na::UnitQuaternion<T>> for Rotor<T> {
    /// Converts a nalgebra unit quaternion to a 3D rotor.
    ///
    /// # Mathematical Correspondence
    ///
    /// A quaternion `q = w + i·i + j·j + k·k` maps to rotor
    /// `R = s + xy·e₁₂ + xz·e₁₃ + yz·e₂₃` where:
    ///
    /// - `s = w` (scalar part)
    /// - `xy = k` (xy-plane ↔ z-axis rotation)
    /// - `xz = -j` (xz-plane ↔ y-axis rotation, note sign)
    /// - `yz = i` (yz-plane ↔ x-axis rotation)
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::{Rotor, Vector};
    /// use nalgebra::{UnitQuaternion, Vector3};
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // Create quaternion for 90° rotation around z-axis
    /// let q = UnitQuaternion::from_axis_angle(
    ///     &nalgebra::Unit::new_normalize(Vector3::z()),
    ///     FRAC_PI_2,
    /// );
    /// let rotor: Rotor<f64> = q.into();
    ///
    /// // Apply rotation: x-axis should become y-axis
    /// let v = Vector::unit_x();
    /// let rotated = rotor.rotate(v);
    /// assert!((rotated.y - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(q: na::UnitQuaternion<T>) -> Self {
        let q = q.quaternion();
        // Inverse mapping: s=w, xy=k, xz=-j, yz=i
        // Use new_unchecked since unit quaternion guarantees unit rotor
        Rotor::new_unchecked(q.w, q.k, -q.j, q.i)
    }
}

// ============================================================================
// Rotor <-> Rotation3
// ============================================================================

impl<T: Float + na::RealField> From<Rotor<T>> for na::Rotation3<T> {
    /// Converts a 3D rotor to a nalgebra 3D rotation matrix.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::{Bivector, Rotor};
    /// use nalgebra::Rotation3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let rotor = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let rot: Rotation3<f64> = rotor.into();
    ///
    /// // Rotation matrix should rotate x to y
    /// let result = rot * nalgebra::Vector3::x();
    /// assert!((result.y - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(rotor: Rotor<T>) -> Self {
        let q: na::UnitQuaternion<T> = rotor.into();
        q.into()
    }
}

impl<T: Float + na::RealField> From<na::Rotation3<T>> for Rotor<T> {
    /// Converts a nalgebra 3D rotation matrix to a rotor.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim3::{Rotor, Vector};
    /// use nalgebra::{Rotation3, Vector3};
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let rot = Rotation3::from_axis_angle(
    ///     &nalgebra::Unit::new_normalize(Vector3::z()),
    ///     FRAC_PI_2,
    /// );
    /// let rotor: Rotor<f64> = rot.into();
    ///
    /// let v = Vector::unit_x();
    /// let rotated = rotor.rotate(v);
    /// assert!((rotated.y - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(rot: na::Rotation3<T>) -> Self {
        let q: na::UnitQuaternion<T> = rot.into();
        q.into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn vector_roundtrip(v in any::<Vector<f64>>()) {
            let na_v: na::Vector3<f64> = v.into();
            let back: Vector<f64> = na_v.into();
            prop_assert!(abs_diff_eq!(v, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivector_dual_roundtrip(b in any::<Bivector<f64>>()) {
            let na_v: na::Vector3<f64> = b.into();
            let back: Bivector<f64> = na_v.into();
            prop_assert!(abs_diff_eq!(b, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivector_matrix_roundtrip(b in any::<Bivector<f64>>()) {
            let m: na::Matrix3<f64> = b.into();
            let back: Bivector<f64> = m.try_into().unwrap();
            prop_assert!(abs_diff_eq!(b, back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivector_matrix_is_antisymmetric(b in any::<Bivector<f64>>()) {
            let m: na::Matrix3<f64> = b.into();

            // Check M + M^T = 0
            let sum = m + m.transpose();
            for i in 0..3 {
                for j in 0..3 {
                    prop_assert!(sum[(i, j)].abs() < ABS_DIFF_EQ_EPS);
                }
            }
        }

        #[test]
        fn rotor_quaternion_roundtrip(r in any::<Rotor<f64>>()) {
            // Normalize for this test
            let r = r.normalize();
            let q: na::UnitQuaternion<f64> = r.into();
            let back: Rotor<f64> = q.into();

            // Rotors have double cover: r and -r represent the same rotation
            // So we test rotation equivalence instead of component equality
            let test_v = Vector::new(1.0, 2.0, 3.0);
            let rotated_orig = r.rotate(test_v);
            let rotated_back = back.rotate(test_v);
            prop_assert!(abs_diff_eq!(rotated_orig, rotated_back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_rotation3_roundtrip(r in any::<Rotor<f64>>()) {
            let r = r.normalize();
            let rot: na::Rotation3<f64> = r.into();
            let back: Rotor<f64> = rot.into();

            let test_v = Vector::new(1.0, 2.0, 3.0);
            let rotated_orig = r.rotate(test_v);
            let rotated_back = back.rotate(test_v);
            prop_assert!(abs_diff_eq!(rotated_orig, rotated_back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_quaternion_rotation_equivalence(
            r in any::<Rotor<f64>>(),
            v in any::<Vector<f64>>(),
        ) {
            let r = r.normalize();
            let na_v: na::Vector3<f64> = v.into();

            // Rotate with clifford rotor
            let rotated_ga = r.rotate(v);

            // Rotate with nalgebra quaternion
            let q: na::UnitQuaternion<f64> = r.into();
            let rotated_na = q * na_v;

            let rotated_back: Vector<f64> = rotated_na.into();
            prop_assert!(abs_diff_eq!(rotated_ga, rotated_back, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn bivector_dual_matches_method(b in any::<Bivector<f64>>()) {
            // Conversion to nalgebra should match the dual() method
            let na_v: na::Vector3<f64> = b.into();
            let dual_v = b.dual();

            prop_assert!(abs_diff_eq!(na_v.x, dual_v.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(na_v.y, dual_v.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(na_v.z, dual_v.z(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn matrix_not_antisymmetric_error() {
        // Non-zero diagonal
        let m = na::Matrix3::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(
            Bivector::<f64>::try_from(m),
            Err(NalgebraConversionError::NotAntisymmetric)
        );

        // Not antisymmetric off-diagonal
        let m = na::Matrix3::new(0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(
            Bivector::<f64>::try_from(m),
            Err(NalgebraConversionError::NotAntisymmetric)
        );
    }
}
