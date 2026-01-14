//! Integration tests for nalgebra interoperability.
//!
//! These tests verify that conversions between clifford and nalgebra types
//! preserve mathematical properties and are consistent across both libraries.

#![cfg(any(feature = "nalgebra-0_33", feature = "nalgebra-0_34"))]

#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use approx::abs_diff_eq;
use proptest::prelude::*;

use clifford::specialized::euclidean::{dim2, dim3};

/// Epsilon for approximate equality in tests.
const ABS_DIFF_EQ_EPS: f64 = 1e-10;

mod dim2_tests {
    use super::*;

    proptest! {
        /// Vector roundtrip: clifford -> nalgebra -> clifford
        #[test]
        fn vector_roundtrip(x in -100.0..100.0, y in -100.0..100.0) {
            let v = dim2::Vector::new(x, y);
            let na_v: na::Vector2<f64> = v.into();
            let back: dim2::Vector<f64> = na_v.into();

            prop_assert!(abs_diff_eq!(v.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(v.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Rotor roundtrip preserves rotation angle
        #[test]
        fn rotor_angle_preserved(angle in -std::f64::consts::PI..std::f64::consts::PI) {
            let rotor = dim2::Rotor::from_angle(angle);
            let rotation: na::Rotation2<f64> = rotor.into();
            let back: dim2::Rotor<f64> = rotation.into();

            prop_assert!(abs_diff_eq!(rotor.angle(), back.angle(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Rotation equivalence: clifford and nalgebra rotations produce same result
        #[test]
                fn rotation_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            vx in -100.0..100.0,
            vy in -100.0..100.0,
        ) {
            let v = dim2::Vector::new(vx, vy);
            let rotor = dim2::Rotor::from_angle(angle);

            // Rotate with clifford
            let rotated_cliff = rotor.rotate(v);

            // Rotate with nalgebra
            let rotation: na::Rotation2<f64> = rotor.into();
            let na_v: na::Vector2<f64> = v.into();
            let rotated_na = rotation * na_v;

            // Convert back and compare
            let rotated_back: dim2::Vector<f64> = rotated_na.into();

            prop_assert!(abs_diff_eq!(rotated_cliff.x(), rotated_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotated_cliff.y(), rotated_back.y(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}

mod dim3_tests {
    use super::*;

    proptest! {
        /// Vector roundtrip: clifford -> nalgebra -> clifford
        #[test]
        fn vector_roundtrip(x in -100.0..100.0, y in -100.0..100.0, z in -100.0..100.0) {
            let v = dim3::Vector::new(x, y, z);
            let na_v: na::Vector3<f64> = v.into();
            let back: dim3::Vector<f64> = na_v.into();

            prop_assert!(abs_diff_eq!(v.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(v.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(v.z(), back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Bivector <-> Vector3 via Hodge dual roundtrip
        #[test]
        fn bivector_dual_roundtrip(xy in -100.0..100.0, xz in -100.0..100.0, yz in -100.0..100.0) {
            let b = dim3::Bivector::new(xy, xz, yz);
            let na_v: na::Vector3<f64> = b.into();
            let back: dim3::Bivector<f64> = na_v.into();

            prop_assert!(abs_diff_eq!(b.xy(), back.xy(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(b.xz(), back.xz(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(b.yz(), back.yz(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Bivector to nalgebra vector matches Bivector::dual()
        #[test]
                fn bivector_to_vector_matches_dual(xy in -100.0..100.0, xz in -100.0..100.0, yz in -100.0..100.0) {
            let b = dim3::Bivector::new(xy, xz, yz);
            let dual = b.dual();
            let na_v: na::Vector3<f64> = b.into();

            prop_assert!(abs_diff_eq!(dual.x(), na_v.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(dual.y(), na_v.y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(dual.z(), na_v.z, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Bivector <-> Matrix3 (antisymmetric) roundtrip
        #[test]
        fn bivector_matrix_roundtrip(xy in -100.0..100.0, xz in -100.0..100.0, yz in -100.0..100.0) {
            let b = dim3::Bivector::new(xy, xz, yz);
            let m: na::Matrix3<f64> = b.into();
            let back: dim3::Bivector<f64> = m.try_into().unwrap();

            prop_assert!(abs_diff_eq!(b.xy(), back.xy(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(b.xz(), back.xz(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(b.yz(), back.yz(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Matrix from bivector is antisymmetric
        #[test]
        fn bivector_matrix_antisymmetric(xy in -100.0..100.0, xz in -100.0..100.0, yz in -100.0..100.0) {
            let b = dim3::Bivector::new(xy, xz, yz);
            let m: na::Matrix3<f64> = b.into();

            // Check M + M^T = 0
            let sum = m + m.transpose();
            for i in 0..3 {
                for j in 0..3 {
                    prop_assert!(abs_diff_eq!(sum[(i, j)], 0.0, epsilon = ABS_DIFF_EQ_EPS));
                }
            }
        }

        /// Rotor <-> UnitQuaternion rotation equivalence
        #[test]
        fn rotor_quaternion_rotation_equivalence(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            axis_x in -1.0f64..1.0,
            axis_y in -1.0f64..1.0,
            axis_z in 0.1f64..1.0, // ensure non-zero z to avoid zero axis
            vx in -100.0f64..100.0,
            vy in -100.0f64..100.0,
            vz in -100.0f64..100.0,
        ) {
            // Normalize axis
            let axis_norm = (axis_x * axis_x + axis_y * axis_y + axis_z * axis_z).sqrt();
            let na_axis = na::Unit::new_normalize(na::Vector3::new(
                axis_x / axis_norm,
                axis_y / axis_norm,
                axis_z / axis_norm,
            ));

            // Create rotor from quaternion
            let q = na::UnitQuaternion::from_axis_angle(&na_axis, angle);
            let rotor: dim3::Rotor<f64> = q.into();

            let v = dim3::Vector::new(vx, vy, vz);

            // Rotate with clifford rotor
            let rotated_cliff = rotor.rotate(v);

            // Rotate with nalgebra quaternion
            let q_back: na::UnitQuaternion<f64> = rotor.into();
            let na_v: na::Vector3<f64> = v.into();
            let rotated_na = q_back * na_v;

            // Convert back and compare
            let rotated_back: dim3::Vector<f64> = rotated_na.into();

            prop_assert!(abs_diff_eq!(rotated_cliff.x(), rotated_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotated_cliff.y(), rotated_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotated_cliff.z(), rotated_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Rotor <-> Rotation3 roundtrip via rotation equivalence
        #[test]
        fn rotor_rotation3_roundtrip(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            axis_x in -1.0f64..1.0,
            axis_y in -1.0f64..1.0,
            axis_z in 0.1f64..1.0,
        ) {
            let axis_norm = (axis_x * axis_x + axis_y * axis_y + axis_z * axis_z).sqrt();
            let na_axis = na::Unit::new_normalize(na::Vector3::new(
                axis_x / axis_norm,
                axis_y / axis_norm,
                axis_z / axis_norm,
            ));

            // Create rotor from quaternion
            let q = na::UnitQuaternion::from_axis_angle(&na_axis, angle);
            let rotor: dim3::Rotor<f64> = q.into();

            let rot: na::Rotation3<f64> = rotor.into();
            let back: dim3::Rotor<f64> = rot.into();

            // Test rotation equivalence (rotors have double cover)
            let test_v = dim3::Vector::new(1.0, 2.0, 3.0);
            let rotated_orig = rotor.rotate(test_v);
            let rotated_back = back.rotate(test_v);

            prop_assert!(abs_diff_eq!(rotated_orig.x(), rotated_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotated_orig.y(), rotated_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rotated_orig.z(), rotated_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn matrix_not_antisymmetric_rejected() {
        // Non-zero diagonal should fail
        let m = na::Matrix3::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(dim3::Bivector::<f64>::try_from(m).is_err());

        // Non-antisymmetric off-diagonal should fail
        let m = na::Matrix3::new(0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(dim3::Bivector::<f64>::try_from(m).is_err());
    }

    #[test]
    fn specific_quaternion_mapping() {
        use std::f64::consts::FRAC_PI_2;

        // 90° rotation around z-axis
        let na_axis = na::Unit::new_normalize(na::Vector3::z());
        let q = na::UnitQuaternion::from_axis_angle(&na_axis, FRAC_PI_2);
        let rotor: dim3::Rotor<f64> = q.into();

        // Convert back and verify roundtrip
        let q_back: na::UnitQuaternion<f64> = rotor.into();
        let q_inner = q_back.quaternion();

        // The quaternion should represent a 90° rotation around z
        // w = cos(45°), k = sin(45°), i = j = 0
        let expected_w = (std::f64::consts::FRAC_PI_4).cos();
        let expected_k = (std::f64::consts::FRAC_PI_4).sin();
        assert!(abs_diff_eq!(
            q_inner.w,
            expected_w,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(
            q_inner.k,
            expected_k,
            epsilon = ABS_DIFF_EQ_EPS
        ));
        assert!(abs_diff_eq!(q_inner.i, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(q_inner.j, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}
