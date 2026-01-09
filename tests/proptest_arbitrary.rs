//! Tests verifying that `Arbitrary` implementations are accessible to crate consumers.
//!
//! These integration tests ensure that the `proptest-support` feature correctly
//! exposes all arbitrary types and implementations for external use.
//!
//! Run with: `cargo test --features proptest-support`

#![cfg(feature = "proptest-support")]

use clifford::algebra::Multivector;
use clifford::algebra::arbitrary::{NonZeroVectorE3, UnitVectorE3, VectorE3};
use clifford::prelude::abs_diff_eq;
use clifford::signature::Euclidean3;
use clifford::specialized::ga2d::Vec2;
use clifford::specialized::ga2d::arbitrary::{NonZeroVec2, UnitRotor2, UnitVec2};
use clifford::specialized::ga3d::arbitrary::{NonZeroVec3, UnitRotor3, UnitVec3};
use clifford::specialized::ga3d::{Bivec3, Vec3};
use proptest::prelude::*;
use proptest::prop_assert;

/// Standard epsilon for absolute difference comparisons in tests.
const ABS_DIFF_EQ_EPS: f64 = 1e-10;

// ============================================================================
// Generic Multivector Arbitrary accessibility
// ============================================================================

proptest! {
    #[test]
    fn multivector_arbitrary_accessible(mv in any::<Multivector<f64, Euclidean3>>()) {
        // Verify we can use the arbitrary multivector
        let _ = mv.scalar_part();
    }

    #[test]
    fn vector_e3_arbitrary_accessible(v in any::<VectorE3>()) {
        // Verify VectorE3 wrapper is accessible and usable
        let _ = v.scalar_part();
    }

    #[test]
    fn non_zero_vector_e3_arbitrary_accessible(v in any::<NonZeroVectorE3>()) {
        // Verify NonZeroVectorE3 wrapper is accessible and has non-zero norm
        prop_assert!(v.norm_squared() > 0.01);
    }

    #[test]
    fn unit_vector_e3_arbitrary_accessible(v in any::<UnitVectorE3>()) {
        // Verify UnitVectorE3 wrapper is accessible and has unit norm
        prop_assert!(abs_diff_eq!(v.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}

// ============================================================================
// GA2D Arbitrary accessibility
// ============================================================================

proptest! {
    #[test]
    fn vec2_arbitrary_accessible(v in any::<Vec2<f64>>()) {
        // Verify Vec2 arbitrary is accessible
        let _ = v.norm();
    }

    #[test]
    fn non_zero_vec2_arbitrary_accessible(v in any::<NonZeroVec2>()) {
        // Verify NonZeroVec2 wrapper is accessible and has non-zero norm
        prop_assert!(v.norm_squared() > 0.01);
    }

    #[test]
    fn unit_vec2_arbitrary_accessible(v in any::<UnitVec2>()) {
        // Verify UnitVec2 wrapper is accessible and has unit norm
        prop_assert!(abs_diff_eq!(v.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn unit_rotor2_arbitrary_accessible(r in any::<UnitRotor2>()) {
        // Verify UnitRotor2 wrapper is accessible and has unit norm
        prop_assert!(abs_diff_eq!(r.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}

// ============================================================================
// GA3D Arbitrary accessibility
// ============================================================================

proptest! {
    #[test]
    fn vec3_arbitrary_accessible(v in any::<Vec3<f64>>()) {
        // Verify Vec3 arbitrary is accessible
        let _ = v.norm();
    }

    #[test]
    fn bivec3_arbitrary_accessible(b in any::<Bivec3<f64>>()) {
        // Verify Bivec3 arbitrary is accessible
        let _ = b.norm();
    }

    #[test]
    fn non_zero_vec3_arbitrary_accessible(v in any::<NonZeroVec3>()) {
        // Verify NonZeroVec3 wrapper is accessible and has non-zero norm
        prop_assert!(v.norm_squared() > 0.01);
    }

    #[test]
    fn unit_vec3_arbitrary_accessible(v in any::<UnitVec3>()) {
        // Verify UnitVec3 wrapper is accessible and has unit norm
        prop_assert!(abs_diff_eq!(v.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn unit_rotor3_arbitrary_accessible(r in any::<UnitRotor3>()) {
        // Verify UnitRotor3 wrapper is accessible and has unit norm
        prop_assert!(abs_diff_eq!(r.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}
