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
use clifford::specialized::euclidean::{dim2, dim3};
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
// dim2 Arbitrary accessibility
// ============================================================================

proptest! {
    #[test]
    fn vec2_arbitrary_accessible(v in any::<dim2::Vector<f64>>()) {
        // Verify Vector arbitrary is accessible
        let _ = v.norm();
    }

    #[test]
    fn non_zero_vec2_arbitrary_accessible(v in any::<dim2::arbitrary::NonZeroVector<f64>>()) {
        // Verify NonZeroVector wrapper is accessible and has non-zero norm
        prop_assert!(v.norm_squared() > 0.01);
    }

    #[test]
    fn unit_vec2_arbitrary_accessible(v in any::<dim2::arbitrary::UnitVector<f64>>()) {
        // Verify UnitVector wrapper is accessible and has unit norm
        prop_assert!(abs_diff_eq!(v.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn unit_rotor2_arbitrary_accessible(r in any::<dim2::arbitrary::UnitRotor<f64>>()) {
        // Verify UnitRotor wrapper is accessible and has unit norm
        prop_assert!(abs_diff_eq!(r.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}

// ============================================================================
// dim3 Arbitrary accessibility
// ============================================================================

proptest! {
    #[test]
    fn vec3_arbitrary_accessible(v in any::<dim3::Vector<f64>>()) {
        // Verify Vector arbitrary is accessible from generated code
        let _ = v.norm();
    }

    #[test]
    fn bivec3_arbitrary_accessible(b in any::<dim3::Bivector<f64>>()) {
        // Verify Bivector arbitrary is accessible from generated code
        let _ = b.norm();
    }

    #[test]
    fn rotor3_arbitrary_accessible(r in any::<dim3::Rotor<f64>>()) {
        // Verify Rotor arbitrary is accessible from generated code
        let _ = r.norm();
    }
}
