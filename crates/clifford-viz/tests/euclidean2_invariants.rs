//! Property-based tests for 2D Euclidean rotor visualization.
//!
//! These tests verify that the rotor operations used in the euclidean2 demo
//! maintain the expected geometric invariants:
//!
//! - Rotations preserve vector magnitude
//! - Rotations change angle by the expected amount
//! - Rotor components match the half-angle formula
//! - Rotations are composable (R2 * R1 = R(a1 + a2))

use std::f64::consts::PI;

use approx::assert_relative_eq;
use clifford::ops::Transform;
use clifford::specialized::euclidean::dim2::{Rotor, Vector};
use clifford_viz::testing::{
    VIZ_EPSILON, approx_eq, magnitude_2d, normalize_angle, rotation_changes_angle_by,
    rotation_preserves_magnitude, rotor_is_normalized, rotor_matches_half_angle,
    viz_proptest_config,
};
use proptest::prelude::*;

proptest! {
    #![proptest_config(viz_proptest_config())]

    /// Rotation preserves the magnitude (length) of vectors.
    ///
    /// This is a fundamental property: rotations are isometries.
    #[test]
    fn rotation_preserves_vector_magnitude(
        angle in -PI..PI,
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let rotor = Rotor::<f64>::from_angle(angle);
        let input = Vector::new(x, y);
        let output = rotor.transform(&input);

        prop_assert!(
            rotation_preserves_magnitude(x, y, output.x(), output.y(), VIZ_EPSILON),
            "Rotation should preserve magnitude: input=({}, {}), output=({}, {})",
            x, y, output.x(), output.y()
        );
    }

    /// Rotation changes the vector angle by exactly the rotor angle.
    #[test]
    fn rotation_changes_angle_correctly(
        angle in -PI..PI,
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        // Skip near-zero vectors where angle is undefined
        prop_assume!(magnitude_2d(x, y) > 0.01);

        let rotor = Rotor::<f64>::from_angle(angle);
        let input = Vector::new(x, y);
        let output = rotor.transform(&input);

        prop_assert!(
            rotation_changes_angle_by(x, y, output.x(), output.y(), angle, 1e-5),
            "Rotation should change angle by {}: input=({}, {}), output=({}, {})",
            angle, x, y, output.x(), output.y()
        );
    }

    /// Rotor is always normalized (unit magnitude in the even subalgebra).
    #[test]
    fn rotor_is_unit_magnitude(angle in -PI..PI) {
        let rotor = Rotor::<f64>::from_angle(angle);

        prop_assert!(
            rotor_is_normalized(rotor.s(), rotor.b(), VIZ_EPSILON),
            "Rotor should have unit norm: s={}, b={}, |R|={}",
            rotor.s(), rotor.b(),
            (rotor.s().powi(2) + rotor.b().powi(2)).sqrt()
        );
    }

    /// Rotor components match the half-angle formula: R = cos(θ/2) + sin(θ/2)e₁₂
    #[test]
    fn rotor_components_match_half_angle_formula(angle in -PI..PI) {
        let rotor = Rotor::<f64>::from_angle(angle);

        prop_assert!(
            rotor_matches_half_angle(rotor.s(), rotor.b(), angle, VIZ_EPSILON),
            "Rotor should match half-angle formula for angle {}: s={} (expected {}), b={} (expected {})",
            angle, rotor.s(), (angle / 2.0).cos(), rotor.b(), (angle / 2.0).sin()
        );
    }

    /// Rotations compose additively: R(a) * R(b) applied to v = R(a+b) applied to v
    #[test]
    fn rotations_compose_additively(
        angle_a in -PI..PI,
        angle_b in -PI..PI,
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        // Skip near-zero vectors
        prop_assume!(magnitude_2d(x, y) > 0.01);

        let rotor_a = Rotor::<f64>::from_angle(angle_a);
        let rotor_b = Rotor::<f64>::from_angle(angle_b);
        let rotor_combined = Rotor::<f64>::from_angle(angle_a + angle_b);

        let input = Vector::new(x, y);

        // Apply rotations sequentially
        let after_a = rotor_a.transform(&input);
        let after_ab = rotor_b.transform(&after_a);

        // Apply combined rotation
        let combined_result = rotor_combined.transform(&input);

        // Results should match
        prop_assert!(
            approx_eq(after_ab.x(), combined_result.x(), 1e-5) &&
            approx_eq(after_ab.y(), combined_result.y(), 1e-5),
            "Sequential rotations should equal combined: R_b(R_a(v))=({}, {}), R_(a+b)(v)=({}, {})",
            after_ab.x(), after_ab.y(), combined_result.x(), combined_result.y()
        );
    }

    /// Identity rotor (angle=0) leaves vectors unchanged.
    #[test]
    fn identity_rotor_preserves_vector(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let rotor = Rotor::<f64>::from_angle(0.0);
        let input = Vector::new(x, y);
        let output = rotor.transform(&input);

        prop_assert!(
            approx_eq(output.x(), x, VIZ_EPSILON) && approx_eq(output.y(), y, VIZ_EPSILON),
            "Identity rotor should preserve vector: input=({}, {}), output=({}, {})",
            x, y, output.x(), output.y()
        );
    }

    /// 90° rotation transforms (1, 0) to (0, 1).
    #[test]
    fn ninety_degree_rotation_is_correct(
        scale in 0.1f64..10.0,
    ) {
        let rotor = Rotor::<f64>::from_angle(PI / 2.0);
        let input = Vector::new(scale, 0.0);
        let output = rotor.transform(&input);

        prop_assert!(
            approx_eq(output.x(), 0.0, 1e-5) && approx_eq(output.y(), scale, 1e-5),
            "90° rotation of ({}, 0) should be (0, {}): got ({}, {})",
            scale, scale, output.x(), output.y()
        );
    }

    /// 180° rotation negates vectors.
    #[test]
    fn one_eighty_rotation_negates_vector(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let rotor = Rotor::<f64>::from_angle(PI);
        let input = Vector::new(x, y);
        let output = rotor.transform(&input);

        prop_assert!(
            approx_eq(output.x(), -x, 1e-5) && approx_eq(output.y(), -y, 1e-5),
            "180° rotation should negate vector: input=({}, {}), output=({}, {})",
            x, y, output.x(), output.y()
        );
    }

    /// 360° rotation returns to original position.
    #[test]
    fn full_rotation_returns_to_original(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let rotor = Rotor::<f64>::from_angle(2.0 * PI);
        let input = Vector::new(x, y);
        let output = rotor.transform(&input);

        prop_assert!(
            approx_eq(output.x(), x, 1e-5) && approx_eq(output.y(), y, 1e-5),
            "360° rotation should return to original: input=({}, {}), output=({}, {})",
            x, y, output.x(), output.y()
        );
    }
}

/// Unit tests for specific cases and edge cases.
#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn rotor_at_zero_is_identity() {
        let rotor = Rotor::<f64>::from_angle(0.0);
        assert_relative_eq!(rotor.s(), 1.0, epsilon = VIZ_EPSILON);
        assert_relative_eq!(rotor.b(), 0.0, epsilon = VIZ_EPSILON);
    }

    #[test]
    fn rotor_at_pi_has_zero_scalar() {
        let rotor = Rotor::<f64>::from_angle(PI);
        // cos(π/2) = 0, -sin(π/2) = -1 (note: bivector is negated in this library)
        assert_relative_eq!(rotor.s(), 0.0, epsilon = VIZ_EPSILON);
        assert_relative_eq!(rotor.b(), -1.0, epsilon = VIZ_EPSILON);
    }

    #[test]
    fn rotor_at_negative_angle_has_opposite_bivector() {
        let rotor_pos = Rotor::<f64>::from_angle(PI / 4.0);
        let rotor_neg = Rotor::<f64>::from_angle(-PI / 4.0);

        // Same scalar (cos is even)
        assert_relative_eq!(rotor_pos.s(), rotor_neg.s(), epsilon = VIZ_EPSILON);
        // Opposite bivector (sin is odd, and library negates it)
        assert_relative_eq!(rotor_pos.b(), -rotor_neg.b(), epsilon = VIZ_EPSILON);
    }

    #[test]
    fn normalize_angle_works_correctly() {
        assert_relative_eq!(normalize_angle(0.0), 0.0, epsilon = VIZ_EPSILON);
        assert_relative_eq!(normalize_angle(PI), PI, epsilon = VIZ_EPSILON);
        assert_relative_eq!(normalize_angle(-PI), -PI, epsilon = VIZ_EPSILON);
        // 3π wraps to π (or -π, both are equivalent at the branch cut)
        let normalized = normalize_angle(3.0 * PI);
        assert!(
            (normalized - PI).abs() < 1e-5 || (normalized + PI).abs() < 1e-5,
            "normalize_angle(3π) should be ±π, got {}",
            normalized
        );
    }
}
