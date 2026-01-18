//! Property-based tests for 2D Projective GA visualization.
//!
//! These tests verify that the projective geometry operations used in the
//! projective2 demo maintain the expected geometric invariants:
//!
//! - Join of two points creates a line incident with both points
//! - Meet of two lines creates a point incident with both lines
//! - Parallel lines meet at an ideal point
//! - Motor transformations preserve distances and incidence
//! - Reflections are self-inverse

use std::f64::consts::PI;

use clifford::ops::{Join, Meet, Transform};
use clifford::specialized::projective::dim2::{Flector, Line, Motor, Point};
use clifford_viz::testing::{VIZ_EPSILON, approx_eq, viz_proptest_config};
use proptest::prelude::*;

/// Epsilon for projective geometry tests (slightly larger due to homogeneous coords).
const PGA_EPSILON: f64 = 1e-5;

/// Check if a point lies on a line (within tolerance).
///
/// For point P and line L, they are incident if L·P ≈ 0.
fn point_on_line(point: &Point<f64>, line: &Line<f64>, epsilon: f64) -> bool {
    // Line equation: nx*x + ny*y + d*w = 0
    let incidence = line.nx() * point.x() + line.ny() * point.y() + line.d() * point.w();
    incidence.abs() < epsilon * (point.w().abs() + 1.0)
}

/// Check if a point is ideal (at infinity).
fn is_ideal_point(point: &Point<f64>, epsilon: f64) -> bool {
    point.w().abs() < epsilon
}

proptest! {
    #![proptest_config(viz_proptest_config())]

    // =========================================================================
    // Join operation invariants
    // =========================================================================

    /// Join of two points creates a line incident with both points.
    ///
    /// This is the fundamental property of the join operation.
    #[test]
    fn join_creates_incident_line(
        x1 in -10.0f64..10.0,
        y1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0,
        y2 in -10.0f64..10.0,
    ) {
        // Skip coincident points (line is undefined)
        let dist_sq = (x2 - x1).powi(2) + (y2 - y1).powi(2);
        prop_assume!(dist_sq > 0.01);

        let p1 = Point::from_cartesian(x1, y1);
        let p2 = Point::from_cartesian(x2, y2);
        let line = p1.join(&p2);

        // Both points should lie on the line
        prop_assert!(
            point_on_line(&p1, &line, PGA_EPSILON),
            "Point P1=({}, {}) should lie on line through P1 and P2. Incidence: {}",
            x1, y1,
            line.nx() * p1.x() + line.ny() * p1.y() + line.d() * p1.w()
        );
        prop_assert!(
            point_on_line(&p2, &line, PGA_EPSILON),
            "Point P2=({}, {}) should lie on line through P1 and P2. Incidence: {}",
            x2, y2,
            line.nx() * p2.x() + line.ny() * p2.y() + line.d() * p2.w()
        );
    }

    /// Join is anticommutative up to sign.
    ///
    /// P1 ∧ P2 = -(P2 ∧ P1) for points.
    #[test]
    fn join_anticommutative(
        x1 in -10.0f64..10.0,
        y1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0,
        y2 in -10.0f64..10.0,
    ) {
        let dist_sq = (x2 - x1).powi(2) + (y2 - y1).powi(2);
        prop_assume!(dist_sq > 0.01);

        let p1 = Point::from_cartesian(x1, y1);
        let p2 = Point::from_cartesian(x2, y2);

        let l12 = p1.join(&p2);
        let l21 = p2.join(&p1);

        // Should be negatives of each other
        prop_assert!(
            approx_eq(l12.d(), -l21.d(), PGA_EPSILON) &&
            approx_eq(l12.nx(), -l21.nx(), PGA_EPSILON) &&
            approx_eq(l12.ny(), -l21.ny(), PGA_EPSILON),
            "Join should be anticommutative: P1∧P2 = -P2∧P1"
        );
    }

    // =========================================================================
    // Meet operation invariants
    // =========================================================================

    /// Meet of two non-parallel lines gives a finite intersection point.
    #[test]
    fn meet_non_parallel_gives_finite_point(
        // First line: ax + by + c = 0
        a1 in -5.0f64..5.0,
        b1 in -5.0f64..5.0,
        c1 in -5.0f64..5.0,
        // Second line (different slope)
        a2 in -5.0f64..5.0,
        b2 in -5.0f64..5.0,
        c2 in -5.0f64..5.0,
    ) {
        // Skip degenerate lines
        prop_assume!(a1.abs() + b1.abs() > 0.1);
        prop_assume!(a2.abs() + b2.abs() > 0.1);

        // Skip parallel lines (cross product of normals ≈ 0)
        let cross = a1 * b2 - b1 * a2;
        prop_assume!(cross.abs() > 0.1);

        let l1 = Line::from_equation(a1, b1, c1);
        let l2 = Line::from_equation(a2, b2, c2);
        let intersection = l1.meet(&l2);

        // Intersection should be a finite point
        prop_assert!(
            !is_ideal_point(&intersection, PGA_EPSILON),
            "Non-parallel lines should meet at a finite point, got w={}",
            intersection.w()
        );

        // Intersection should lie on both lines
        prop_assert!(
            point_on_line(&intersection, &l1, PGA_EPSILON),
            "Intersection should lie on L1"
        );
        prop_assert!(
            point_on_line(&intersection, &l2, PGA_EPSILON),
            "Intersection should lie on L2"
        );
    }

    /// Meet of parallel lines gives an ideal point.
    #[test]
    fn meet_parallel_gives_ideal_point(
        a in -5.0f64..5.0,
        b in -5.0f64..5.0,
        c1 in -5.0f64..5.0,
        c2 in -5.0f64..5.0,
    ) {
        // Skip degenerate lines
        prop_assume!(a.abs() + b.abs() > 0.1);

        // Skip coincident lines (same c value would give the same line)
        prop_assume!((c2 - c1).abs() > 0.1);

        // Two parallel lines with same normal but different offsets
        let l1 = Line::from_equation(a, b, c1);
        let l2 = Line::from_equation(a, b, c2);
        let intersection = l1.meet(&l2);

        // Intersection should be an ideal point (w ≈ 0)
        prop_assert!(
            is_ideal_point(&intersection, PGA_EPSILON),
            "Parallel lines should meet at ideal point, got w={}",
            intersection.w()
        );
    }

    // =========================================================================
    // Motor transformation invariants
    // =========================================================================

    /// Motor rotation preserves distance from origin.
    #[test]
    fn motor_rotation_preserves_distance(
        angle in -PI..PI,
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let motor = Motor::from_rotation(angle);
        let point = Point::from_cartesian(x, y);
        let rotated = motor.transform(&point);

        let original_dist = (x * x + y * y).sqrt();
        let rotated_coords = rotated.to_cartesian().expect("rotated point should be finite");
        let rotated_dist = (rotated_coords.0.powi(2) + rotated_coords.1.powi(2)).sqrt();

        prop_assert!(
            approx_eq(original_dist, rotated_dist, PGA_EPSILON),
            "Rotation should preserve distance: original={}, rotated={}",
            original_dist, rotated_dist
        );
    }

    /// Motor rotation by 90° rotates (1, 0) to (0, 1).
    #[test]
    fn motor_rotation_90_degrees(scale in 0.1f64..10.0) {
        let motor = Motor::from_rotation(PI / 2.0);
        let point = Point::from_cartesian(scale, 0.0);
        let rotated = motor.transform(&point);

        let (rx, ry) = rotated.to_cartesian().expect("rotated point should be finite");

        prop_assert!(
            approx_eq(rx, 0.0, PGA_EPSILON) && approx_eq(ry, scale, PGA_EPSILON),
            "90° rotation of ({}, 0) should be (0, {}), got ({}, {})",
            scale, scale, rx, ry
        );
    }

    /// Motor rotation by 180° negates coordinates.
    #[test]
    fn motor_rotation_180_degrees(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let motor = Motor::from_rotation(PI);
        let point = Point::from_cartesian(x, y);
        let rotated = motor.transform(&point);

        let (rx, ry) = rotated.to_cartesian().expect("rotated point should be finite");

        prop_assert!(
            approx_eq(rx, -x, PGA_EPSILON) && approx_eq(ry, -y, PGA_EPSILON),
            "180° rotation of ({}, {}) should be ({}, {}), got ({}, {})",
            x, y, -x, -y, rx, ry
        );
    }

    /// Motor translation moves points by the translation vector.
    #[test]
    fn motor_translation_moves_point(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
        dx in -5.0f64..5.0,
        dy in -5.0f64..5.0,
    ) {
        let motor = Motor::from_translation(dx, dy);
        let point = Point::from_cartesian(x, y);
        let translated = motor.transform(&point);

        let (tx, ty) = translated.to_cartesian().expect("translated point should be finite");

        prop_assert!(
            approx_eq(tx, x + dx, PGA_EPSILON) && approx_eq(ty, y + dy, PGA_EPSILON),
            "Translation of ({}, {}) by ({}, {}) should be ({}, {}), got ({}, {})",
            x, y, dx, dy, x + dx, y + dy, tx, ty
        );
    }

    /// Motor identity preserves points.
    #[test]
    fn motor_identity_preserves_point(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let motor = Motor::identity();
        let point = Point::from_cartesian(x, y);
        let result = motor.transform(&point);

        let (rx, ry) = result.to_cartesian().expect("result should be finite");

        prop_assert!(
            approx_eq(rx, x, PGA_EPSILON) && approx_eq(ry, y, PGA_EPSILON),
            "Identity motor should preserve ({}, {}), got ({}, {})",
            x, y, rx, ry
        );
    }

    /// Motor transformation preserves incidence (points on a line stay on the transformed line).
    #[test]
    fn motor_preserves_incidence(
        x1 in -5.0f64..5.0,
        y1 in -5.0f64..5.0,
        x2 in -5.0f64..5.0,
        y2 in -5.0f64..5.0,
        angle in -PI..PI,
        dx in -3.0f64..3.0,
        dy in -3.0f64..3.0,
    ) {
        // Skip coincident points
        let dist_sq = (x2 - x1).powi(2) + (y2 - y1).powi(2);
        prop_assume!(dist_sq > 0.01);

        let p1 = Point::from_cartesian(x1, y1);
        let p2 = Point::from_cartesian(x2, y2);
        let line = p1.join(&p2);

        // Apply rotation then translation
        let rot_motor = Motor::from_rotation(angle);
        let trans_motor = Motor::from_translation(dx, dy);

        let p1_rot = rot_motor.transform(&p1);
        let p2_rot = rot_motor.transform(&p2);
        let line_rot = rot_motor.transform(&line);

        let p1_trans = trans_motor.transform(&p1_rot);
        let p2_trans = trans_motor.transform(&p2_rot);
        let line_trans = trans_motor.transform(&line_rot);

        // Transformed points should still lie on transformed line
        prop_assert!(
            point_on_line(&p1_trans, &line_trans, PGA_EPSILON * 10.0),
            "P1 should remain on transformed line"
        );
        prop_assert!(
            point_on_line(&p2_trans, &line_trans, PGA_EPSILON * 10.0),
            "P2 should remain on transformed line"
        );
    }

    /// Motor composition: rotating twice equals rotating by sum of angles.
    #[test]
    fn motor_rotation_composition(
        angle1 in -PI..PI,
        angle2 in -PI..PI,
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        prop_assume!(x.abs() + y.abs() > 0.1);

        let m1 = Motor::from_rotation(angle1);
        let m2 = Motor::from_rotation(angle2);
        let m_combined = Motor::from_rotation(angle1 + angle2);

        let point = Point::from_cartesian(x, y);

        // Apply sequentially
        let after_m1 = m1.transform(&point);
        let after_m1_m2 = m2.transform(&after_m1);

        // Apply combined
        let after_combined = m_combined.transform(&point);

        let (seq_x, seq_y) = after_m1_m2.to_cartesian().expect("should be finite");
        let (com_x, com_y) = after_combined.to_cartesian().expect("should be finite");

        prop_assert!(
            approx_eq(seq_x, com_x, PGA_EPSILON * 10.0) &&
            approx_eq(seq_y, com_y, PGA_EPSILON * 10.0),
            "R(a)R(b) should equal R(a+b): sequential=({}, {}), combined=({}, {})",
            seq_x, seq_y, com_x, com_y
        );
    }

    // =========================================================================
    // Flector (reflection) invariants
    // =========================================================================

    /// Reflection is self-inverse: reflecting twice returns to original.
    #[test]
    fn flector_is_self_inverse(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
        line_angle in -PI..PI,
    ) {
        // Create a reflection through a line at given angle
        let flector = Flector::reflect_through_angle(line_angle);
        let point = Point::from_cartesian(x, y);

        let reflected_once = flector.transform(&point);
        let reflected_twice = flector.transform(&reflected_once);

        let (rx, ry) = reflected_twice.to_cartesian().expect("should be finite");

        prop_assert!(
            approx_eq(rx, x, PGA_EPSILON) && approx_eq(ry, y, PGA_EPSILON),
            "Double reflection should return to original: ({}, {}) → ({}, {})",
            x, y, rx, ry
        );
    }

    /// Reflection through x-axis negates y coordinate.
    #[test]
    fn flector_x_axis_negates_y(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let flector = Flector::reflect_x_axis();
        let point = Point::from_cartesian(x, y);
        let reflected = flector.transform(&point);

        let (rx, ry) = reflected.to_cartesian().expect("should be finite");

        prop_assert!(
            approx_eq(rx, x, PGA_EPSILON) && approx_eq(ry, -y, PGA_EPSILON),
            "X-axis reflection of ({}, {}) should be ({}, {}), got ({}, {})",
            x, y, x, -y, rx, ry
        );
    }

    /// Reflection through y-axis negates x coordinate.
    #[test]
    fn flector_y_axis_negates_x(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let flector = Flector::reflect_y_axis();
        let point = Point::from_cartesian(x, y);
        let reflected = flector.transform(&point);

        let (rx, ry) = reflected.to_cartesian().expect("should be finite");

        prop_assert!(
            approx_eq(rx, -x, PGA_EPSILON) && approx_eq(ry, y, PGA_EPSILON),
            "Y-axis reflection of ({}, {}) should be ({}, {}), got ({}, {})",
            x, y, -x, y, rx, ry
        );
    }

    /// Reflection preserves distance from origin (for reflections through origin).
    #[test]
    fn flector_through_origin_preserves_distance(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
        angle in -PI..PI,
    ) {
        let flector = Flector::reflect_through_angle(angle);
        let point = Point::from_cartesian(x, y);
        let reflected = flector.transform(&point);

        let original_dist = (x * x + y * y).sqrt();
        let (rx, ry) = reflected.to_cartesian().expect("should be finite");
        let reflected_dist = (rx * rx + ry * ry).sqrt();

        prop_assert!(
            approx_eq(original_dist, reflected_dist, PGA_EPSILON),
            "Reflection should preserve distance from origin: {} vs {}",
            original_dist, reflected_dist
        );
    }

    // =========================================================================
    // Distance and geometric property invariants
    // =========================================================================

    /// Distance is symmetric: d(A, B) = d(B, A).
    #[test]
    fn distance_is_symmetric(
        x1 in -10.0f64..10.0,
        y1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0,
        y2 in -10.0f64..10.0,
    ) {
        let p1 = Point::from_cartesian(x1, y1);
        let p2 = Point::from_cartesian(x2, y2);

        let d12 = p1.distance(&p2);
        let d21 = p2.distance(&p1);

        prop_assert!(
            approx_eq(d12, d21, VIZ_EPSILON),
            "Distance should be symmetric: d(P1, P2)={}, d(P2, P1)={}",
            d12, d21
        );
    }

    /// Distance to self is zero.
    #[test]
    fn distance_to_self_is_zero(
        x in -10.0f64..10.0,
        y in -10.0f64..10.0,
    ) {
        let p = Point::from_cartesian(x, y);
        let dist = p.distance(&p);

        prop_assert!(
            approx_eq(dist, 0.0, VIZ_EPSILON),
            "Distance to self should be zero, got {}",
            dist
        );
    }

    /// Midpoint is equidistant from both endpoints.
    #[test]
    fn midpoint_is_equidistant(
        x1 in -10.0f64..10.0,
        y1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0,
        y2 in -10.0f64..10.0,
    ) {
        let p1 = Point::from_cartesian(x1, y1);
        let p2 = Point::from_cartesian(x2, y2);
        let mid = p1.midpoint(&p2);

        let d1 = mid.distance(&p1);
        let d2 = mid.distance(&p2);

        prop_assert!(
            approx_eq(d1, d2, PGA_EPSILON),
            "Midpoint should be equidistant: d(mid, P1)={}, d(mid, P2)={}",
            d1, d2
        );
    }
}

/// Unit tests for specific edge cases.
#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn axes_meet_at_origin() {
        let x_axis = Line::x_axis();
        let y_axis = Line::y_axis();
        let origin = x_axis.meet(&y_axis);

        let (x, y) = origin.to_cartesian().expect("should be finite point");
        assert!(approx_eq(x, 0.0, VIZ_EPSILON), "x should be 0, got {}", x);
        assert!(approx_eq(y, 0.0, VIZ_EPSILON), "y should be 0, got {}", y);
    }

    #[test]
    fn line_through_origin_has_zero_d() {
        let p1 = Point::origin();
        let p2 = Point::from_cartesian(1.0, 1.0);
        let line = p1.join(&p2);

        assert!(
            approx_eq(line.d(), 0.0, VIZ_EPSILON),
            "Line through origin should have d=0, got {}",
            line.d()
        );
    }

    #[test]
    fn motor_360_rotation_returns_to_original() {
        let motor = Motor::from_rotation(2.0 * PI);
        let point = Point::from_cartesian(3.0, 4.0);
        let rotated = motor.transform(&point);

        let (rx, ry) = rotated.to_cartesian().expect("should be finite");
        assert!(
            approx_eq(rx, 3.0, PGA_EPSILON) && approx_eq(ry, 4.0, PGA_EPSILON),
            "360° rotation should return to original: got ({}, {})",
            rx,
            ry
        );
    }

    #[test]
    fn point_on_line_check_works() {
        let p1 = Point::from_cartesian(0.0, 0.0);
        let p2 = Point::from_cartesian(1.0, 1.0);
        let line = p1.join(&p2);

        // Point on the line
        let p_on = Point::from_cartesian(0.5, 0.5);
        assert!(point_on_line(&p_on, &line, PGA_EPSILON));

        // Point not on the line
        let p_off = Point::from_cartesian(0.5, 1.0);
        assert!(!point_on_line(&p_off, &line, PGA_EPSILON));
    }
}
