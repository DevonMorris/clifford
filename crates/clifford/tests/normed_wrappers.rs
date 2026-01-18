//! Test that Normed trait implementations work with wrapper types.

use approx::relative_eq;
use clifford::norm::Normed;
use clifford::specialized::euclidean::dim3::Vector;
use clifford::wrappers::Unit;

const EPS: f64 = 1e-10;

#[test]
fn vector_implements_normed_trait() {
    let v = Vector::new(3.0, 4.0, 0.0);

    // Test norm_squared via trait
    let ns = Normed::norm_squared(&v);
    assert!(relative_eq!(ns, 25.0, epsilon = EPS, max_relative = EPS));

    // Test norm via trait
    let n = Normed::norm(&v);
    assert!(relative_eq!(n, 5.0, epsilon = EPS, max_relative = EPS));
}

#[test]
fn unit_wrapper_works_with_vector() {
    let v = Vector::new(3.0, 4.0, 0.0);

    // Create a Unit<Vector> using try_new
    let unit: Unit<Vector<f64>> = Unit::try_new(v).expect("should normalize");

    // Verify norm is 1
    assert!(relative_eq!(
        unit.norm(),
        1.0,
        epsilon = EPS,
        max_relative = EPS
    ));

    // Verify components (use as_inner() for inherent Vector methods)
    assert!(relative_eq!(
        unit.as_inner().x(),
        0.6,
        epsilon = EPS,
        max_relative = EPS
    ));
    assert!(relative_eq!(
        unit.as_inner().y(),
        0.8,
        epsilon = EPS,
        max_relative = EPS
    ));
    assert!(relative_eq!(
        unit.as_inner().z(),
        0.0,
        epsilon = EPS,
        max_relative = EPS
    ));
}

#[test]
fn unit_new_normalize_works() {
    let v = Vector::new(1.0, 1.0, 1.0);

    // Create via new_normalize
    let unit = Unit::new_normalize(v);

    // Verify norm is 1
    assert!(relative_eq!(
        unit.norm(),
        1.0,
        epsilon = EPS,
        max_relative = EPS
    ));
}

#[test]
fn unit_try_new_returns_none_for_zero() {
    let zero = Vector::new(0.0, 0.0, 0.0);
    assert!(Unit::try_new(zero).is_none());
}

// ============================================================================
// DegenerateNormed tests for PGA types
// ============================================================================

use clifford::norm::DegenerateNormed;
use clifford::specialized::projective::dim3::{Motor, Point};
use clifford::wrappers::Bulk;

#[test]
fn motor_implements_degenerate_normed() {
    // A pure rotation motor (no translation)
    // Per RGA convention: rotation uses velocity bivectors (tx, ry, rz) and ps for identity
    // For x-axis rotation by pi/2: tx = sin(pi/4), ps = cos(pi/4)
    let motor = Motor::new_unchecked(
        0.0,            // s = 0 (no translation coupling)
        0.0,            // tz (e12) - moment z
        0.0,            // ty (e13) - moment y
        0.5_f64.sqrt(), // tx (e14) = sin(pi/4) - velocity x (rotation around x)
        0.0,            // rx (e23) - moment x
        0.0,            // ry (e24) - velocity y
        0.0,            // rz (e34) - velocity z
        0.5_f64.sqrt(), // ps = cos(pi/4) - identity
    );

    // Bulk norm should be 0 for pure rotation (no translation/moment)
    // Bulk = sqrt(s² + tz² + ty² + rx²) = 0
    let bulk = motor.bulk_norm();
    assert!(bulk.abs() < EPS);

    // Weight norm should be 1 for a unit rotation
    // Weight = sqrt(tx² + ry² + rz² + ps²) = sqrt(0.5 + 0 + 0 + 0.5) = 1
    let weight = motor.weight_norm();
    assert!(relative_eq!(weight, 1.0, epsilon = EPS, max_relative = EPS));
}

#[test]
fn bulk_wrapper_works_with_motor() {
    let motor = Motor::new_unchecked(2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    // Create Bulk<Motor> - should normalize by bulk norm
    let bulk: Bulk<Motor<f64>> = Bulk::try_new(motor).expect("should bulk-normalize");

    // Bulk norm should be 1
    assert!(relative_eq!(
        bulk.bulk_norm(),
        1.0,
        epsilon = EPS,
        max_relative = EPS
    ));
}

#[test]
fn point_has_bulk_and_weight_norm() {
    // Point::new_unchecked(e1, e2, e3, e0)
    let point = Point::new_unchecked(1.0, 2.0, 3.0, 4.0);

    // Bulk is from e1, e2, e3 components
    let bulk_sq = point.bulk_norm_squared();
    assert!(relative_eq!(
        bulk_sq,
        1.0 + 4.0 + 9.0, // e1² + e2² + e3² = 1² + 2² + 3²
        epsilon = EPS,
        max_relative = EPS
    ));

    // Weight is from e0 component
    let weight_sq = point.weight_norm_squared();
    assert!(relative_eq!(
        weight_sq,
        16.0,
        epsilon = EPS,
        max_relative = EPS
    )); // e0² = 4²
}
