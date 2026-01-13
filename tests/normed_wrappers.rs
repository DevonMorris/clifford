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

    // Verify components
    assert!(relative_eq!(
        unit.x(),
        0.6,
        epsilon = EPS,
        max_relative = EPS
    ));
    assert!(relative_eq!(
        unit.y(),
        0.8,
        epsilon = EPS,
        max_relative = EPS
    ));
    assert!(relative_eq!(
        unit.z(),
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
    let motor = Motor::new_unchecked(
        0.5_f64.sqrt(), // s = cos(pi/4)
        0.0,
        0.0,
        0.5_f64.sqrt(), // e12 = sin(pi/4)
        0.0,
        0.0,
        0.0,
        0.0,
    );

    // Bulk norm should be 1 for a unit rotation
    let bulk = motor.bulk_norm();
    assert!(relative_eq!(bulk, 1.0, epsilon = EPS, max_relative = EPS));

    // Weight norm should be 0 for pure rotation (no translation)
    let weight = motor.weight_norm();
    assert!(weight.abs() < EPS);
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
