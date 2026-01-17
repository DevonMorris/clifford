//! Geometric invariant testing utilities.
//!
//! Provides helper functions for testing that geometric operations preserve
//! expected invariants. These are used in property-based tests to verify
//! correctness without relying on golden images.

/// Default epsilon for floating-point comparisons in visualization tests.
pub const VIZ_EPSILON: f64 = 1e-6;

/// Default max_relative for floating-point comparisons.
pub const VIZ_MAX_RELATIVE: f64 = 1e-5;

/// Check if two f64 values are approximately equal.
#[must_use]
pub fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
    (a - b).abs() < epsilon
}

/// Check if two f32 values are approximately equal.
#[must_use]
pub fn approx_eq_f32(a: f32, b: f32, epsilon: f32) -> bool {
    (a - b).abs() < epsilon
}

/// Check if a value is approximately zero.
#[must_use]
pub fn approx_zero(a: f64, epsilon: f64) -> bool {
    a.abs() < epsilon
}

/// Compute the distance between two 2D points.
#[must_use]
pub fn distance_2d(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt()
}

/// Compute the distance from origin for a 2D point.
#[must_use]
pub fn magnitude_2d(x: f64, y: f64) -> f64 {
    (x * x + y * y).sqrt()
}

/// Compute the angle of a 2D vector from the positive x-axis.
#[must_use]
pub fn angle_2d(x: f64, y: f64) -> f64 {
    y.atan2(x)
}

/// Compute the angle between two 2D vectors.
#[must_use]
pub fn angle_between_2d(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let dot = x1 * x2 + y1 * y2;
    let mag1 = magnitude_2d(x1, y1);
    let mag2 = magnitude_2d(x2, y2);
    if mag1 < VIZ_EPSILON || mag2 < VIZ_EPSILON {
        return 0.0;
    }
    (dot / (mag1 * mag2)).clamp(-1.0, 1.0).acos()
}

/// Normalize the angle to be in the range [-PI, PI].
#[must_use]
pub fn normalize_angle(angle: f64) -> f64 {
    let mut a = angle % std::f64::consts::TAU;
    if a > std::f64::consts::PI {
        a -= std::f64::consts::TAU;
    } else if a < -std::f64::consts::PI {
        a += std::f64::consts::TAU;
    }
    a
}

/// Check that a rotation preserves the magnitude of a vector.
///
/// This is a fundamental invariant: rotations should not change vector length.
#[must_use]
pub fn rotation_preserves_magnitude(
    original_x: f64,
    original_y: f64,
    rotated_x: f64,
    rotated_y: f64,
    epsilon: f64,
) -> bool {
    let original_mag = magnitude_2d(original_x, original_y);
    let rotated_mag = magnitude_2d(rotated_x, rotated_y);
    approx_eq(original_mag, rotated_mag, epsilon)
}

/// Check that a rotation changes the angle by the expected amount.
///
/// Note: This handles angle wrapping correctly.
#[must_use]
pub fn rotation_changes_angle_by(
    original_x: f64,
    original_y: f64,
    rotated_x: f64,
    rotated_y: f64,
    expected_rotation: f64,
    epsilon: f64,
) -> bool {
    // Skip near-zero vectors (angle is undefined)
    if magnitude_2d(original_x, original_y) < epsilon {
        return true;
    }

    let original_angle = angle_2d(original_x, original_y);
    let rotated_angle = angle_2d(rotated_x, rotated_y);
    let actual_rotation = normalize_angle(rotated_angle - original_angle);
    let expected_normalized = normalize_angle(expected_rotation);

    approx_eq(actual_rotation, expected_normalized, epsilon)
}

/// Check that two rotations composed equal the sum of their angles.
///
/// For 2D rotations: R(a) * R(b) = R(a + b)
#[must_use]
pub fn rotations_compose_additively(
    angle_a: f64,
    angle_b: f64,
    composed_angle: f64,
    epsilon: f64,
) -> bool {
    let expected = normalize_angle(angle_a + angle_b);
    let actual = normalize_angle(composed_angle);
    approx_eq(expected, actual, epsilon)
}

/// Check that a rotor has unit norm (is normalized).
///
/// For a 2D rotor R = s + b*e12, the norm is sqrt(s^2 + b^2).
#[must_use]
pub fn rotor_is_normalized(scalar: f64, bivector: f64, epsilon: f64) -> bool {
    let norm = (scalar * scalar + bivector * bivector).sqrt();
    approx_eq(norm, 1.0, epsilon)
}

/// Check that the rotor components match the half-angle formula.
///
/// Note: The clifford library uses the convention `R = cos(θ/2) - sin(θ/2)e₁₂`
/// (negated bivector) so that the sandwich product `RvR†` rotates in the
/// positive direction.
#[must_use]
pub fn rotor_matches_half_angle(scalar: f64, bivector: f64, angle: f64, epsilon: f64) -> bool {
    let expected_scalar = (angle / 2.0).cos();
    // Note: bivector is negated in this convention
    let expected_bivector = -(angle / 2.0).sin();
    approx_eq(scalar, expected_scalar, epsilon) && approx_eq(bivector, expected_bivector, epsilon)
}
