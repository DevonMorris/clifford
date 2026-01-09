//! Operator implementations for 2D geometric algebra types.

use super::types::{Bivec2, Rotor2, Vec2};
use crate::scalar::Float;
use std::ops::{Add, BitXor, Mul, Neg, Sub};

// ============================================================================
// Vec2 operations
// ============================================================================

impl<T: Float> Neg for Vec2<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y)
    }
}

impl<T: Float> Add for Vec2<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y)
    }
}

impl<T: Float> Sub for Vec2<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.x - other.x, self.y - other.y)
    }
}

impl<T: Float> Mul<T> for Vec2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self::new(self.x * scalar, self.y * scalar)
    }
}

impl Mul<Vec2<f32>> for f32 {
    type Output = Vec2<f32>;

    #[inline]
    fn mul(self, v: Vec2<f32>) -> Self::Output {
        Vec2::new(self * v.x, self * v.y)
    }
}

impl Mul<Vec2<f64>> for f64 {
    type Output = Vec2<f64>;

    #[inline]
    fn mul(self, v: Vec2<f64>) -> Self::Output {
        Vec2::new(self * v.x, self * v.y)
    }
}

/// Wedge product operator: `a ^ b`.
impl<T: Float> BitXor for Vec2<T> {
    type Output = Bivec2<T>;

    #[inline]
    fn bitxor(self, other: Self) -> Self::Output {
        self.wedge(other)
    }
}

// ============================================================================
// Bivec2 operations
// ============================================================================

impl<T: Float> Neg for Bivec2<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl<T: Float> Add for Bivec2<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0)
    }
}

impl<T: Float> Sub for Bivec2<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self(self.0 - other.0)
    }
}

impl<T: Float> Mul<T> for Bivec2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self(self.0 * scalar)
    }
}

// ============================================================================
// Rotor2 operations
// ============================================================================

impl<T: Float> Neg for Rotor2<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.s, -self.xy)
    }
}

impl<T: Float> Add for Rotor2<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.s + other.s, self.xy + other.xy)
    }
}

impl<T: Float> Sub for Rotor2<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.s - other.s, self.xy - other.xy)
    }
}

/// Rotor multiplication (composition).
impl<T: Float> Mul for Rotor2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        self.compose(other)
    }
}

impl<T: Float> Mul<T> for Rotor2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self::new(self.s * scalar, self.xy * scalar)
    }
}

/// Rotor applied to vector.
impl<T: Float> Mul<Vec2<T>> for Rotor2<T> {
    type Output = Vec2<T>;

    #[inline]
    fn mul(self, v: Vec2<T>) -> Self::Output {
        self.rotate(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    /// Arbitrary Vec2 for property testing.
    fn arb_vec2() -> impl Strategy<Value = Vec2<f64>> {
        (-100.0..100.0, -100.0..100.0).prop_map(|(x, y)| Vec2::new(x, y))
    }

    /// Arbitrary non-zero Vec2.
    fn arb_nonzero_vec2() -> impl Strategy<Value = Vec2<f64>> {
        arb_vec2().prop_filter("non-zero vector", |v| v.norm_squared() > 1e-10)
    }

    /// Arbitrary unit Vec2.
    fn arb_unit_vec2() -> impl Strategy<Value = Vec2<f64>> {
        arb_nonzero_vec2().prop_map(|v| v.normalized())
    }

    /// Arbitrary unit Rotor2.
    fn arb_unit_rotor2() -> impl Strategy<Value = Rotor2<f64>> {
        (0.0..2.0 * PI).prop_map(|angle| Rotor2::from_angle(angle))
    }

    // ========================================================================
    // Vec2 tests
    // ========================================================================

    proptest! {
        #[test]
        fn vec2_add_commutative(a in arb_vec2(), b in arb_vec2()) {
            let ab = a + b;
            let ba = b + a;
            prop_assert!((ab.x - ba.x).abs() < 1e-10);
            prop_assert!((ab.y - ba.y).abs() < 1e-10);
        }

        #[test]
        fn vec2_dot_commutative(a in arb_vec2(), b in arb_vec2()) {
            prop_assert!((a.dot(b) - b.dot(a)).abs() < 1e-10);
        }

        #[test]
        fn vec2_wedge_anticommutative(a in arb_vec2(), b in arb_vec2()) {
            let ab = a.wedge(b);
            let ba = b.wedge(a);
            prop_assert!((ab.0 + ba.0).abs() < 1e-10);
        }

        #[test]
        fn vec2_normalized_has_unit_length(v in arb_nonzero_vec2()) {
            let n = v.normalized();
            prop_assert!((n.norm() - 1.0).abs() < 1e-10);
        }

        #[test]
        fn vec2_perp_is_perpendicular(v in arb_nonzero_vec2()) {
            let p = v.perp();
            prop_assert!(v.dot(p).abs() < 1e-10);
        }
    }

    // ========================================================================
    // Rotor tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor2_preserves_norm(r in arb_unit_rotor2(), v in arb_vec2()) {
            let rotated = r.rotate(v);
            prop_assert!((v.norm() - rotated.norm()).abs() < 1e-9);
        }

        #[test]
        fn rotor2_composition(
            r1 in arb_unit_rotor2(),
            r2 in arb_unit_rotor2(),
            v in arb_vec2()
        ) {
            let sequential = r2.rotate(r1.rotate(v));
            let composed = r2.compose(r1).rotate(v);
            prop_assert!((sequential.x - composed.x).abs() < 1e-8);
            prop_assert!((sequential.y - composed.y).abs() < 1e-8);
        }

        #[test]
        fn rotor2_inverse(r in arb_unit_rotor2(), v in arb_vec2()) {
            let roundtrip = r.inverse().rotate(r.rotate(v));
            prop_assert!((roundtrip.x - v.x).abs() < 1e-8);
            prop_assert!((roundtrip.y - v.y).abs() < 1e-8);
        }

        #[test]
        fn rotor2_from_vectors(a in arb_unit_vec2(), b in arb_unit_vec2()) {
            let r = Rotor2::from_vectors(a, b);
            let rotated = r.rotate(a);
            prop_assert!((rotated.x - b.x).abs() < 1e-8);
            prop_assert!((rotated.y - b.y).abs() < 1e-8);
        }
    }

    // ========================================================================
    // Specific rotation tests
    // ========================================================================

    #[test]
    fn rotor2_90_deg_rotation() {
        let rotor = Rotor2::from_angle(FRAC_PI_2);
        let v = Vec2::unit_x();
        let rotated = rotor.rotate(v);

        assert!(rotated.x.abs() < 1e-10);
        assert!((rotated.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rotor2_180_deg_rotation() {
        let rotor = Rotor2::from_angle(PI);
        let v = Vec2::unit_x();
        let rotated = rotor.rotate(v);

        assert!((rotated.x + 1.0).abs() < 1e-10);
        assert!(rotated.y.abs() < 1e-10);
    }

    #[test]
    fn rotor2_45_plus_45_equals_90() {
        let r45 = Rotor2::from_angle(FRAC_PI_4);
        let r90 = r45.compose(r45);
        let v = Vec2::unit_x();
        let rotated = r90.rotate(v);

        assert!(rotated.x.abs() < 1e-10);
        assert!((rotated.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rotor2_angle_roundtrip() {
        let angle = 1.234;
        let rotor = Rotor2::from_angle(angle);
        assert!((rotor.angle() - angle).abs() < 1e-10);
    }
}
