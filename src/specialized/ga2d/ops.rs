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
    use crate::specialized::ga2d::arbitrary::{NonZeroVec2, UnitRotor2, UnitVec2};
    use approx::abs_diff_eq;
    use proptest::prelude::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    // ========================================================================
    // Vec2 tests
    // ========================================================================

    proptest! {
        #[test]
        fn vec2_add_commutative(a in any::<Vec2<f64>>(), b in any::<Vec2<f64>>()) {
            let ab = a + b;
            let ba = b + a;
            prop_assert!(abs_diff_eq!(ab, ba, epsilon = 1e-10));
        }

        #[test]
        fn vec2_dot_commutative(a in any::<Vec2<f64>>(), b in any::<Vec2<f64>>()) {
            prop_assert!(abs_diff_eq!(a.dot(b), b.dot(a), epsilon = 1e-10));
        }

        #[test]
        fn vec2_wedge_anticommutative(a in any::<Vec2<f64>>(), b in any::<Vec2<f64>>()) {
            let ab = a.wedge(b);
            let ba = b.wedge(a);
            prop_assert!(abs_diff_eq!(ab.0, -ba.0, epsilon = 1e-10));
        }

        #[test]
        fn vec2_normalized_has_unit_length(v in any::<NonZeroVec2>()) {
            let n = v.normalized();
            prop_assert!(abs_diff_eq!(n.norm(), 1.0, epsilon = 1e-10));
        }

        #[test]
        fn vec2_perp_is_perpendicular(v in any::<NonZeroVec2>()) {
            let p = v.perp();
            prop_assert!(abs_diff_eq!(v.dot(p), 0.0, epsilon = 1e-10));
        }
    }

    // ========================================================================
    // Rotor tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor2_preserves_norm(r in any::<UnitRotor2>(), v in any::<Vec2<f64>>()) {
            let rotated = r.rotate(v);
            prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = 1e-9));
        }

        #[test]
        fn rotor2_composition(
            r1 in any::<UnitRotor2>(),
            r2 in any::<UnitRotor2>(),
            v in any::<Vec2<f64>>()
        ) {
            let sequential = r2.rotate(r1.rotate(v));
            let composed = r2.compose(*r1).rotate(v);
            prop_assert!(abs_diff_eq!(sequential, composed, epsilon = 1e-8));
        }

        #[test]
        fn rotor2_inverse(r in any::<UnitRotor2>(), v in any::<Vec2<f64>>()) {
            let roundtrip = r.inverse().rotate(r.rotate(v));
            prop_assert!(abs_diff_eq!(roundtrip, v, epsilon = 1e-8));
        }

        #[test]
        fn rotor2_from_vectors(a in any::<UnitVec2>(), b in any::<UnitVec2>()) {
            let r = Rotor2::from_vectors(*a, *b);
            let rotated = r.rotate(*a);
            prop_assert!(abs_diff_eq!(rotated, *b, epsilon = 1e-8));
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

        assert!(abs_diff_eq!(rotated, Vec2::unit_y(), epsilon = 1e-10));
    }

    #[test]
    fn rotor2_180_deg_rotation() {
        let rotor = Rotor2::from_angle(PI);
        let v = Vec2::unit_x();
        let rotated = rotor.rotate(v);

        assert!(abs_diff_eq!(rotated, -Vec2::unit_x(), epsilon = 1e-10));
    }

    #[test]
    fn rotor2_45_plus_45_equals_90() {
        let r45 = Rotor2::from_angle(FRAC_PI_4);
        let r90 = r45.compose(r45);
        let v = Vec2::unit_x();
        let rotated = r90.rotate(v);

        assert!(abs_diff_eq!(rotated, Vec2::unit_y(), epsilon = 1e-10));
    }

    #[test]
    fn rotor2_angle_roundtrip() {
        let angle = 1.234;
        let rotor = Rotor2::from_angle(angle);
        assert!(abs_diff_eq!(rotor.angle(), angle, epsilon = 1e-10));
    }
}
