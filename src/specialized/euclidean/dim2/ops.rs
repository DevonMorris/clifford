//! Operator implementations for 2D geometric algebra types.

use super::types::{Bivector, Rotor, Vector};
use crate::scalar::Float;
use std::ops::{Add, BitXor, Mul, Neg, Sub};

// ============================================================================
// Vector operations
// ============================================================================

impl<T: Float> Neg for Vector<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y)
    }
}

impl<T: Float> Add for Vector<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y)
    }
}

impl<T: Float> Sub for Vector<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.x - other.x, self.y - other.y)
    }
}

impl<T: Float> Mul<T> for Vector<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self::new(self.x * scalar, self.y * scalar)
    }
}

impl Mul<Vector<f32>> for f32 {
    type Output = Vector<f32>;

    #[inline]
    fn mul(self, v: Vector<f32>) -> Self::Output {
        Vector::new(self * v.x, self * v.y)
    }
}

impl Mul<Vector<f64>> for f64 {
    type Output = Vector<f64>;

    #[inline]
    fn mul(self, v: Vector<f64>) -> Self::Output {
        Vector::new(self * v.x, self * v.y)
    }
}

/// Wedge product operator: `a ^ b`.
impl<T: Float> BitXor for Vector<T> {
    type Output = Bivector<T>;

    #[inline]
    fn bitxor(self, other: Self) -> Self::Output {
        self.wedge(other)
    }
}

// ============================================================================
// Bivector operations
// ============================================================================

impl<T: Float> Neg for Bivector<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl<T: Float> Add for Bivector<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0)
    }
}

impl<T: Float> Sub for Bivector<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self(self.0 - other.0)
    }
}

impl<T: Float> Mul<T> for Bivector<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self(self.0 * scalar)
    }
}

// ============================================================================
// Rotor operations
// ============================================================================

impl<T: Float> Neg for Rotor<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.s, -self.xy)
    }
}

impl<T: Float> Add for Rotor<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.s + other.s, self.xy + other.xy)
    }
}

impl<T: Float> Sub for Rotor<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.s - other.s, self.xy - other.xy)
    }
}

/// Rotor multiplication (composition).
impl<T: Float> Mul for Rotor<T> {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        self.compose(other)
    }
}

impl<T: Float> Mul<T> for Rotor<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self::new(self.s * scalar, self.xy * scalar)
    }
}

/// Rotor applied to vector.
impl<T: Float> Mul<Vector<T>> for Rotor<T> {
    type Output = Vector<T>;

    #[inline]
    fn mul(self, v: Vector<T>) -> Self::Output {
        self.rotate(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::specialized::euclidean::dim2::arbitrary::{NonZeroVector, UnitRotor, UnitVector};
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    // ========================================================================
    // Vector tests
    // ========================================================================

    proptest! {
        #[test]
        fn vec2_add_commutative(a in any::<Vector<f64>>(), b in any::<Vector<f64>>()) {
            let ab = a + b;
            let ba = b + a;
            prop_assert!(abs_diff_eq!(ab, ba, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec2_dot_commutative(a in any::<Vector<f64>>(), b in any::<Vector<f64>>()) {
            prop_assert!(abs_diff_eq!(a.dot(b), b.dot(a), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec2_wedge_anticommutative(a in any::<Vector<f64>>(), b in any::<Vector<f64>>()) {
            let ab = a.wedge(b);
            let ba = b.wedge(a);
            prop_assert!(abs_diff_eq!(ab.0, -ba.0, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec2_normalized_has_unit_length(v in any::<NonZeroVector<f64>>()) {
            let n = v.normalized();
            prop_assert!(abs_diff_eq!(n.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec2_perp_is_perpendicular(v in any::<NonZeroVector<f64>>()) {
            let p = v.perp();
            prop_assert!(abs_diff_eq!(v.dot(p), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Rotor tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor2_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
            let rotated = r.rotate(v);
            prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor2_composition(
            r1 in any::<UnitRotor<f64>>(),
            r2 in any::<UnitRotor<f64>>(),
            v in any::<Vector<f64>>()
        ) {
            let sequential = r2.rotate(r1.rotate(v));
            let composed = r2.compose(*r1).rotate(v);
            prop_assert!(abs_diff_eq!(sequential, composed, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor2_inverse(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
            let roundtrip = r.inverse().rotate(r.rotate(v));
            prop_assert!(abs_diff_eq!(roundtrip, v, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor2_from_vectors(a in any::<UnitVector<f64>>(), b in any::<UnitVector<f64>>()) {
            let r = Rotor::from_vectors(*a, *b);
            let rotated = r.rotate(*a);
            prop_assert!(abs_diff_eq!(rotated, *b, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Specific rotation tests
    // ========================================================================

    #[test]
    fn rotor2_90_deg_rotation() {
        let rotor = Rotor::from_angle(FRAC_PI_2);
        let v = Vector::unit_x();
        let rotated = rotor.rotate(v);

        assert!(abs_diff_eq!(
            rotated,
            Vector::unit_y(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn rotor2_180_deg_rotation() {
        let rotor = Rotor::from_angle(PI);
        let v = Vector::unit_x();
        let rotated = rotor.rotate(v);

        assert!(abs_diff_eq!(
            rotated,
            -Vector::unit_x(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn rotor2_45_plus_45_equals_90() {
        let r45 = Rotor::from_angle(FRAC_PI_4);
        let r90 = r45.compose(r45);
        let v = Vector::unit_x();
        let rotated = r90.rotate(v);

        assert!(abs_diff_eq!(
            rotated,
            Vector::unit_y(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn rotor2_angle_roundtrip() {
        let angle = 1.234;
        let rotor = Rotor::from_angle(angle);
        assert!(abs_diff_eq!(
            rotor.angle(),
            angle,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }
}
