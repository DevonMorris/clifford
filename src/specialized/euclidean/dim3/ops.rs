//! Operator implementations for 3D geometric algebra types.

use super::types::{Bivec3, Even3, Rotor3, Trivec3, Vec3};
use crate::scalar::Float;
use std::ops::{Add, BitXor, Mul, Neg, Sub};

// ============================================================================
// Vec3 operations
// ============================================================================

impl<T: Float> Neg for Vec3<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl<T: Float> Add for Vec3<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl<T: Float> Sub for Vec3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl<T: Float> Mul<T> for Vec3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Mul<Vec3<f32>> for f32 {
    type Output = Vec3<f32>;

    #[inline]
    fn mul(self, v: Vec3<f32>) -> Self::Output {
        Vec3::new(self * v.x, self * v.y, self * v.z)
    }
}

impl Mul<Vec3<f64>> for f64 {
    type Output = Vec3<f64>;

    #[inline]
    fn mul(self, v: Vec3<f64>) -> Self::Output {
        Vec3::new(self * v.x, self * v.y, self * v.z)
    }
}

/// Wedge product operator for vectors: `a ^ b = a ∧ b`.
impl<T: Float> BitXor for Vec3<T> {
    type Output = Bivec3<T>;

    #[inline]
    fn bitxor(self, other: Self) -> Self::Output {
        self.wedge(other)
    }
}

// ============================================================================
// Bivec3 operations
// ============================================================================

impl<T: Float> Neg for Bivec3<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.xy, -self.xz, -self.yz)
    }
}

impl<T: Float> Add for Bivec3<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.xy + other.xy, self.xz + other.xz, self.yz + other.yz)
    }
}

impl<T: Float> Sub for Bivec3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.xy - other.xy, self.xz - other.xz, self.yz - other.yz)
    }
}

impl<T: Float> Mul<T> for Bivec3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self::new(self.xy * scalar, self.xz * scalar, self.yz * scalar)
    }
}

impl Mul<Bivec3<f32>> for f32 {
    type Output = Bivec3<f32>;

    #[inline]
    fn mul(self, b: Bivec3<f32>) -> Self::Output {
        Bivec3::new(self * b.xy, self * b.xz, self * b.yz)
    }
}

impl Mul<Bivec3<f64>> for f64 {
    type Output = Bivec3<f64>;

    #[inline]
    fn mul(self, b: Bivec3<f64>) -> Self::Output {
        Bivec3::new(self * b.xy, self * b.xz, self * b.yz)
    }
}

// ============================================================================
// Trivec3 operations
// ============================================================================

impl<T: Float> Neg for Trivec3<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl<T: Float> Add for Trivec3<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0)
    }
}

impl<T: Float> Sub for Trivec3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self(self.0 - other.0)
    }
}

impl<T: Float> Mul<T> for Trivec3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self(self.0 * scalar)
    }
}

// ============================================================================
// Rotor3 operations
// ============================================================================

impl<T: Float> Neg for Rotor3<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.s, -self.b)
    }
}

impl<T: Float> Add for Rotor3<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.s + other.s, self.b + other.b)
    }
}

impl<T: Float> Sub for Rotor3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.s - other.s, self.b - other.b)
    }
}

/// Rotor multiplication (composition).
impl<T: Float> Mul for Rotor3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        self.compose(other)
    }
}

impl<T: Float> Mul<T> for Rotor3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Self::new(self.s * scalar, self.b * scalar)
    }
}

/// Rotor applied to vector.
impl<T: Float> Mul<Vec3<T>> for Rotor3<T> {
    type Output = Vec3<T>;

    #[inline]
    fn mul(self, v: Vec3<T>) -> Self::Output {
        self.rotate(v)
    }
}

// ============================================================================
// Even3 operations
// ============================================================================

impl<T: Float> Neg for Even3<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.s, -self.b)
    }
}

impl<T: Float> Add for Even3<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.s + other.s, self.b + other.b)
    }
}

impl<T: Float> Sub for Even3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.s - other.s, self.b - other.b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::specialized::euclidean::dim3::arbitrary::{NonZeroVec3, UnitRotor3, UnitVec3};
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    // ========================================================================
    // Vec3 tests
    // ========================================================================

    proptest! {
        #[test]
        fn vec3_add_commutative(a in any::<Vec3<f64>>(), b in any::<Vec3<f64>>()) {
            let ab = a + b;
            let ba = b + a;
            prop_assert!(abs_diff_eq!(ab, ba, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec3_dot_commutative(a in any::<Vec3<f64>>(), b in any::<Vec3<f64>>()) {
            prop_assert!(abs_diff_eq!(a.dot(b), b.dot(a), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec3_wedge_anticommutative(a in any::<Vec3<f64>>(), b in any::<Vec3<f64>>()) {
            let ab = a.wedge(b);
            let ba = b.wedge(a);
            prop_assert!(abs_diff_eq!(ab, -ba, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec3_cross_anticommutative(a in any::<Vec3<f64>>(), b in any::<Vec3<f64>>()) {
            let ab = a.cross(b);
            let ba = b.cross(a);
            prop_assert!(abs_diff_eq!(ab, -ba, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn vec3_normalized_has_unit_length(v in any::<NonZeroVec3<f64>>()) {
            let n = v.normalized();
            prop_assert!(abs_diff_eq!(n.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Rotor tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor_preserves_norm(r in any::<UnitRotor3<f64>>(), v in any::<Vec3<f64>>()) {
            let rotated = r.rotate(v);
            prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_composition(
            r1 in any::<UnitRotor3<f64>>(),
            r2 in any::<UnitRotor3<f64>>(),
            v in any::<Vec3<f64>>()
        ) {
            let sequential = r2.rotate(r1.rotate(v));
            let composed = r2.compose(*r1).rotate(v);
            prop_assert!(abs_diff_eq!(sequential, composed, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_inverse(r in any::<UnitRotor3<f64>>(), v in any::<Vec3<f64>>()) {
            let roundtrip = r.inverse().rotate(r.rotate(v));
            prop_assert!(abs_diff_eq!(roundtrip, v, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_identity_is_noop(v in any::<Vec3<f64>>()) {
            let rotated = Rotor3::identity().rotate(v);
            prop_assert!(abs_diff_eq!(rotated, v, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn rotor_from_vectors(a in any::<UnitVec3<f64>>(), b in any::<UnitVec3<f64>>()) {
            let r = Rotor3::from_vectors(*a, *b);
            let rotated = r.rotate(*a);
            prop_assert!(abs_diff_eq!(rotated, *b, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Specific rotation tests
    // ========================================================================

    #[test]
    fn rotor_90_deg_xy_rotation() {
        let rotor = Rotor3::from_angle_plane(FRAC_PI_2, Bivec3::unit_xy());
        let v = Vec3::unit_x();
        let rotated = rotor.rotate(v);

        assert!(abs_diff_eq!(
            rotated,
            Vec3::unit_y(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn rotor_180_deg_rotation() {
        let rotor = Rotor3::from_angle_plane(PI, Bivec3::unit_xy());
        let v = Vec3::unit_x();
        let rotated = rotor.rotate(v);

        assert!(abs_diff_eq!(
            rotated,
            -Vec3::unit_x(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn rotor_45_plus_45_equals_90() {
        let r45 = Rotor3::from_angle_plane(FRAC_PI_4, Bivec3::unit_xy());
        let r90 = r45.compose(r45);
        let v = Vec3::unit_x();
        let rotated = r90.rotate(v);

        assert!(abs_diff_eq!(
            rotated,
            Vec3::unit_y(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn rotor_slerp_endpoints() {
        let r1 = Rotor3::identity();
        let r2 = Rotor3::from_angle_plane(FRAC_PI_2, Bivec3::unit_xy());

        let at_0 = r1.slerp(r2, 0.0);
        let at_1 = r1.slerp(r2, 1.0);

        assert!(abs_diff_eq!(at_0.s, r1.s, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(at_1.s, r2.s, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn rotor_slerp_midpoint() {
        let r1 = Rotor3::identity();
        let r2 = Rotor3::from_angle_plane(FRAC_PI_2, Bivec3::unit_xy());

        let mid = r1.slerp(r2, 0.5);
        let v = Vec3::unit_x();
        let rotated = mid.rotate(v);

        // 45° rotation should give (√2/2, √2/2, 0)
        let expected = std::f64::consts::FRAC_1_SQRT_2;
        assert!(abs_diff_eq!(rotated.x, expected, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rotated.y, expected, epsilon = ABS_DIFF_EQ_EPS));
    }
}
