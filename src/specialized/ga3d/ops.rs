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
    use crate::specialized::ga3d::arbitrary::{NonZeroVec3, UnitRotor3, UnitVec3};
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
            prop_assert!((ab.x - ba.x).abs() < 1e-10);
            prop_assert!((ab.y - ba.y).abs() < 1e-10);
            prop_assert!((ab.z - ba.z).abs() < 1e-10);
        }

        #[test]
        fn vec3_dot_commutative(a in any::<Vec3<f64>>(), b in any::<Vec3<f64>>()) {
            prop_assert!((a.dot(b) - b.dot(a)).abs() < 1e-10);
        }

        #[test]
        fn vec3_wedge_anticommutative(a in any::<Vec3<f64>>(), b in any::<Vec3<f64>>()) {
            let ab = a.wedge(b);
            let ba = b.wedge(a);
            prop_assert!((ab.xy + ba.xy).abs() < 1e-10);
            prop_assert!((ab.xz + ba.xz).abs() < 1e-10);
            prop_assert!((ab.yz + ba.yz).abs() < 1e-10);
        }

        #[test]
        fn vec3_cross_anticommutative(a in any::<Vec3<f64>>(), b in any::<Vec3<f64>>()) {
            let ab = a.cross(b);
            let ba = b.cross(a);
            prop_assert!((ab.x + ba.x).abs() < 1e-10);
            prop_assert!((ab.y + ba.y).abs() < 1e-10);
            prop_assert!((ab.z + ba.z).abs() < 1e-10);
        }

        #[test]
        fn vec3_normalized_has_unit_length(v in any::<NonZeroVec3>()) {
            let n = v.0.normalized();
            prop_assert!((n.norm() - 1.0).abs() < 1e-10);
        }
    }

    // ========================================================================
    // Rotor tests
    // ========================================================================

    proptest! {
        #[test]
        fn rotor_preserves_norm(r in any::<UnitRotor3>(), v in any::<Vec3<f64>>()) {
            let rotated = r.0.rotate(v);
            prop_assert!((v.norm() - rotated.norm()).abs() < 1e-9);
        }

        #[test]
        fn rotor_composition(
            r1 in any::<UnitRotor3>(),
            r2 in any::<UnitRotor3>(),
            v in any::<Vec3<f64>>()
        ) {
            let sequential = r2.0.rotate(r1.0.rotate(v));
            let composed = r2.0.compose(r1.0).rotate(v);
            prop_assert!((sequential.x - composed.x).abs() < 1e-8);
            prop_assert!((sequential.y - composed.y).abs() < 1e-8);
            prop_assert!((sequential.z - composed.z).abs() < 1e-8);
        }

        #[test]
        fn rotor_inverse(r in any::<UnitRotor3>(), v in any::<Vec3<f64>>()) {
            let roundtrip = r.0.inverse().rotate(r.0.rotate(v));
            prop_assert!((roundtrip.x - v.x).abs() < 1e-8);
            prop_assert!((roundtrip.y - v.y).abs() < 1e-8);
            prop_assert!((roundtrip.z - v.z).abs() < 1e-8);
        }

        #[test]
        fn rotor_identity_is_noop(v in any::<Vec3<f64>>()) {
            let rotated = Rotor3::identity().rotate(v);
            prop_assert!((rotated.x - v.x).abs() < 1e-10);
            prop_assert!((rotated.y - v.y).abs() < 1e-10);
            prop_assert!((rotated.z - v.z).abs() < 1e-10);
        }

        #[test]
        fn rotor_from_vectors(a in any::<UnitVec3>(), b in any::<UnitVec3>()) {
            let r = Rotor3::from_vectors(a.0, b.0);
            let rotated = r.rotate(a.0);
            prop_assert!((rotated.x - b.0.x).abs() < 1e-8);
            prop_assert!((rotated.y - b.0.y).abs() < 1e-8);
            prop_assert!((rotated.z - b.0.z).abs() < 1e-8);
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

        assert!(rotated.x.abs() < 1e-10);
        assert!((rotated.y - 1.0).abs() < 1e-10);
        assert!(rotated.z.abs() < 1e-10);
    }

    #[test]
    fn rotor_180_deg_rotation() {
        let rotor = Rotor3::from_angle_plane(PI, Bivec3::unit_xy());
        let v = Vec3::unit_x();
        let rotated = rotor.rotate(v);

        assert!((rotated.x + 1.0).abs() < 1e-10);
        assert!(rotated.y.abs() < 1e-10);
        assert!(rotated.z.abs() < 1e-10);
    }

    #[test]
    fn rotor_45_plus_45_equals_90() {
        let r45 = Rotor3::from_angle_plane(FRAC_PI_4, Bivec3::unit_xy());
        let r90 = r45.compose(r45);
        let v = Vec3::unit_x();
        let rotated = r90.rotate(v);

        assert!(rotated.x.abs() < 1e-10);
        assert!((rotated.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rotor_slerp_endpoints() {
        let r1 = Rotor3::identity();
        let r2 = Rotor3::from_angle_plane(FRAC_PI_2, Bivec3::unit_xy());

        let at_0 = r1.slerp(r2, 0.0);
        let at_1 = r1.slerp(r2, 1.0);

        assert!((at_0.s - r1.s).abs() < 1e-10);
        assert!((at_1.s - r2.s).abs() < 1e-10);
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
        assert!((rotated.x - expected).abs() < 1e-10);
        assert!((rotated.y - expected).abs() < 1e-10);
    }
}
