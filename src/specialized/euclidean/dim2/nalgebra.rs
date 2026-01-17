//! nalgebra interoperability for 2D types.
//!
//! This module provides bidirectional conversions between clifford's 2D types
//! and nalgebra's equivalent types.
//!
//! Enable with feature `nalgebra-0_33` or `nalgebra-0_34`.
//!
//! # Conversions
//!
//! | clifford | nalgebra |
//! |----------|----------|
//! | [`Vector<T>`] | [`na::Vector2<T>`] |
//! | [`Rotor<T>`] | [`na::Rotation2<T>`] |

#[cfg(feature = "nalgebra-0_32")]
use nalgebra_0_32 as na;
#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use crate::scalar::Float;

use super::{Rotor, Vector};

// ============================================================================
// Vector <-> Vector2
// ============================================================================

impl<T: Float + na::Scalar> From<na::Vector2<T>> for Vector<T> {
    /// Converts a nalgebra 2D vector to a clifford 2D vector.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Vector;
    /// use nalgebra::Vector2;
    ///
    /// let na_v = Vector2::new(1.0, 2.0);
    /// let v: Vector<f64> = na_v.into();
    /// assert_eq!(v.x(), 1.0);
    /// assert_eq!(v.y(), 2.0);
    /// ```
    #[inline]
    fn from(v: na::Vector2<T>) -> Self {
        Vector::new(v.x, v.y)
    }
}

impl<T: Float + na::Scalar> From<Vector<T>> for na::Vector2<T> {
    /// Converts a clifford 2D vector to a nalgebra 2D vector.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Vector;
    /// use nalgebra::Vector2;
    ///
    /// let v = Vector::new(1.0, 2.0);
    /// let na_v: Vector2<f64> = v.into();
    /// assert_eq!(na_v.x, 1.0);
    /// assert_eq!(na_v.y, 2.0);
    /// ```
    #[inline]
    fn from(v: Vector<T>) -> Self {
        na::Vector2::new(v.x(), v.y())
    }
}

// ============================================================================
// Rotor <-> Rotation2
// ============================================================================

impl<T: Float + na::RealField> From<na::Rotation2<T>> for Rotor<T> {
    /// Converts a nalgebra 2D rotation to a clifford 2D rotor.
    ///
    /// # Mathematical Correspondence
    ///
    /// A 2D rotation by angle θ corresponds to a rotor:
    /// - `R = cos(θ/2) + sin(θ/2)·e₁₂`
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use nalgebra::Rotation2;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let rotation = Rotation2::new(FRAC_PI_2); // 90° rotation
    /// let rotor: Rotor<f64> = rotation.into();
    /// ```
    #[inline]
    fn from(rotation: na::Rotation2<T>) -> Self {
        Rotor::from_angle(rotation.angle())
    }
}

impl<T: Float + na::RealField> From<Rotor<T>> for na::Rotation2<T> {
    /// Converts a clifford 2D rotor to a nalgebra 2D rotation.
    ///
    /// # Mathematical Correspondence
    ///
    /// A rotor `R = s + xy·e₁₂ = cos(θ/2) + sin(θ/2)·e₁₂` corresponds to
    /// a nalgebra rotation by angle `θ`.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::euclidean::dim2::Rotor;
    /// use nalgebra::Rotation2;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let rotor = Rotor::from_angle(FRAC_PI_2); // 90° rotation
    /// let rotation: Rotation2<f64> = rotor.into();
    /// ```
    #[inline]
    fn from(rotor: Rotor<T>) -> Self {
        na::Rotation2::new(rotor.angle())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::Transform;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use crate::wrappers::Unit;
    use approx::relative_eq;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn vector_roundtrip(v in any::<Vector<f64>>()) {
            let na_v: na::Vector2<f64> = v.into();
            let back: Vector<f64> = na_v.into();
            prop_assert!(relative_eq!(v, back, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn rotor_angle_roundtrip(angle in -std::f64::consts::PI..std::f64::consts::PI) {
            let rotor = Rotor::from_angle(angle);
            let rotation: na::Rotation2<f64> = rotor.into();
            let back: Rotor<f64> = rotation.into();

            // Verify angle is preserved
            prop_assert!(relative_eq!(rotor.angle(), back.angle(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn rotor_rotation_equivalence(
            r in any::<Unit<Rotor<f64>>>(),
            v in any::<Vector<f64>>(),
        ) {
            let na_v: na::Vector2<f64> = v.into();

            // Rotate with clifford rotor
            let rotated_ga = r.as_inner().transform(&v);

            // Rotate with nalgebra Rotation2
            let rotation: na::Rotation2<f64> = r.into_inner().into();
            let rotated_na = rotation * na_v;

            let rotated_back: Vector<f64> = rotated_na.into();
            prop_assert!(relative_eq!(rotated_ga, rotated_back, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }
    }
}
