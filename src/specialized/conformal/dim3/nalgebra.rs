//! nalgebra interoperability for 3D CGA types.
//!
//! This module provides bidirectional conversions between clifford's 3D CGA types
//! and nalgebra's equivalent types.
//!
//! Enable with feature `nalgebra-0_32`, `nalgebra-0_33`, or `nalgebra-0_34`.
//!
//! # Conversions
//!
//! | clifford | nalgebra | Notes |
//! |----------|----------|-------|
//! | [`Point<T>`] | [`na::Point3<T>`] | Euclidean coordinates |
//! | [`Sphere<T>`] | `(na::Point3<T>, T)` | Center and radius |

#[cfg(feature = "nalgebra-0_32")]
use nalgebra_0_32 as na;
#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use crate::scalar::Float;

use super::{Point, Sphere};

// ============================================================================
// Point <-> Point3
// ============================================================================

impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T> {
    /// Converts a nalgebra 3D point to a CGA conformal point.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::conformal::dim3::Point;
    /// use nalgebra::Point3;
    ///
    /// let na_p = Point3::new(1.0, 2.0, 3.0);
    /// let p: Point<f64> = na_p.into();
    /// assert!((p.x() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(p: na::Point3<T>) -> Self {
        Point::new(p.x, p.y, p.z)
    }
}

impl<T: Float + na::Scalar> From<Point<T>> for na::Point3<T> {
    /// Converts a CGA conformal point to a nalgebra 3D point.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::conformal::dim3::Point;
    /// use nalgebra::Point3;
    ///
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// let na_p: Point3<f64> = p.into();
    /// assert!((na_p.x - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(p: Point<T>) -> Self {
        na::Point3::new(p.x(), p.y(), p.z())
    }
}

// ============================================================================
// Sphere <-> (Point3, T)
// ============================================================================

impl<T: Float + na::Scalar> From<Sphere<T>> for (na::Point3<T>, T) {
    /// Converts a CGA sphere to (center, radius) tuple.
    ///
    /// # Panics
    ///
    /// Panics if the sphere is imaginary (negative squared radius).
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::conformal::dim3::Sphere;
    /// use nalgebra::Point3;
    ///
    /// let sphere = Sphere::from_center_radius(1.0, 2.0, 3.0, 5.0);
    /// let (center, radius): (Point3<f64>, f64) = sphere.into();
    /// assert!((center.x - 1.0).abs() < 1e-10);
    /// assert!((radius - 5.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from(s: Sphere<T>) -> Self {
        let center = s.center();
        (
            na::Point3::new(center.x(), center.y(), center.z()),
            s.radius(),
        )
    }
}

impl<T: Float + na::Scalar> From<(na::Point3<T>, T)> for Sphere<T> {
    /// Creates a CGA sphere from (center, radius) tuple.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::specialized::conformal::dim3::Sphere;
    /// use nalgebra::Point3;
    ///
    /// let center = Point3::new(1.0, 2.0, 3.0);
    /// let radius = 5.0;
    /// let sphere: Sphere<f64> = (center, radius).into();
    /// assert!((sphere.radius() - 5.0).abs() < 1e-10);
    /// ```
    #[inline]
    fn from((center, radius): (na::Point3<T>, T)) -> Self {
        Sphere::from_center_radius(center.x, center.y, center.z, radius)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn point_roundtrip(x in -100.0..100.0, y in -100.0..100.0, z in -100.0..100.0) {
            let p = Point::new(x, y, z);
            let na_p: na::Point3<f64> = p.into();
            let back: Point<f64> = na_p.into();

            prop_assert!(abs_diff_eq!(p.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), back.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn sphere_roundtrip(
            cx in -100.0..100.0,
            cy in -100.0..100.0,
            cz in -100.0..100.0,
            r in 0.01f64..100.0,
        ) {
            let sphere = Sphere::from_center_radius(cx, cy, cz, r);
            let (na_center, na_radius): (na::Point3<f64>, f64) = sphere.into();
            let back: Sphere<f64> = (na_center, na_radius).into();

            prop_assert!(abs_diff_eq!(sphere.cx(), back.cx(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(sphere.cy(), back.cy(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(sphere.cz(), back.cz(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(sphere.radius(), back.radius(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn point_distance_matches_nalgebra(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);

            let na_p1 = na::Point3::new(x1, y1, z1);
            let na_p2 = na::Point3::new(x2, y2, z2);

            let cga_dist = p1.distance(&p2);
            let na_dist = na::distance(&na_p1, &na_p2);

            prop_assert!(abs_diff_eq!(cga_dist, na_dist, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
