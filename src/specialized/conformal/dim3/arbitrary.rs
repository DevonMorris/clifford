//! Proptest `Arbitrary` implementations for 3D CGA types.
//!
//! This module provides [`Arbitrary`] implementations for property-based testing
//! with [`mod@proptest`]. It is available when either running tests or when the
//! `proptest-support` feature is enabled.
//!
//! # Example
//!
//! ```
//! use clifford::specialized::conformal::dim3::Point;
//! use proptest::prelude::*;
//!
//! proptest! {
//!     #[test]
//!     fn point_is_always_null(p in any::<Point<f64>>()) {
//!         prop_assert!(p.is_null(1e-9));
//!     }
//! }
//! ```

use super::{Circle, Dilator, Dipole, FlatPoint, Line, Plane, Point, Rotor, Sphere, Translator};
use crate::scalar::Float;
use crate::specialized::euclidean::dim3::Bivector;
use core::fmt::Debug;
use proptest::arbitrary::Arbitrary;
use proptest::prelude::*;
use proptest::strategy::BoxedStrategy;

// ============================================================================
// Point Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Point<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| Point::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}

// ============================================================================
// FlatPoint Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for FlatPoint<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| FlatPoint::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}

// ============================================================================
// Sphere Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Sphere<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            0.01f64..100.0,
        )
            .prop_map(|(cx, cy, cz, r)| {
                Sphere::from_center_radius(
                    T::from_f64(cx),
                    T::from_f64(cy),
                    T::from_f64(cz),
                    T::from_f64(r),
                )
            })
            .boxed()
    }
}

// ============================================================================
// Plane Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Plane<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate a point on the plane and a normal direction
        (
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            -1.0f64..1.0,
            -1.0f64..1.0,
            -1.0f64..1.0,
        )
            .prop_filter_map("non-zero normal", |(px, py, pz, nx, ny, nz)| {
                let norm_sq = nx * nx + ny * ny + nz * nz;
                if norm_sq < 0.01 {
                    return None;
                }
                let norm = norm_sq.sqrt();
                let nx = nx / norm;
                let ny = ny / norm;
                let nz = nz / norm;
                // d = -(n · p)
                let d = -(nx * px + ny * py + nz * pz);
                Some(Plane::from_normal_distance(
                    T::from_f64(nx),
                    T::from_f64(ny),
                    T::from_f64(nz),
                    T::from_f64(d),
                ))
            })
            .boxed()
    }
}

// ============================================================================
// Circle Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Circle<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            // Center
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            // Radius (positive)
            0.01f64..100.0,
            // Normal direction
            -1.0f64..1.0,
            -1.0f64..1.0,
            -1.0f64..1.0,
        )
            .prop_filter_map("non-zero normal", |(cx, cy, cz, r, nx, ny, nz)| {
                let norm_sq = nx * nx + ny * ny + nz * nz;
                if norm_sq < 0.01 {
                    return None;
                }
                let norm = norm_sq.sqrt();
                Some(Circle::from_center_radius_normal(
                    &crate::specialized::euclidean::dim3::Vector::new(
                        T::from_f64(cx),
                        T::from_f64(cy),
                        T::from_f64(cz),
                    ),
                    T::from_f64(r),
                    &crate::specialized::euclidean::dim3::Vector::new(
                        T::from_f64(nx / norm),
                        T::from_f64(ny / norm),
                        T::from_f64(nz / norm),
                    ),
                ))
            })
            .boxed()
    }
}

// ============================================================================
// Line Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Line<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            // Point on line
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            // Direction
            -1.0f64..1.0,
            -1.0f64..1.0,
            -1.0f64..1.0,
        )
            .prop_filter_map("non-zero direction", |(px, py, pz, dx, dy, dz)| {
                let norm_sq = dx * dx + dy * dy + dz * dz;
                if norm_sq < 0.01 {
                    return None;
                }
                Some(Line::from_point_direction(
                    &Point::new(T::from_f64(px), T::from_f64(py), T::from_f64(pz)),
                    &crate::specialized::euclidean::dim3::Vector::new(
                        T::from_f64(dx),
                        T::from_f64(dy),
                        T::from_f64(dz),
                    ),
                ))
            })
            .boxed()
    }
}

// ============================================================================
// Dipole Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Dipole<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            // First point
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
            // Second point
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
        )
            .prop_filter_map("distinct points", |(x1, y1, z1, x2, y2, z2)| {
                let dist_sq = (x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2);
                if dist_sq < 0.01 {
                    return None;
                }
                Some(Dipole::from_two_points(
                    &Point::new(T::from_f64(x1), T::from_f64(y1), T::from_f64(z1)),
                    &Point::new(T::from_f64(x2), T::from_f64(y2), T::from_f64(z2)),
                ))
            })
            .boxed()
    }
}

// ============================================================================
// Translator Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Translator<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(tx, ty, tz)| {
                Translator::new(T::from_f64(tx), T::from_f64(ty), T::from_f64(tz))
            })
            .boxed()
    }
}

// ============================================================================
// Rotor Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Rotor<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            // Angle
            -std::f64::consts::PI..std::f64::consts::PI,
            // Axis direction (will be normalized)
            -1.0f64..1.0,
            -1.0f64..1.0,
            -1.0f64..1.0,
        )
            .prop_filter_map("non-zero axis", |(angle, ax, ay, az)| {
                let norm_sq = ax * ax + ay * ay + az * az;
                if norm_sq < 0.01 {
                    return None;
                }
                let norm = norm_sq.sqrt();
                let plane = Bivector::new(
                    T::from_f64(az / norm),
                    T::from_f64(-ay / norm),
                    T::from_f64(ax / norm),
                );
                Some(Rotor::from_angle_plane(T::from_f64(angle), plane))
            })
            .boxed()
    }
}

// ============================================================================
// Dilator Arbitrary implementation
// ============================================================================

impl<T: Float + Debug> Arbitrary for Dilator<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            // Scale factor (positive)
            0.1f64..10.0,
            // Center
            -100.0f64..100.0,
            -100.0f64..100.0,
            -100.0f64..100.0,
        )
            .prop_map(|(scale, cx, cy, cz)| {
                Dilator::from_center_scale(
                    &Point::new(T::from_f64(cx), T::from_f64(cy), T::from_f64(cz)),
                    T::from_f64(scale),
                )
            })
            .boxed()
    }
}
