//! Clifford - A Rust library for Geometric Algebra (Clifford Algebra).
//!
//! This library provides types and operations for working with geometric algebra,
//! a mathematical framework that unifies and extends linear algebra, complex numbers,
//! quaternions, and more.
//!
//! # Overview
//!
//! Geometric Algebra (GA), also known as Clifford Algebra, provides a unified
//! mathematical language for geometry. It extends traditional vector algebra with
//! the geometric product, which combines the dot product and wedge product into
//! a single operation.
//!
//! # Module Structure
//!
//! - [`algebra`] - Core algebraic types ([`Multivector`](algebra::Multivector))
//! - [`scalar`] - Floating-point scalar type abstraction
//! - [`signature`] - Metric signatures defining the algebra
//! - [`basis`] - Basis blade representation and utilities
//! - [`specialized`] - Optimized types for 2D and 3D Euclidean geometry
//!
//! # Features
//!
//! - `proptest-support` - Enable [`proptest`](https://docs.rs/proptest) strategies
//!   for property-based testing
//! - `serde` - Enable serialization/deserialization via [`serde`](https://docs.rs/serde)
//! - `nalgebra-0_33` - Enable conversions to/from [`nalgebra`](https://docs.rs/nalgebra) 0.33.x types
//! - `nalgebra-0_34` - Enable conversions to/from [`nalgebra`](https://docs.rs/nalgebra) 0.34.x types
//!
//! Note: `nalgebra-0_33` and `nalgebra-0_34` are mutually exclusive.
//!
//! # Getting Started
//!
//! ## Generic Multivector API
//!
//! The generic API works with any metric signature:
//!
//! ```
//! use clifford::algebra::Multivector;
//! use clifford::signature::Euclidean3;
//!
//! // Create two vectors
//! let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 0.0]);
//! let b: Multivector<f64, Euclidean3> = Multivector::vector(&[0.0, 1.0, 0.0]);
//!
//! // Geometric product combines dot and wedge products
//! let ab = &a * &b;
//!
//! // Dot product is the scalar part: a·b = 1*0 + 2*1 + 0*0 = 2
//! assert!((ab.scalar_part() - 2.0).abs() < 1e-10);
//! ```
//!
//! ## Specialized 2D/3D Types
//!
//! For common 2D and 3D Euclidean geometry, use the optimized specialized types:
//!
//! ```
//! use clifford::specialized::euclidean::{dim2, dim3};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // 2D: Rotate a vector by 90 degrees
//! let v = dim2::Vector::new(1.0, 0.0);
//! let rotor = dim2::Rotor::from_angle(FRAC_PI_2);
//! let rotated = rotor.rotate(v);
//! assert!((rotated.x).abs() < 1e-10);
//! assert!((rotated.y - 1.0).abs() < 1e-10);
//!
//! // 3D: Rotate around an axis
//! let v = dim3::Vector::new(1.0, 0.0, 0.0);
//! let plane = dim3::Bivector::unit_xy(); // rotation in xy-plane (around z-axis)
//! let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_2, plane);
//! let rotated = rotor.rotate(v);
//! assert!((rotated.x).abs() < 1e-10);
//! assert!((rotated.y - 1.0).abs() < 1e-10);
//! assert!((rotated.z).abs() < 1e-10);
//! ```

// Ensure nalgebra feature flags are mutually exclusive
#[cfg(any(
    all(feature = "nalgebra-0_32", feature = "nalgebra-0_33"),
    all(feature = "nalgebra-0_32", feature = "nalgebra-0_34"),
    all(feature = "nalgebra-0_33", feature = "nalgebra-0_34"),
))]
compile_error!(
    "Features `nalgebra-0_32`, `nalgebra-0_33`, and `nalgebra-0_34` are mutually exclusive. Enable only one."
);

pub mod algebra;
pub mod basis;
pub mod prelude;
pub mod scalar;
pub mod signature;
pub mod specialized;

/// Test utilities available only during testing.
#[cfg(test)]
pub(crate) mod test_utils {
    /// Standard epsilon for absolute difference comparisons in tests.
    ///
    /// Use this constant instead of magic numbers like `1e-10` or `1e-9`.
    /// This value is chosen to be strict enough for most operations while
    /// allowing for reasonable floating-point accumulation in compound operations.
    pub const ABS_DIFF_EQ_EPS: f64 = 1e-10;
}
