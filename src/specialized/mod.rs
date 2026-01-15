//! Specialized types for common geometric algebras.
//!
//! This module provides optimized, ergonomic types organized by algebra type
//! and dimension. Currently supported:
//!
//! - **[`euclidean`]**: Euclidean geometric algebra (EGA) for 2D and 3D
//! - **[`elliptic`]**: Elliptic projective geometry for 2D
//! - **[`projective`]**: Projective geometric algebra (PGA) for 2D and 3D
//! - **[`complex`]**: Complex numbers Cl(0,1,0)
//! - **[`dual`]**: Dual numbers Cl(0,0,1) for automatic differentiation
//! - **[`dualquat`]**: Dual quaternions Cl(0,2,1) for rigid transformations
//! - **[`hyperbolic`]**: Hyperbolic algebras: numbers Cl(1,0,0) and plane Cl(2,1,0)
//! - **[`minkowski`]**: Minkowski spacetime algebras (indefinite signature)
//! - **[`quaternion`]**: Quaternions Cl(0,2,0)
//!
//! # Type Organization
//!
//! ```text
//! specialized/
//!   euclidean/
//!     dim2/   - 2D Euclidean GA: Vector, Bivector, Rotor
//!     dim3/   - 3D Euclidean GA: Vector, Bivector, Trivector, Rotor
//!   elliptic/
//!     dim2/   - 2D Elliptic: Point, Line, Rotor
//!   projective/
//!     dim2/   - 2D PGA: Point, Line, Motor
//!     dim3/   - 3D PGA: Point, Line, Plane, Motor
//!   complex/    - Complex numbers: Scalar, ImagUnit, Complex
//!   dual/       - Dual numbers: Scalar, DualUnit, Dual
//!   dualquat/   - Dual quaternions: Scalar, Vector, Bivector, DualQuaternion
//!   hyperbolic/ - Hyperbolic numbers: Scalar, HypUnit, Hyperbolic
//!     dim2/   - 2D Hyperbolic plane: Point, Line, Rotor
//!   minkowski/
//!     dim2/   - 2D Minkowski: Vector, Bivector, Eventor, Spacetime
//!     dim3/   - 3D Minkowski: Vector, Bivector, Trivector, Eventor
//!   quaternion/ - Quaternions: Scalar, Imaginary, Bivector, Quaternion
//! ```
//!
//! # Benefits of Specialized Types
//!
//! - **Named fields**: Access components as `x`, `y`, `z` instead of indices
//! - **Type safety**: Vectors, bivectors, and rotors are distinct types
//! - **Ergonomics**: Common operations with intuitive method names
//!
//! # Choosing Generic vs Specialized
//!
//! Use **specialized types** when:
//! - Working exclusively in 2D or 3D Euclidean space
//! - You want named field access (`v.x`, `v.y`, `v.z`)
//! - Type safety between grades is important
//!
//! Use **generic [`Multivector`](crate::algebra::Multivector)** when:
//! - Working with arbitrary signatures (PGA, CGA, etc.)
//! - You need to mix grades freely
//! - Writing dimension-generic algorithms
//!
//! # Example
//!
//! ```ignore
//! use clifford::ops::Transform;
//! use clifford::specialized::euclidean::dim3::{Vector, Bivector, Rotor};
//!
//! // Create a vector
//! let v = Vector::new(1.0, 0.0, 0.0);
//!
//! // Create a 90° rotation in the xy-plane
//! let rotor = Rotor::from_angle_plane(
//!     std::f64::consts::FRAC_PI_2,
//!     Bivector::unit_xy(),
//! );
//!
//! // Apply rotation
//! let rotated = rotor.transform(&v);
//! assert!((rotated.x()).abs() < 1e-10);
//! assert!((rotated.y() - 1.0).abs() < 1e-10);
//! ```

pub mod complex;
pub mod dual;
pub mod dualquat;
pub mod elliptic;
pub mod euclidean;
pub mod hyperbolic;
pub mod minkowski;
pub mod projective;
pub mod quaternion;

#[cfg(feature = "rerun-0_28")]
pub mod visualization;
