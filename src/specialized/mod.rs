//! Specialized types for common geometric algebras.
//!
//! This module provides optimized, ergonomic types organized by algebra type
//! and dimension. Currently supported:
//!
//! - **[`euclidean`]**: Euclidean geometric algebra (EGA) for 2D and 3D
//!
//! # Type Organization
//!
//! ```text
//! specialized/
//!   euclidean/
//!     dim2/   - 2D Euclidean GA: Vec2, Bivec2, Rotor2
//!     dim3/   - 3D Euclidean GA: Vec3, Bivec3, Trivec3, Rotor3
//!   projective/   (future)
//!   conformal/    (future)
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
//! ```
//! use clifford::specialized::euclidean::dim3::{Vec3, Bivec3, Rotor3};
//!
//! // Create a vector
//! let v = Vec3::new(1.0, 0.0, 0.0);
//!
//! // Create a 90° rotation in the xy-plane
//! let rotor = Rotor3::from_angle_plane(
//!     std::f64::consts::FRAC_PI_2,
//!     Bivec3::unit_xy(),
//! );
//!
//! // Apply rotation
//! let rotated = rotor.rotate(v);
//! assert!((rotated.x).abs() < 1e-10);
//! assert!((rotated.y - 1.0).abs() < 1e-10);
//! ```

pub mod euclidean;
