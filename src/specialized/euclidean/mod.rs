//! Euclidean Geometric Algebra (EGA) specialized types.
//!
//! Euclidean GA models standard Euclidean space where all basis vectors
//! square to +1. This is the most common GA for graphics, physics, and robotics.
//!
//! # Modules
//!
//! - [`dim2`]: 2D Euclidean GA `Cl(2,0,0)` - vectors, bivectors, rotors
//! - [`dim3`]: 3D Euclidean GA `Cl(3,0,0)` - vectors, bivectors, trivectors, rotors
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim3::{Vector, Rotor, Bivector};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // 90° rotation around z-axis
//! let rotor = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
//! let v = Vector::new(1.0, 0.0, 0.0);
//! let rotated = rotor.rotate(v);
//!
//! assert!((rotated.y - 1.0).abs() < 1e-10);
//! ```

pub mod dim2;
pub mod dim3;
