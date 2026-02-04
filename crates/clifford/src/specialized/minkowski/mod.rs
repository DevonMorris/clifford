//! Minkowski spacetime algebras.
//!
//! This module provides specialized types for Minkowski (spacetime) geometric algebras,
//! which have mixed signature (both positive and negative squaring basis vectors).
//! These algebras are fundamental for relativistic physics and spacetime geometry.
//!
//! # Available Algebras
//!
//! - **[`dim2`]**: Minkowski plane Cl(1,1,0) - 2D spacetime (1 space + 1 time)
//! - **[`dim3`]**: Minkowski spacetime Cl(1,2,0) - 3D spacetime (1 time + 2 space)
//!
//! # What is Minkowski Space?
//!
//! Minkowski space is the geometric setting of special relativity. Unlike Euclidean
//! space where distance is always positive, Minkowski space has a metric that can
//! be positive, negative, or zero:
//!
//! ```text
//! Euclidean:  ds² = dx² + dy² + dz²        (always >= 0)
//! Minkowski:  ds² = -dt² + dx² + dy² + dz²  (can be <0, =0, or >0)
//! ```
//!
//! # Metric Signature
//!
//! We use the "mostly plus" convention (-,+,+,+):
//! - **Time basis** `e_t` squares to **-1** (negative signature)
//! - **Space bases** `e_x, e_y, ...` square to **+1** (positive signature)
//!
//! # Causal Structure
//!
//! Minkowski algebras have an **indefinite metric**, meaning vectors are classified by
//! their squared norm:
//!
//! | Type | Condition | Physical Meaning |
//! |------|-----------|------------------|
//! | **Timelike** | v² < 0 | Slower than light (massive particles) |
//! | **Lightlike/Null** | v² = 0 | Speed of light (photons) |
//! | **Spacelike** | v² > 0 | Faster than light (no causal connection) |
//!
//! Null vectors form the "light cone" that separates causal (timelike) from
//! non-causal (spacelike) regions.
//!
//! # Applications
//!
//! - **Special Relativity**: Lorentz boosts, time dilation, length contraction
//! - **Particle Physics**: Representing 4-momentum, energy-momentum conservation
//! - **Electromagnetism**: The electromagnetic field tensor
//! - **Spacetime Intervals**: Computing proper time and proper length
//!
//! # Example
//!
//! ```ignore
//! use clifford::specialized::minkowski::dim2::{Vector, Bivector};
//!
//! // A timelike vector (inside the light cone)
//! let v = Vector::new(2.0, 1.0);  // (t=2, x=1)
//! // v² = -t² + x² = -4 + 1 = -3 < 0 (timelike)
//!
//! // A null vector (on the light cone)
//! let light = Vector::new(1.0, 1.0);  // (t=1, x=1)
//! // light² = -1 + 1 = 0 (null/lightlike)
//! ```
//!
//! # References
//!
//! - [Minkowski space](https://en.wikipedia.org/wiki/Minkowski_space)
//! - [Spacetime algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)

pub mod dim2;
pub mod dim3;
