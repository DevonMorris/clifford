//! Minkowski Plane - Cl(1,1,0)
//!
//! The Minkowski plane is a 2-dimensional spacetime algebra with:
//! - One spacelike basis vector `e₁` where `e₁² = +1`
//! - One timelike basis vector `e₂` where `e₂² = -1`
//!
//! # Structure
//!
//! | Grade | Blades | Description |
//! |-------|--------|-------------|
//! | 0 | 1 | Scalar |
//! | 1 | e₁, e₂ | Vectors (spacelike and timelike) |
//! | 2 | e₁₂ | Bivector (pseudoscalar) |
//!
//! # Indefinite Metric
//!
//! The mixed signature gives three types of vectors:
//! - **Spacelike**: v² > 0 (e.g., pure e₁ component)
//! - **Timelike**: v² < 0 (e.g., pure e₂ component)
//! - **Null/Lightlike**: v² = 0 (e.g., e₁ + e₂)
//!
//! # Norm
//!
//! The Minkowski plane uses **Clifford conjugate** for norm:
//! ```text
//! conjugate(s + xe₁ + te₂ + b·e₁₂) = s - xe₁ - te₂ - b·e₁₂
//! a * conjugate(a) = s² - x² + t² - b²
//! ```
//!
//! Note: This norm is **indefinite** - it can be positive, negative, or zero.
//!
//! # Even Subalgebra
//!
//! The even subalgebra (grades 0 and 2) is isomorphic to hyperbolic numbers:
//! - `(e₁₂)² = e₁·e₂·e₁·e₂ = +1`
//!
//! This gives split-complex structure in the even part.
//!
//! # References
//!
//! - [Spacetime algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
//! - [Minkowski space](https://en.wikipedia.org/wiki/Minkowski_space)

mod generated;

pub use generated::types::*;
