//! 3D Conformal Geometric Algebra - Cl(4,1,0)
//!
//! The 3D Conformal Geometric Algebra is a 32-dimensional algebra that embeds
//! 3D Euclidean space into a 5D conformal space. This enables elegant representation
//! of spheres, circles, planes, and conformal transformations.
//!
//! # Structure
//!
//! | Grade | Count | Geometric Object |
//! |-------|-------|-----------------|
//! | 0 | 1 | Scalar |
//! | 1 | 5 | RoundPoint (null vectors = actual points) |
//! | 2 | 10 | Dipole (point pairs, flat points) |
//! | 3 | 10 | Circle (circles, lines) |
//! | 4 | 5 | Sphere (spheres, planes) |
//! | 5 | 1 | Pseudoscalar |
//!
//! # Basis Convention
//!
//! Following the convention from conformalgeometricalgebra.org:
//!
//! - `e₁, e₂, e₃`: Euclidean basis vectors (spatial directions)
//! - `e₄ = ½(e₋ - e₊)`: Origin representation
//! - `e₅ = e₋ + e₊`: Point at infinity
//!
//! Where `e₋² = -1` and `e₊² = +1` are the conformal basis vectors.
//!
//! # Null Vectors and Points
//!
//! In CGA, actual geometric points are represented as **null vectors** (vectors
//! with zero norm). A 3D Euclidean point `(x, y, z)` is embedded as:
//!
//! ```text
//! P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
//! ```
//!
//! Where `e₀` is the origin and `e∞` is the point at infinity.
//!
//! # Key Properties
//!
//! - **Inner product of points**: `P₁ · P₂ = -½|p₁ - p₂|²`
//! - **Sphere from 4 points**: `S = P₁ ∧ P₂ ∧ P₃ ∧ P₄`
//! - **Circle from 3 points**: `C = P₁ ∧ P₂ ∧ P₃`
//! - **Intersection**: Computed via wedge/antiwedge products
//!
//! # Applications
//!
//! - **Computer graphics**: Ray tracing, collision detection
//! - **Robotics**: Inverse kinematics, motion planning
//! - **Physics**: Electromagnetism, quantum mechanics
//! - **Geometry**: Sphere packing, circle fitting
//!
//! # Reference
//!
//! See <https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page>
//! for detailed documentation on the types and operations.

mod generated;

pub use generated::types::*;
