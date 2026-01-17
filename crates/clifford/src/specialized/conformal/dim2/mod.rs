//! 2D Conformal Geometric Algebra - Cl(3,1,0)
//!
//! The 2D Conformal Geometric Algebra is a 16-dimensional algebra that embeds
//! 2D Euclidean space into a 4D conformal space. This enables elegant representation
//! of circles, lines, and conformal transformations in the plane.
//!
//! # Structure
//!
//! | Grade | Count | Geometric Object |
//! |-------|-------|-----------------|
//! | 0 | 1 | Scalar |
//! | 1 | 4 | RoundPoint (null vectors = actual 2D points) |
//! | 2 | 6 | PointPair (two points, flat points) |
//! | 3 | 4 | Circle (circles, lines) |
//! | 4 | 1 | Pseudoscalar |
//!
//! # Basis Convention
//!
//! Following the convention from conformalgeometricalgebra.org:
//!
//! - `e₁, e₂`: Euclidean basis vectors (spatial directions)
//! - `e₃ = ½(e₋ - e₊)`: Origin representation
//! - `e₄ = e₋ + e₊`: Point at infinity
//!
//! Where `e₋² = -1` and `e₊² = +1` are the conformal basis vectors.
//!
//! # Null Vectors and Points
//!
//! In CGA, actual geometric 2D points are represented as **null vectors** (vectors
//! with zero norm). A 2D Euclidean point `(x, y)` is embedded as:
//!
//! ```text
//! P = x·e₁ + y·e₂ + e₃ + ½(x² + y²)·e₄
//! ```
//!
//! Where `e₃` is the origin and `e₄` is the point at infinity.
//!
//! # Key Properties
//!
//! - **Inner product of points**: `P₁ · P₂ = -½|p₁ - p₂|²`
//! - **Circle from 3 points**: `C = P₁ ∧ P₂ ∧ P₃`
//! - **Line from 2 points**: `L = P₁ ∧ P₂ ∧ e∞` (circle through infinity)
//! - **Intersection**: Computed via wedge/antiwedge products
//!
//! # Applications
//!
//! - **2D graphics**: Vector graphics, font rendering, collision detection
//! - **Robotics**: 2D motion planning, linkage kinematics
//! - **Geometry**: Circle packing, Apollonian gaskets
//! - **Physics**: 2D wave optics, electromagnetic problems
//!
//! # Reference
//!
//! See <https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page>
//! for detailed documentation on the types and operations.

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions
mod extensions;

// Re-export generated types
pub use generated::types::*;
