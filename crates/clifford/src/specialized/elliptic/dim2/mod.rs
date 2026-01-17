//! 2D Elliptic Projective Geometry - Cl(3,0,0)
//!
//! The elliptic plane is the real projective plane RP² with an elliptic
//! metric. It can be visualized as the unit sphere with antipodal points
//! identified.
//!
//! # Structure
//!
//! | Grade | Blades | Description |
//! |-------|--------|-------------|
//! | 0 | 1 | Scalar |
//! | 1 | e₁, e₂, e₃ | Points (homogeneous coordinates) |
//! | 2 | e₁₂, e₁₃, e₂₃ | Lines (great circles) |
//! | 3 | e₁₂₃ | Pseudoscalar |
//!
//! # Key Properties
//!
//! - **No parallel lines**: Any two distinct lines intersect at exactly one point
//! - **Finite extent**: The plane "wraps around" like a sphere
//! - **Point-line duality**: Perfect symmetry between points and lines
//!
//! # Geometric Operations
//!
//! | Operation | Input | Output | Meaning |
//! |-----------|-------|--------|---------|
//! | `p1.wedge(p2)` | Point, Point | Line | Line through two points |
//! | `l1.wedge(l2)` | Line, Line | Point | Intersection of two lines |
//! | `r.transform(p)` | Rotor, Point | Point | Rotate point |
//!
//! # Relationship to Euclidean 3D
//!
//! This algebra has the same signature as 3D Euclidean GA (Cl(3,0,0)).
//! The difference is purely interpretive:
//! - Euclidean: grade-1 = vectors (directions), grade-2 = bivectors (areas)
//! - Elliptic: grade-1 = points, grade-2 = lines
//!
//! # Norm
//!
//! Non-degenerate algebra with all positive signature.
//! Uses reversal for the norm: `|x|² = x * x̃`.
//!
//! # Applications
//!
//! - **Spherical geometry**: Great circle navigation, astronomy
//! - **Computer graphics**: Environment maps, omnidirectional cameras
//! - **Robotics**: Orientation representation
//! - **Projective geometry**: Incidence relations, cross-ratios

mod generated;

pub use generated::types::*;
