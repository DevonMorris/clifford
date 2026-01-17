//! Euclidean Geometric Algebra (EGA) specialized types.
//!
//! Euclidean GA models standard Euclidean space where all basis vectors
//! square to +1. This is the most common GA for graphics, physics, and robotics.
//!
//! # What "Euclidean" Means
//!
//! In the metric signature `Cl(p, q, r)`:
//! - `p` = number of basis vectors that square to **+1**
//! - `q` = number that square to **-1** (used in relativity)
//! - `r` = number that square to **0** (used in projective GA)
//!
//! **Euclidean GA has signature `Cl(n, 0, 0)`**: all basis vectors square to +1.
//! This matches the familiar dot product: `eᵢ · eᵢ = +1`.
//!
//! # Familiar Operations Work as Expected
//!
//! | Linear Algebra | Euclidean GA | Notes |
//! |----------------|--------------|-------|
//! | `a · b` (dot) | `a.dot(&b)` | Same formula: `Σ aᵢbᵢ` |
//! | `‖a‖` (norm) | `a.norm()` | Same: `√(a·a)` |
//! | `a × b` (cross) | `(a.wedge(&b)).dual()` | Only in 3D; see below |
//! | Rotation matrix | `Rotor` | Half the storage, no drift |
//!
//! # Wedge Product vs Cross Product
//!
//! The **wedge product** `a ∧ b` generalizes the cross product to any dimension:
//!
//! | Property | Cross Product | Wedge Product |
//! |----------|---------------|---------------|
//! | Dimensions | 3D only | Any dimension |
//! | Result type | Vector | Bivector (oriented plane) |
//! | `a × b` magnitude | `|a||b|sin θ` | Same: `|a||b|sin θ` |
//! | `b × a` | `-a × b` | `-(a ∧ b)` (same anti-symmetry) |
//!
//! **Connection**: In 3D, `a × b = *(a ∧ b)` where `*` is the Hodge dual.
//! The cross product is the dual of the wedge product.
//!
//! # Types by Dimension
//!
//! ## 2D: `Cl(2,0,0)`
//!
//! | Type | Grade | Basis | Meaning |
//! |------|-------|-------|---------|
//! | `Scalar` | 0 | `1` | Magnitude |
//! | `Vector` | 1 | `e₁, e₂` | Point or direction |
//! | `Bivector` | 2 | `e₁₂` | Oriented area (like signed angle) |
//!
//! ## 3D: `Cl(3,0,0)`
//!
//! | Type | Grade | Basis | Meaning |
//! |------|-------|-------|---------|
//! | `Scalar` | 0 | `1` | Magnitude |
//! | `Vector` | 1 | `e₁, e₂, e₃` | Point or direction |
//! | `Bivector` | 2 | `e₁₂, e₂₃, e₁₃` | Oriented plane (rotation axis dual) |
//! | `Trivector` | 3 | `e₁₂₃` | Oriented volume (pseudoscalar) |
//!
//! # Rotors: Better Than Rotation Matrices
//!
//! A **rotor** represents a rotation using the geometric product. Benefits:
//!
//! - **Compact**: 3 components (2D) or 4 components (3D) vs 4 or 9 for matrices
//! - **No drift**: Renormalizing a rotor is cheap (`R / |R|`)
//! - **Composable**: `R₂ * R₁` rotates by R₁ then R₂
//! - **Interpolatable**: SLERP works naturally
//!
//! In 3D, a rotor has the same structure as a quaternion (see [`dim3`] for the mapping).
//!
//! # Modules
//!
//! - [`dim2`]: 2D Euclidean GA `Cl(2,0,0)` - vectors, bivectors, rotors
//! - [`dim3`]: 3D Euclidean GA `Cl(3,0,0)` - vectors, bivectors, trivectors, rotors
//!
//! # Example
//!
//! ```ignore
//! use clifford::ops::Transform;
//! use clifford::specialized::euclidean::dim3::{Vector, Rotor, Bivector};
//! use std::f64::consts::FRAC_PI_2;
//!
//! // 90° rotation around z-axis
//! let rotor = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
//! let v = Vector::new(1.0, 0.0, 0.0);
//! let rotated = rotor.transform(&v);
//!
//! assert!((rotated.y() - 1.0).abs() < 1e-10);
//! ```

pub mod dim2;
pub mod dim3;
