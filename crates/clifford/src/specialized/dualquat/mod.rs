//! Dual Quaternions - Cl(0,2,1)
//!
//! Dual quaternions extend quaternions with a dual (nilpotent) component,
//! enabling unified representation of rigid body transformations
//! (rotation + translation) in a single algebraic element.
//!
//! # Structure
//!
//! | Grade | Blades | Description |
//! |-------|--------|-------------|
//! | 0 | 1 | Scalar |
//! | 1 | e₁, e₂, e₃ | Vectors (quaternion i,j and dual unit) |
//! | 2 | e₁₂, e₁₃, e₂₃ | Bivectors (quaternion k and dual i,j) |
//! | 3 | e₁₂₃ | Trivector (dual k, pseudoscalar) |
//!
//! # Basis Properties
//!
//! - `e₁² = e₂² = -1` (quaternion basis)
//! - `e₃² = 0` (dual/nilpotent unit)
//! - `e₁₂` is quaternion k (squares to -1)
//! - `e₁₂₃` is the pseudoscalar (squares to 0)
//!
//! # Dual Quaternion Form
//!
//! A dual quaternion `q = q_r + ε·q_d` consists of:
//! - `q_r`: Real quaternion part (rotation)
//! - `q_d`: Dual quaternion part (translation encoded)
//! - `ε = e₃`: Dual unit with ε² = 0
//!
//! # Applications
//!
//! - **Rigid body transformations**: Rotation and translation combined
//! - **Screw motion**: Natural representation of screw axes
//! - **Skeletal animation**: Efficient bone transforms with linear blending
//! - **Robotics**: Forward and inverse kinematics
//!
//! # Norm
//!
//! Dual quaternions use **Clifford conjugate** for their norm.
//! Since r > 0 (degenerate), types implement `DegenerateNormed`.
//!
//! Unit dual quaternions satisfy: `q * conjugate(q) = 1`
//!
//! # References
//!
//! - [Dual quaternion](https://en.wikipedia.org/wiki/Dual_quaternion)
//! - [Clifford algebra](https://en.wikipedia.org/wiki/Clifford_algebra)

mod generated;

pub use generated::types::*;
