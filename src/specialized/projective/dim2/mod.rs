//! 2D Projective Geometric Algebra types.
//!
//! This module provides specialized types for 2D PGA `Cl(2,0,1)`.
//!
//! # Types by Grade
//!
//! | Grade | Type | Components | Geometric Meaning |
//! |-------|------|------------|-------------------|
//! | 1 | [`Point`] | e₁, e₂, e₀ | Point (homogeneous) |
//! | 2 | [`Line`] | 3 components | Line |
//! | 1+3 | [`Motor`] | 4 components | Rigid transformation (rotation + translation) |
//! | 0+2 | [`Flector`] | 4 components | Reflection transformation |
//!
//! # Points as Homogeneous Coordinates
//!
//! A 2D PGA point is a homogeneous coordinate:
//!
//! ```text
//! P = x·e₁ + y·e₂ + w·e₀
//! ```
//!
//! | Type | Condition | Cartesian |
//! |------|-----------|-----------|
//! | Finite point | `w ≠ 0` | `(x/w, y/w)` |
//! | Point at infinity | `w = 0` | Direction `(x, y)` |
//!
//! # Motors in 2D PGA
//!
//! In 2D PGA (odd total dimension n=3), the complement operation flips grade parity.
//! This means the **odd subalgebra** `[1, 3]` is closed under the antiproduct,
//! making it the natural choice for motors (rigid transformations).
//!
//! Motors transform objects via the **antisandwich product**: `X' = M̃ ⊙ X ⊙ M`

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions (hand-written)
mod extensions;

// Re-export generated types at module level
pub use generated::types::*;
