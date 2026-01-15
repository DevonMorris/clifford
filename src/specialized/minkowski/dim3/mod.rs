//! 3D Minkowski Spacetime - Cl(3,1,0)
//!
//! The spacetime algebra (STA) is a 16-dimensional algebra that provides
//! the natural geometric framework for special relativity. It has 3 spatial
//! dimensions and 1 temporal dimension.
//!
//! # Structure
//!
//! | Grade | Blades | Squares | Description |
//! |-------|--------|---------|-------------|
//! | 0 | 1 | +1 | Scalar |
//! | 1 | e₁, e₂, e₃, e₄ | +1, +1, +1, -1 | Spacetime vectors |
//! | 2 | e₁₂, e₁₃, e₁₄, e₂₃, e₂₄, e₃₄ | -1, -1, +1, -1, +1, +1 | Bivectors |
//! | 3 | e₁₂₃, e₁₂₄, e₁₃₄, e₂₃₄ | -1, +1, +1, +1 | Trivectors |
//! | 4 | e₁₂₃₄ | -1 | Pseudoscalar |
//!
//! # Causal Structure
//!
//! The indefinite metric creates a causal structure for vectors:
//!
//! | Condition | Type | Physical Interpretation |
//! |-----------|------|------------------------|
//! | v² < 0 | Timelike | Massive particles, observers |
//! | v² = 0 | Null/Lightlike | Light rays, photons |
//! | v² > 0 | Spacelike | Spatial separations |
//!
//! # Bivector Structure
//!
//! The six bivectors split into two groups:
//!
//! - **Spatial rotations** (xy, xz, yz): Square to -1, generate SO(3) rotations
//! - **Spacetime boosts** (xt, yt, zt): Square to +1, generate Lorentz boosts
//!
//! The even subalgebra (Eventor) is isomorphic to the group Spin(3,1),
//! the double cover of the proper orthochronous Lorentz group SO⁺(3,1).
//!
//! # Physical Interpretation
//!
//! | Element | Physical Meaning |
//! |---------|-----------------|
//! | Vector | Events, 4-velocities, 4-momenta |
//! | Spatial bivector | Angular momentum, magnetic field |
//! | Spacetime bivector | Boost velocity, electric field |
//! | Electromagnetic F | F = E + iB (unified field tensor) |
//!
//! # Applications
//!
//! - **Special relativity**: Lorentz transformations, proper time
//! - **Electromagnetism**: Maxwell's equations in compact form
//! - **Particle physics**: Dirac equation, spinors
//! - **Cosmology**: Relativistic mechanics

mod generated;

pub use generated::types::*;
