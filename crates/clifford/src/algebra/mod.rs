//! Core algebraic types for geometric algebra.
//!
//! This module provides the [`Multivector`] type, the fundamental element of
//! geometric algebra that can represent scalars, vectors, bivectors, and
//! arbitrary combinations thereof.
//!
//! # What is a Grade?
//!
//! In linear algebra, you work with vectors (1D directed quantities). Geometric
//! algebra generalizes this to **k-dimensional** directed quantities called **blades**:
//!
//! | Grade | Name | Geometric Meaning | LA Analogue |
//! |-------|------|-------------------|-------------|
//! | 0 | Scalar | Magnitude only | Real number |
//! | 1 | Vector | Directed line segment | Vector in Rⁿ |
//! | 2 | Bivector | Oriented plane segment | *(none)* |
//! | 3 | Trivector | Oriented volume | *(determinant, roughly)* |
//! | n | Pseudoscalar | Oriented n-volume | *(n×n determinant)* |
//!
//! **Key insight**: The cross product in 3D actually produces a bivector (an
//! oriented plane), not a vector. The "vector" result is the **dual** of that plane.
//!
//! # Basis Blade Notation
//!
//! Basis blades are written as products of basis vectors:
//!
//! ```text
//! Grade 0: 1           (the scalar, just "1")
//! Grade 1: e₁, e₂, e₃  (basis vectors)
//! Grade 2: e₁₂, e₂₃, e₁₃  (basis bivectors = e₁∧e₂, e₂∧e₃, e₁∧e₃)
//! Grade 3: e₁₂₃        (pseudoscalar = e₁∧e₂∧e₃)
//! ```
//!
//! The subscript `e₁₂` means "e₁ wedge e₂" — an oriented plane containing both vectors.
//! Swapping order flips sign: `e₂₁ = -e₁₂`.
//!
//! # The Multivector
//!
//! A multivector is a linear combination of basis blades from different grades:
//!
//! ```text
//! M = s + v₁e₁ + v₂e₂ + v₃e₃ + b₁e₁₂ + b₂e₂₃ + b₃e₁₃ + pe₁₂₃
//!     ─   ─────────────────   ───────────────────────   ──────
//!     │         │                      │                  │
//!   scalar   vector              bivector          pseudoscalar
//!   (grade 0) (grade 1)          (grade 2)          (grade 3)
//! ```
//!
//! # Why Mix Grades?
//!
//! The geometric product of two vectors yields mixed grades:
//!
//! ```text
//! ab = a·b + a∧b
//!      ───   ───
//!       │     │
//!    scalar  bivector
//!    (how aligned)  (plane they span)
//! ```
//!
//! - **Parallel vectors**: `ab = |a||b|` (pure scalar, bivector = 0)
//! - **Perpendicular vectors**: `ab = a∧b` (pure bivector, scalar = 0)
//! - **General case**: Both parts are non-zero
//!
//! This is why GA unifies the dot and cross products — they're two parts of one operation!
//!
//! # The Geometric Product
//!
//! Unlike the dot product (scalar result) or cross product (vector result), the
//! **geometric product** returns a multivector. For vectors:
//!
//! | Vector Relationship | Dot Product | Wedge Product | Geometric Product |
//! |---------------------|-------------|---------------|-------------------|
//! | Parallel (a ∥ b) | `|a||b|` | 0 | `|a||b|` (scalar) |
//! | Perpendicular (a ⊥ b) | 0 | `a∧b` | `a∧b` (bivector) |
//! | General | `|a||b|cos θ` | `|a||b|sin θ B̂` | `|a||b|(cos θ + sin θ B̂)` |
//!
//! where `B̂` is the unit bivector of the plane containing a and b.
//!
//! **Special property**: `a² = a·a = |a|²` (a vector times itself is its squared length)
//!
//! # Operations
//!
//! The main operations on multivectors are:
//!
//! - **Geometric product** (`*`): The fundamental product combining inner and outer
//! - **Inner product**: Grade-lowering contraction (generalizes dot product)
//! - **Outer/wedge product** (`∧`): Grade-raising (generalizes cross product)
//! - **Addition/Subtraction** (`+`, `-`): Component-wise
//! - **Scalar multiplication** (`* scalar`): Scale all components
//! - **Reverse** (`.reverse()`): Reverses blade orientation
//! - **Grade projection** (`.grade(k)`): Extract the grade-k part
//! - **Inverse** (`.inverse()`): Multiplicative inverse (when it exists)
//!
//! # When to Use Generic Multivector
//!
//! This generic [`Multivector`] type works with any metric signature. Use it when:
//!
//! - Learning GA concepts (see all the components explicitly)
//! - Working with exotic algebras (Minkowski space, higher dimensions)
//! - Maximum flexibility is needed
//!
//! For production 2D/3D work, prefer the optimized types in [`crate::specialized`]:
//!
//! - [`crate::specialized::euclidean`]: Optimized `Vector`, `Bivector`, `Rotor`
//! - [`crate::specialized::projective`]: PGA types for rigid transforms
//!
//! The specialized types have the same mathematical semantics but better performance.
//!
//! # Example
//!
//! ```
//! use clifford::algebra::Multivector;
//! use clifford::signature::Euclidean3;
//!
//! // Create two vectors
//! let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 0.0, 0.0]);
//! let b: Multivector<f64, Euclidean3> = Multivector::vector(&[0.0, 1.0, 0.0]);
//!
//! // Geometric product: ab = a·b + a∧b = 0 + e₁₂
//! let ab = &a * &b;
//!
//! // The result is a pure bivector (oriented plane)
//! assert!(ab.scalar_part().abs() < 1e-10); // No scalar part (perpendicular)
//! ```

mod multivector;

#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;

pub use multivector::Multivector;
