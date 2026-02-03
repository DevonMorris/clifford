#![doc = include_str!("../README.md")]
//!
//! # For Linear Algebra Users
//!
//! If you know linear algebra, you're already halfway to understanding Geometric Algebra (GA).
//! This section maps concepts you know to their GA equivalents.
//!
//! ## From Vectors to Blades
//!
//! In LA, you work with vectors. In GA, vectors are just one type of **blade**:
//!
//! | Grade | Name | LA Analogue | Geometric Meaning |
//! |-------|------|-------------|-------------------|
//! | 0 | Scalar | Real number | Magnitude, no direction |
//! | 1 | Vector | Vector in R^n | Directed line segment |
//! | 2 | Bivector | (no direct equivalent) | Oriented plane segment |
//! | 3 | Trivector | (no direct equivalent) | Oriented volume |
//! | n | Pseudoscalar | Determinant (sort of) | Oriented n-volume |
//!
//! **Key insight**: The cross product `a × b` in 3D actually produces a bivector
//! (oriented plane), not a vector. The "vector" you get is the **dual** of that plane.
//!
//! ## From Products to the Geometric Product
//!
//! In LA, you have separate operations: dot product, cross product, matrix multiplication.
//! In GA, the **geometric product** unifies them:
//!
//! ```text
//! a * b = a·b + a∧b
//!       = (inner product) + (outer product)
//!       = (scalar part) + (bivector part)
//! ```
//!
//! | LA Operation | GA Equivalent | Result |
//! |--------------|---------------|--------|
//! | Dot product `a·b` | Inner product | Scalar (grade 0) |
//! | Cross product `a×b` | `*(a∧b)` (dual of wedge) | Vector (grade 1) |
//! | Matrix multiply | Geometric product | Mixed grades |
//!
//! ## From Rotations to Rotors
//!
//! This is where GA really shines. Rotations are represented by **rotors**:
//!
//! | Dimension | LA Representation | GA Representation |
//! |-----------|-------------------|-------------------|
//! | 2D | 2×2 rotation matrix | Rotor (scalar + bivector) |
//! | 2D | Complex number `e^{iθ}` | Rotor `cos(θ/2) + sin(θ/2)e₁₂` |
//! | 3D | 3×3 rotation matrix | Rotor (scalar + bivector) |
//! | 3D | Quaternion | Rotor `s + xy·e₁₂ + xz·e₁₃ + yz·e₂₃` |
//! | nD | SO(n) matrix | Rotor (even-grade multivector) |
//!
//! **Quaternions ARE rotors!** The imaginary units i, j, k are bivectors e₂₃, e₃₁, e₁₂.
//!
//! ## Module Guide
//!
//! - [`algebra`]: Generic [`Multivector`](algebra::Multivector) for any metric signature.
//!   Use when you need maximum flexibility or exotic algebras.
//!
//! - [`specialized::euclidean`]: Optimized types for 2D/3D Euclidean geometry.
//!   Use for standard vector math, rotations, reflections.
//!
//! - [`specialized::projective`]: Projective GA (PGA) for rigid body transforms.
//!   Use when you need unified rotation + translation (like 4×4 matrices, but better).
//!
//! - [`basis`]: Low-level blade representation. Rarely needed directly.
//!
//! - [`signature`]: Metric signatures defining the algebra. Use [`prelude`] instead.
//!
//! ## Quick Decision Guide
//!
//! | Task | Recommended Module |
//! |------|-------------------|
//! | 2D/3D rotations | [`specialized::euclidean`] |
//! | Rigid body transforms | [`specialized::projective`] |
//! | Robotics kinematics | [`specialized::projective`] |
//! | Graphics transforms | [`specialized::projective`] |
//! | Physics simulations | [`algebra`] with appropriate signature |
//! | Learning GA | Start with [`specialized::euclidean::dim3`] |
//!
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/DevonMorris/clifford/main/assets/clifford.png",
    html_favicon_url = "https://raw.githubusercontent.com/DevonMorris/clifford/main/assets/clifford.png"
)]

// Ensure nalgebra feature flags are mutually exclusive
#[cfg(any(
    all(feature = "nalgebra-0_32", feature = "nalgebra-0_33"),
    all(feature = "nalgebra-0_32", feature = "nalgebra-0_34"),
    all(feature = "nalgebra-0_33", feature = "nalgebra-0_34"),
))]
compile_error!(
    "Features `nalgebra-0_32`, `nalgebra-0_33`, and `nalgebra-0_34` are mutually exclusive. Enable only one."
);

pub mod algebra;
pub mod basis;
pub mod norm;
pub mod ops;
pub mod prelude;
pub mod scalar;
pub mod signature;
pub mod specialized;
pub mod wrappers;

// ============================================================================
// Error Types
// ============================================================================

/// Error returned when a geometric constraint is not satisfied.
///
/// This error is returned by `new_checked()` constructors on constrained types
/// (like `Motor`, `Bivector`, `Flector` in PGA) when the provided coefficients
/// violate the geometric constraint.
///
/// # Example
///
/// ```ignore
/// use clifford::generated::projective3::Motor;
///
/// // Valid motor - constraint satisfied
/// let m = Motor::new(1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03);
///
/// // Check with explicit coefficients
/// let result = Motor::new_checked(1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 999.0, 1e-10);
/// assert!(result.is_err()); // e0123 = 999.0 violates constraint
/// ```
#[derive(Debug, Clone)]
pub struct ConstraintError {
    /// The type that failed the constraint check.
    pub type_name: &'static str,
    /// The constraint expression that was violated.
    pub constraint: &'static str,
    /// The residual value (how far from zero).
    pub residual: f64,
}

impl ConstraintError {
    /// Creates a new constraint error.
    pub fn new(type_name: &'static str, constraint: &'static str, residual: f64) -> Self {
        Self {
            type_name,
            constraint,
            residual,
        }
    }
}

impl std::fmt::Display for ConstraintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} constraint violated: {} (residual: {:.2e})",
            self.type_name, self.constraint, self.residual
        )
    }
}

impl std::error::Error for ConstraintError {}

/// Test utilities available only during testing.
#[cfg(test)]
pub(crate) mod test_utils {
    /// Standard epsilon for relative comparisons in tests.
    ///
    /// Use this constant instead of magic numbers like `1e-10` or `1e-9`.
    /// Relative comparisons are more robust than absolute comparisons because
    /// they scale with the magnitude of values being compared.
    ///
    /// # Example
    /// ```ignore
    /// use approx::relative_eq;
    /// use crate::test_utils::RELATIVE_EQ_EPS;
    ///
    /// prop_assert!(relative_eq!(a, b, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    /// ```
    pub const RELATIVE_EQ_EPS: f64 = 1e-10;
}
