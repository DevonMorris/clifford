#![doc = include_str!("../README.md")]
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
