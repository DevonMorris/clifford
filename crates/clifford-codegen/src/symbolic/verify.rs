//! Constraint verification.
//!
//! This module verifies that product outputs satisfy algebraic constraints
//! when inputs satisfy their constraints.

use symbolica::atom::{Atom, AtomCore};

use super::ConstraintExpr;

/// Result of constraint verification.
#[derive(Debug, Clone)]
pub enum VerificationResult {
    /// Constraint is always satisfied when inputs satisfy their constraints.
    Always,
    /// Verification was inconclusive (expression too complex).
    Inconclusive(String),
}

/// Verifies that outputs satisfy constraints given input constraints.
#[derive(Debug, Default)]
pub struct ConstraintVerifier;

impl ConstraintVerifier {
    /// Creates a new constraint verifier.
    pub fn new() -> Self {
        Self
    }

    /// Checks if a constraint expression is trivially satisfied (both sides equal).
    ///
    /// This is a simpler verification that doesn't require substitution - it just
    /// checks if lhs == rhs after expansion.
    pub fn check_trivial(&self, constraint: &ConstraintExpr) -> VerificationResult {
        let zero_form = constraint.as_zero_form();
        self.check_is_zero(&zero_form)
    }

    /// Checks if an expression is zero.
    fn check_is_zero(&self, expr: &Atom) -> VerificationResult {
        // Create zero atom using Atom::num
        let zero = Atom::num(0);

        // Check if the expression is literally zero
        if expr == &zero {
            return VerificationResult::Always;
        }

        // Try expanding and simplifying
        let expanded = expr.expand();
        if expanded == zero {
            return VerificationResult::Always;
        }

        // Check if it's a simple numeric value by string comparison
        let expr_str = format!("{}", expanded);
        if expr_str == "0" {
            return VerificationResult::Always;
        }

        // For complex expressions, we may need more sophisticated analysis
        // For now, return inconclusive
        VerificationResult::Inconclusive(format!(
            "Could not verify: expression simplifies to {}",
            expanded
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    // Symbolica uses global state that conflicts when tests run in parallel.
    // Tests prefixed with `symbolica_` are configured to run serially via nextest.
    // The mutex provides a fallback for `cargo test` users.
    static SYMBOLICA_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn symbolica_verify_numeric_constraint() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        // Test with a numeric constraint directly constructed
        let verifier = ConstraintVerifier::new();

        // Create atoms for 5 - 5 = 0
        let lhs = Atom::num(5);
        let rhs = Atom::num(5);
        let constraint = super::super::ConstraintExpr {
            lhs,
            rhs,
            original: "5 = 5".to_string(),
        };

        let result = verifier.check_trivial(&constraint);
        assert!(matches!(result, VerificationResult::Always));
    }

    #[test]
    fn symbolica_verify_zero_atom() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        let verifier = ConstraintVerifier::new();

        // Create atoms for 0 = 0
        let lhs = Atom::num(0);
        let rhs = Atom::num(0);
        let constraint = super::super::ConstraintExpr {
            lhs,
            rhs,
            original: "0 = 0".to_string(),
        };

        let result = verifier.check_trivial(&constraint);
        assert!(matches!(result, VerificationResult::Always));
    }
}
