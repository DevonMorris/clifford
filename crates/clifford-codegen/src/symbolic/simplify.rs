//! Expression simplification utilities.
//!
//! This module provides tools for simplifying symbolic expressions before
//! converting them to Rust code.

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;

/// Utilities for simplifying Symbolica expressions.
///
/// This struct provides methods for simplifying mathematical expressions,
/// which reduces the number of terms and operations in generated code.
pub struct ExpressionSimplifier;

impl ExpressionSimplifier {
    /// Creates a new expression simplifier.
    pub fn new() -> Self {
        Self
    }

    /// Simplifies an expression by expanding and collecting like terms.
    ///
    /// This is the main simplification entry point. It:
    /// 1. Expands all products and powers
    /// 2. Collects like terms
    /// 3. Simplifies numeric coefficients
    ///
    /// # Example
    ///
    /// ```ignore
    /// let simplifier = ExpressionSimplifier::new();
    /// let expr = Atom::parse("(a + b) * (a + b)").unwrap();
    /// let simplified = simplifier.simplify(&expr);
    /// // Result: a*a + 2*a*b + b*b
    /// ```
    pub fn simplify(&self, expr: &Atom) -> Atom {
        expr.expand()
    }

    /// Checks if an expression is zero.
    ///
    /// Returns true if the expression simplifies to zero.
    pub fn is_zero(&self, expr: &Atom) -> bool {
        let expanded = expr.expand();
        match expanded.as_atom_view() {
            AtomView::Num(n) => {
                let coeff = n.get_coeff_view();
                matches!(coeff, CoefficientView::Natural(0, _, 0, _))
            }
            _ => false,
        }
    }

    /// Counts the number of terms in an expression.
    ///
    /// This is useful for comparing simplified vs unsimplified expressions.
    pub fn term_count(&self, expr: &Atom) -> usize {
        let expanded = expr.expand();
        match expanded.as_atom_view() {
            AtomView::Add(a) => a.iter().count(),
            AtomView::Num(n) => {
                let coeff = n.get_coeff_view();
                if matches!(coeff, CoefficientView::Natural(0, _, 0, _)) {
                    0
                } else {
                    1
                }
            }
            _ => 1,
        }
    }

    /// Counts the total number of operations in an expression.
    ///
    /// This counts additions, multiplications, and other operations.
    pub fn operation_count(&self, expr: &Atom) -> usize {
        self.count_ops_view(expr.as_atom_view())
    }

    /// Counts operations recursively for an AtomView.
    fn count_ops_view(&self, view: AtomView<'_>) -> usize {
        match view {
            AtomView::Num(_) | AtomView::Var(_) => 0,
            AtomView::Add(a) => {
                let children: usize = a.iter().map(|c| self.count_ops_view(c)).sum();
                let additions = a.iter().count().saturating_sub(1);
                children + additions
            }
            AtomView::Mul(m) => {
                let children: usize = m.iter().map(|c| self.count_ops_view(c)).sum();
                let multiplications = m.iter().count().saturating_sub(1);
                children + multiplications
            }
            AtomView::Pow(p) => {
                let (base, exp) = p.get_base_exp();
                self.count_ops_view(base) + self.count_ops_view(exp) + 1
            }
            AtomView::Fun(f) => {
                let children: usize = f.iter().map(|c| self.count_ops_view(c)).sum();
                children + 1
            }
        }
    }
}

impl Default for ExpressionSimplifier {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore = "Symbolica global state conflicts"]
    fn simplify_expands_product() {
        let simplifier = ExpressionSimplifier::new();

        // Create (1 + 1) which should simplify to 2
        let one = Atom::num(1);
        let sum = &one + &one;
        let simplified = simplifier.simplify(&sum);

        // Should be the number 2
        assert!(!simplifier.is_zero(&simplified));
    }

    #[test]
    #[ignore = "Symbolica global state conflicts"]
    fn is_zero_detects_zero() {
        let simplifier = ExpressionSimplifier::new();

        let zero = Atom::num(0);
        assert!(simplifier.is_zero(&zero));

        let one = Atom::num(1);
        assert!(!simplifier.is_zero(&one));
    }

    #[test]
    #[ignore = "Symbolica global state conflicts"]
    fn term_count_for_sum() {
        let simplifier = ExpressionSimplifier::new();

        // a + b has 2 terms
        let a = Atom::num(1);
        let b = Atom::num(2);
        let sum = &a + &b;

        // After simplification, 1 + 2 = 3, which is 1 term
        assert_eq!(simplifier.term_count(&sum), 1);
    }
}
