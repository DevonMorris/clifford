//! Constraint-based expression simplification.
//!
//! This module applies type constraints to simplify symbolic expressions.
//! For example, if a Rotor satisfies `s*s + xy*xy + xz*xz + yz*yz = 1`,
//! any occurrence of that pattern in an expression can be replaced with `1`.

use symbolica::atom::{Atom, AtomCore};

use crate::spec::TypeSpec;

/// Simplifies expressions using type constraints.
///
/// This simplifier recognizes patterns from type constraints and substitutes
/// them with the constraint's constant value.
///
/// # Example
///
/// For a Rotor with constraint `s*s + xy*xy + xz*xz + yz*yz = 1`:
/// - Pattern: `a_s*a_s + a_xy*a_xy + a_xz*a_xz + a_yz*a_yz`
/// - Substituted with: `1`
pub struct ConstraintSimplifier {
    /// Map from pattern atom to replacement value.
    /// Key: the expanded constraint LHS as Atom
    /// Value: the RHS constant as Atom
    substitutions: Vec<(Atom, Atom)>,
}

impl ConstraintSimplifier {
    /// Creates a new constraint simplifier for the given types.
    ///
    /// # Arguments
    ///
    /// * `types` - The types being operated on
    /// * `prefixes` - The variable prefixes used (e.g., "a", "b")
    ///
    /// Note: Constraints are now auto-derived during code generation rather than
    /// stored in TypeSpec. This simplifier currently returns an empty instance.
    /// Constraint-based simplification can be reimplemented using ConstraintDeriver
    /// if needed.
    #[allow(unused_variables)]
    pub fn new(types: &[&TypeSpec], prefixes: &[&str]) -> Self {
        // Constraints are now auto-derived, not stored in TypeSpec.
        // Return empty simplifier for now.
        Self {
            substitutions: Vec::new(),
        }
    }

    /// Applies constraint substitutions to an expression.
    ///
    /// Looks for subexpressions matching constraint patterns and replaces them.
    pub fn apply(&self, expr: &Atom) -> Atom {
        if self.substitutions.is_empty() {
            return expr.clone();
        }

        // Expand the expression to canonical form for matching
        let mut result = expr.expand();

        // Try each substitution
        for (pattern, value) in &self.substitutions {
            result = Self::substitute_pattern(&result, pattern, value);
        }

        result
    }

    /// Substitutes a pattern in an expression.
    ///
    /// This is a simple approach that checks if the expression contains
    /// the pattern as a subexpression.
    fn substitute_pattern(expr: &Atom, pattern: &Atom, value: &Atom) -> Atom {
        // Check if expr equals pattern directly
        if Self::atoms_equal(expr, pattern) {
            return value.clone();
        }

        // Check if expr contains pattern as a subexpression in a sum
        // expr = pattern + rest => expr = value + rest
        if let Some(result) = Self::try_substitute_in_sum(expr, pattern, value) {
            return result;
        }

        // For more complex cases, we could use Symbolica's pattern matching
        // but for now, return expr unchanged
        expr.clone()
    }

    /// Checks if two atoms are equal.
    ///
    /// Uses Symbolica's native `PartialEq` implementation on expanded forms.
    fn atoms_equal(a: &Atom, b: &Atom) -> bool {
        // Expand both to canonical form for comparison
        let a_expanded = a.expand();
        let b_expanded = b.expand();

        // Use Symbolica's PartialEq implementation directly
        a_expanded == b_expanded
    }

    /// Tries to substitute a pattern within a sum expression.
    ///
    /// If expr = pattern + rest, returns value + rest.
    fn try_substitute_in_sum(expr: &Atom, pattern: &Atom, value: &Atom) -> Option<Atom> {
        // Compute expr - pattern
        let difference = expr - pattern;
        let simplified_diff = difference.expand();
        let expanded_expr = expr.expand();

        // Count terms in the expanded expressions
        let expr_terms = Self::count_terms(&expanded_expr);
        let diff_terms = Self::count_terms(&simplified_diff);
        let pattern_terms = Self::count_terms(&pattern.expand());

        // If the pattern was present, subtracting it should reduce terms
        // The difference should have fewer terms than the original
        // (approximately: expr_terms - pattern_terms + possible cancellations)
        if diff_terms < expr_terms || Self::is_simpler(&simplified_diff, &expanded_expr) {
            // Check that we removed approximately the right number of terms
            // Allow for some cancellation effects
            if expr_terms.saturating_sub(diff_terms) > 0 || diff_terms + pattern_terms > expr_terms
            {
                // Pattern was found and removed, add value back
                let result = &simplified_diff + value;
                return Some(result.expand());
            }
        }

        None
    }

    /// Counts the number of terms in an atom.
    ///
    /// For Add expressions, returns the number of addends.
    /// For other atoms, returns 1.
    fn count_terms(atom: &Atom) -> usize {
        atom.as_add_view().map(|add| add.get_nargs()).unwrap_or(1)
    }

    /// Checks if expr_a is simpler than expr_b.
    ///
    /// Uses term count as the complexity metric.
    fn is_simpler(a: &Atom, b: &Atom) -> bool {
        Self::count_terms(a) < Self::count_terms(b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::borrow::Cow;
    use std::sync::Mutex;
    use symbolica::atom::DefaultNamespace;
    use symbolica::parser::ParseSettings;

    // Symbolica uses global state that conflicts when tests run in parallel.
    // Tests prefixed with `symbolica_` are configured to run serially via nextest.
    // The mutex provides a fallback for `cargo test` users.
    static SYMBOLICA_LOCK: Mutex<()> = Mutex::new(());

    fn parse_atom(s: &str) -> Atom {
        let input = DefaultNamespace {
            namespace: Cow::Borrowed(env!("CARGO_CRATE_NAME")),
            data: s,
            file: Cow::Borrowed(file!()),
            line: line!() as usize,
        };
        Atom::parse(input, ParseSettings::symbolica()).unwrap()
    }

    #[test]
    fn symbolica_empty_simplifier_returns_unchanged() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        let simplifier = ConstraintSimplifier::new(&[], &[]);
        let expr = parse_atom("a + b");
        let result = simplifier.apply(&expr);
        assert_eq!(result, expr);
    }
}
