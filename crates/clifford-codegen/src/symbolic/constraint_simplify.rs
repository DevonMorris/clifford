//! Constraint-based expression simplification.
//!
//! This module applies type constraints to simplify symbolic expressions.
//! For example, if a Rotor satisfies `s*s + xy*xy + xz*xz + yz*yz = 1`,
//! any occurrence of that pattern in an expression can be replaced with `1`.

use std::borrow::Cow;

use symbolica::atom::{Atom, AtomCore, DefaultNamespace};
use symbolica::parser::ParseSettings;

use crate::spec::{TypeSpec, UserConstraint};

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
    /// * `types` - The types being operated on (with their constraints)
    /// * `prefixes` - The variable prefixes used (e.g., "a", "b")
    pub fn new(types: &[&TypeSpec], prefixes: &[&str]) -> Self {
        let mut substitutions = Vec::new();

        for (ty, prefix) in types.iter().zip(prefixes.iter()) {
            for constraint in &ty.constraints {
                if let Some((pattern, value)) = Self::build_substitution(ty, constraint, prefix) {
                    substitutions.push((pattern, value));
                }
            }
        }

        Self { substitutions }
    }

    /// Builds a substitution rule from a constraint.
    ///
    /// Returns (pattern_atom, replacement_atom) if the constraint can be used
    /// for substitution.
    fn build_substitution(
        ty: &TypeSpec,
        constraint: &UserConstraint,
        prefix: &str,
    ) -> Option<(Atom, Atom)> {
        // Parse constraint expression: "s*s + xy*xy + xz*xz + yz*yz = 1"
        let parts: Vec<&str> = constraint.expression.split('=').collect();
        if parts.len() != 2 {
            return None;
        }

        let lhs = parts[0].trim();
        let rhs = parts[1].trim();

        // Parse RHS as a constant
        let rhs_value: i64 = rhs.parse().ok()?;

        // Build the pattern with prefixed variables
        let pattern_str = Self::prefix_expression(lhs, &ty.fields, prefix);

        // Parse as Symbolica atoms
        let pattern = Self::parse_atom(&pattern_str)?;
        let value = Atom::num(rhs_value);

        // Expand the pattern to canonical form
        let expanded_pattern = pattern.expand();

        Some((expanded_pattern, value))
    }

    /// Prefixes all field names in an expression.
    ///
    /// Converts "s*s + xy*xy" to "a_s*a_s + a_xy*a_xy" with prefix "a".
    fn prefix_expression(expr: &str, fields: &[crate::spec::FieldSpec], prefix: &str) -> String {
        let field_names: Vec<&str> = fields.iter().map(|f| f.name.as_str()).collect();

        let mut result = expr.to_string();

        // Sort by length descending to avoid partial replacements
        let mut sorted_names = field_names.clone();
        sorted_names.sort_by_key(|b| std::cmp::Reverse(b.len()));

        for name in sorted_names {
            // Replace field name with prefixed version
            // Use word boundaries to avoid partial matches
            let prefixed = format!("{}_{}", prefix, name);

            // Simple replacement - assumes field names are distinct tokens
            // This works because field names don't contain operators
            result = Self::replace_identifier(&result, name, &prefixed);
        }

        result
    }

    /// Replaces an identifier in an expression, respecting word boundaries.
    fn replace_identifier(expr: &str, from: &str, to: &str) -> String {
        let mut result = String::new();
        let mut current_word = String::new();

        for c in expr.chars() {
            if c.is_alphanumeric() || c == '_' {
                current_word.push(c);
            } else {
                if !current_word.is_empty() {
                    if current_word == from {
                        result.push_str(to);
                    } else {
                        result.push_str(&current_word);
                    }
                    current_word.clear();
                }
                result.push(c);
            }
        }

        // Don't forget the last word
        if !current_word.is_empty() {
            if current_word == from {
                result.push_str(to);
            } else {
                result.push_str(&current_word);
            }
        }

        result
    }

    /// Parses a string into a Symbolica Atom.
    fn parse_atom(s: &str) -> Option<Atom> {
        let input = DefaultNamespace {
            namespace: Cow::Borrowed(env!("CARGO_CRATE_NAME")),
            data: s,
            file: Cow::Borrowed(file!()),
            line: line!() as usize,
        };
        Atom::parse(input, ParseSettings::symbolica()).ok()
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
    use crate::spec::FieldSpec;

    fn make_rotor_type() -> TypeSpec {
        TypeSpec {
            name: "Rotor".to_string(),
            grades: vec![0, 2],
            description: None,
            fields: vec![
                FieldSpec {
                    name: "s".to_string(),
                    blade_index: 0,
                    grade: 0,
                },
                FieldSpec {
                    name: "xy".to_string(),
                    blade_index: 3,
                    grade: 2,
                },
                FieldSpec {
                    name: "xz".to_string(),
                    blade_index: 5,
                    grade: 2,
                },
                FieldSpec {
                    name: "yz".to_string(),
                    blade_index: 6,
                    grade: 2,
                },
            ],
            alias_of: None,
            constraints: vec![crate::spec::UserConstraint {
                name: "unit".to_string(),
                description: None,
                expression: "s*s + xy*xy + xz*xz + yz*yz = 1".to_string(),
                solve_for: Some("s".to_string()),
                sign: crate::spec::SignConvention::Positive,
                enforce: None,
                has_domain_restriction: true,
            }],
            versor: None,
        }
    }

    #[test]
    #[ignore = "Symbolica global state conflicts"]
    fn prefix_expression_works() {
        let rotor = make_rotor_type();
        let result = ConstraintSimplifier::prefix_expression(
            "s*s + xy*xy + xz*xz + yz*yz",
            &rotor.fields,
            "a",
        );
        assert!(result.contains("a_s"));
        assert!(result.contains("a_xy"));
        assert!(result.contains("a_xz"));
        assert!(result.contains("a_yz"));
    }

    #[test]
    #[ignore = "Symbolica global state conflicts"]
    fn replace_identifier_works() {
        let result = ConstraintSimplifier::replace_identifier("s*s + s", "s", "a_s");
        assert_eq!(result, "a_s*a_s + a_s");

        // Shouldn't replace partial matches
        let result = ConstraintSimplifier::replace_identifier("xy*xyz", "xy", "a_xy");
        assert_eq!(result, "a_xy*xyz");
    }
}
