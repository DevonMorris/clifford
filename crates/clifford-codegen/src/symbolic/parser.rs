//! Constraint expression parsing.
//!
//! This module parses constraint strings like `s * s + xy * xy = 1` into
//! Symbolica expressions for symbolic manipulation.

use std::borrow::Cow;

use symbolica::atom::{Atom, DefaultNamespace};
use symbolica::parser::ParseSettings;
use thiserror::Error;

/// Errors that can occur during constraint parsing.
#[derive(Debug, Error)]
pub enum ParseError {
    /// Missing equality operator in constraint.
    #[error("constraint must contain '=' operator: {0}")]
    MissingEquality(String),

    /// Failed to parse expression.
    #[error("failed to parse expression '{expr}': {reason}")]
    InvalidExpression {
        /// The expression that failed to parse.
        expr: String,
        /// The reason for the failure.
        reason: String,
    },
}

/// A parsed constraint expression of the form `lhs = rhs`.
#[derive(Debug, Clone)]
pub struct ConstraintExpr {
    /// Left-hand side of the constraint.
    pub lhs: Atom,
    /// Right-hand side of the constraint.
    pub rhs: Atom,
    /// The original expression string.
    pub original: String,
}

impl ConstraintExpr {
    /// Returns the constraint as `lhs - rhs = 0` form.
    pub fn as_zero_form(&self) -> Atom {
        &self.lhs - &self.rhs
    }
}

/// Parser for constraint expressions.
///
/// Converts constraint strings into Symbolica atoms for symbolic manipulation.
#[derive(Debug, Default)]
pub struct ConstraintParser;

impl ConstraintParser {
    /// Creates a new constraint parser.
    pub fn new() -> Self {
        Self
    }

    /// Parses a constraint expression.
    ///
    /// The expression must be of the form `expr = value`, where:
    /// - `expr` is a mathematical expression using field names
    /// - `value` is a numeric constant or expression
    ///
    /// # Arguments
    ///
    /// * `expr` - The constraint expression string (e.g., "s * s + xy * xy = 1")
    /// * `field_names` - The field names that appear in the expression
    ///
    /// # Returns
    ///
    /// A parsed constraint expression.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use clifford_codegen::symbolic::ConstraintParser;
    ///
    /// let parser = ConstraintParser::new();
    /// let constraint = parser.parse("s * s + xy * xy = 1", &["s", "xy"])?;
    /// ```
    pub fn parse(&self, expr: &str, field_names: &[&str]) -> Result<ConstraintExpr, ParseError> {
        // Split on '=' to get lhs and rhs
        let parts: Vec<&str> = expr.split('=').collect();
        if parts.len() != 2 {
            return Err(ParseError::MissingEquality(expr.to_string()));
        }

        let lhs_str = parts[0].trim();
        let rhs_str = parts[1].trim();

        // Rename fields to valid Symbolica identifiers if needed
        let lhs_normalized = self.normalize_fields(lhs_str, field_names);
        let rhs_normalized = self.normalize_fields(rhs_str, field_names);

        // Parse using Symbolica
        let lhs = self.parse_expression(&lhs_normalized, expr)?;
        let rhs = self.parse_expression(&rhs_normalized, expr)?;

        Ok(ConstraintExpr {
            lhs,
            rhs,
            original: expr.to_string(),
        })
    }

    /// Normalizes field names in an expression.
    ///
    /// Some field names might conflict with Symbolica reserved words or
    /// need special handling.
    fn normalize_fields(&self, expr: &str, _field_names: &[&str]) -> String {
        // For now, field names are used as-is since Symbolica handles
        // arbitrary identifiers well
        expr.to_string()
    }

    /// Parses a mathematical expression into a Symbolica Atom.
    fn parse_expression(&self, expr: &str, original: &str) -> Result<Atom, ParseError> {
        // Create a DefaultNamespace for the input using the crate's namespace
        let input = DefaultNamespace {
            namespace: Cow::Borrowed(env!("CARGO_CRATE_NAME")),
            data: expr,
            file: Cow::Borrowed(file!()),
            line: line!() as usize,
        };

        // Parse using Symbolica
        Atom::parse(input, ParseSettings::symbolica()).map_err(|e| ParseError::InvalidExpression {
            expr: original.to_string(),
            reason: e.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Note: Symbolica's global state can cause crashes when running multiple
    // parsing tests in parallel. These tests are marked as ignored but can
    // be run individually with `cargo test --package clifford-codegen <test_name>`.

    #[test]
    #[ignore = "Symbolica global state conflicts"]
    fn parse_constraint_basic() {
        let parser = ConstraintParser::new();
        // Use numeric constraint to avoid symbol conflicts
        let result = parser.parse("1 + 2 = 3", &[]);
        assert!(result.is_ok());
        let constraint = result.unwrap();
        assert_eq!(constraint.original, "1 + 2 = 3");
    }

    #[test]
    fn reject_missing_equality() {
        // This test doesn't use Symbolica parsing at all - just string splitting
        let parser = ConstraintParser::new();
        let result = parser.parse("no equals sign", &[]);
        assert!(matches!(result, Err(ParseError::MissingEquality(_))));
    }
}
