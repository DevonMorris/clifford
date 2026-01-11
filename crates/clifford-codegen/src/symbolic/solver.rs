//! Constraint solver for computing dependent coefficients.
//!
//! This module solves linear constraint equations for a specified variable,
//! generating the Rust expression to compute that variable from the others.
//!
//! The constraints we handle are always linear in form:
//! `coeff1*var1*var2 + coeff2*var3*var4 + ... = 0`
//!
//! Where one of the variables is the "solve_for" target.

use thiserror::Error;

/// Errors that can occur during constraint solving.
#[derive(Debug, Error)]
pub enum SolveError {
    /// Failed to parse the constraint expression.
    #[error("failed to parse constraint: {0}")]
    ParseError(String),

    /// Missing equality in constraint.
    #[error("constraint must contain '=' operator: {0}")]
    MissingEquality(String),

    /// The solve_for variable doesn't appear in the constraint.
    #[error("variable '{0}' does not appear in constraint")]
    VariableNotFound(String),
}

/// Result of solving a constraint for a variable.
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// The variable being solved for.
    pub variable: String,
    /// The Rust expression for the numerator.
    pub numerator: String,
    /// The Rust expression for the divisor (if division is needed).
    pub divisor: Option<String>,
    /// The original constraint.
    pub constraint: String,
}

/// A term in a linear constraint expression.
#[derive(Debug, Clone)]
struct Term {
    /// Numeric coefficient (e.g., 2, -2).
    coefficient: i32,
    /// Variables in the term (e.g., ["s", "e0123"]).
    variables: Vec<String>,
}

/// Solver for linear constraint equations.
#[derive(Debug, Default)]
pub struct ConstraintSolver;

impl ConstraintSolver {
    /// Creates a new constraint solver.
    pub fn new() -> Self {
        Self
    }

    /// Solves a constraint for the specified variable.
    ///
    /// Given a constraint like `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`,
    /// solving for `e0123` yields:
    /// - numerator: `e12 * e03 - e13 * e02 + e23 * e01`
    /// - divisor: `s`
    ///
    /// # Arguments
    ///
    /// * `constraint` - The constraint expression string
    /// * `solve_for` - The variable to solve for
    ///
    /// # Returns
    ///
    /// A `SolveResult` containing the numerator and optional divisor expressions.
    pub fn solve(&self, constraint: &str, solve_for: &str) -> Result<SolveResult, SolveError> {
        // Parse constraint into terms
        let terms = self.parse_constraint(constraint)?;

        // Separate terms containing solve_for from the rest
        let mut solve_for_terms: Vec<Term> = Vec::new();
        let mut other_terms: Vec<Term> = Vec::new();

        for term in terms {
            if term.variables.contains(&solve_for.to_string()) {
                solve_for_terms.push(term);
            } else {
                other_terms.push(term);
            }
        }

        if solve_for_terms.is_empty() {
            return Err(SolveError::VariableNotFound(solve_for.to_string()));
        }

        // For linear constraints, there should be exactly one term with solve_for
        // (or multiple terms with the same other factors)
        // Extract the coefficient of solve_for
        let (divisor_parts, total_coeff) = self.extract_coefficient(&solve_for_terms, solve_for);

        // The solution is: solve_for = -other_terms / coefficient
        // Negate other_terms
        let negated_other: Vec<Term> = other_terms
            .into_iter()
            .map(|mut t| {
                t.coefficient = -t.coefficient;
                t
            })
            .collect();

        // Build numerator expression
        let numerator = self.build_rust_expression(&negated_other, total_coeff);

        // Build divisor expression (if not a constant)
        let divisor = if divisor_parts.is_empty() {
            None
        } else {
            Some(divisor_parts.join(" * "))
        };

        Ok(SolveResult {
            variable: solve_for.to_string(),
            numerator,
            divisor,
            constraint: constraint.to_string(),
        })
    }

    /// Parse a constraint string into terms.
    ///
    /// Input: "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
    /// Output: Vec of Terms
    fn parse_constraint(&self, constraint: &str) -> Result<Vec<Term>, SolveError> {
        // Split on '='
        let parts: Vec<&str> = constraint.split('=').collect();
        if parts.len() != 2 {
            return Err(SolveError::MissingEquality(constraint.to_string()));
        }

        let lhs = parts[0].trim();
        let rhs = parts[1].trim();

        // For constraints in form "expr = 0", we just parse the lhs
        // For constraints in form "expr = value", we'd need to handle differently
        if rhs != "0" {
            return Err(SolveError::ParseError(format!(
                "expected '= 0' but got '= {}'",
                rhs
            )));
        }

        self.parse_expression(lhs)
    }

    /// Parse an expression into terms.
    fn parse_expression(&self, expr: &str) -> Result<Vec<Term>, SolveError> {
        let mut terms = Vec::new();
        let mut current_term = String::new();
        let mut sign = 1;

        // Normalize: remove spaces around operators
        let normalized = expr.replace(" ", "");

        let chars: Vec<char> = normalized.chars().collect();
        let mut i = 0;

        while i < chars.len() {
            let c = chars[i];

            match c {
                '+' => {
                    if !current_term.is_empty() {
                        terms.push(self.parse_term(&current_term, sign)?);
                        current_term.clear();
                    }
                    sign = 1;
                }
                '-' => {
                    if !current_term.is_empty() {
                        terms.push(self.parse_term(&current_term, sign)?);
                        current_term.clear();
                    }
                    sign = -1;
                }
                _ => {
                    current_term.push(c);
                }
            }
            i += 1;
        }

        // Don't forget the last term
        if !current_term.is_empty() {
            terms.push(self.parse_term(&current_term, sign)?);
        }

        Ok(terms)
    }

    /// Parse a single term like "2*s*e0123" into a Term.
    fn parse_term(&self, term: &str, sign: i32) -> Result<Term, SolveError> {
        let factors: Vec<&str> = term.split('*').collect();

        let mut coefficient = sign;
        let mut variables = Vec::new();

        for factor in factors {
            let factor = factor.trim();
            if factor.is_empty() {
                continue;
            }

            // Check if it's a number
            if let Ok(num) = factor.parse::<i32>() {
                coefficient *= num;
            } else {
                variables.push(factor.to_string());
            }
        }

        Ok(Term {
            coefficient,
            variables,
        })
    }

    /// Extract the coefficient (other variables) of the solve_for variable.
    ///
    /// For terms like [2*s*e0123], extracting e0123's coefficient gives:
    /// - divisor_parts: ["s"]
    /// - total_coeff: 2
    fn extract_coefficient(&self, terms: &[Term], solve_for: &str) -> (Vec<String>, i32) {
        // For now, assume there's only one term containing solve_for
        // (our constraints are simple enough for this)
        if let Some(term) = terms.first() {
            let other_vars: Vec<String> = term
                .variables
                .iter()
                .filter(|v| *v != solve_for)
                .cloned()
                .collect();

            (other_vars, term.coefficient)
        } else {
            (Vec::new(), 1)
        }
    }

    /// Build a Rust expression from terms.
    ///
    /// Takes into account the coefficient divisor to simplify expressions.
    fn build_rust_expression(&self, terms: &[Term], divisor_coeff: i32) -> String {
        if terms.is_empty() {
            return "T::zero()".to_string();
        }

        let parts: Vec<String> = terms
            .iter()
            .map(|term| {
                // Simplify coefficient if divisor matches
                let simplified_coeff =
                    if divisor_coeff != 0 && term.coefficient % divisor_coeff == 0 {
                        term.coefficient / divisor_coeff
                    } else {
                        term.coefficient
                    };

                let vars_expr = if term.variables.is_empty() {
                    format!("T::from_i8({})", simplified_coeff)
                } else {
                    term.variables.join(" * ")
                };

                match simplified_coeff {
                    1 => vars_expr,
                    -1 => format!("-{}", vars_expr),
                    _ if term.variables.is_empty() => format!("T::from_i8({})", simplified_coeff),
                    _ => format!("T::from_i8({}) * {}", simplified_coeff, vars_expr),
                }
            })
            .collect();

        // Join with proper operators
        let mut result = String::new();
        for (i, part) in parts.iter().enumerate() {
            if i == 0 {
                result.push_str(part);
            } else if part.starts_with('-') {
                result.push_str(" - ");
                result.push_str(&part[1..]);
            } else {
                result.push_str(" + ");
                result.push_str(part);
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_term() {
        let solver = ConstraintSolver::new();

        let term = solver.parse_term("2*s*e0123", 1).unwrap();
        assert_eq!(term.coefficient, 2);
        assert_eq!(term.variables, vec!["s", "e0123"]);

        let term = solver.parse_term("e12*e03", -1).unwrap();
        assert_eq!(term.coefficient, -1);
        assert_eq!(term.variables, vec!["e12", "e03"]);
    }

    #[test]
    fn test_parse_expression() {
        let solver = ConstraintSolver::new();

        let terms = solver
            .parse_expression("2*s*e0123 - 2*e12*e03 + 2*e13*e02")
            .unwrap();
        assert_eq!(terms.len(), 3);
        assert_eq!(terms[0].coefficient, 2);
        assert_eq!(terms[1].coefficient, -2);
        assert_eq!(terms[2].coefficient, 2);
    }

    #[test]
    fn test_solve_motor_constraint() {
        let solver = ConstraintSolver::new();

        let result = solver
            .solve("2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0", "e0123")
            .unwrap();

        assert_eq!(result.variable, "e0123");
        assert_eq!(result.divisor, Some("s".to_string()));
        // Numerator should be: e12*e03 - e13*e02 + e23*e01
        assert!(result.numerator.contains("e12 * e03"));
        assert!(result.numerator.contains("e13 * e02"));
        assert!(result.numerator.contains("e23 * e01"));
    }

    #[test]
    fn test_solve_bivector_constraint() {
        let solver = ConstraintSolver::new();

        let result = solver
            .solve("-2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0", "e03")
            .unwrap();

        assert_eq!(result.variable, "e03");
        assert_eq!(result.divisor, Some("e12".to_string()));
    }

    #[test]
    fn test_variable_not_found() {
        let solver = ConstraintSolver::new();

        let result = solver.solve("2*s*e0123 = 0", "nonexistent");
        assert!(matches!(result, Err(SolveError::VariableNotFound(_))));
    }
}
