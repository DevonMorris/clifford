//! Constraint solver for computing dependent coefficients.
//!
//! This module solves constraint equations for a specified variable,
//! generating the Rust expression to compute that variable from the others.
//!
//! Supports two types of constraints:
//!
//! 1. **Linear**: `coeff1*var1*var2 + coeff2*var3*var4 + ... = 0`
//!    Solution is direct division.
//!
//! 2. **Quadratic**: `var*var + other_terms = constant`
//!    Solution requires square root, with domain restrictions.

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

    /// Variable appears in a form that can't be solved.
    #[error("variable '{0}' appears in a form that cannot be algebraically solved")]
    UnsolvableForm(String),
}

/// The type of solution for a constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolutionType {
    /// Linear solution: var = expr / divisor
    Linear,
    /// Quadratic solution: var = sqrt(expr), has domain restrictions
    Quadratic,
}

/// Result of solving a constraint for a variable.
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// The variable being solved for.
    pub variable: String,
    /// The Rust expression for the numerator (or sqrt argument for quadratic).
    pub numerator: String,
    /// The Rust expression for the divisor (if division is needed).
    pub divisor: Option<String>,
    /// The type of solution (affects code generation).
    pub solution_type: SolutionType,
    /// Whether to use positive or negative root (for quadratic).
    pub positive_root: bool,
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
    /// Handles both linear and quadratic constraints:
    ///
    /// **Linear**: `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`
    /// solving for `e0123` yields: `e0123 = (e12*e03 - e13*e02 + e23*e01) / s`
    ///
    /// **Quadratic**: `s*s + e12*e12 + e13*e13 + e23*e23 = 1`
    /// solving for `s` yields: `s = sqrt(1 - e12*e12 - e13*e13 - e23*e23)`
    ///
    /// # Arguments
    ///
    /// * `constraint` - The constraint expression string
    /// * `solve_for` - The variable to solve for
    /// * `positive_root` - Whether to use positive root for quadratic (default true)
    ///
    /// # Returns
    ///
    /// A `SolveResult` containing the solution expression and type.
    pub fn solve(&self, constraint: &str, solve_for: &str) -> Result<SolveResult, SolveError> {
        self.solve_with_sign(constraint, solve_for, true)
    }

    /// Solves a constraint with explicit sign convention for quadratic solutions.
    pub fn solve_with_sign(
        &self,
        constraint: &str,
        solve_for: &str,
        positive_root: bool,
    ) -> Result<SolveResult, SolveError> {
        // Parse constraint into terms and RHS constant
        let (terms, rhs_constant) = self.parse_constraint(constraint)?;

        // Check if variable appears squared (quadratic)
        let is_quadratic = self.is_quadratic_in_variable(&terms, solve_for);

        if is_quadratic {
            self.solve_quadratic(&terms, rhs_constant, solve_for, positive_root, constraint)
        } else {
            self.solve_linear(&terms, solve_for, constraint)
        }
    }

    /// Checks if a variable appears squared in any term.
    fn is_quadratic_in_variable(&self, terms: &[Term], var: &str) -> bool {
        for term in terms {
            let count = term.variables.iter().filter(|v| *v == var).count();
            if count >= 2 {
                return true;
            }
        }
        false
    }

    /// Solves a linear constraint.
    fn solve_linear(
        &self,
        terms: &[Term],
        solve_for: &str,
        constraint: &str,
    ) -> Result<SolveResult, SolveError> {
        // Separate terms containing solve_for from the rest
        let mut solve_for_terms: Vec<Term> = Vec::new();
        let mut other_terms: Vec<Term> = Vec::new();

        for term in terms {
            if term.variables.contains(&solve_for.to_string()) {
                solve_for_terms.push(term.clone());
            } else {
                other_terms.push(term.clone());
            }
        }

        if solve_for_terms.is_empty() {
            return Err(SolveError::VariableNotFound(solve_for.to_string()));
        }

        // Extract the coefficient of solve_for
        let (divisor_parts, total_coeff) = self.extract_coefficient(&solve_for_terms, solve_for);

        // The solution is: solve_for = -other_terms / coefficient
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
            solution_type: SolutionType::Linear,
            positive_root: true,
            constraint: constraint.to_string(),
        })
    }

    /// Solves a quadratic constraint of form: coeff*var*var + other_terms = constant
    fn solve_quadratic(
        &self,
        terms: &[Term],
        rhs_constant: i32,
        solve_for: &str,
        positive_root: bool,
        constraint: &str,
    ) -> Result<SolveResult, SolveError> {
        // Find the term with var*var
        let mut squared_coeff = 0i32;
        let mut other_terms: Vec<Term> = Vec::new();

        for term in terms {
            let var_count = term.variables.iter().filter(|v| *v == solve_for).count();

            if var_count == 2 && term.variables.len() == 2 {
                // Pure var*var term
                squared_coeff += term.coefficient;
            } else if var_count == 0 {
                other_terms.push(term.clone());
            } else {
                // Mixed term (var appears but not squared alone) - can't solve simply
                return Err(SolveError::UnsolvableForm(solve_for.to_string()));
            }
        }

        if squared_coeff == 0 {
            return Err(SolveError::VariableNotFound(solve_for.to_string()));
        }

        // Solution: var = sqrt((constant - other_terms) / squared_coeff)
        // For typical unit norm constraints where squared_coeff = 1:
        // var = sqrt(constant - other_terms)

        // Build the sqrt argument: (constant - other_terms) / squared_coeff
        // First negate other_terms
        let negated_other: Vec<Term> = other_terms
            .into_iter()
            .map(|mut t| {
                t.coefficient = -t.coefficient;
                t
            })
            .collect();

        // Build expression for constant + negated_other
        let mut sqrt_arg = if rhs_constant != 0 {
            if squared_coeff == 1 {
                format!("T::from_i8({})", rhs_constant)
            } else {
                format!(
                    "T::from_i8({}) / T::from_i8({})",
                    rhs_constant, squared_coeff
                )
            }
        } else {
            String::new()
        };

        // Add negated other terms
        for term in &negated_other {
            let term_expr = self.term_to_rust_expression(term, squared_coeff);
            if sqrt_arg.is_empty() {
                sqrt_arg = term_expr;
            } else if let Some(stripped) = term_expr.strip_prefix('-') {
                sqrt_arg = format!("{} - {}", sqrt_arg, stripped);
            } else {
                sqrt_arg = format!("{} + {}", sqrt_arg, term_expr);
            }
        }

        if sqrt_arg.is_empty() {
            sqrt_arg = "T::zero()".to_string();
        }

        Ok(SolveResult {
            variable: solve_for.to_string(),
            numerator: sqrt_arg,
            divisor: None,
            solution_type: SolutionType::Quadratic,
            positive_root,
            constraint: constraint.to_string(),
        })
    }

    /// Converts a single term to a Rust expression.
    fn term_to_rust_expression(&self, term: &Term, divisor_coeff: i32) -> String {
        let simplified_coeff = if divisor_coeff != 0 && term.coefficient % divisor_coeff == 0 {
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
    }

    /// Parse a constraint string into terms and RHS constant.
    ///
    /// Input: "2*s*e0123 - 2*e12*e03 = 0" or "s*s + b*b = 1"
    /// Output: (Vec of Terms, RHS constant)
    fn parse_constraint(&self, constraint: &str) -> Result<(Vec<Term>, i32), SolveError> {
        // Split on '='
        let parts: Vec<&str> = constraint.split('=').collect();
        if parts.len() != 2 {
            return Err(SolveError::MissingEquality(constraint.to_string()));
        }

        let lhs = parts[0].trim();
        let rhs = parts[1].trim();

        // Parse RHS as a constant
        let rhs_constant: i32 = rhs.parse().map_err(|_| {
            SolveError::ParseError(format!("RHS must be an integer constant, got '{}'", rhs))
        })?;

        let terms = self.parse_expression(lhs)?;
        Ok((terms, rhs_constant))
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
            } else if let Some(stripped) = part.strip_prefix('-') {
                result.push_str(" - ");
                result.push_str(stripped);
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
        assert_eq!(result.solution_type, SolutionType::Linear);
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
        assert_eq!(result.solution_type, SolutionType::Linear);
    }

    #[test]
    fn test_solve_quadratic_unit_norm() {
        let solver = ConstraintSolver::new();

        // Unit norm constraint: s*s + e12*e12 + e13*e13 + e23*e23 = 1
        let result = solver
            .solve("s*s + e12*e12 + e13*e13 + e23*e23 = 1", "s")
            .unwrap();

        assert_eq!(result.variable, "s");
        assert_eq!(result.solution_type, SolutionType::Quadratic);
        assert!(result.positive_root);
        assert!(result.divisor.is_none());
        // sqrt argument should be: 1 - e12*e12 - e13*e13 - e23*e23
        assert!(result.numerator.contains("T::from_i8(1)"));
        assert!(result.numerator.contains("e12 * e12"));
        assert!(result.numerator.contains("e13 * e13"));
        assert!(result.numerator.contains("e23 * e23"));
    }

    #[test]
    fn test_solve_quadratic_negative_root() {
        let solver = ConstraintSolver::new();

        let result = solver.solve_with_sign("a*a + b*b = 1", "a", false).unwrap();

        assert_eq!(result.variable, "a");
        assert_eq!(result.solution_type, SolutionType::Quadratic);
        assert!(!result.positive_root);
    }

    #[test]
    fn test_solve_simple_quadratic() {
        let solver = ConstraintSolver::new();

        // Simple case: x*x = 1
        let result = solver.solve("x*x = 1", "x").unwrap();

        assert_eq!(result.variable, "x");
        assert_eq!(result.solution_type, SolutionType::Quadratic);
        assert_eq!(result.numerator, "T::from_i8(1)");
    }

    #[test]
    fn test_variable_not_found() {
        let solver = ConstraintSolver::new();

        let result = solver.solve("2*s*e0123 = 0", "nonexistent");
        assert!(matches!(result, Err(SolveError::VariableNotFound(_))));
    }
}
