//! Field constraint derivation for geometric entities.
//!
//! When a grade combination doesn't automatically satisfy geometric constraints,
//! this module derives the field constraint expression that must hold.
//!
//! There are two types of constraints:
//!
//! 1. **Geometric Product Constraint**: `u * ũ = scalar`
//!    For example, PGA bivectors (lines) have the constraint:
//!    `e01*e23 + e02*e31 + e03*e12 = 0` (direction · moment = 0)
//!
//! 2. **Antiproduct Constraint**: `u ⊟ ũ̃ = antiscalar`
//!    This is the dual constraint that must also be satisfied.

use std::collections::HashMap;

use crate::algebra::{
    Algebra, Blade, ProductTable, antireverse_sign, blades_of_grades, grade, reverse_sign,
};

/// Derives the field constraint expression for a grade combination.
///
/// If the grade combination automatically satisfies geometric constraints,
/// returns `None`. Otherwise, returns the constraint expression that must
/// equal zero for the constraint to hold.
///
/// # Arguments
///
/// * `grades` - The grades present in the type
/// * `algebra` - The algebra definition
///
/// # Returns
///
/// `None` if no constraint needed, or `Some(expression)` where the expression
/// must equal zero. The expression uses blade names from the algebra.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::derive_field_constraint;
/// use clifford_codegen::algebra::Algebra;
///
/// let algebra = Algebra::pga(3);
///
/// // PGA bivectors need a constraint
/// let constraint = derive_field_constraint(&[2], &algebra);
/// assert!(constraint.is_some());
/// // Returns something like "e01*e23 + e02*e31 + e03*e12 = 0"
/// ```
pub fn derive_field_constraint(grades: &[usize], algebra: &Algebra) -> Option<String> {
    let table = ProductTable::new(algebra);
    let blades = blades_of_grades(algebra.dim(), grades);

    // Collect non-scalar contributions from u * reverse(u)
    // For each non-scalar output blade, collect the terms that contribute to it
    let mut non_scalar_terms: HashMap<usize, Vec<(usize, usize, i32)>> = HashMap::new();

    for (i, &a) in blades.iter().enumerate() {
        for &b in &blades[i..] {
            let ga = grade(a);
            let gb = grade(b);
            let rev_a = reverse_sign(ga);
            let rev_b = reverse_sign(gb);

            let (sign_ab, result_ab) = table.geometric(a, b);

            if grade(result_ab) == 0 {
                continue; // Scalar output, no constraint needed
            }

            if a == b {
                // Same blade producing non-scalar - this is a fundamental problem
                // that can't be solved with field constraints
                if sign_ab != 0 {
                    return None; // Can't derive a useful constraint
                }
            } else {
                // Different blades producing non-scalar
                let (sign_ba, _) = table.geometric(b, a);

                // Total coefficient for result blade when terms don't cancel
                let coef =
                    i32::from(sign_ab) * i32::from(rev_b) + i32::from(sign_ba) * i32::from(rev_a);

                if coef != 0 {
                    // This pair contributes to a non-scalar result
                    non_scalar_terms
                        .entry(result_ab)
                        .or_default()
                        .push((a, b, coef));
                }
            }
        }
    }

    if non_scalar_terms.is_empty() {
        return None; // No constraint needed
    }

    // Build constraint expression
    // For each non-scalar output blade, the sum of contributions must be zero
    // We combine all into one constraint expression
    let mut terms: Vec<String> = Vec::new();

    for contributions in non_scalar_terms.values() {
        for &(a, b, coef) in contributions {
            let blade_a = Blade::from_index(a);
            let blade_b = Blade::from_index(b);
            let name_a = algebra.blade_name(blade_a);
            let name_b = algebra.blade_name(blade_b);

            let term = if coef == 1 {
                format!("{}*{}", name_a, name_b)
            } else if coef == -1 {
                format!("-{}*{}", name_a, name_b)
            } else {
                // Both positive and negative coefficients with magnitude > 1
                format!("{}*{}*{}", coef, name_a, name_b)
            };
            terms.push(term);
        }
    }

    if terms.is_empty() {
        return None;
    }

    // Join terms with + (handling leading -)
    let mut expr = String::new();
    for (i, term) in terms.iter().enumerate() {
        if i == 0 {
            expr.push_str(term);
        } else if let Some(stripped) = term.strip_prefix('-') {
            expr.push_str(" - ");
            expr.push_str(stripped);
        } else {
            expr.push_str(" + ");
            expr.push_str(term);
        }
    }

    Some(format!("{} = 0", expr))
}

/// Derives the antiproduct field constraint expression for a grade combination.
///
/// The antiproduct constraint requires that `u ⊟ ũ̃` produces only an antiscalar
/// (grade n). This function derives the field constraint that must hold for this
/// to be true.
///
/// The antiproduct is computed as: `a ⊟ b = dual(dual(a) * dual(b))`
///
/// # Arguments
///
/// * `grades` - The grades present in the type
/// * `algebra` - The algebra definition
///
/// # Returns
///
/// `None` if no constraint needed, or `Some(expression)` where the expression
/// must equal zero.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::derive_antiproduct_constraint;
/// use clifford_codegen::algebra::Algebra;
///
/// let algebra = Algebra::pga(3);
///
/// // PGA motors need an antiproduct constraint
/// let constraint = derive_antiproduct_constraint(&[0, 2, 4], &algebra);
/// // Returns the constraint expression if one is needed
/// ```
pub fn derive_antiproduct_constraint(grades: &[usize], algebra: &Algebra) -> Option<String> {
    let table = ProductTable::new(algebra);
    let dim = algebra.dim();
    let blades = blades_of_grades(dim, grades);
    let pseudoscalar = (1 << dim) - 1; // All bits set = e123...n

    // Collect non-antiscalar contributions from u ⊟ antireverse(u)
    // For each non-antiscalar output blade, collect the terms that contribute to it
    let mut non_antiscalar_terms: HashMap<usize, Vec<(usize, usize, i32)>> = HashMap::new();

    for (i, &a) in blades.iter().enumerate() {
        for &b in &blades[i..] {
            let ga = grade(a);
            let gb = grade(b);
            let antirev_a = antireverse_sign(ga, dim);
            let antirev_b = antireverse_sign(gb, dim);

            // Compute dual blades
            let dual_a = pseudoscalar ^ a;
            let dual_b = pseudoscalar ^ b;

            // Product of duals
            let (sign_ab, dual_result_ab) = table.geometric(dual_a, dual_b);
            // Dual back to get antiproduct result
            let result_ab = pseudoscalar ^ dual_result_ab;

            if grade(result_ab) == dim {
                continue; // Antiscalar output (grade n), no constraint needed
            }

            if a == b {
                // Same blade producing non-antiscalar - fundamental problem
                if sign_ab != 0 {
                    return None; // Can't derive a useful constraint
                }
            } else {
                // Different blades producing non-antiscalar
                let (sign_ba, _) = table.geometric(dual_b, dual_a);

                // Total coefficient for result blade when terms don't cancel
                let coef = i32::from(sign_ab) * i32::from(antirev_b)
                    + i32::from(sign_ba) * i32::from(antirev_a);

                if coef != 0 {
                    // This pair contributes to a non-antiscalar result
                    non_antiscalar_terms
                        .entry(result_ab)
                        .or_default()
                        .push((a, b, coef));
                }
            }
        }
    }

    if non_antiscalar_terms.is_empty() {
        return None; // No constraint needed
    }

    // Build constraint expression
    let mut terms: Vec<String> = Vec::new();

    for contributions in non_antiscalar_terms.values() {
        for &(a, b, coef) in contributions {
            let blade_a = Blade::from_index(a);
            let blade_b = Blade::from_index(b);
            let name_a = algebra.blade_name(blade_a);
            let name_b = algebra.blade_name(blade_b);

            let term = if coef == 1 {
                format!("{}*{}", name_a, name_b)
            } else if coef == -1 {
                format!("-{}*{}", name_a, name_b)
            } else {
                format!("{}*{}*{}", coef, name_a, name_b)
            };
            terms.push(term);
        }
    }

    if terms.is_empty() {
        return None;
    }

    // Join terms with + (handling leading -)
    let mut expr = String::new();
    for (i, term) in terms.iter().enumerate() {
        if i == 0 {
            expr.push_str(&term);
        } else if let Some(stripped) = term.strip_prefix('-') {
            expr.push_str(" - ");
            expr.push_str(stripped);
        } else {
            expr.push_str(" + ");
            expr.push_str(&term);
        }
    }

    Some(format!("{} = 0", expr))
}

/// Checks if a grade combination can satisfy constraints with field constraints.
///
/// Some grade combinations fundamentally cannot satisfy geometric constraints
/// (e.g., a blade whose square is non-scalar). This function distinguishes between:
/// - Combinations that automatically satisfy constraints (returns `(true, None)`)
/// - Combinations that can satisfy with field constraints (returns `(true, Some(constraint))`)
/// - Combinations that cannot satisfy constraints (returns `(false, None)`)
///
/// # Arguments
///
/// * `grades` - The grades present in the type
/// * `algebra` - The algebra definition
///
/// # Returns
///
/// `(satisfiable, constraint)` where:
/// - `satisfiable` is true if constraints can be satisfied (with or without field constraints)
/// - `constraint` is the field constraint expression, if one is needed
pub fn can_satisfy_constraints(grades: &[usize], algebra: &Algebra) -> (bool, Option<String>) {
    let table = ProductTable::new(algebra);
    let blades = blades_of_grades(algebra.dim(), grades);

    // Check: are there any same-blade non-scalar products?
    // This indicates a fundamental problem that field constraints can't solve
    // (a blade whose square is non-scalar violates the constraint regardless of coefficient)
    for &a in &blades {
        let (sign_aa, result_aa) = table.geometric(a, a);
        if sign_aa != 0 && grade(result_aa) != 0 {
            return (false, None); // Fundamental violation
        }
    }

    // All other combinations can potentially satisfy constraints with field constraints
    // Derive the constraint expression (if any is needed)
    let constraint = derive_field_constraint(grades, algebra);
    (true, constraint)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euclidean_3d_no_constraints() {
        let algebra = Algebra::euclidean(3);

        // Single grades need no constraints
        assert_eq!(derive_field_constraint(&[0], &algebra), None);
        assert_eq!(derive_field_constraint(&[1], &algebra), None);
        assert_eq!(derive_field_constraint(&[2], &algebra), None);
        assert_eq!(derive_field_constraint(&[3], &algebra), None);

        // Even/odd subalgebras need no constraints
        assert_eq!(derive_field_constraint(&[0, 2], &algebra), None);
        assert_eq!(derive_field_constraint(&[1, 3], &algebra), None);
    }

    #[test]
    fn pga_3d_bivector_constraint() {
        let algebra = Algebra::pga(3);

        // Bivectors need a constraint
        let constraint = derive_field_constraint(&[2], &algebra);
        assert!(constraint.is_some());

        let expr = constraint.unwrap();
        // Should contain products of bivector blades
        assert!(expr.contains("="));
    }

    #[test]
    fn pga_3d_satisfiability() {
        let algebra = Algebra::pga(3);

        // Scalar - satisfiable, no constraint
        let (sat, constr) = can_satisfy_constraints(&[0], &algebra);
        assert!(sat);
        assert!(constr.is_none());

        // Bivector - satisfiable with constraint
        let (sat, constr) = can_satisfy_constraints(&[2], &algebra);
        assert!(sat);
        assert!(constr.is_some());

        // Scalar + Vector - satisfiable with constraint (constraint requires s * v_i = 0)
        let (sat, constr) = can_satisfy_constraints(&[0, 1], &algebra);
        assert!(sat);
        assert!(constr.is_some());

        // Vector + Trivector (Flector) - satisfiable with constraint
        let (sat, constr) = can_satisfy_constraints(&[1, 3], &algebra);
        assert!(sat);
        assert!(constr.is_some());
    }
}
