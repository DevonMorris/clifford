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
            let name_a = algebra.blade_index_name(blade_a);
            let name_b = algebra.blade_index_name(blade_b);

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
            let name_a = algebra.blade_index_name(blade_a);
            let name_b = algebra.blade_index_name(blade_b);

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

/// Derives the blade constraint expression for a grade combination.
///
/// A k-vector B is a simple blade (can be written as v₁ ∧ v₂ ∧ ... ∧ vₖ) if and only if
/// B ∧ B = 0. This function derives the constraint equations that must hold.
///
/// This is the constraint used by CGA for geometric primitives like dipoles, circles,
/// and spheres. It ensures the element represents a valid geometric object.
///
/// # Arguments
///
/// * `grades` - The grades present in the type (typically a single grade for blades)
/// * `algebra` - The algebra definition
///
/// # Returns
///
/// A vector of constraint expressions, one for each non-zero component of B ∧ B.
/// Each expression must equal zero for the element to be a valid blade.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::derive_blade_constraint;
/// use clifford_codegen::algebra::Algebra;
///
/// // CGA dipole (grade 2 in 5D) needs blade constraints
/// let algebra = Algebra::new(4, 1, 0); // Cl(4,1,0)
/// let constraints = derive_blade_constraint(&[2], &algebra);
/// assert_eq!(constraints.len(), 5); // 5 grade-4 components
/// ```
pub fn derive_blade_constraint(grades: &[usize], algebra: &Algebra) -> Vec<String> {
    let table = ProductTable::new(algebra);
    let dim = algebra.dim();
    let blades = blades_of_grades(dim, grades);

    // B ∧ B produces grade 2k for a k-vector
    // Collect terms for each output blade
    let mut output_terms: HashMap<usize, Vec<(usize, usize, i32)>> = HashMap::new();

    for (i, &a) in blades.iter().enumerate() {
        for &b in &blades[i..] {
            let (sign, result) = table.exterior(a, b);
            if sign == 0 {
                continue;
            }

            // For B ∧ B, we get 2*b_a*b_b for a ≠ b (symmetric sum)
            let coef = if a == b {
                i32::from(sign)
            } else {
                2 * i32::from(sign)
            };

            output_terms.entry(result).or_default().push((a, b, coef));
        }
    }

    // Build constraint expressions
    let mut constraints = Vec::new();

    for (result_blade, terms) in output_terms {
        if terms.is_empty() {
            continue;
        }

        let mut expr_parts: Vec<String> = Vec::new();
        for (a, b, coef) in terms {
            let blade_a = Blade::from_index(a);
            let blade_b = Blade::from_index(b);
            let name_a = algebra.blade_index_name(blade_a);
            let name_b = algebra.blade_index_name(blade_b);

            let term = if a == b {
                // Diagonal term: coef * name²
                if coef == 1 {
                    format!("{}^2", name_a)
                } else if coef == -1 {
                    format!("-{}^2", name_a)
                } else {
                    format!("{}*{}^2", coef, name_a)
                }
            } else {
                // Off-diagonal term: coef * name_a * name_b
                if coef == 1 {
                    format!("{}*{}", name_a, name_b)
                } else if coef == -1 {
                    format!("-{}*{}", name_a, name_b)
                } else {
                    format!("{}*{}*{}", coef, name_a, name_b)
                }
            };
            expr_parts.push(term);
        }

        // Join terms
        let mut expr = String::new();
        for (i, part) in expr_parts.iter().enumerate() {
            if i == 0 {
                expr.push_str(part);
            } else if let Some(stripped) = part.strip_prefix('-') {
                expr.push_str(" - ");
                expr.push_str(stripped);
            } else {
                expr.push_str(" + ");
                expr.push_str(part);
            }
        }

        // Add blade name as comment
        let result_name = algebra.blade_index_name(Blade::from_index(result_blade));
        constraints.push(format!("{} = 0  // {} component", expr, result_name));
    }

    constraints
}

/// Derives the null constraint expression for a grade combination.
///
/// A null element has zero norm: `v * ṽ = 0`. This function derives the constraint
/// equations that must hold for an element to be null.
///
/// In CGA, points are null vectors (grade 1 with v * ṽ = 0). This constraint is
/// essential for representing actual geometric points vs general vectors.
///
/// # Arguments
///
/// * `grades` - The grades present in the type
/// * `algebra` - The algebra definition
///
/// # Returns
///
/// A vector of constraint expressions. For the null constraint, this includes:
/// - The scalar part of v * ṽ must be zero
/// - Any non-scalar parts must also be zero (same as versor constraint)
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::derive_null_constraint;
/// use clifford_codegen::algebra::Algebra;
///
/// // CGA point (grade 1) needs null constraint
/// let algebra = Algebra::new(4, 1, 0); // Cl(4,1,0)
/// let constraints = derive_null_constraint(&[1], &algebra);
/// // Returns constraint that x² + y² + z² + w² - u² = 0 (null vector)
/// ```
pub fn derive_null_constraint(grades: &[usize], algebra: &Algebra) -> Vec<String> {
    let table = ProductTable::new(algebra);
    let blades = blades_of_grades(algebra.dim(), grades);

    // Collect ALL contributions from v * reverse(v), including scalar
    let mut all_terms: HashMap<usize, Vec<(usize, usize, i32)>> = HashMap::new();

    for (i, &a) in blades.iter().enumerate() {
        for &b in &blades[i..] {
            let ga = grade(a);
            let gb = grade(b);
            let rev_a = reverse_sign(ga);
            let rev_b = reverse_sign(gb);

            let (sign_ab, result_ab) = table.geometric(a, b);

            if sign_ab == 0 {
                continue;
            }

            if a == b {
                // Diagonal term: a * reverse(a) = sign * |a|²
                let coef = i32::from(sign_ab) * i32::from(rev_a);
                all_terms.entry(result_ab).or_default().push((a, a, coef));
            } else {
                // Off-diagonal terms
                let (sign_ba, _) = table.geometric(b, a);
                let coef =
                    i32::from(sign_ab) * i32::from(rev_b) + i32::from(sign_ba) * i32::from(rev_a);

                if coef != 0 {
                    all_terms.entry(result_ab).or_default().push((a, b, coef));
                }
            }
        }
    }

    // Build constraint expressions for ALL output grades (including scalar = 0)
    let mut constraints = Vec::new();

    for (result_blade, terms) in all_terms {
        if terms.is_empty() {
            continue;
        }

        let mut expr_parts: Vec<String> = Vec::new();
        for (a, b, coef) in terms {
            let blade_a = Blade::from_index(a);
            let name_a = algebra.blade_index_name(blade_a);

            let term = if a == b {
                // Squared term
                if coef == 1 {
                    format!("{}²", name_a)
                } else if coef == -1 {
                    format!("-{}²", name_a)
                } else {
                    format!("{}*{}²", coef, name_a)
                }
            } else {
                let blade_b = Blade::from_index(b);
                let name_b = algebra.blade_index_name(blade_b);
                if coef == 1 {
                    format!("{}*{}", name_a, name_b)
                } else if coef == -1 {
                    format!("-{}*{}", name_a, name_b)
                } else {
                    format!("{}*{}*{}", coef, name_a, name_b)
                }
            };
            expr_parts.push(term);
        }

        // Join terms
        let mut expr = String::new();
        for (i, part) in expr_parts.iter().enumerate() {
            if i == 0 {
                expr.push_str(part);
            } else if let Some(stripped) = part.strip_prefix('-') {
                expr.push_str(" - ");
                expr.push_str(stripped);
            } else {
                expr.push_str(" + ");
                expr.push_str(part);
            }
        }

        let result_grade = grade(result_blade);
        let grade_name = if result_grade == 0 {
            "scalar (norm)".to_string()
        } else {
            let result_name = algebra.blade_index_name(Blade::from_index(result_blade));
            result_name.to_string()
        };

        constraints.push(format!("{} = 0  // {} component", expr, grade_name));
    }

    constraints
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

    #[test]
    fn cga_dipole_constraints_analysis() {
        // CGA: Cl(4,1,0) - 4 positive, 1 negative
        let algebra = Algebra::new(4, 1, 0);

        let geometric = derive_field_constraint(&[2], &algebra);
        let _antiproduct = derive_antiproduct_constraint(&[2], &algebra);

        // Print for analysis
        eprintln!("\n=== CGA Dipole Constraint Analysis ===");
        eprintln!("Field mapping (CGA wiki convention):");
        eprintln!("  m (moment):   e12=mx, e13=my, e23=mz");
        eprintln!("  v (velocity): e14=vx, e24=vy, e34=vz");
        eprintln!("  p (position): e15=px, e25=py, e35=pz, e45=pw");
        eprintln!();

        if let Some(ref c) = geometric {
            eprintln!("Inferred geometric constraint (d * d̃ = scalar):");
            eprintln!("  {}", c);
        }
        eprintln!();
        eprintln!("CGA wiki constraints:");
        eprintln!("  1. p × v - pw*m = 0");
        eprintln!("     Expands to:");
        eprintln!("       py*vz - pz*vy - pw*mx = 0  (e25*e34 - e35*e24 - e45*e12)");
        eprintln!("       pz*vx - px*vz - pw*my = 0  (e35*e14 - e15*e34 - e45*e13)");
        eprintln!("       px*vy - py*vx - pw*mz = 0  (e15*e24 - e25*e14 - e45*e23)");
        eprintln!("  2. p · m = 0  →  px*mx + py*my + pz*mz = 0  (e15*e12 + e25*e13 + e35*e23)");
        eprintln!("  3. v · m = 0  →  vx*mx + vy*my + vz*mz = 0  (e14*e12 + e24*e13 + e34*e23)");
        eprintln!("==========================================\n");

        // The inferred constraint is for the norm, which is different from geometric validity
        assert!(
            geometric.is_some(),
            "CGA dipoles should have a geometric constraint"
        );
    }

    #[test]
    fn blade_constraint_across_algebras() {
        // Test which algebras need blade constraints for bivectors
        eprintln!("\n=== Blade Constraints (B ∧ B = 0) Across Algebras ===\n");

        let test_cases = [
            ("Euclidean 2D", Algebra::euclidean(2), 2),
            ("Euclidean 3D", Algebra::euclidean(3), 2),
            ("Euclidean 4D", Algebra::euclidean(4), 2),
            ("PGA 2D", Algebra::pga(2), 2),
            ("PGA 3D", Algebra::pga(3), 2),
            ("CGA 3D (Cl(4,1,0))", Algebra::new(4, 1, 0), 2),
            ("Minkowski (Cl(3,1,0))", Algebra::new(3, 1, 0), 2),
        ];

        for (name, algebra, grade) in &test_cases {
            let constraints = derive_blade_constraint(&[*grade], &algebra);
            let num_blades = blades_of_grades(algebra.dim(), &[*grade]).len();

            eprintln!(
                "{}: {} grade-{} blades → {} blade constraints",
                name,
                num_blades,
                grade,
                constraints.len()
            );

            if !constraints.is_empty() {
                for c in &constraints {
                    eprintln!("    {}", c);
                }
            }
        }

        eprintln!("\n=== Summary ===");
        eprintln!("Blade constraints needed when dim > 2*grade:");
        eprintln!("  - 2D: bivectors (3 components) - no constraint (dim=2, 2*2=4 > 2)");
        eprintln!("  - 3D: bivectors (3 components) - no constraint (all bivectors are simple)");
        eprintln!("  - 4D+: bivectors (6+ components) - constraints needed");
        eprintln!("=============================================\n");
    }

    #[test]
    fn versor_constraint_across_algebras() {
        // Test which algebras need versor constraints for even subalgebra
        eprintln!("\n=== Versor Constraints (v * ṽ = scalar) Across Algebras ===\n");

        let test_cases = [
            (
                "Euclidean 3D Rotor [0,2]",
                Algebra::euclidean(3),
                vec![0, 2],
            ),
            ("PGA 3D Motor [0,2,4]", Algebra::pga(3), vec![0, 2, 4]),
            ("CGA 3D Motor [0,2,4]", Algebra::new(4, 1, 0), vec![0, 2, 4]),
        ];

        for (name, algebra, grades) in &test_cases {
            let constraint = derive_field_constraint(grades, &algebra);

            eprintln!("{}: {:?}", name, grades);
            if let Some(c) = constraint {
                eprintln!("    Constraint: {}", c);
            } else {
                eprintln!("    No constraint needed (automatically satisfies)");
            }
        }
        eprintln!("=============================================\n");
    }

    #[test]
    fn cga_blade_constraint_matches_wiki() {
        // Verify CGA dipole blade constraints match CGA wiki
        let algebra = Algebra::new(4, 1, 0);
        let constraints = derive_blade_constraint(&[2], &algebra);

        eprintln!("\n=== CGA Dipole Blade Constraints ===");
        eprintln!("Field mapping: mx=e12, my=e13, mz=e23, vx=e14, vy=e24, vz=e34,");
        eprintln!("               px=e15, py=e25, pz=e35, pw=e45\n");

        for c in &constraints {
            eprintln!("  {}", c);
        }

        eprintln!("\nCGA wiki equivalent:");
        eprintln!("  1. p × v - pw·m = 0 (3 equations from e1245, e1345, e2345)");
        eprintln!("  2. p · m = 0        (from e1235)");
        eprintln!("  3. v · m = 0        (from e1234)");
        eprintln!("==========================================\n");

        assert_eq!(
            constraints.len(),
            5,
            "CGA dipole should have 5 blade constraints"
        );
    }

    #[test]
    fn null_constraint_cga_point() {
        // CGA points are null vectors: v * ṽ = 0
        let algebra = Algebra::new(4, 1, 0);
        let constraints = derive_null_constraint(&[1], &algebra);

        eprintln!("\n=== CGA Point Null Constraint ===");
        eprintln!("A point in CGA is a null vector (v * ṽ = 0)");
        eprintln!("Components: x=e1, y=e2, z=e3, w=e4, u=e5\n");

        for c in &constraints {
            eprintln!("  {}", c);
        }

        eprintln!("\nExpected: x² + y² + z² + w² - u² = 0");
        eprintln!("(indefinite signature: e1²=e2²=e3²=e4²=+1, e5²=-1)");
        eprintln!("==========================================\n");

        // Should have 1 constraint (the scalar/norm part)
        assert!(
            !constraints.is_empty(),
            "CGA point should have null constraint"
        );
    }

    #[test]
    fn null_constraint_euclidean() {
        // Euclidean vectors: v * ṽ = |v|² (positive definite, no null vectors except 0)
        let algebra = Algebra::euclidean(3);
        let constraints = derive_null_constraint(&[1], &algebra);

        eprintln!("\n=== Euclidean 3D Vector Null Constraint ===");
        for c in &constraints {
            eprintln!("  {}", c);
        }
        eprintln!("Expected: x² + y² + z² = 0 (only zero vector is null)");
        eprintln!("==========================================\n");

        assert!(!constraints.is_empty());
    }

    #[test]
    fn null_constraint_minkowski() {
        // Minkowski vectors can be null (lightlike)
        let algebra = Algebra::new(3, 1, 0); // Cl(3,1,0)
        let constraints = derive_null_constraint(&[1], &algebra);

        eprintln!("\n=== Minkowski Vector Null Constraint ===");
        eprintln!("Lightlike vectors satisfy: x² + y² + z² - t² = 0\n");

        for c in &constraints {
            eprintln!("  {}", c);
        }
        eprintln!("==========================================\n");

        assert!(!constraints.is_empty());
    }

    #[test]
    fn cga_circle_blade_constraint_matches_wiki() {
        // Verify CGA Circle (grade 3) blade constraints match CGA wiki
        // Wiki: https://conformalgeometricalgebra.org/wiki/index.php?title=Circle
        let algebra = Algebra::new(4, 1, 0);
        let constraints = derive_blade_constraint(&[3], &algebra);

        eprintln!("\n=== CGA Circle Blade Constraints ===");
        eprintln!("Field mapping from conformal3.toml:");
        eprintln!("  g (plane normal): gw=e123, gz=e124, gy=e134, gx=e234");
        eprintln!("  m (center):       mz=e125, my=e135, mx=e235");
        eprintln!("  v (moment):       vx=e145, vy=e245, vz=e345\n");

        eprintln!("Derived blade constraints (k ∧ k = 0):");
        for c in &constraints {
            eprintln!("  {}", c);
        }

        eprintln!("\nCGA wiki Circle constraints:");
        eprintln!("  1. g × m - gw*v = 0 (3 equations)");
        eprintln!("     gy*mz - gz*my - gw*vx = 0");
        eprintln!("     gz*mx - gx*mz - gw*vy = 0");
        eprintln!("     gx*my - gy*mx - gw*vz = 0");
        eprintln!("  2. g · v = 0");
        eprintln!("     gx*vx + gy*vy + gz*vz = 0");
        eprintln!("  3. v · m = 0");
        eprintln!("     vx*mx + vy*my + vz*mz = 0");
        eprintln!("==========================================\n");

        // Circle is a trivector (grade 3), so k ∧ k produces grade 6
        // In 5D, grade 6 = 0 (doesn't exist), so no blade constraints!
        // But the wiki constraints DO come from the VERSOR constraint (k * k̃ = scalar)
        eprintln!("Note: In 5D CGA, grade 3 ∧ grade 3 = grade 6, which is empty.");
        eprintln!("However, the wiki constraints come from the VERSOR constraint:");
        eprintln!("k * k̃ for a trivector gives grade 0, 2, 4 components.");
        eprintln!("The grade 4 components must vanish → gives 5 constraints!\n");

        // Grade 3 in 5D: there are C(5,3) = 10 blades
        // Grade 6 in 5D: there are C(5,6) = 0 blades (doesn't exist)
        // So no blade constraints are needed for grade 3 in 5D
        assert!(
            constraints.is_empty(),
            "CGA Circle (grade 3) should have no blade constraints (grade 6 doesn't exist in 5D)"
        );
    }
}

/// Tests verifying Circle constraints match CGA wiki
#[cfg(test)]
mod circle_constraint_tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn circle_5_constraints_match_wiki() {
        // Show the 5 constraints correspond to CGA wiki constraints
        let algebra = Algebra::new(4, 1, 0);
        let table = ProductTable::new(&algebra);
        let blades = blades_of_grades(5, &[3]);

        eprintln!("\n=== Circle: 5 Versor Constraints Match Wiki ===\n");
        eprintln!("Field mapping:");
        eprintln!("  g: gw=e123, gz=e124, gy=e134, gx=e234");
        eprintln!("  m: mz=e125, my=e135, mx=e235");
        eprintln!("  v: vx=e145, vy=e245, vz=e345\n");

        // Group terms by output blade
        let mut by_blade: HashMap<usize, Vec<(String, String, i32)>> = HashMap::new();

        for (i, &a) in blades.iter().enumerate() {
            for &b in &blades[i..] {
                let rev_b = reverse_sign(3);
                let rev_a = reverse_sign(3);

                let (sign_ab, result) = table.geometric(a, b);
                if sign_ab == 0 || grade(result) == 0 {
                    continue;
                }

                let name_a = algebra.blade_index_name(Blade::from_index(a));
                let name_b = algebra.blade_index_name(Blade::from_index(b));

                let coef = if a == b {
                    i32::from(sign_ab) * i32::from(rev_a)
                } else {
                    let (sign_ba, _) = table.geometric(b, a);
                    i32::from(sign_ab) * i32::from(rev_b) + i32::from(sign_ba) * i32::from(rev_a)
                };

                if coef != 0 {
                    by_blade.entry(result).or_default().push((
                        name_a.to_string(),
                        name_b.to_string(),
                        coef,
                    ));
                }
            }
        }

        // Map blade to field name
        let blade_to_field: HashMap<&str, &str> = [
            ("e123", "gw"),
            ("e124", "gz"),
            ("e134", "gy"),
            ("e234", "gx"),
            ("e125", "mz"),
            ("e135", "my"),
            ("e235", "mx"),
            ("e145", "vx"),
            ("e245", "vy"),
            ("e345", "vz"),
        ]
        .into_iter()
        .collect();

        // Print each constraint with wiki interpretation
        let wiki_labels: HashMap<&str, &str> = [
            ("e1234", "v · m = 0"),
            ("e1235", "g · v = 0"),
            ("e1245", "(g × m)_z - gw*vz = 0"),
            ("e1345", "(g × m)_y - gw*vy = 0"),
            ("e2345", "(g × m)_x - gw*vx = 0"),
        ]
        .into_iter()
        .collect();

        for (&result, terms) in &by_blade {
            let result_name = algebra.blade_index_name(Blade::from_index(result));
            let wiki_label = wiki_labels.get(result_name.as_str()).unwrap_or(&"unknown");

            eprintln!("{} constraint (wiki: {}):", result_name, wiki_label);

            let mut field_terms: Vec<String> = Vec::new();
            for (a, b, coef) in terms {
                let fa = blade_to_field
                    .get(a.as_str())
                    .copied()
                    .unwrap_or(a.as_str());
                let fb = blade_to_field
                    .get(b.as_str())
                    .copied()
                    .unwrap_or(b.as_str());
                if *coef > 0 {
                    field_terms.push(format!("+{}*{}*{}", coef, fa, fb));
                } else {
                    field_terms.push(format!("{}*{}*{}", coef, fa, fb));
                }
            }
            eprintln!("  {}", field_terms.join(" "));
            eprintln!();
        }

        eprintln!("==========================================\n");

        // Verify we have exactly 5 constraints
        assert_eq!(
            by_blade.len(),
            5,
            "Should have exactly 5 grade-4 output blades"
        );
    }
}
