//! Geometric constraint checking for blade combinations.
//!
//! According to [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint),
//! geometric entities in GA satisfy two fundamental constraints:
//!
//! 1. **Geometric Product Constraint**: `u ⊡ ũ = u • u`
//!    - The geometric product of an element with its reverse yields a scalar
//!
//! 2. **Geometric Antiproduct Constraint**: `u ⊟ ũ̃ = u ⊙ u`
//!    - The antiproduct of an element with its antireverse yields an antiscalar
//!
//! This module provides functions to check whether a given grade combination
//! satisfies these constraints.
//!
//! # Cancellation of Cross-Terms
//!
//! When checking constraints, we must account for cancellation between symmetric
//! pairs. For blades a and b, the contributions to `u * ũ` are:
//! - From (a, b): `uₐ * uᵦ * rev_sign(grade_b) * sign(a, b)`
//! - From (b, a): `uᵦ * uₐ * rev_sign(grade_a) * sign(b, a)`
//!
//! For non-scalar results to cancel, we need:
//! `rev_sign(grade_b) * sign(a, b) + rev_sign(grade_a) * sign(b, a) = 0`

use super::{Algebra, ProductTable, blades_of_grades, grade, reverse_sign};

/// Checks if a grade combination satisfies the geometric constraint.
///
/// The geometric constraint requires that `u ⊡ ũ` produces only a scalar
/// (grade 0). This is checked by verifying that for all blade pairs in the
/// type, the symmetric contributions either produce scalars or cancel out.
///
/// # Arguments
///
/// * `grades` - The grades present in the type (e.g., `[0, 2]` for a rotor)
/// * `algebra` - The algebra definition
/// * `table` - Precomputed product table
///
/// # Returns
///
/// `true` if the grade combination satisfies the geometric constraint.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::{Algebra, ProductTable, satisfies_geometric_constraint};
///
/// let algebra = Algebra::euclidean(3);
/// let table = ProductTable::new(&algebra);
///
/// // Vectors satisfy the geometric constraint
/// assert!(satisfies_geometric_constraint(&[1], &algebra, &table));
///
/// // Rotors (even subalgebra) satisfy the constraint
/// assert!(satisfies_geometric_constraint(&[0, 2], &algebra, &table));
/// ```
pub fn satisfies_geometric_constraint(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
) -> bool {
    let blades = blades_of_grades(algebra.dim(), grades);

    // Check that u * reverse(u) produces only grade 0
    // We check unordered pairs to account for cancellation between (a,b) and (b,a)
    for (i, &a) in blades.iter().enumerate() {
        for &b in &blades[i..] {
            // b >= a to avoid double-counting
            let ga = grade(a);
            let gb = grade(b);
            let rev_a = reverse_sign(ga);
            let rev_b = reverse_sign(gb);

            let (sign_ab, result_ab) = table.geometric(a, b);

            if a == b {
                // Same blade: a * reverse(a) must be scalar
                if sign_ab != 0 && grade(result_ab) != 0 {
                    return false;
                }
            } else {
                // Different blades: check if cross-terms cancel
                let (sign_ba, _result_ba) = table.geometric(b, a);

                // If result is non-scalar, the symmetric contributions must cancel
                if grade(result_ab) != 0 {
                    // Contribution from (a, b): sign_ab * rev_b
                    // Contribution from (b, a): sign_ba * rev_a
                    // These must sum to zero for cancellation
                    let total = i32::from(sign_ab) * i32::from(rev_b)
                        + i32::from(sign_ba) * i32::from(rev_a);
                    if total != 0 {
                        return false;
                    }
                }
            }
        }
    }

    true
}

/// Checks if a grade combination satisfies the antiproduct constraint.
///
/// The antiproduct constraint requires that `u ⊟ ũ̃` produces only an
/// antiscalar (grade n, the pseudoscalar). This is checked using the
/// duality relationship: `a ⊟ b = dual(dual(a) * dual(b))`.
///
/// Like the geometric constraint, we must account for cancellation between
/// symmetric pairs when checking this constraint.
///
/// # Arguments
///
/// * `grades` - The grades present in the type
/// * `algebra` - The algebra definition
/// * `table` - Precomputed product table
///
/// # Returns
///
/// `true` if the grade combination satisfies the antiproduct constraint.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::{Algebra, ProductTable, satisfies_antiproduct_constraint};
///
/// let algebra = Algebra::euclidean(3);
/// let table = ProductTable::new(&algebra);
///
/// // Vectors satisfy the antiproduct constraint
/// assert!(satisfies_antiproduct_constraint(&[1], &algebra, &table));
///
/// // Rotors satisfy the constraint
/// assert!(satisfies_antiproduct_constraint(&[0, 2], &algebra, &table));
/// ```
pub fn satisfies_antiproduct_constraint(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
) -> bool {
    let dim = algebra.dim();
    let blades = blades_of_grades(dim, grades);
    let pseudoscalar = (1 << dim) - 1; // All bits set = e123...n

    // For antiproduct: a ⊟ b = dual(dual(a) * dual(b))
    // The dual of a k-blade is the (n-k)-blade obtained by XOR with pseudoscalar
    // We need the result to be grade n (antiscalar = pseudoscalar)
    //
    // Check unordered pairs to account for cancellation between (a,b) and (b,a)

    for (i, &a) in blades.iter().enumerate() {
        for &b in &blades[i..] {
            // b >= a to avoid double-counting
            let ga = grade(a);
            let gb = grade(b);

            // Compute dual blade indices (complement with pseudoscalar)
            let dual_a = pseudoscalar ^ a;
            let dual_b = pseudoscalar ^ b;

            // Antireverse signs: reverse sign of the dual grades
            let antirev_a = reverse_sign(dim - ga);
            let antirev_b = reverse_sign(dim - gb);

            // Product of duals
            let (sign_ab, dual_result_ab) = table.geometric(dual_a, dual_b);
            let result_ab = pseudoscalar ^ dual_result_ab;

            if a == b {
                // Same blade: a ⊟ antireverse(a) must be antiscalar
                if sign_ab != 0 && grade(result_ab) != dim {
                    return false;
                }
            } else {
                // Different blades: check if cross-terms cancel
                let (sign_ba, _dual_result_ba) = table.geometric(dual_b, dual_a);

                // If result is not antiscalar, the symmetric contributions must cancel
                if grade(result_ab) != dim {
                    // Contribution from (a, b): sign_ab * antirev_b
                    // Contribution from (b, a): sign_ba * antirev_a
                    // These must sum to zero for cancellation
                    let total = i32::from(sign_ab) * i32::from(antirev_b)
                        + i32::from(sign_ba) * i32::from(antirev_a);
                    if total != 0 {
                        return false;
                    }
                }
            }
        }
    }

    true
}

/// Checks if a grade combination satisfies both geometric constraints.
///
/// A valid geometric entity must satisfy:
/// 1. `u ⊡ ũ = u • u` (produces only scalar)
/// 2. `u ⊟ ũ̃ = u ⊙ u` (produces only antiscalar)
///
/// # Arguments
///
/// * `grades` - The grades present in the type
/// * `algebra` - The algebra definition
/// * `table` - Precomputed product table
///
/// # Returns
///
/// `true` if the grade combination satisfies both constraints.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::{Algebra, ProductTable, satisfies_all_constraints};
///
/// let algebra = Algebra::euclidean(3);
/// let table = ProductTable::new(&algebra);
///
/// // Single grades satisfy both constraints
/// assert!(satisfies_all_constraints(&[0], &algebra, &table));       // Scalar
/// assert!(satisfies_all_constraints(&[1], &algebra, &table));       // Vector
/// assert!(satisfies_all_constraints(&[2], &algebra, &table));       // Bivector
/// assert!(satisfies_all_constraints(&[3], &algebra, &table));       // Trivector
///
/// // Even and odd subalgebras satisfy constraints
/// assert!(satisfies_all_constraints(&[0, 2], &algebra, &table));    // Rotor
/// assert!(satisfies_all_constraints(&[1, 3], &algebra, &table));    // Odd
///
/// // Full multivector does NOT satisfy: scalar*vector produces non-canceling terms
/// assert!(!satisfies_all_constraints(&[0, 1, 2, 3], &algebra, &table));
/// ```
pub fn satisfies_all_constraints(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
) -> bool {
    satisfies_geometric_constraint(grades, algebra, table)
        && satisfies_antiproduct_constraint(grades, algebra, table)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euclidean_3d_single_grades() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // All single grades satisfy both constraints
        assert!(satisfies_all_constraints(&[0], &algebra, &table)); // Scalar
        assert!(satisfies_all_constraints(&[1], &algebra, &table)); // Vector
        assert!(satisfies_all_constraints(&[2], &algebra, &table)); // Bivector
        assert!(satisfies_all_constraints(&[3], &algebra, &table)); // Trivector
    }

    #[test]
    fn euclidean_3d_even_odd() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Even subalgebra (rotor) satisfies both constraints
        assert!(satisfies_all_constraints(&[0, 2], &algebra, &table));

        // Odd subalgebra [1, 3] satisfies both constraints
        assert!(satisfies_all_constraints(&[1, 3], &algebra, &table));

        // Full multivector does NOT satisfy constraints:
        // scalar * vector produces non-canceling vector terms
        assert!(!satisfies_all_constraints(&[0, 1, 2, 3], &algebra, &table));
    }

    #[test]
    fn euclidean_2d() {
        let algebra = Algebra::euclidean(2);
        let table = ProductTable::new(&algebra);

        // Single grades satisfy constraints
        assert!(satisfies_all_constraints(&[0], &algebra, &table)); // Scalar
        assert!(satisfies_all_constraints(&[1], &algebra, &table)); // Vector
        assert!(satisfies_all_constraints(&[2], &algebra, &table)); // Bivector

        // Even subalgebra [0, 2] satisfies constraints
        assert!(satisfies_all_constraints(&[0, 2], &algebra, &table));

        // Full multivector [0, 1, 2] does NOT satisfy:
        // mixing scalar with vector creates non-canceling terms
        assert!(!satisfies_all_constraints(&[0, 1, 2], &algebra, &table));
    }

    #[test]
    fn pga_3d_constraints() {
        // PGA for 3D is a 4D algebra with one degenerate basis (e0² = 0).
        // In 4D, the constraint situation is more complex:
        // - Disjoint k-blades commute when (-1)^(k²) = 1
        // - Cross-terms between grades can produce non-canceling results

        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Scalar and pseudoscalar always satisfy (only one blade each)
        assert!(satisfies_all_constraints(&[0], &algebra, &table)); // Scalar
        assert!(satisfies_all_constraints(&[4], &algebra, &table)); // Pseudoscalar

        // In 4D, single vectors satisfy the geometric constraint
        // (same-grade blades of grade 1 anticommute)
        assert!(satisfies_all_constraints(&[1], &algebra, &table)); // Vectors

        // Single trivectors in 4D: similar to vectors
        assert!(satisfies_all_constraints(&[3], &algebra, &table)); // Trivectors

        // Bivectors in 4D: disjoint bivectors commute (k²=4 is even)
        // This means e01 and e23 produce non-canceling e0123 terms
        assert!(!satisfies_all_constraints(&[2], &algebra, &table)); // Bivectors fail

        // In 4D, combining grade 1 and 3 creates cross-terms that don't cancel
        // because of the reverse sign interaction
        assert!(!satisfies_all_constraints(&[1, 3], &algebra, &table)); // Flector fails

        // Full even subalgebra [0, 2, 4] also fails due to bivector issues
        assert!(!satisfies_all_constraints(&[0, 2, 4], &algebra, &table)); // Motor fails

        // Note: In actual PGA usage, motors and flectors are defined by specific
        // properties beyond just grade constraints. The geometric constraint
        // is one characterization that helps identify certain structures.
    }

    #[test]
    fn geometric_constraint_only() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Check geometric constraint individually
        assert!(satisfies_geometric_constraint(&[1], &algebra, &table));
        assert!(satisfies_geometric_constraint(&[0, 2], &algebra, &table));
    }

    #[test]
    fn antiproduct_constraint_only() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Check antiproduct constraint individually
        assert!(satisfies_antiproduct_constraint(&[1], &algebra, &table));
        assert!(satisfies_antiproduct_constraint(&[0, 2], &algebra, &table));
    }
}
