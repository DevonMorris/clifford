//! Entity discovery for geometric algebra code generation.
//!
//! This module automatically discovers valid geometric entities from an algebra
//! signature by analyzing geometric constraints.
//!
//! # Overview
//!
//! Given an algebra signature (p, q, r), entity discovery:
//! 1. Enumerates all possible grade combinations
//! 2. Analyzes which grade pairs conflict (create non-canceling terms)
//! 3. Determines what additional constraints are needed for satisfaction
//! 4. Names known entities using heuristics
//! 5. Can generate a TOML template for further customization
//!
//! # Constraint Analysis
//!
//! Not all grade combinations automatically satisfy the geometric constraint
//! `u * ũ = scalar`. For combinations that don't, we analyze:
//! - Which pairs of grades create conflicting (non-canceling) terms
//! - What constraints would be needed to satisfy the constraint
//!
//! For example, [0, 1] (scalar + vector) has a conflict because:
//! - scalar * vector = vector, and these don't cancel
//! - Constraint needed: either scalar = 0 OR vector = 0
//!
//! # Example
//!
//! ```
//! use clifford_codegen::discovery::{discover_entities, analyze_constraints, DiscoveredEntity};
//! use clifford_codegen::algebra::Algebra;
//!
//! let algebra = Algebra::euclidean(3);
//!
//! // Analyze what constraints would be needed for the full multivector
//! let conflicts = analyze_constraints(&[0, 1, 2, 3], &algebra);
//! // Returns: conflicting grade pairs like (0, 1), (0, 3), (1, 2), (2, 3)
//!
//! // Discover entities that automatically satisfy constraints
//! let entities = discover_entities(&algebra);
//! assert!(entities.iter().any(|e| e.name == "Entity_0_2")); // [0, 2] - no conflicts
//! ```

mod constraints;
mod entity;
mod naming;
pub mod products;
mod template;

pub use constraints::{
    can_satisfy_constraints, derive_antiproduct_constraint, derive_blade_constraint,
    derive_field_constraint, derive_null_constraint,
};
pub use entity::DiscoveredEntity;
pub use naming::suggest_name;
pub use products::{
    BladeProductResult, EntityBladeSet, ProductResult, ProductTable2D, ProductType,
    infer_all_products, infer_all_products_blades, infer_output_blades, infer_output_grades,
    infer_product, infer_product_blades,
};
pub use template::generate_toml_template;

use crate::algebra::{
    Algebra, ProductTable, blades_of_grades, grade, reverse_sign, satisfies_all_constraints,
};

/// Represents a conflict between two grades in the geometric constraint.
///
/// When grades `grade_a` and `grade_b` have a conflict, it means that
/// elements with non-zero components in both grades will produce
/// non-canceling cross-terms in `u * ũ`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GradeConflict {
    /// First grade in the conflict.
    pub grade_a: usize,
    /// Second grade in the conflict.
    pub grade_b: usize,
    /// Description of why these grades conflict.
    pub reason: String,
}

/// Analyzes which grade pairs conflict in the geometric constraint.
///
/// For a grade combination, this identifies pairs of grades where the
/// cross-terms in `u * ũ` don't cancel, producing non-scalar results.
///
/// # Returns
///
/// A list of conflicting grade pairs. If empty, the grade combination
/// automatically satisfies the geometric constraint.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::analyze_constraints;
/// use clifford_codegen::algebra::Algebra;
///
/// let algebra = Algebra::euclidean(3);
///
/// // Rotor [0, 2] has no conflicts - cross-terms cancel
/// let conflicts = analyze_constraints(&[0, 2], &algebra);
/// assert!(conflicts.is_empty());
///
/// // Full multivector [0, 1, 2, 3] has conflicts
/// let conflicts = analyze_constraints(&[0, 1, 2, 3], &algebra);
/// assert!(!conflicts.is_empty());
/// ```
pub fn analyze_constraints(grades: &[usize], algebra: &Algebra) -> Vec<GradeConflict> {
    let table = ProductTable::new(algebra);
    let mut conflicts = Vec::new();

    // Check each pair of grades
    for (i, &ga) in grades.iter().enumerate() {
        for &gb in &grades[i..] {
            if ga == gb {
                // Same grade - check if self-products are scalar
                if !same_grade_satisfies(ga, algebra, &table) {
                    conflicts.push(GradeConflict {
                        grade_a: ga,
                        grade_b: gb,
                        reason: format!("Grade {} blades produce non-scalar self-products", ga),
                    });
                }
            } else {
                // Different grades - check if cross-terms cancel
                if !cross_grades_cancel(ga, gb, algebra, &table) {
                    conflicts.push(GradeConflict {
                        grade_a: ga,
                        grade_b: gb,
                        reason: format!(
                            "Grades {} and {} produce non-canceling cross-terms",
                            ga, gb
                        ),
                    });
                }
            }
        }
    }

    conflicts
}

/// Checks if blades of the same grade produce only scalar self-products.
fn same_grade_satisfies(g: usize, algebra: &Algebra, table: &ProductTable) -> bool {
    let blades = blades_of_grades(algebra.dim(), &[g]);
    let rev_g = reverse_sign(g);

    for (i, &a) in blades.iter().enumerate() {
        // Check a * reverse(a)
        let (sign_aa, result_aa) = table.geometric(a, a);
        if sign_aa != 0 && grade(result_aa) != 0 {
            return false;
        }

        // Check cross-products between same-grade blades
        for &b in &blades[i + 1..] {
            let (sign_ab, result_ab) = table.geometric(a, b);
            let (sign_ba, _) = table.geometric(b, a);

            if grade(result_ab) != 0 {
                let total =
                    i32::from(sign_ab) * i32::from(rev_g) + i32::from(sign_ba) * i32::from(rev_g);
                if total != 0 {
                    return false;
                }
            }
        }
    }

    true
}

/// Checks if cross-terms between two different grades cancel.
fn cross_grades_cancel(ga: usize, gb: usize, algebra: &Algebra, table: &ProductTable) -> bool {
    let blades_a = blades_of_grades(algebra.dim(), &[ga]);
    let blades_b = blades_of_grades(algebra.dim(), &[gb]);
    let rev_a = reverse_sign(ga);
    let rev_b = reverse_sign(gb);

    for &a in &blades_a {
        for &b in &blades_b {
            let (sign_ab, result_ab) = table.geometric(a, b);
            let (sign_ba, _) = table.geometric(b, a);

            // If result is non-scalar, check if contributions cancel
            if grade(result_ab) != 0 {
                let total =
                    i32::from(sign_ab) * i32::from(rev_b) + i32::from(sign_ba) * i32::from(rev_a);
                if total != 0 {
                    return false;
                }
            }
        }
    }

    true
}

/// Suggests constraints needed to satisfy the geometric constraint.
///
/// Given a list of grade conflicts, suggests what constraints would
/// make the grade combination satisfy `u * ũ = scalar`.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::{analyze_constraints, suggest_required_constraints};
/// use clifford_codegen::algebra::Algebra;
///
/// let algebra = Algebra::euclidean(3);
/// let conflicts = analyze_constraints(&[0, 1, 2, 3], &algebra);
/// let suggestions = suggest_required_constraints(&conflicts);
///
/// // Will suggest something like: "At most one of {0,2} or {1,3} can be non-zero"
/// ```
pub fn suggest_required_constraints(conflicts: &[GradeConflict]) -> Vec<String> {
    if conflicts.is_empty() {
        return vec!["No additional constraints needed".to_string()];
    }

    let mut suggestions = Vec::new();

    // Group conflicts by grade
    let mut conflicting_grades: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for conflict in conflicts {
        conflicting_grades.insert(conflict.grade_a);
        conflicting_grades.insert(conflict.grade_b);
    }

    // Find maximal non-conflicting subsets
    // This is a simplified heuristic - could be made more sophisticated
    if conflicting_grades.len() > 2 {
        suggestions.push(format!(
            "Grades {:?} have mutual conflicts. Consider using subsets that satisfy constraints.",
            conflicting_grades.iter().collect::<Vec<_>>()
        ));
    } else {
        for conflict in conflicts {
            suggestions.push(format!(
                "Either grade {} or grade {} must be zero (constraint: {})",
                conflict.grade_a, conflict.grade_b, conflict.reason
            ));
        }
    }

    suggestions
}

/// Enumerates all non-empty grade combinations for an n-dimensional algebra.
///
/// For an n-dimensional algebra with grades 0 through n, this generates
/// all 2^(n+1) - 1 non-empty subsets of grades.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::enumerate_grade_combinations;
///
/// // 2D algebra has grades 0, 1, 2
/// let combinations = enumerate_grade_combinations(2);
/// assert_eq!(combinations.len(), 7); // 2^3 - 1 = 7 non-empty subsets
/// ```
pub fn enumerate_grade_combinations(dim: usize) -> Vec<Vec<usize>> {
    let num_grades = dim + 1; // grades 0 through dim

    // 2^(num_grades) - 1 non-empty subsets
    (1..(1 << num_grades))
        .map(|mask| (0..num_grades).filter(|&g| (mask >> g) & 1 == 1).collect())
        .collect()
}

/// Discovers all valid grade combinations that satisfy geometric constraints.
///
/// Returns grade combinations that satisfy both:
/// 1. Geometric product constraint: `u * ũ` produces only scalar
/// 2. Antiproduct constraint: `u ⊟ ũ̃` produces only antiscalar
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::discover_valid_combinations;
/// use clifford_codegen::algebra::Algebra;
///
/// let algebra = Algebra::euclidean(3);
/// let valid = discover_valid_combinations(&algebra);
///
/// // All single grades satisfy constraints
/// assert!(valid.contains(&vec![0])); // Scalar
/// assert!(valid.contains(&vec![1])); // Vector
/// assert!(valid.contains(&vec![2])); // Bivector
/// assert!(valid.contains(&vec![3])); // Trivector
///
/// // Even and odd subalgebras satisfy
/// assert!(valid.contains(&vec![0, 2])); // Rotor
/// assert!(valid.contains(&vec![1, 3])); // Odd
/// ```
pub fn discover_valid_combinations(algebra: &Algebra) -> Vec<Vec<usize>> {
    let table = ProductTable::new(algebra);

    enumerate_grade_combinations(algebra.dim())
        .into_iter()
        .filter(|grades| satisfies_all_constraints(grades, algebra, &table))
        .collect()
}

/// Discovers the minimal closed set of geometric entities.
///
/// This returns the standard set of types that form a closed algebra:
/// - All single-grade types (scalar, vector, bivector, etc.)
/// - Even subalgebra (grades 0, 2, 4, ...) - motors/rotors
/// - Odd subalgebra (grades 1, 3, 5, ...) - flectors
///
/// This minimal set is sufficient for most geometric algebra applications
/// and matches the standard types used in PGA, CGA, etc.
///
/// For grade combinations that need field constraints (like PGA bivectors),
/// the constraint expression is included in the entity.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::discover_entities;
/// use clifford_codegen::algebra::Algebra;
///
/// let algebra = Algebra::pga(3);
/// let entities = discover_entities(&algebra);
///
/// // Returns 7 entities for PGA:
/// // [0], [1], [2], [3], [4], [0,2,4], [1,3]
/// assert_eq!(entities.len(), 7);
/// ```
pub fn discover_entities(algebra: &Algebra) -> Vec<DiscoveredEntity> {
    let dim = algebra.dim();

    // Build the minimal set of grade combinations
    let mut grade_sets: Vec<Vec<usize>> = Vec::new();

    // Single grades (scalar, vector, bivector, ..., pseudoscalar)
    for g in 0..=dim {
        grade_sets.push(vec![g]);
    }

    // Even subalgebra (grades 0, 2, 4, ...) - only if more than one grade
    let even: Vec<usize> = (0..=dim).step_by(2).collect();
    if even.len() > 1 {
        grade_sets.push(even);
    }

    // Odd subalgebra (grades 1, 3, 5, ...) - only if more than one grade
    let odd: Vec<usize> = (1..=dim).step_by(2).collect();
    if odd.len() > 1 {
        grade_sets.push(odd);
    }

    // Convert to entities with constraints
    grade_sets
        .into_iter()
        .map(|grades| {
            let geometric_constraint = derive_field_constraint(&grades, algebra);
            let antiproduct_constraint = derive_antiproduct_constraint(&grades, algebra);
            let name = suggest_name(&grades, algebra.dim());

            DiscoveredEntity {
                name,
                grades,
                geometric_constraint,
                antiproduct_constraint,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumerate_2d() {
        let combinations = enumerate_grade_combinations(2);
        // 2^3 - 1 = 7 non-empty subsets
        assert_eq!(combinations.len(), 7);

        // Check all single grades present
        assert!(combinations.contains(&vec![0]));
        assert!(combinations.contains(&vec![1]));
        assert!(combinations.contains(&vec![2]));
    }

    #[test]
    fn enumerate_3d() {
        let combinations = enumerate_grade_combinations(3);
        // 2^4 - 1 = 15 non-empty subsets
        assert_eq!(combinations.len(), 15);
    }

    #[test]
    fn discover_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let valid = discover_valid_combinations(&algebra);

        // Single grades
        assert!(valid.contains(&vec![0])); // Scalar
        assert!(valid.contains(&vec![1])); // Vector
        assert!(valid.contains(&vec![2])); // Bivector
        assert!(valid.contains(&vec![3])); // Trivector

        // Even and odd subalgebras
        assert!(valid.contains(&vec![0, 2])); // Rotor
        assert!(valid.contains(&vec![1, 3])); // Odd

        // Full multivector does NOT satisfy constraints
        assert!(!valid.contains(&vec![0, 1, 2, 3]));
    }

    #[test]
    fn discover_euclidean_2d() {
        let algebra = Algebra::euclidean(2);
        let valid = discover_valid_combinations(&algebra);

        // Single grades
        assert!(valid.contains(&vec![0])); // Scalar
        assert!(valid.contains(&vec![1])); // Vector
        assert!(valid.contains(&vec![2])); // Bivector

        // Even subalgebra
        assert!(valid.contains(&vec![0, 2])); // Rotor
    }

    #[test]
    fn discover_entities_minimal_set() {
        // Euclidean 3D: 4 single grades + even [0,2] + odd [1,3] = 6 entities
        let algebra = Algebra::euclidean(3);
        let entities = discover_entities(&algebra);
        assert_eq!(entities.len(), 6);

        // Check all expected entities are present
        assert!(entities.iter().any(|e| e.grades == vec![0]));
        assert!(entities.iter().any(|e| e.grades == vec![1]));
        assert!(entities.iter().any(|e| e.grades == vec![2]));
        assert!(entities.iter().any(|e| e.grades == vec![3]));
        assert!(entities.iter().any(|e| e.grades == vec![0, 2]));
        assert!(entities.iter().any(|e| e.grades == vec![1, 3]));
    }

    #[test]
    fn discover_pga_entities_minimal_set() {
        // PGA 3D (4D algebra): 5 single grades + even [0,2,4] + odd [1,3] = 7 entities
        let algebra = Algebra::pga(3);
        let entities = discover_entities(&algebra);
        assert_eq!(entities.len(), 7);

        // Check bivector has constraint(s)
        let bivector = entities.iter().find(|e| e.grades == vec![2]).unwrap();
        assert!(
            bivector.has_constraints(),
            "Bivector should have constraints"
        );

        // Check motor has constraint(s) - should have 2 constraints for 6 DOF
        let motor = entities.iter().find(|e| e.grades == vec![0, 2, 4]).unwrap();
        assert!(motor.has_constraints(), "Motor should have constraints");
        // Motor should have 2 constraints (8 coeffs - 2 constraints = 6 DOF)
        assert_eq!(
            motor.constraint_count(),
            2,
            "Motor should have exactly 2 constraints for 6 DOF"
        );

        // Check flector has constraint(s)
        let flector = entities.iter().find(|e| e.grades == vec![1, 3]).unwrap();
        assert!(flector.has_constraints(), "Flector should have constraints");
    }

    #[test]
    fn analyze_rotor_no_conflicts() {
        let algebra = Algebra::euclidean(3);
        let conflicts = analyze_constraints(&[0, 2], &algebra);
        assert!(
            conflicts.is_empty(),
            "Rotor [0, 2] should have no conflicts"
        );
    }

    #[test]
    fn analyze_full_has_conflicts() {
        let algebra = Algebra::euclidean(3);
        let conflicts = analyze_constraints(&[0, 1, 2, 3], &algebra);
        assert!(
            !conflicts.is_empty(),
            "Full multivector should have conflicts"
        );

        // Should find conflicts between grades that don't cancel
        // e.g., (0, 1), (0, 3), (1, 2), (2, 3)
        assert!(conflicts.iter().any(|c| c.grade_a == 0 && c.grade_b == 1));
    }

    #[test]
    fn analyze_scalar_vector_conflict() {
        let algebra = Algebra::euclidean(3);
        let conflicts = analyze_constraints(&[0, 1], &algebra);
        assert_eq!(conflicts.len(), 1);
        assert_eq!(conflicts[0].grade_a, 0);
        assert_eq!(conflicts[0].grade_b, 1);
    }

    #[test]
    fn discover_pga_3d() {
        let algebra = Algebra::pga(3); // 4D
        let valid = discover_valid_combinations(&algebra);

        // Scalar and pseudoscalar always satisfy
        assert!(valid.contains(&vec![0]));
        assert!(valid.contains(&vec![4]));

        // Vectors and trivectors satisfy
        assert!(valid.contains(&vec![1]));
        assert!(valid.contains(&vec![3]));

        // Bivectors do NOT satisfy in 4D (disjoint bivectors commute)
        assert!(!valid.contains(&vec![2]));
    }
}
