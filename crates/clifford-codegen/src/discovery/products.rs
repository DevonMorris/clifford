//! Product output inference for discovered entities.
//!
//! This module infers the output type for products between discovered entities.
//! Given two entities and a product type, it computes which grades appear in
//! the output and matches them to known entities.
//!
//! When the raw output doesn't match a known entity exactly, constraint
//! propagation finds the largest known entity that is a subset of the output.

use crate::algebra::{
    Algebra, ProductTable, blades_of_grades, geometric_grades, grade, inner_grade,
    left_contraction_grade, outer_grade,
};
use std::collections::BTreeSet;

/// Represents a product type for inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProductType {
    /// Geometric product: `a * b`
    Geometric,
    /// Outer (wedge) product: `a ∧ b`
    Outer,
    /// Inner product: `a · b`
    Inner,
    /// Left contraction: `a ⌋ b`
    LeftContraction,
}

impl ProductType {
    /// Returns all standard product types.
    pub fn all() -> &'static [ProductType] {
        &[
            ProductType::Geometric,
            ProductType::Outer,
            ProductType::Inner,
            ProductType::LeftContraction,
        ]
    }

    /// Returns the name for TOML output.
    pub fn toml_name(&self) -> &'static str {
        match self {
            ProductType::Geometric => "geometric",
            ProductType::Outer => "outer",
            ProductType::Inner => "inner",
            ProductType::LeftContraction => "left_contraction",
        }
    }
}

/// Result of product inference between two entities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductResult {
    /// Grades present in the output after applying constraints.
    /// This is the constrained/simplified output that satisfies geometric constraints.
    pub output_grades: Vec<usize>,
    /// Name of matching entity, if one exists.
    pub matching_entity: Option<String>,
    /// Whether the output is always zero (after constraint application).
    pub is_zero: bool,
}

/// Infers the output grades for a product between two grade sets.
///
/// # Arguments
///
/// * `lhs_grades` - Grades in the left operand
/// * `rhs_grades` - Grades in the right operand
/// * `product_type` - The type of product
/// * `algebra` - The algebra (for dimension and actual computation)
///
/// # Returns
///
/// The set of grades that can appear in the output.
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::products::{ProductType, infer_output_grades};
/// use clifford_codegen::algebra::Algebra;
///
/// let algebra = Algebra::euclidean(3);
///
/// // Vector * Vector produces scalar + bivector (grades 0, 2)
/// let output = infer_output_grades(&[1], &[1], ProductType::Geometric, &algebra);
/// assert_eq!(output, vec![0, 2]);
///
/// // Vector ∧ Vector produces bivector (grade 2)
/// let output = infer_output_grades(&[1], &[1], ProductType::Outer, &algebra);
/// assert_eq!(output, vec![2]);
/// ```
pub fn infer_output_grades(
    lhs_grades: &[usize],
    rhs_grades: &[usize],
    product_type: ProductType,
    algebra: &Algebra,
) -> Vec<usize> {
    let dim = algebra.dim();
    let mut output_set = BTreeSet::new();

    for &ga in lhs_grades {
        for &gb in rhs_grades {
            match product_type {
                ProductType::Geometric => {
                    for g in geometric_grades(ga, gb, dim) {
                        output_set.insert(g);
                    }
                }
                ProductType::Outer => {
                    if let Some(g) = outer_grade(ga, gb, dim) {
                        output_set.insert(g);
                    }
                }
                ProductType::Inner => {
                    output_set.insert(inner_grade(ga, gb));
                }
                ProductType::LeftContraction => {
                    if let Some(g) = left_contraction_grade(ga, gb) {
                        output_set.insert(g);
                    }
                }
            }
        }
    }

    output_set.into_iter().collect()
}

/// Infers the output grades using actual product computation.
///
/// This is more precise than `infer_output_grades` because it considers
/// the metric signature. Some products may be zero due to the metric
/// even though the grades would theoretically be present.
///
/// # Arguments
///
/// * `lhs_grades` - Grades in the left operand
/// * `rhs_grades` - Grades in the right operand
/// * `product_type` - The type of product
/// * `algebra` - The algebra (for dimension and product computation)
/// * `table` - Precomputed product table
///
/// # Returns
///
/// The set of grades that actually appear in non-zero output.
pub fn infer_output_grades_precise(
    lhs_grades: &[usize],
    rhs_grades: &[usize],
    product_type: ProductType,
    algebra: &Algebra,
    table: &ProductTable,
) -> Vec<usize> {
    let lhs_blades = blades_of_grades(algebra.dim(), lhs_grades);
    let rhs_blades = blades_of_grades(algebra.dim(), rhs_grades);
    let mut output_set = BTreeSet::new();

    for &a in &lhs_blades {
        for &b in &rhs_blades {
            let (sign, result) = table.geometric(a, b);
            if sign == 0 {
                continue;
            }

            let result_grade = grade(result);
            let ga = grade(a);
            let gb = grade(b);

            // Check if this product should be included based on product type
            let include = match product_type {
                ProductType::Geometric => true,
                ProductType::Outer => {
                    // Outer product: only grade ga + gb terms
                    result_grade == ga + gb
                }
                ProductType::Inner => {
                    // Inner product: only grade |ga - gb| terms
                    result_grade == ga.abs_diff(gb)
                }
                ProductType::LeftContraction => {
                    // Left contraction: only grade gb - ga terms (when ga <= gb)
                    ga <= gb && result_grade == gb - ga
                }
            };

            if include {
                output_set.insert(result_grade);
            }
        }
    }

    output_set.into_iter().collect()
}

/// Infers the product output and matches to known entities.
///
/// This function:
/// 1. Computes raw output grades from the product
/// 2. Matches directly to known entities (exact match)
/// 3. If no exact match, finds the best matching known entity that is a subset
///    of the raw output (constraint propagation from input types)
///
/// # Arguments
///
/// * `lhs_grades` - Grades in the left operand
/// * `rhs_grades` - Grades in the right operand
/// * `product_type` - The type of product
/// * `known_entities` - Map from grade set to entity name
/// * `algebra` - The algebra
/// * `table` - Precomputed product table
///
/// # Returns
///
/// The inferred product result.
pub fn infer_product(
    lhs_grades: &[usize],
    rhs_grades: &[usize],
    product_type: ProductType,
    known_entities: &[(Vec<usize>, String)],
    algebra: &Algebra,
    table: &ProductTable,
) -> ProductResult {
    let raw_output = infer_output_grades_precise(lhs_grades, rhs_grades, product_type, algebra, table);

    if raw_output.is_empty() {
        return ProductResult {
            output_grades: vec![],
            matching_entity: None,
            is_zero: true,
        };
    }

    let raw_set: BTreeSet<usize> = raw_output.iter().copied().collect();

    // First try exact match
    if let Some((_, name)) = known_entities.iter().find(|(grades, _)| grades == &raw_output) {
        return ProductResult {
            output_grades: raw_output,
            matching_entity: Some(name.clone()),
            is_zero: false,
        };
    }

    // No exact match - find the best matching known entity that is a subset of raw output
    // This implements constraint propagation: the input types' constraints limit the output
    // to grades that form a valid entity
    let best_match = find_best_matching_subset(&raw_set, known_entities);

    if let Some((grades, name)) = best_match {
        ProductResult {
            output_grades: grades,
            matching_entity: Some(name),
            is_zero: false,
        }
    } else {
        // No matching entity - return raw output without entity match
        ProductResult {
            output_grades: raw_output,
            matching_entity: None,
            is_zero: false,
        }
    }
}

/// Finds the largest known entity that is a subset of the given grades.
///
/// This implements constraint propagation: when the raw product output
/// doesn't match any known entity exactly, we find the largest known
/// entity whose grades are all present in the output.
fn find_best_matching_subset(
    raw_grades: &BTreeSet<usize>,
    known_entities: &[(Vec<usize>, String)],
) -> Option<(Vec<usize>, String)> {
    // Filter to entities that are subsets of raw_grades
    let matching: Vec<_> = known_entities
        .iter()
        .filter(|(grades, _)| {
            grades.iter().all(|g| raw_grades.contains(g))
        })
        .collect();

    // Find the one with the most grades (largest subset)
    matching
        .into_iter()
        .max_by_key(|(grades, _)| grades.len())
        .map(|(grades, name)| (grades.clone(), name.clone()))
}

/// Represents a complete product table between entities.
#[derive(Debug, Clone)]
pub struct ProductTable2D {
    /// Product type.
    pub product_type: ProductType,
    /// Entries as (lhs_name, rhs_name, result).
    pub entries: Vec<(String, String, ProductResult)>,
}

/// Infers all products between a set of entities.
///
/// # Arguments
///
/// * `entities` - List of (name, grades) pairs
/// * `product_type` - The type of product
/// * `algebra` - The algebra
///
/// # Returns
///
/// A complete product table.
pub fn infer_all_products(
    entities: &[(String, Vec<usize>)],
    product_type: ProductType,
    algebra: &Algebra,
) -> ProductTable2D {
    let table = ProductTable::new(algebra);
    let known_entities: Vec<(Vec<usize>, String)> = entities
        .iter()
        .map(|(name, grades)| (grades.clone(), name.clone()))
        .collect();

    let mut entries = Vec::new();

    for (lhs_name, lhs_grades) in entities {
        for (rhs_name, rhs_grades) in entities {
            let result = infer_product(
                lhs_grades,
                rhs_grades,
                product_type,
                &known_entities,
                algebra,
                &table,
            );
            entries.push((lhs_name.clone(), rhs_name.clone(), result));
        }
    }

    ProductTable2D {
        product_type,
        entries,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometric_vector_vector() {
        let algebra = Algebra::euclidean(3);

        // Vector * Vector = Scalar + Bivector
        let output = infer_output_grades(&[1], &[1], ProductType::Geometric, &algebra);
        assert_eq!(output, vec![0, 2]);
    }

    #[test]
    fn geometric_vector_vector_precise() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Vector * Vector = Scalar + Bivector (same result with precise computation)
        let output = infer_output_grades_precise(&[1], &[1], ProductType::Geometric, &algebra, &table);
        assert_eq!(output, vec![0, 2]);
    }

    #[test]
    fn outer_vector_vector() {
        let algebra = Algebra::euclidean(3);

        // Vector ∧ Vector = Bivector
        let output = infer_output_grades(&[1], &[1], ProductType::Outer, &algebra);
        assert_eq!(output, vec![2]);
    }

    #[test]
    fn outer_vector_bivector() {
        let algebra = Algebra::euclidean(3);

        // Vector ∧ Bivector = Trivector
        let output = infer_output_grades(&[1], &[2], ProductType::Outer, &algebra);
        assert_eq!(output, vec![3]);
    }

    #[test]
    fn outer_bivector_bivector() {
        let algebra = Algebra::euclidean(3);

        // Bivector ∧ Bivector = 0 in 3D (grade 4 > dim 3)
        let output = infer_output_grades(&[2], &[2], ProductType::Outer, &algebra);
        assert!(output.is_empty());
    }

    #[test]
    fn inner_vector_vector() {
        let algebra = Algebra::euclidean(3);

        // Vector · Vector = Scalar
        let output = infer_output_grades(&[1], &[1], ProductType::Inner, &algebra);
        assert_eq!(output, vec![0]);
    }

    #[test]
    fn inner_bivector_vector() {
        let algebra = Algebra::euclidean(3);

        // Bivector · Vector = Vector
        let output = infer_output_grades(&[2], &[1], ProductType::Inner, &algebra);
        assert_eq!(output, vec![1]);
    }

    #[test]
    fn left_contraction_vector_bivector() {
        let algebra = Algebra::euclidean(3);

        // Vector ⌋ Bivector = Vector
        let output = infer_output_grades(&[1], &[2], ProductType::LeftContraction, &algebra);
        assert_eq!(output, vec![1]);
    }

    #[test]
    fn left_contraction_bivector_vector() {
        let algebra = Algebra::euclidean(3);

        // Bivector ⌋ Vector = 0 (grade 2 > grade 1)
        let output = infer_output_grades(&[2], &[1], ProductType::LeftContraction, &algebra);
        assert!(output.is_empty());
    }

    #[test]
    fn rotor_geometric_products() {
        let algebra = Algebra::euclidean(3);

        // Even * Even = Even (grades 0, 2)
        let output = infer_output_grades(&[0, 2], &[0, 2], ProductType::Geometric, &algebra);
        assert_eq!(output, vec![0, 2]);

        // Even * Vector = Odd (grades 1, 3)
        let output = infer_output_grades(&[0, 2], &[1], ProductType::Geometric, &algebra);
        assert_eq!(output, vec![1, 3]);
    }

    #[test]
    fn infer_product_with_matching() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        let entities = vec![
            (vec![0], "Entity_0".to_string()),
            (vec![1], "Entity_1".to_string()),
            (vec![2], "Entity_2".to_string()),
            (vec![0, 2], "Entity_0_2".to_string()),
        ];

        // Vector * Vector should match Entity_0_2 (grades 0, 2)
        let result = infer_product(&[1], &[1], ProductType::Geometric, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Entity_0_2".to_string()));
        assert!(!result.is_zero);

        // Vector ∧ Vector should match Entity_2 (grade 2)
        let result = infer_product(&[1], &[1], ProductType::Outer, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![2]);
        assert_eq!(result.matching_entity, Some("Entity_2".to_string()));
        assert!(!result.is_zero);

        // Vector · Vector should match Entity_0 (grade 0)
        let result = infer_product(&[1], &[1], ProductType::Inner, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![0]);
        assert_eq!(result.matching_entity, Some("Entity_0".to_string()));
        assert!(!result.is_zero);
    }

    #[test]
    fn infer_product_constraint_propagation() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Only define Entity_0 and Entity_1 (no Entity_0_2)
        let entities = vec![
            (vec![0], "Entity_0".to_string()),
            (vec![1], "Entity_1".to_string()),
        ];

        // Vector * Vector produces raw [0, 2], but Entity_0_2 doesn't exist
        // Constraint propagation finds Entity_0 as the best matching subset
        let result = infer_product(&[1], &[1], ProductType::Geometric, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![0]);
        assert_eq!(result.matching_entity, Some("Entity_0".to_string()));
        assert!(!result.is_zero);
    }

    #[test]
    fn infer_product_no_matching_entity() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Define only Entity_3 (trivector) - no subset of [0, 2]
        let entities = vec![(vec![3], "Entity_3".to_string())];

        // Vector * Vector produces raw [0, 2], no subset matches
        let result = infer_product(&[1], &[1], ProductType::Geometric, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![0, 2]);
        assert!(result.matching_entity.is_none());
        assert!(!result.is_zero);
    }

    #[test]
    fn infer_product_zero() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        let entities = vec![(vec![2], "Entity_2".to_string())];

        // Bivector ∧ Bivector = 0 in 3D
        let result = infer_product(&[2], &[2], ProductType::Outer, &entities, &algebra, &table);
        assert!(result.output_grades.is_empty());
        assert!(result.matching_entity.is_none());
        assert!(result.is_zero);
    }

    #[test]
    fn infer_all_products_basic() {
        let algebra = Algebra::euclidean(3);

        // Include Entity_0_2 so Vector * Vector can match exactly
        let entities = vec![
            ("Entity_0".to_string(), vec![0]),
            ("Entity_1".to_string(), vec![1]),
            ("Entity_0_2".to_string(), vec![0, 2]),
        ];

        let table = infer_all_products(&entities, ProductType::Geometric, &algebra);
        assert_eq!(table.product_type, ProductType::Geometric);
        assert_eq!(table.entries.len(), 9); // 3x3 = 9 combinations

        // Find Entity_1 * Entity_1 entry
        let v_times_v = table.entries.iter().find(|(lhs, rhs, _)| lhs == "Entity_1" && rhs == "Entity_1");
        assert!(v_times_v.is_some());
        let (_, _, result) = v_times_v.unwrap();
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Entity_0_2".to_string()));
    }

    #[test]
    fn pga_null_basis() {
        let algebra = Algebra::pga(3); // 4D with degenerate basis

        // In PGA, the null basis e4 squares to 0
        // This affects products involving grade 4 or higher
        let table = ProductTable::new(&algebra);

        // Single vector still produces scalar + bivector
        let output = infer_output_grades_precise(&[1], &[1], ProductType::Geometric, &algebra, &table);
        assert!(output.contains(&0));
        assert!(output.contains(&2));
    }

    #[test]
    fn constraint_propagation_finds_best_subset() {
        // Test that find_best_matching_subset finds the largest matching entity
        let known_entities = vec![
            (vec![0], "Entity_0".to_string()),
            (vec![1], "Entity_1".to_string()),
            (vec![0, 2], "Entity_0_2".to_string()),
            (vec![1, 3], "Entity_1_3".to_string()),
        ];

        // Raw output [0, 1, 2, 3] - should match the largest subset
        let raw_set: BTreeSet<usize> = [0, 1, 2, 3].iter().copied().collect();
        let result = find_best_matching_subset(&raw_set, &known_entities);

        // Should match Entity_0_2 or Entity_1_3 (both have 2 grades)
        assert!(result.is_some());
        let (grades, _) = result.unwrap();
        assert_eq!(grades.len(), 2);

        // Raw output [0, 2] - should match Entity_0_2 exactly
        let raw_set: BTreeSet<usize> = [0, 2].iter().copied().collect();
        let result = find_best_matching_subset(&raw_set, &known_entities);
        assert!(result.is_some());
        let (grades, name) = result.unwrap();
        assert_eq!(grades, vec![0, 2]);
        assert_eq!(name, "Entity_0_2");
    }

    #[test]
    fn product_with_constraint_application() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Define entities including Even (rotor) and Odd
        let entities = vec![
            (vec![0], "Entity_0".to_string()),
            (vec![1], "Entity_1".to_string()),
            (vec![2], "Entity_2".to_string()),
            (vec![3], "Entity_3".to_string()),
            (vec![0, 2], "Entity_0_2".to_string()),
            (vec![1, 3], "Entity_1_3".to_string()),
        ];

        // Even * Vector = Odd [1, 3]
        let result = infer_product(&[0, 2], &[1], ProductType::Geometric, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![1, 3]);
        assert_eq!(result.matching_entity, Some("Entity_1_3".to_string()));

        // Odd * Vector = Even [0, 2]
        let result = infer_product(&[1, 3], &[1], ProductType::Geometric, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Entity_0_2".to_string()));

        // Even * Even = Even
        let result = infer_product(&[0, 2], &[0, 2], ProductType::Geometric, &entities, &algebra, &table);
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Entity_0_2".to_string()));
    }
}
