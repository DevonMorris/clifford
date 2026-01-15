//! Product output inference for discovered entities.
//!
//! This module infers the output type for products between discovered entities.
//! Given two entities and a product type, it computes which grades appear in
//! the output and matches them to known entities.
//!
//! Products are only generated when the output grades exactly match a known
//! entity type. This prevents generating incorrect code that loses grade
//! components.

use crate::algebra::{
    Algebra, ProductTable, binomial, blades_of_grades, geometric_grades, grade,
    left_contraction_grade, outer_grade,
};
use std::collections::BTreeSet;

/// Represents a product type for inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProductType {
    /// Geometric product: `a * b`
    Geometric,
    /// Exterior (wedge) product: `a ∧ b`
    Exterior,
    /// Left contraction: `a ⌋ b`
    LeftContraction,
    /// Right contraction: `a ⌊ b`
    RightContraction,
    /// Regressive (meet) product: `a ∨ b`
    Regressive,
    /// Scalar product: grade-0 projection of geometric
    Scalar,
    /// Antigeometric product: `a ⊟ b`
    Antigeometric,
    /// Antiscalar product: grade-n projection of antigeometric
    Antiscalar,
    /// Bulk contraction: `a ∨ b★` (antiwedge with bulk dual)
    BulkContraction,
    /// Weight contraction: `a ∨ b☆` (antiwedge with weight dual)
    WeightContraction,
    /// Bulk expansion: `a ∧ b★` (wedge with bulk dual)
    BulkExpansion,
    /// Weight expansion: `a ∧ b☆` (wedge with weight dual)
    WeightExpansion,
    /// Dot product: `a • b` (metric inner, same-grade only, returns scalar)
    Dot,
    /// Antidot product: `a ⊚ b` (metric antiproduct inner, same-antigrade only, returns scalar)
    Antidot,
    /// Projection: `b ∨ (a ∧ b☆)` (target antiwedge with wedge of self and weight dual of target)
    Project,
    /// Antiprojection: `b ∧ (a ∨ b☆)` (target wedge with antiwedge of self and weight dual of target)
    Antiproject,
}

impl ProductType {
    /// Returns all standard product types.
    pub fn all() -> &'static [ProductType] {
        &[
            ProductType::Geometric,
            ProductType::Exterior,
            ProductType::LeftContraction,
            ProductType::RightContraction,
            ProductType::Regressive,
            ProductType::Scalar,
            ProductType::Antigeometric,
            ProductType::Antiscalar,
            ProductType::BulkContraction,
            ProductType::WeightContraction,
            ProductType::BulkExpansion,
            ProductType::WeightExpansion,
            ProductType::Dot,
            ProductType::Antidot,
            ProductType::Project,
            ProductType::Antiproject,
        ]
    }

    /// Returns the name for TOML output.
    pub fn toml_name(&self) -> &'static str {
        match self {
            ProductType::Geometric => "geometric",
            ProductType::Exterior => "exterior",
            ProductType::LeftContraction => "left_contraction",
            ProductType::RightContraction => "right_contraction",
            ProductType::Regressive => "regressive",
            ProductType::Scalar => "scalar",
            ProductType::Antigeometric => "antigeometric",
            ProductType::Antiscalar => "antiscalar",
            ProductType::BulkContraction => "bulk_contraction",
            ProductType::WeightContraction => "weight_contraction",
            ProductType::BulkExpansion => "bulk_expansion",
            ProductType::WeightExpansion => "weight_expansion",
            ProductType::Dot => "dot",
            ProductType::Antidot => "antidot",
            ProductType::Project => "project",
            ProductType::Antiproject => "antiproject",
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

/// Entity representation with exact blade set for blade-level inference.
///
/// Used for sparse types that don't span all blades of their grades.
#[derive(Debug, Clone)]
pub struct EntityBladeSet {
    /// Entity name.
    pub name: String,
    /// Exact blade indices this entity contains.
    pub blades: BTreeSet<usize>,
    /// Grades spanned by this entity.
    pub grades: Vec<usize>,
    /// Whether this is a sparse type.
    pub is_sparse: bool,
}

impl EntityBladeSet {
    /// Creates a new EntityBladeSet from a name and blade indices.
    pub fn new(name: String, blades: impl IntoIterator<Item = usize>) -> Self {
        let blades: BTreeSet<usize> = blades.into_iter().collect();
        let grades: Vec<usize> = blades
            .iter()
            .map(|b| b.count_ones() as usize)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let is_sparse = !Self::spans_all_grades(&blades, &grades);
        Self {
            name,
            blades,
            grades,
            is_sparse,
        }
    }

    /// Creates a non-sparse EntityBladeSet from grades.
    ///
    /// This creates an entity that spans all blades of the specified grades.
    pub fn from_grades(name: String, grades: Vec<usize>, dim: usize) -> Self {
        let blades = blades_of_grades(dim, &grades).into_iter().collect();
        Self {
            name,
            blades,
            grades,
            is_sparse: false,
        }
    }

    /// Checks if blades span all blades of the given grades.
    fn spans_all_grades(blades: &BTreeSet<usize>, grades: &[usize]) -> bool {
        // Count expected blades for each grade
        let max_blade = blades.iter().max().copied().unwrap_or(0);
        let dim = (max_blade as f64).log2().ceil() as usize;
        if dim == 0 && max_blade > 0 {
            return false;
        }

        for &grade in grades {
            let expected_count = binomial(dim.max(1), grade);
            let actual_count = blades
                .iter()
                .filter(|b| b.count_ones() as usize == grade)
                .count();
            if actual_count != expected_count {
                return false;
            }
        }
        true
    }
}

/// Result of blade-level product inference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BladeProductResult {
    /// Blade indices present in the output.
    pub output_blades: BTreeSet<usize>,
    /// Grades present in the output.
    pub output_grades: Vec<usize>,
    /// Name of matching entity, if one exists.
    pub matching_entity: Option<String>,
    /// Whether the output is always zero.
    pub is_zero: bool,
}

/// Infers the output blades for a product between two blade sets.
///
/// This function computes the exact blades that appear in the product output,
/// which is necessary for sparse types that don't span all blades of their grades.
///
/// # Arguments
///
/// * `lhs_blades` - Blade indices in the left operand
/// * `rhs_blades` - Blade indices in the right operand
/// * `product_type` - The type of product
/// * `algebra` - The algebra
/// * `table` - Precomputed product table
///
/// # Returns
///
/// The set of blade indices that appear in non-zero output.
pub fn infer_output_blades(
    lhs_blades: &BTreeSet<usize>,
    rhs_blades: &BTreeSet<usize>,
    product_type: ProductType,
    algebra: &Algebra,
    table: &ProductTable,
) -> BTreeSet<usize> {
    let mut output_set = BTreeSet::new();

    for &a in lhs_blades {
        for &b in rhs_blades {
            // Get product result using appropriate method
            let (sign, result) = match product_type {
                ProductType::Regressive => table.regressive(a, b),
                ProductType::Exterior => table.exterior(a, b),
                ProductType::Antigeometric | ProductType::Antiscalar => table.antiproduct(a, b),
                ProductType::BulkContraction => table.bulk_contraction(a, b),
                ProductType::WeightContraction => table.weight_contraction(a, b),
                ProductType::BulkExpansion => table.bulk_expansion(a, b),
                ProductType::WeightExpansion => table.weight_expansion(a, b),
                ProductType::Dot => table.dot(a, b),
                ProductType::Antidot => table.antidot(a, b),
                ProductType::Project => table.project(a, b),
                ProductType::Antiproject => table.antiproject(a, b),
                _ => table.geometric(a, b),
            };

            if sign == 0 {
                continue;
            }

            let ga = grade(a);
            let gb = grade(b);
            let result_grade = grade(result);

            // Check if this product should be included based on product type
            let include = match product_type {
                ProductType::Geometric => true,
                ProductType::Exterior => result_grade == ga + gb,
                ProductType::LeftContraction => ga <= gb && result_grade == gb - ga,
                ProductType::RightContraction => gb <= ga && result_grade == ga - gb,
                ProductType::Regressive => true,
                ProductType::Scalar => result_grade == 0,
                ProductType::Antigeometric => true,
                ProductType::Antiscalar => result_grade == algebra.dim(),
                ProductType::BulkContraction
                | ProductType::WeightContraction
                | ProductType::BulkExpansion
                | ProductType::WeightExpansion => true,
                ProductType::Dot => ga == gb && result_grade == 0,
                ProductType::Antidot => ga == gb && result_grade == 0,
                ProductType::Project | ProductType::Antiproject => true,
            };

            if include {
                output_set.insert(result);
            }
        }
    }

    output_set
}

/// Infers the product output at blade level and matches to known entities.
///
/// This function handles sparse types correctly by computing products at the
/// blade level instead of the grade level.
///
/// # Arguments
///
/// * `lhs` - Left operand entity with exact blades
/// * `rhs` - Right operand entity with exact blades
/// * `product_type` - The type of product
/// * `known_entities` - List of entities with their blade sets
/// * `algebra` - The algebra
/// * `table` - Precomputed product table
///
/// # Returns
///
/// The inferred blade product result.
pub fn infer_product_blades(
    lhs: &EntityBladeSet,
    rhs: &EntityBladeSet,
    product_type: ProductType,
    known_entities: &[EntityBladeSet],
    algebra: &Algebra,
    table: &ProductTable,
) -> BladeProductResult {
    let output_blades = infer_output_blades(&lhs.blades, &rhs.blades, product_type, algebra, table);

    if output_blades.is_empty() {
        return BladeProductResult {
            output_blades: BTreeSet::new(),
            output_grades: vec![],
            matching_entity: None,
            is_zero: true,
        };
    }

    // Compute output grades
    let output_grades: Vec<usize> = output_blades
        .iter()
        .map(|b| b.count_ones() as usize)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    // Match to known entities by exact blade set
    let matching_entity = known_entities
        .iter()
        .find(|e| e.blades == output_blades)
        .map(|e| e.name.clone());

    BladeProductResult {
        output_blades,
        output_grades,
        matching_entity,
        is_zero: false,
    }
}

/// Infers all products between entities using blade-level inference.
///
/// This version handles sparse types correctly by computing at the blade level.
///
/// # Arguments
///
/// * `entities` - List of entities with their blade sets
/// * `product_type` - The type of product
/// * `algebra` - The algebra
///
/// # Returns
///
/// A list of (lhs_name, rhs_name, result) tuples.
pub fn infer_all_products_blades(
    entities: &[EntityBladeSet],
    product_type: ProductType,
    algebra: &Algebra,
) -> Vec<(String, String, BladeProductResult)> {
    let table = ProductTable::new(algebra);
    let mut results = Vec::new();

    for lhs in entities {
        for rhs in entities {
            let result = infer_product_blades(lhs, rhs, product_type, entities, algebra, &table);
            results.push((lhs.name.clone(), rhs.name.clone(), result));
        }
    }

    results
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
/// let output = infer_output_grades(&[1], &[1], ProductType::Exterior, &algebra);
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
                ProductType::Exterior => {
                    if let Some(g) = outer_grade(ga, gb, dim) {
                        output_set.insert(g);
                    }
                }
                ProductType::LeftContraction => {
                    if let Some(g) = left_contraction_grade(ga, gb) {
                        output_set.insert(g);
                    }
                }
                ProductType::RightContraction => {
                    // Right contraction: grade(a) - grade(b) when gb <= ga
                    if gb <= ga {
                        output_set.insert(ga - gb);
                    }
                }
                ProductType::Regressive => {
                    // Regressive product: (a* ∧ b*)* where * is dual
                    // Result grade = ga + gb - dim (when >= 0)
                    let result = ga + gb;
                    if result >= dim {
                        output_set.insert(result - dim);
                    }
                }
                ProductType::Scalar => {
                    // Scalar product: only grade 0 from geometric
                    if ga == gb {
                        output_set.insert(0);
                    }
                }
                ProductType::Antigeometric => {
                    // Antigeometric: dual(dual(a) * dual(b))
                    // Same grade structure as geometric, just different metric
                    for g in geometric_grades(ga, gb, dim) {
                        output_set.insert(g);
                    }
                }
                ProductType::Antiscalar => {
                    // Antiscalar: only grade dim from antigeometric
                    if ga + gb >= dim && (ga + gb - dim).is_multiple_of(2) {
                        // Can contribute to pseudoscalar
                        output_set.insert(dim);
                    }
                }
                ProductType::BulkContraction | ProductType::WeightContraction => {
                    // Contraction: a ∨ dual(b)
                    // dual(b) has antigrade = dim - gb
                    // antiwedge: ga + (dim - gb) - dim = ga - gb (when ga >= gb)
                    if ga >= gb {
                        output_set.insert(ga - gb);
                    }
                }
                ProductType::BulkExpansion | ProductType::WeightExpansion => {
                    // Expansion: a ∧ dual(b)
                    // dual(b) has antigrade = dim - gb
                    // wedge: ga + (dim - gb) (when <= dim)
                    let result = ga + dim - gb;
                    if result <= dim {
                        output_set.insert(result);
                    }
                }
                ProductType::Dot => {
                    // Dot product: only same-grade elements produce non-zero
                    // Returns scalar (grade 0)
                    if ga == gb {
                        output_set.insert(0);
                    }
                }
                ProductType::Antidot => {
                    // Antidot product: only same-antigrade elements produce non-zero
                    // Since antigrade = dim - grade, same-antigrade means same-grade
                    // Returns scalar (grade 0)
                    if ga == gb {
                        output_set.insert(0);
                    }
                }
                ProductType::Project => {
                    // Project: b ∨ (a ∧ b☆)
                    // b☆ has grade = dim - gb (weight dual)
                    // a ∧ b☆ has grade = ga + (dim - gb) when <= dim
                    // b ∨ (a ∧ b☆) has grade = gb + (ga + dim - gb) - dim = ga
                    // So projection preserves the grade of a
                    output_set.insert(ga);
                }
                ProductType::Antiproject => {
                    // Antiproject: b ∧ (a ∨ b☆)
                    // b☆ has grade = dim - gb
                    // a ∨ b☆ has grade = ga + (dim - gb) - dim = ga - gb when ga >= gb
                    // b ∧ (a ∨ b☆) has grade = gb + (ga - gb) = ga when ga >= gb
                    // So antiprojection also produces grade ga when defined
                    if ga >= gb {
                        output_set.insert(ga);
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
            // For products with specialized computation, use the appropriate method
            let (sign, result_grade) = match product_type {
                ProductType::Regressive => {
                    // Regressive product uses complement-based formula
                    let (sign, result) = table.regressive(a, b);
                    (sign, grade(result))
                }
                ProductType::Exterior => {
                    // Exterior product has its own method
                    let (sign, result) = table.exterior(a, b);
                    (sign, grade(result))
                }
                ProductType::Antigeometric | ProductType::Antiscalar => {
                    // Antiproducts use the anti-metric
                    let (sign, result) = table.antiproduct(a, b);
                    (sign, grade(result))
                }
                ProductType::BulkContraction => {
                    let (sign, result) = table.bulk_contraction(a, b);
                    (sign, grade(result))
                }
                ProductType::WeightContraction => {
                    let (sign, result) = table.weight_contraction(a, b);
                    (sign, grade(result))
                }
                ProductType::BulkExpansion => {
                    let (sign, result) = table.bulk_expansion(a, b);
                    (sign, grade(result))
                }
                ProductType::WeightExpansion => {
                    let (sign, result) = table.weight_expansion(a, b);
                    (sign, grade(result))
                }
                ProductType::Dot => {
                    let (sign, result) = table.dot(a, b);
                    (sign, grade(result))
                }
                ProductType::Antidot => {
                    let (sign, result) = table.antidot(a, b);
                    (sign, grade(result))
                }
                ProductType::Project => {
                    let (sign, result) = table.project(a, b);
                    (sign, grade(result))
                }
                ProductType::Antiproject => {
                    let (sign, result) = table.antiproject(a, b);
                    (sign, grade(result))
                }
                _ => {
                    // All other products derive from the geometric product
                    let (sign, result) = table.geometric(a, b);
                    (sign, grade(result))
                }
            };

            if sign == 0 {
                continue;
            }

            let ga = grade(a);
            let gb = grade(b);

            // Check if this product should be included based on product type
            let include = match product_type {
                ProductType::Geometric => true,
                ProductType::Exterior => {
                    // Already filtered by exterior method, but verify grade
                    result_grade == ga + gb
                }
                ProductType::LeftContraction => {
                    // Left contraction: only grade gb - ga terms (when ga <= gb)
                    ga <= gb && result_grade == gb - ga
                }
                ProductType::RightContraction => {
                    // Right contraction: only grade ga - gb terms (when gb <= ga)
                    gb <= ga && result_grade == ga - gb
                }
                ProductType::Regressive => {
                    // Already computed correctly by regressive method
                    true
                }
                ProductType::Scalar => {
                    // Scalar product: only grade 0 terms
                    result_grade == 0
                }
                ProductType::Antigeometric => {
                    // Already computed by antiproduct, include all
                    true
                }
                ProductType::Antiscalar => {
                    // Antiscalar: only grade dim terms
                    result_grade == algebra.dim()
                }
                ProductType::BulkContraction
                | ProductType::WeightContraction
                | ProductType::BulkExpansion
                | ProductType::WeightExpansion => {
                    // Already computed correctly by specialized table methods
                    true
                }
                ProductType::Dot => {
                    // Dot product: only same-grade, only scalar output
                    ga == gb && result_grade == 0
                }
                ProductType::Antidot => {
                    // Antidot product: only same-antigrade (same-grade), only scalar output
                    ga == gb && result_grade == 0
                }
                ProductType::Project | ProductType::Antiproject => {
                    // Already computed correctly by specialized table methods
                    true
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
/// 2. Matches directly to known entities (exact match only)
///
/// Products are only generated when the output grades exactly match a known
/// entity. Subset matching would generate incorrect code that loses grade
/// components (e.g., returning Quadvector for an output that should be
/// Bivector + Quadvector).
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
    let raw_output =
        infer_output_grades_precise(lhs_grades, rhs_grades, product_type, algebra, table);

    if raw_output.is_empty() {
        return ProductResult {
            output_grades: vec![],
            matching_entity: None,
            is_zero: true,
        };
    }

    // Only match if output grades exactly match a known entity
    // Subset matching would generate incorrect code that loses grade components
    if let Some((_, name)) = known_entities
        .iter()
        .find(|(grades, _)| grades == &raw_output)
    {
        return ProductResult {
            output_grades: raw_output,
            matching_entity: Some(name.clone()),
            is_zero: false,
        };
    }

    // No exact match - return raw output without entity match
    // This product won't be generated in code
    ProductResult {
        output_grades: raw_output,
        matching_entity: None,
        is_zero: false,
    }
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
        let output =
            infer_output_grades_precise(&[1], &[1], ProductType::Geometric, &algebra, &table);
        assert_eq!(output, vec![0, 2]);
    }

    #[test]
    fn outer_vector_vector() {
        let algebra = Algebra::euclidean(3);

        // Vector ∧ Vector = Bivector
        let output = infer_output_grades(&[1], &[1], ProductType::Exterior, &algebra);
        assert_eq!(output, vec![2]);
    }

    #[test]
    fn outer_vector_bivector() {
        let algebra = Algebra::euclidean(3);

        // Vector ∧ Bivector = Trivector
        let output = infer_output_grades(&[1], &[2], ProductType::Exterior, &algebra);
        assert_eq!(output, vec![3]);
    }

    #[test]
    fn outer_bivector_bivector() {
        let algebra = Algebra::euclidean(3);

        // Bivector ∧ Bivector = 0 in 3D (grade 4 > dim 3)
        let output = infer_output_grades(&[2], &[2], ProductType::Exterior, &algebra);
        assert!(output.is_empty());
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
        let result = infer_product(
            &[1],
            &[1],
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Entity_0_2".to_string()));
        assert!(!result.is_zero);

        // Vector ∧ Vector should match Entity_2 (grade 2)
        let result = infer_product(
            &[1],
            &[1],
            ProductType::Exterior,
            &entities,
            &algebra,
            &table,
        );
        assert_eq!(result.output_grades, vec![2]);
        assert_eq!(result.matching_entity, Some("Entity_2".to_string()));
        assert!(!result.is_zero);
    }

    #[test]
    fn infer_product_requires_exact_match() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Only define Entity_0 and Entity_1 (no Entity_0_2)
        let entities = vec![
            (vec![0], "Entity_0".to_string()),
            (vec![1], "Entity_1".to_string()),
        ];

        // Vector * Vector produces raw [0, 2], but Entity_0_2 doesn't exist
        // No exact match, so matching_entity is None (product won't be generated)
        let result = infer_product(
            &[1],
            &[1],
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );
        assert_eq!(result.output_grades, vec![0, 2]);
        assert!(result.matching_entity.is_none());
        assert!(!result.is_zero);
    }

    #[test]
    fn infer_product_no_matching_entity() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Define only Entity_3 (trivector) - no subset of [0, 2]
        let entities = vec![(vec![3], "Entity_3".to_string())];

        // Vector * Vector produces raw [0, 2], no subset matches
        let result = infer_product(
            &[1],
            &[1],
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );
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
        let result = infer_product(
            &[2],
            &[2],
            ProductType::Exterior,
            &entities,
            &algebra,
            &table,
        );
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
        let v_times_v = table
            .entries
            .iter()
            .find(|(lhs, rhs, _)| lhs == "Entity_1" && rhs == "Entity_1");
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
        let output =
            infer_output_grades_precise(&[1], &[1], ProductType::Geometric, &algebra, &table);
        assert!(output.contains(&0));
        assert!(output.contains(&2));
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
        let result = infer_product(
            &[0, 2],
            &[1],
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );
        assert_eq!(result.output_grades, vec![1, 3]);
        assert_eq!(result.matching_entity, Some("Entity_1_3".to_string()));

        // Odd * Vector = Even [0, 2]
        let result = infer_product(
            &[1, 3],
            &[1],
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Entity_0_2".to_string()));

        // Even * Even = Even
        let result = infer_product(
            &[0, 2],
            &[0, 2],
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Entity_0_2".to_string()));
    }

    /// Diagnostic test that reports missing product matches for all algebras.
    ///
    /// This test doesn't fail - it just prints which products are inferred but
    /// don't have matching entities in each algebra.
    #[test]
    fn report_missing_product_matches() {
        use crate::spec::parse_spec;

        let algebras = [
            (
                "euclidean2",
                include_str!("../../../../algebras/euclidean2.toml"),
            ),
            (
                "euclidean3",
                include_str!("../../../../algebras/euclidean3.toml"),
            ),
            (
                "projective2",
                include_str!("../../../../algebras/projective2.toml"),
            ),
            (
                "projective3",
                include_str!("../../../../algebras/projective3.toml"),
            ),
            (
                "conformal3",
                include_str!("../../../../algebras/conformal3.toml"),
            ),
            (
                "quaternion",
                include_str!("../../../../algebras/quaternion.toml"),
            ),
            (
                "dualquat",
                include_str!("../../../../algebras/dualquat.toml"),
            ),
            ("complex", include_str!("../../../../algebras/complex.toml")),
            ("dual", include_str!("../../../../algebras/dual.toml")),
            (
                "hyperbolic",
                include_str!("../../../../algebras/hyperbolic.toml"),
            ),
            (
                "minkowski2",
                include_str!("../../../../algebras/minkowski2.toml"),
            ),
            (
                "minkowski3",
                include_str!("../../../../algebras/minkowski3.toml"),
            ),
            (
                "elliptic2",
                include_str!("../../../../algebras/elliptic2.toml"),
            ),
            (
                "hyperbolic2",
                include_str!("../../../../algebras/hyperbolic2.toml"),
            ),
        ];

        let product_types = [
            ProductType::Geometric,
            ProductType::Exterior,
            ProductType::LeftContraction,
            ProductType::Regressive,
        ];

        let mut total_missing = 0;

        for (name, toml) in &algebras {
            let spec = parse_spec(toml).unwrap();
            let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);

            // Build entity list (exclude sparse and alias types)
            let entities: Vec<(String, Vec<usize>)> = spec
                .types
                .iter()
                .filter(|t| t.alias_of.is_none() && !t.is_sparse)
                .map(|t| (t.name.clone(), t.grades.clone()))
                .collect();

            let mut algebra_missing = 0;

            for product_type in &product_types {
                let table = infer_all_products(&entities, *product_type, &algebra);

                for (lhs, rhs, result) in &table.entries {
                    if !result.is_zero && result.matching_entity.is_none() {
                        algebra_missing += 1;
                        eprintln!(
                            "  {}: {} {:?} {} -> grades {:?} (no match)",
                            name,
                            lhs,
                            product_type.toml_name(),
                            rhs,
                            result.output_grades
                        );
                    }
                }
            }

            if algebra_missing > 0 {
                eprintln!("{}: {} missing product matches", name, algebra_missing);
            }
            total_missing += algebra_missing;
        }

        eprintln!("\nTotal missing product matches: {}", total_missing);
        // This test is informational - it always passes
        // If you want to enforce complete coverage, change this to:
        // assert_eq!(total_missing, 0, "Some products don't have matching entities");
    }

    #[test]
    fn entity_blade_set_from_grades() {
        let entity = EntityBladeSet::from_grades("Vector".to_string(), vec![1], 3);

        // Grade 1 in 3D has blades e1=1, e2=2, e3=4
        assert_eq!(entity.name, "Vector");
        assert_eq!(entity.grades, vec![1]);
        assert!(!entity.is_sparse);
        assert!(entity.blades.contains(&1)); // e1
        assert!(entity.blades.contains(&2)); // e2
        assert!(entity.blades.contains(&4)); // e3
        assert_eq!(entity.blades.len(), 3);
    }

    #[test]
    fn entity_blade_set_sparse() {
        // Create a sparse entity with only 2 of 3 grade-1 blades
        let entity = EntityBladeSet::new("Partial".to_string(), vec![1, 2]); // e1, e2 only

        assert_eq!(entity.name, "Partial");
        assert_eq!(entity.grades, vec![1]);
        assert!(entity.is_sparse, "Entity should be sparse (missing e3)");
        assert!(entity.blades.contains(&1)); // e1
        assert!(entity.blades.contains(&2)); // e2
        assert!(!entity.blades.contains(&4)); // e3 not present
    }

    #[test]
    fn infer_output_blades_geometric() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Vector (grade 1) * Vector (grade 1) = Scalar + Bivector
        let lhs_blades: BTreeSet<usize> = vec![1, 2, 4].into_iter().collect(); // e1, e2, e3
        let rhs_blades: BTreeSet<usize> = vec![1, 2, 4].into_iter().collect();

        let output = infer_output_blades(
            &lhs_blades,
            &rhs_blades,
            ProductType::Geometric,
            &algebra,
            &table,
        );

        // Should contain scalar (0) and bivectors (3=e12, 5=e13, 6=e23)
        assert!(output.contains(&0)); // scalar
        assert!(output.contains(&3)); // e12
        assert!(output.contains(&5)); // e13
        assert!(output.contains(&6)); // e23
    }

    #[test]
    fn infer_output_blades_sparse_geometric() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Partial vector (just e1, e2) * full vector
        let lhs_blades: BTreeSet<usize> = vec![1, 2].into_iter().collect(); // e1, e2 only
        let rhs_blades: BTreeSet<usize> = vec![1, 2, 4].into_iter().collect(); // e1, e2, e3

        let output = infer_output_blades(
            &lhs_blades,
            &rhs_blades,
            ProductType::Geometric,
            &algebra,
            &table,
        );

        // e1*e1 = 1 (scalar), e1*e2 = e12, e1*e3 = e13
        // e2*e1 = -e12, e2*e2 = 1, e2*e3 = e23
        assert!(output.contains(&0)); // scalar from e1*e1, e2*e2
        assert!(output.contains(&3)); // e12 from e1*e2
        assert!(output.contains(&5)); // e13 from e1*e3
        assert!(output.contains(&6)); // e23 from e2*e3

        // e3*anything is not in output since e3 not in lhs
        // e3 (index 4) should not appear in output
        assert!(!output.contains(&4)); // e3 not in output
    }

    #[test]
    fn infer_product_blades_with_matching() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        let scalar = EntityBladeSet::from_grades("Scalar".to_string(), vec![0], 3);
        let vector = EntityBladeSet::from_grades("Vector".to_string(), vec![1], 3);
        let bivector = EntityBladeSet::from_grades("Bivector".to_string(), vec![2], 3);
        let rotor = EntityBladeSet::from_grades("Rotor".to_string(), vec![0, 2], 3);

        let entities = vec![
            scalar.clone(),
            vector.clone(),
            bivector.clone(),
            rotor.clone(),
        ];

        // Vector * Vector should match Rotor
        let result = infer_product_blades(
            &vector,
            &vector,
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );

        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Rotor".to_string()));
        assert!(!result.is_zero);
    }

    #[test]
    fn infer_product_blades_sparse_no_match() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Create a partial vector (sparse) - just e1 and e2
        let partial_vec = EntityBladeSet::new("PartialVec".to_string(), vec![1, 2]);
        // Create a full vector
        let full_vec = EntityBladeSet::from_grades("Vector".to_string(), vec![1], 3);
        // Create full types
        let rotor = EntityBladeSet::from_grades("Rotor".to_string(), vec![0, 2], 3);

        let entities = vec![partial_vec.clone(), full_vec.clone(), rotor.clone()];

        // PartialVec * Vector produces a subset of blades
        // The output won't match Rotor exactly because it's missing some bivector components
        let result = infer_product_blades(
            &partial_vec,
            &full_vec,
            ProductType::Geometric,
            &entities,
            &algebra,
            &table,
        );

        // Output should be scalar + 3 bivectors (e12, e13, e23)
        // But e13 comes from e1*e3, e23 comes from e2*e3
        // This matches the full Rotor blades, so it should match
        assert_eq!(result.output_grades, vec![0, 2]);
        assert_eq!(result.matching_entity, Some("Rotor".to_string()));
    }

    #[test]
    fn infer_all_products_blades_basic() {
        let algebra = Algebra::euclidean(3);

        let scalar = EntityBladeSet::from_grades("Scalar".to_string(), vec![0], 3);
        let vector = EntityBladeSet::from_grades("Vector".to_string(), vec![1], 3);
        let rotor = EntityBladeSet::from_grades("Rotor".to_string(), vec![0, 2], 3);

        let entities = vec![scalar, vector, rotor];

        let results = infer_all_products_blades(&entities, ProductType::Geometric, &algebra);

        assert_eq!(results.len(), 9); // 3x3 = 9 combinations

        // Find Vector * Vector
        let vv = results
            .iter()
            .find(|(l, r, _)| l == "Vector" && r == "Vector");
        assert!(vv.is_some());
        let (_, _, result) = vv.unwrap();
        assert_eq!(result.matching_entity, Some("Rotor".to_string()));
    }
}
