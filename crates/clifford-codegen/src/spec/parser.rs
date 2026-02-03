//! Parser for algebra specifications.
//!
//! Converts raw TOML structures to validated IR types.

use std::collections::{HashMap, HashSet};

use crate::algebra::{Algebra, binomial, versor_parity};
use crate::discovery::{
    EntityBladeSet, ProductType, infer_all_products, infer_all_products_blades,
};

use super::error::{MissingProduct, ParseError};
use super::ir::{
    AlgebraSpec, BasisVector, FieldSpec, InvolutionKind, NormSpec, ProductEntry, ProductsSpec,
    SignatureSpec, TypeSpec, VersorSpec,
};
use super::raw::{RawAlgebraSpec, RawNormSpec, RawSignature, RawTypeSpec};

/// Maximum supported dimension.
const MAX_DIM: usize = 6;

/// Parses a TOML specification into the IR.
///
/// # Arguments
///
/// * `toml_content` - The TOML specification as a string.
///
/// # Returns
///
/// A validated `AlgebraSpec` or an error describing what went wrong.
///
/// # Example
///
/// ```
/// use clifford_codegen::spec::parse_spec;
///
/// let spec = parse_spec(r#"
/// [algebra]
/// name = "euclidean2"
/// complete = false
///
/// [signature]
/// positive = ["e1", "e2"]
///
/// [types.Vector]
/// grades = [1]
/// field_map = [
///   { name = "x", blade = "e1" },
///   { name = "y", blade = "e2" }
/// ]
/// "#).unwrap();
///
/// assert_eq!(spec.name, "euclidean2");
/// assert_eq!(spec.signature.dim(), 2);
/// ```
pub fn parse_spec(toml_content: &str) -> Result<AlgebraSpec, ParseError> {
    let raw: RawAlgebraSpec = toml::from_str(toml_content)?;

    // Build signature
    let signature = parse_signature(&raw.signature)?;

    // Build norm configuration
    let norm = parse_norm(&raw.norm)?;

    // Build blade name map
    let blade_names = parse_blade_names(&raw.blades, &signature)?;

    // Build types
    let types = parse_types(&raw.types, &signature, &blade_names)?;

    // Auto-infer products from types (products section in TOML is ignored)
    let products = infer_products_from_types(&types, &signature);

    // Validate the complete specification
    validate_spec(&types)?;

    // Check algebra completeness if requested
    let complete = raw.algebra.complete;
    if complete {
        let missing = check_algebra_completeness(&types, &signature);
        if !missing.is_empty() {
            return Err(format_completeness_error(&raw.algebra.name, missing));
        }
    }

    Ok(AlgebraSpec {
        name: raw.algebra.name,
        module_path: raw.algebra.module_path,
        description: raw.algebra.description,
        signature,
        norm,
        blade_names,
        types,
        products,
        complete,
    })
}

/// Parses the signature section.
///
/// Basis vectors can be specified in any order across the positive/negative/zero
/// arrays. The index of each basis is determined by the number in its name
/// (e.g., "e1" → index 0, "e2" → index 1), NOT by its position in the array.
///
/// This allows physics conventions like Minkowski with space (e1) negative and
/// time (e2) positive:
/// ```toml
/// [signature]
/// positive = ["e2"]  # e2² = +1 (time)
/// negative = ["e1"]  # e1² = -1 (space)
/// ```
fn parse_signature(raw: &RawSignature) -> Result<SignatureSpec, ParseError> {
    let p = raw.positive.len();
    let q = raw.negative.len();
    let r = raw.zero.len();
    let dim = p + q + r;

    if dim == 0 {
        return Err(ParseError::EmptySignature);
    }
    if dim > MAX_DIM {
        return Err(ParseError::DimensionTooLarge(dim));
    }

    // Check for duplicate names
    let mut names = HashSet::new();
    for name in raw
        .positive
        .iter()
        .chain(raw.negative.iter())
        .chain(raw.zero.iter())
    {
        if !names.insert(name) {
            return Err(ParseError::DuplicateBasisName(name.clone()));
        }
    }

    // Determine indexing mode: numeric (e1, e2, etc.) or positional (ep, em, etc.)
    // If ALL names are numeric, use numeric indexing; otherwise use positional
    let all_names: Vec<&String> = raw
        .positive
        .iter()
        .chain(raw.negative.iter())
        .chain(raw.zero.iter())
        .collect();

    let use_numeric_indexing = all_names
        .iter()
        .all(|name| try_parse_basis_index(name, dim).is_some());

    // Build basis vectors
    let mut basis = Vec::with_capacity(dim);

    if use_numeric_indexing {
        // Numeric indexing: extract index from name (e.g., "e1" → index 0)
        // This allows arbitrary metric assignment regardless of array order
        for name in &raw.positive {
            let index = try_parse_basis_index(name, dim).unwrap();
            basis.push(BasisVector {
                name: name.clone(),
                index,
                metric: 1,
            });
        }
        for name in &raw.negative {
            let index = try_parse_basis_index(name, dim).unwrap();
            basis.push(BasisVector {
                name: name.clone(),
                index,
                metric: -1,
            });
        }
        for name in &raw.zero {
            let index = try_parse_basis_index(name, dim).unwrap();
            basis.push(BasisVector {
                name: name.clone(),
                index,
                metric: 0,
            });
        }

        // Sort by index for consistent ordering
        basis.sort_by_key(|b| b.index);

        // Validate that indices are contiguous 0..dim
        for (expected, bv) in basis.iter().enumerate() {
            if bv.index != expected {
                return Err(ParseError::NonContiguousBasisIndices {
                    expected,
                    found: bv.index,
                    name: bv.name.clone(),
                });
            }
        }
    } else {
        // Positional indexing: assign indices in order (positive, then negative, then zero)
        // Used for non-numeric basis names like "ep", "em" in CGA
        let mut index = 0;
        for name in &raw.positive {
            basis.push(BasisVector {
                name: name.clone(),
                index,
                metric: 1,
            });
            index += 1;
        }
        for name in &raw.negative {
            basis.push(BasisVector {
                name: name.clone(),
                index,
                metric: -1,
            });
            index += 1;
        }
        for name in &raw.zero {
            basis.push(BasisVector {
                name: name.clone(),
                index,
                metric: 0,
            });
            index += 1;
        }
    }

    Ok(SignatureSpec { basis, p, q, r })
}

/// Tries to parse a basis vector name like "e1" or "e2" into its 0-based index.
///
/// Returns `Some(index)` if the name follows a numeric pattern, `None` otherwise.
/// This is used to determine whether to use numeric or positional indexing for signatures.
///
/// # Naming Conventions
///
/// - Standard bases: "e1" → 0, "e2" → 1, "e3" → 2, etc. (1-indexed names, 0-indexed result)
/// - PGA convention: "e0" → last index (dim-1) for backward compatibility
fn try_parse_basis_index(name: &str, dim: usize) -> Option<usize> {
    if !name.starts_with('e') {
        return None;
    }

    let digits = &name[1..];
    if digits.is_empty() {
        return None;
    }

    let num: usize = digits.parse().ok()?;

    // Special case: e0 maps to the last index (PGA convention for degenerate basis)
    if num == 0 {
        return Some(dim - 1);
    }

    if num > dim {
        return None;
    }

    // Convert to 0-based (e1 → 0, e2 → 1, etc.)
    Some(num - 1)
}

/// Parses the norm configuration section.
///
/// The `[norm]` section specifies which involution the algebra uses for its
/// canonical norm computation. This affects how `Involute` is generated.
///
/// Options for `primary_involution`:
/// - `"reverse"` (default): Uses reverse involution, `(-1)^(k(k-1)/2)` for grade k
/// - `"grade_involution"`: Uses grade involution, `(-1)^k` for grade k
/// - `"clifford_conjugate"`: Uses Clifford conjugate, `(-1)^(k(k+1)/2)` for grade k
fn parse_norm(raw: &RawNormSpec) -> Result<NormSpec, ParseError> {
    let primary_involution = match raw.primary_involution.as_deref() {
        None | Some("reverse") => InvolutionKind::Reverse,
        Some("grade_involution") => InvolutionKind::GradeInvolution,
        Some("clifford_conjugate") => InvolutionKind::CliffordConjugate,
        Some(other) => {
            return Err(ParseError::InvalidValue {
                field: "norm.primary_involution".to_string(),
                value: other.to_string(),
                expected: "\"reverse\", \"grade_involution\", or \"clifford_conjugate\""
                    .to_string(),
            });
        }
    };

    Ok(NormSpec { primary_involution })
}

/// Parses the blade names section.
fn parse_blade_names(
    raw: &HashMap<String, String>,
    sig: &SignatureSpec,
) -> Result<HashMap<usize, String>, ParseError> {
    let dim = sig.dim();
    let mut blade_names = HashMap::new();

    for (blade_name, field_name) in raw {
        let index = parse_blade_index(blade_name, dim)?;
        blade_names.insert(index, field_name.clone());
    }

    Ok(blade_names)
}

/// Result of parsing a blade name, including sign correction.
#[derive(Debug, Clone, Copy)]
pub struct BladeParseResult {
    /// Canonical blade index (bitmask with bits in standard order).
    pub index: usize,
    /// Sign relative to canonical ordering (+1 or -1).
    /// e.g., e20 = -e02, so sign = -1.
    pub sign: i8,
    /// Grade of the blade (number of basis vectors).
    pub grade: usize,
}

/// Parses a blade name like "e12", "e20", "e312", or "s" into its canonical index and sign.
///
/// Non-canonical orderings (e.g., "e20" instead of "e02") are supported and will
/// compute the appropriate sign correction based on permutation parity.
///
/// # Special cases
///
/// - "s" → index=0, sign=+1, grade=0 (scalar blade)
///
/// # Examples
///
/// - "s" → index=0, sign=+1, grade=0 (scalar)
/// - "e1" → index=1 (0b1), sign=+1, grade=1
/// - "e12" → index=3 (0b11), sign=+1 (already canonical)
/// - "e21" → index=3 (0b11), sign=-1 (one swap needed)
/// - "e20" → index=5 (0b101 for bases 0,2), sign=-1 (one swap: 20→02)
/// - "e312" → index=7 (0b111), sign=+1 (312→132→123 = 2 swaps = even)
fn parse_blade_with_sign(name: &str, dim: usize) -> Result<BladeParseResult, ParseError> {
    // Special case: "s" for scalar (grade 0)
    if name == "s" {
        return Ok(BladeParseResult {
            index: 0,
            sign: 1,
            grade: 0,
        });
    }

    // Must start with 'e'
    if !name.starts_with('e') {
        return Err(ParseError::InvalidBladeName {
            name: name.to_string(),
        });
    }

    let digits = &name[1..];
    if digits.is_empty() {
        return Err(ParseError::InvalidBladeName {
            name: name.to_string(),
        });
    }

    // Parse indices in the order given
    let mut indices: Vec<usize> = Vec::new();
    for c in digits.chars() {
        let digit = c.to_digit(10).ok_or_else(|| ParseError::InvalidBladeName {
            name: name.to_string(),
        })? as usize;

        if digit == 0 || digit > dim {
            return Err(ParseError::BladeIndexOutOfBounds {
                name: name.to_string(),
                index: digit,
                dim,
            });
        }

        // Convert to 0-based
        let idx = digit - 1;

        // Check for duplicates
        if indices.contains(&idx) {
            return Err(ParseError::InvalidBladeName {
                name: name.to_string(),
            });
        }
        indices.push(idx);
    }

    // Compute canonical bitmask
    let mut index = 0usize;
    for &idx in &indices {
        index |= 1 << idx;
    }

    // Compute sign from permutation parity
    // Count inversions: pairs (i,j) where i < j but indices[i] > indices[j]
    let sign = permutation_sign(&indices);

    let grade = indices.len();

    Ok(BladeParseResult { index, sign, grade })
}

/// Computes the sign of a permutation by counting inversions.
///
/// Returns +1 for even permutations, -1 for odd permutations.
fn permutation_sign(indices: &[usize]) -> i8 {
    let mut inversions = 0;
    for i in 0..indices.len() {
        for j in (i + 1)..indices.len() {
            if indices[i] > indices[j] {
                inversions += 1;
            }
        }
    }
    if inversions % 2 == 0 { 1 } else { -1 }
}

/// Parses a blade name like "e12" or "e123" into its canonical index (legacy compatibility).
///
/// This is a wrapper around `parse_blade_with_sign` that only returns the index,
/// ignoring the sign. Use `parse_blade_with_sign` for full sign support.
fn parse_blade_index(name: &str, dim: usize) -> Result<usize, ParseError> {
    Ok(parse_blade_with_sign(name, dim)?.index)
}

/// Parses the types section.
fn parse_types(
    raw: &HashMap<String, RawTypeSpec>,
    sig: &SignatureSpec,
    blade_names: &HashMap<usize, String>,
) -> Result<Vec<TypeSpec>, ParseError> {
    let dim = sig.dim();
    let mut types = Vec::with_capacity(raw.len());
    let mut type_names = HashSet::new();

    for (name, raw_type) in raw {
        if !type_names.insert(name.clone()) {
            return Err(ParseError::DuplicateTypeName(name.clone()));
        }

        let type_spec = parse_type(name, raw_type, dim, sig, blade_names)?;
        types.push(type_spec);
    }

    // Sort by name for deterministic output
    types.sort_by(|a, b| a.name.cmp(&b.name));

    Ok(types)
}

/// Parses a single type definition.
fn parse_type(
    name: &str,
    raw: &RawTypeSpec,
    dim: usize,
    _sig: &SignatureSpec,
    _blade_names: &HashMap<usize, String>,
) -> Result<TypeSpec, ParseError> {
    // Validate grades
    for &grade in &raw.grades {
        if grade > dim {
            return Err(ParseError::InvalidGrade {
                type_name: name.to_string(),
                grade,
                max: dim,
            });
        }
    }

    // field_map is required for all non-alias types
    if raw.alias_of.is_none() && raw.field_map.is_empty() {
        return Err(ParseError::MissingFieldMap {
            type_name: name.to_string(),
        });
    }

    // Build fields from field_map
    let fields = build_fields_from_field_map(&raw.field_map, &raw.grades, dim, name)?;

    // Check for duplicate field names
    let mut field_names = HashSet::new();
    for field in &fields {
        if !field_names.insert(&field.name) {
            return Err(ParseError::DuplicateFieldName {
                type_name: name.to_string(),
                field: field.name.clone(),
            });
        }
    }

    // Check self-alias
    if let Some(alias) = &raw.alias_of {
        if alias == name {
            return Err(ParseError::SelfAlias {
                type_name: name.to_string(),
            });
        }
    }

    // Auto-identify versors based on grade parity
    // A type is a versor if all grades have the same parity (all even or all odd)
    // This includes single-grade types: vectors are reflectors, circles in CGA are inversors
    let versor = if versor_parity(&raw.grades).is_some() {
        Some(VersorSpec {
            // is_unit will be determined by inferred constraints in Phase 4
            is_unit: false,
            // sandwich_targets are auto-inferred during code generation
            sandwich_targets: Vec::new(),
        })
    } else {
        None
    };

    // Determine if this is a sparse type (not all blades of the grade are present)
    let expected_blade_count: usize = raw.grades.iter().map(|&g| binomial(dim, g)).sum();
    let is_sparse = fields.len() < expected_blade_count;

    Ok(TypeSpec {
        name: name.to_string(),
        grades: raw.grades.clone(),
        description: raw.description.clone(),
        fields,
        alias_of: raw.alias_of.clone(),
        versor,
        is_sparse,
        inverse_sandwich_targets: raw.inverse_sandwich_targets.clone(),
    })
}

/// Builds field specs from explicit field_map entries.
///
/// Each entry in field_map specifies a field name and its corresponding blade.
/// Blade names can use non-canonical ordering (e.g., "e20" instead of "e02"),
/// and the appropriate sign correction will be computed.
fn build_fields_from_field_map(
    field_map: &[super::raw::RawFieldMapping],
    grades: &[usize],
    dim: usize,
    type_name: &str,
) -> Result<Vec<FieldSpec>, ParseError> {
    let mut fields = Vec::with_capacity(field_map.len());

    for mapping in field_map {
        // Parse blade name with sign detection
        let blade_result = parse_blade_with_sign(&mapping.blade, dim)?;

        // Validate blade grade matches specified grades
        if !grades.contains(&blade_result.grade) {
            return Err(ParseError::FieldMapGradeMismatch {
                type_name: type_name.to_string(),
                field: mapping.name.clone(),
                blade: mapping.blade.clone(),
                blade_grade: blade_result.grade,
                grades: grades.to_vec(),
            });
        }

        fields.push(FieldSpec {
            name: mapping.name.clone(),
            blade_index: blade_result.index,
            grade: blade_result.grade,
            sign: blade_result.sign,
        });
    }

    Ok(fields)
}

/// Validates that fields in a TypeSpec have canonical blade ordering.
///
/// Fields must be ordered by grade first (ascending), then by blade index
/// within each grade (ascending). This ensures consistent behavior in
/// product computations and conversions.
///
/// # Returns
///
/// `true` if the fields are in canonical order, `false` otherwise.
#[cfg(test)]
pub fn validate_canonical_field_order(ty: &TypeSpec) -> bool {
    if ty.fields.is_empty() {
        return true;
    }

    let mut prev_grade = 0;
    let mut prev_blade_index = 0;

    for (i, field) in ty.fields.iter().enumerate() {
        // Fields must be ordered by grade first
        if field.grade < prev_grade {
            return false;
        }

        // Within a grade, blade indices must be ascending
        if field.grade == prev_grade && i > 0 && field.blade_index <= prev_blade_index {
            return false;
        }

        prev_grade = field.grade;
        prev_blade_index = field.blade_index;
    }

    true
}

/// Infers products automatically from types.
///
/// Products are always auto-inferred from the defined types.
/// All standard product types are generated: geometric, exterior, inner, left/right contraction,
/// regressive, scalar, antigeometric, and antiscalar.
///
/// Sparse types are included in product inference using blade-level computation (PRD-45).
/// This ensures that sparse types like Line and Plane can participate in products.
fn infer_products_from_types(types: &[TypeSpec], signature: &SignatureSpec) -> ProductsSpec {
    // Build algebra for product computation
    let algebra = Algebra::from_metrics(signature.metrics_by_index());
    let dim = signature.dim();

    // Check if any types are sparse
    let has_sparse = types.iter().any(|t| t.is_sparse && t.alias_of.is_none());

    if has_sparse {
        // Use blade-level inference for sparse type support
        infer_products_blade_level(types, &algebra, dim)
    } else {
        // Use grade-level inference for better performance
        infer_products_grade_level(types, &algebra)
    }
}

/// Grade-level product inference (fast path for non-sparse algebras).
fn infer_products_grade_level(types: &[TypeSpec], algebra: &Algebra) -> ProductsSpec {
    // Build entity list for inference (exclude aliases)
    let entities: Vec<(String, Vec<usize>)> = types
        .iter()
        .filter(|t| t.alias_of.is_none())
        .map(|t| (t.name.clone(), t.grades.clone()))
        .collect();

    // Convert inferred products to ProductEntry format
    let convert_entries = |table: crate::discovery::ProductTable2D| -> Vec<ProductEntry> {
        table
            .entries
            .into_iter()
            .filter(|(_, _, result)| !result.is_zero && result.matching_entity.is_some())
            .map(|(lhs, rhs, result)| {
                let output = result.matching_entity.unwrap();
                ProductEntry {
                    lhs,
                    rhs,
                    output: output.clone(),
                    output_constrained: false,
                }
            })
            .collect()
    };

    ProductsSpec {
        geometric: convert_entries(infer_all_products(
            &entities,
            ProductType::Geometric,
            algebra,
        )),
        wedge: convert_entries(infer_all_products(
            &entities,
            ProductType::Exterior,
            algebra,
        )),
        left_contraction: convert_entries(infer_all_products(
            &entities,
            ProductType::LeftContraction,
            algebra,
        )),
        right_contraction: convert_entries(infer_all_products(
            &entities,
            ProductType::RightContraction,
            algebra,
        )),
        antiwedge: convert_entries(infer_all_products(
            &entities,
            ProductType::Regressive,
            algebra,
        )),
        scalar: convert_entries(infer_all_products(&entities, ProductType::Scalar, algebra)),
        antigeometric: convert_entries(infer_all_products(
            &entities,
            ProductType::Antigeometric,
            algebra,
        )),
        antiscalar: convert_entries(infer_all_products(
            &entities,
            ProductType::Antiscalar,
            algebra,
        )),
        bulk_contraction: convert_entries(infer_all_products(
            &entities,
            ProductType::BulkContraction,
            algebra,
        )),
        weight_contraction: convert_entries(infer_all_products(
            &entities,
            ProductType::WeightContraction,
            algebra,
        )),
        bulk_expansion: convert_entries(infer_all_products(
            &entities,
            ProductType::BulkExpansion,
            algebra,
        )),
        weight_expansion: convert_entries(infer_all_products(
            &entities,
            ProductType::WeightExpansion,
            algebra,
        )),
        dot: convert_entries(infer_all_products(&entities, ProductType::Dot, algebra)),
        antidot: convert_entries(infer_all_products(&entities, ProductType::Antidot, algebra)),
        project: convert_entries(infer_all_products(&entities, ProductType::Project, algebra)),
        antiproject: convert_entries(infer_all_products(
            &entities,
            ProductType::Antiproject,
            algebra,
        )),
    }
}

/// Blade-level product inference (supports sparse types).
fn infer_products_blade_level(types: &[TypeSpec], algebra: &Algebra, dim: usize) -> ProductsSpec {
    // Build EntityBladeSet for each type
    let entities: Vec<EntityBladeSet> = types
        .iter()
        .filter(|t| t.alias_of.is_none())
        .map(|t| {
            if t.is_sparse {
                // Use exact blade indices from fields
                let blades = t.fields.iter().map(|f| f.blade_index);
                EntityBladeSet::new(t.name.clone(), blades)
            } else {
                // Use all blades of the grades
                EntityBladeSet::from_grades(t.name.clone(), t.grades.clone(), dim)
            }
        })
        .collect();

    // Convert blade-level results to ProductEntry format
    let convert_entries =
        |results: Vec<(String, String, crate::discovery::BladeProductResult)>| -> Vec<ProductEntry> {
            results
                .into_iter()
                .filter(|(_, _, result)| !result.is_zero && result.matching_entity.is_some())
                .map(|(lhs, rhs, result)| {
                    let output = result.matching_entity.unwrap();
                    ProductEntry {
                        lhs,
                        rhs,
                        output: output.clone(),
                        output_constrained: false,
                    }
                })
                .collect()
        };

    ProductsSpec {
        geometric: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Geometric,
            algebra,
        )),
        wedge: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Exterior,
            algebra,
        )),
        left_contraction: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::LeftContraction,
            algebra,
        )),
        right_contraction: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::RightContraction,
            algebra,
        )),
        antiwedge: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Regressive,
            algebra,
        )),
        scalar: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Scalar,
            algebra,
        )),
        antigeometric: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Antigeometric,
            algebra,
        )),
        antiscalar: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Antiscalar,
            algebra,
        )),
        bulk_contraction: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::BulkContraction,
            algebra,
        )),
        weight_contraction: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::WeightContraction,
            algebra,
        )),
        bulk_expansion: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::BulkExpansion,
            algebra,
        )),
        weight_expansion: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::WeightExpansion,
            algebra,
        )),
        dot: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Dot,
            algebra,
        )),
        antidot: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Antidot,
            algebra,
        )),
        project: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Project,
            algebra,
        )),
        antiproject: convert_entries(infer_all_products_blades(
            &entities,
            ProductType::Antiproject,
            algebra,
        )),
    }
}

/// Validates the complete specification.
fn validate_spec(types: &[TypeSpec]) -> Result<(), ParseError> {
    let type_names: HashSet<_> = types.iter().map(|t| t.name.as_str()).collect();

    // Check alias references
    for ty in types {
        if let Some(alias) = &ty.alias_of {
            if !type_names.contains(alias.as_str()) {
                return Err(ParseError::UnknownType(alias.clone()));
            }
        }
    }

    // Check for alias cycles (simple check - no transitive cycles)
    for ty in types {
        if let Some(alias) = &ty.alias_of {
            let target = types.iter().find(|t| t.name == *alias);
            if let Some(target_type) = target {
                if target_type.alias_of.as_ref() == Some(&ty.name) {
                    return Err(ParseError::AliasCycle {
                        type_name: ty.name.clone(),
                    });
                }
            }
        }
    }

    Ok(())
}

/// Checks if all products between defined types have matching output types.
///
/// Returns a list of products that don't have matching output types.
/// This is used when `complete = true` to enforce algebra completeness.
fn check_algebra_completeness(
    types: &[TypeSpec],
    signature: &SignatureSpec,
) -> Vec<MissingProduct> {
    let algebra = Algebra::from_metrics(signature.metrics_by_index());

    // Build entity list (exclude sparse and alias types for now)
    // TODO: PRD-45 will add blade-level inference for sparse types
    let entities: Vec<(String, Vec<usize>)> = types
        .iter()
        .filter(|t| t.alias_of.is_none() && !t.is_sparse)
        .map(|t| (t.name.clone(), t.grades.clone()))
        .collect();

    let mut missing = Vec::new();

    // Check all product types
    let product_types = [
        ProductType::Geometric,
        ProductType::Exterior,
        ProductType::LeftContraction,
        ProductType::RightContraction,
        ProductType::Regressive,
    ];

    for product_type in &product_types {
        let table = infer_all_products(&entities, *product_type, &algebra);

        for (lhs, rhs, result) in table.entries {
            if !result.is_zero && result.matching_entity.is_none() {
                missing.push(MissingProduct {
                    lhs,
                    rhs,
                    product_type: product_type.toml_name().to_string(),
                    output_grades: result.output_grades,
                });
            }
        }
    }

    missing
}

/// Formats a completeness error with detailed information about missing products.
fn format_completeness_error(name: &str, missing: Vec<MissingProduct>) -> ParseError {
    // Group by output grades
    let mut by_grades: std::collections::HashMap<Vec<usize>, Vec<String>> =
        std::collections::HashMap::new();

    for m in &missing {
        let entry = by_grades.entry(m.output_grades.clone()).or_default();
        entry.push(format!("{} {} {}", m.lhs, m.product_type, m.rhs));
    }

    // Format details
    let mut details = String::new();
    for (grades, products) in &by_grades {
        details.push_str(&format!("  Missing type for grades {:?}:\n", grades));
        for (i, product) in products.iter().enumerate() {
            if i < 3 {
                details.push_str(&format!("    - {}\n", product));
            } else if i == 3 {
                details.push_str(&format!("    ... and {} more\n", products.len() - 3));
                break;
            }
        }
    }

    ParseError::IncompleteAlgebra {
        name: name.to_string(),
        count: missing.len(),
        details,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_minimal_spec() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]
            "#,
        )
        .unwrap();

        assert_eq!(spec.name, "test");
        assert_eq!(spec.signature.p, 2);
        assert_eq!(spec.signature.q, 0);
        assert_eq!(spec.signature.r, 0);
        assert_eq!(spec.signature.dim(), 2);
    }

    #[test]
    fn parse_with_types() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "euclidean2"

            [signature]
            positive = ["e1", "e2"]

            [types.Scalar]
            grades = [0]
            field_map = [{ name = "s", blade = "s" }]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "y", blade = "e2" }
            ]

            [types.Bivector]
            grades = [2]
            field_map = [{ name = "b", blade = "e12" }]

            [types.Rotor]
            grades = [0, 2]
            field_map = [
                { name = "s", blade = "s" },
                { name = "xy", blade = "e12" }
            ]
            "#,
        )
        .unwrap();

        assert_eq!(spec.types.len(), 4);

        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();
        assert_eq!(vector.grades, vec![1]);
        assert_eq!(vector.fields.len(), 2);

        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.grades, vec![0, 2]);
        assert_eq!(rotor.fields.len(), 2);
    }

    #[test]
    fn parse_blade_names_section() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [blades]
            e1 = "x"
            e2 = "y"
            e12 = "xy"
            "#,
        )
        .unwrap();

        assert_eq!(spec.blade_names.get(&1), Some(&"x".to_string()));
        assert_eq!(spec.blade_names.get(&2), Some(&"y".to_string()));
        assert_eq!(spec.blade_names.get(&3), Some(&"xy".to_string()));
    }

    #[test]
    fn reject_empty_signature() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            "#,
        );

        assert!(matches!(result, Err(ParseError::EmptySignature)));
    }

    #[test]
    fn reject_dimension_too_large() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3", "e4", "e5", "e6", "e7"]
            "#,
        );

        assert!(matches!(result, Err(ParseError::DimensionTooLarge(7))));
    }

    #[test]
    fn reject_duplicate_basis_name() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e1"]
            "#,
        );

        assert!(matches!(result, Err(ParseError::DuplicateBasisName(_))));
    }

    #[test]
    fn reject_invalid_grade() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Bad]
            grades = [5]
            "#,
        );

        assert!(matches!(
            result,
            Err(ParseError::InvalidGrade {
                grade: 5,
                max: 2,
                ..
            })
        ));
    }

    #[test]
    fn sparse_type_with_subset_of_blades() {
        // With field_map, providing fewer blades than the grade has is valid
        // and marks the type as sparse
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "y", blade = "e2" }
            ]
            "#,
        )
        .unwrap();

        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();
        assert!(
            vector.is_sparse,
            "Type with subset of blades should be sparse"
        );
        assert_eq!(vector.fields.len(), 2);
    }

    #[test]
    fn reject_duplicate_field_name() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "x", blade = "e2" }
            ]
            "#,
        );

        assert!(matches!(result, Err(ParseError::DuplicateFieldName { .. })));
    }

    #[test]
    fn reject_self_alias() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Rotor]
            grades = [0, 2]
            alias_of = "Rotor"
            "#,
        );

        assert!(matches!(result, Err(ParseError::SelfAlias { .. })));
    }

    #[test]
    fn parse_pga_signature() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "pga3"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]
            "#,
        )
        .unwrap();

        assert_eq!(spec.signature.p, 3);
        assert_eq!(spec.signature.q, 0);
        assert_eq!(spec.signature.r, 1);
        assert_eq!(spec.signature.dim(), 4);
    }

    #[test]
    fn parse_cga_signature() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "cga3"

            [signature]
            positive = ["e1", "e2", "e3", "ep"]
            negative = ["em"]
            "#,
        )
        .unwrap();

        assert_eq!(spec.signature.p, 4);
        assert_eq!(spec.signature.q, 1);
        assert_eq!(spec.signature.r, 0);
        assert_eq!(spec.signature.dim(), 5);
    }

    #[test]
    fn products_are_inferred_euclidean3() {
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();

        // All product types should be inferred
        assert!(
            !spec.products.geometric.is_empty(),
            "Geometric products should be inferred"
        );
        assert!(
            !spec.products.wedge.is_empty(),
            "Wedge products should be inferred"
        );
        assert!(
            !spec.products.left_contraction.is_empty(),
            "Left contraction products should be inferred"
        );
    }

    #[test]
    fn blade_indices_are_canonical_euclidean2() {
        let spec = parse_spec(include_str!("../../../../algebras/euclidean2.toml")).unwrap();

        for ty in &spec.types {
            assert!(
                super::validate_canonical_field_order(ty),
                "Type {} in euclidean2.toml has non-canonical field ordering:\n{:?}",
                ty.name,
                ty.fields
            );
        }
    }

    #[test]
    fn blade_indices_are_canonical_euclidean3() {
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();

        for ty in &spec.types {
            assert!(
                super::validate_canonical_field_order(ty),
                "Type {} in euclidean3.toml has non-canonical field ordering:\n{:?}",
                ty.name,
                ty.fields
            );
        }
    }

    #[test]
    fn validate_canonical_order_function() {
        // Test the validation function with synthetic data
        use super::super::ir::FieldSpec;

        // Valid canonical ordering for grade 2 in 3D
        let valid_type = super::super::ir::TypeSpec {
            name: "Bivector".to_string(),
            grades: vec![2],
            description: None,
            fields: vec![
                FieldSpec {
                    name: "xy".to_string(),
                    blade_index: 3,
                    grade: 2,
                    sign: 1,
                },
                FieldSpec {
                    name: "xz".to_string(),
                    blade_index: 5,
                    grade: 2,
                    sign: 1,
                },
                FieldSpec {
                    name: "yz".to_string(),
                    blade_index: 6,
                    grade: 2,
                    sign: 1,
                },
            ],
            alias_of: None,
            versor: None,
            is_sparse: false,
            inverse_sandwich_targets: vec![],
        };
        assert!(super::validate_canonical_field_order(&valid_type));

        // Invalid: wrong order within grade
        let invalid_type = super::super::ir::TypeSpec {
            name: "Bivector".to_string(),
            grades: vec![2],
            description: None,
            fields: vec![
                FieldSpec {
                    name: "yz".to_string(),
                    blade_index: 6,
                    grade: 2,
                    sign: 1,
                }, // Wrong!
                FieldSpec {
                    name: "xz".to_string(),
                    blade_index: 5,
                    grade: 2,
                    sign: 1,
                }, // Wrong!
                FieldSpec {
                    name: "xy".to_string(),
                    blade_index: 3,
                    grade: 2,
                    sign: 1,
                }, // Wrong!
            ],
            alias_of: None,
            versor: None,
            is_sparse: false,
            inverse_sandwich_targets: vec![],
        };
        assert!(!super::validate_canonical_field_order(&invalid_type));

        // Valid multi-grade type
        let valid_rotor = super::super::ir::TypeSpec {
            name: "Rotor".to_string(),
            grades: vec![0, 2],
            description: None,
            fields: vec![
                FieldSpec {
                    name: "s".to_string(),
                    blade_index: 0,
                    grade: 0,
                    sign: 1,
                },
                FieldSpec {
                    name: "xy".to_string(),
                    blade_index: 3,
                    grade: 2,
                    sign: 1,
                },
                FieldSpec {
                    name: "xz".to_string(),
                    blade_index: 5,
                    grade: 2,
                    sign: 1,
                },
                FieldSpec {
                    name: "yz".to_string(),
                    blade_index: 6,
                    grade: 2,
                    sign: 1,
                },
            ],
            alias_of: None,
            versor: None,
            is_sparse: false,
            inverse_sandwich_targets: vec![],
        };
        assert!(super::validate_canonical_field_order(&valid_rotor));
    }

    #[test]
    fn parse_norm_default() {
        use super::super::ir::InvolutionKind;

        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]
            "#,
        )
        .unwrap();

        // Default should be Reverse
        assert_eq!(spec.norm.primary_involution, InvolutionKind::Reverse);
    }

    #[test]
    fn parse_norm_reverse_explicit() {
        use super::super::ir::InvolutionKind;

        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [norm]
            primary_involution = "reverse"
            "#,
        )
        .unwrap();

        assert_eq!(spec.norm.primary_involution, InvolutionKind::Reverse);
    }

    #[test]
    fn parse_norm_grade_involution() {
        use super::super::ir::InvolutionKind;

        let spec = parse_spec(
            r#"
            [algebra]
            name = "hyperbolic"

            [signature]
            positive = ["e1"]

            [norm]
            primary_involution = "grade_involution"
            "#,
        )
        .unwrap();

        assert_eq!(
            spec.norm.primary_involution,
            InvolutionKind::GradeInvolution
        );
    }

    #[test]
    fn parse_norm_clifford_conjugate() {
        use super::super::ir::InvolutionKind;

        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1"]

            [norm]
            primary_involution = "clifford_conjugate"
            "#,
        )
        .unwrap();

        assert_eq!(
            spec.norm.primary_involution,
            InvolutionKind::CliffordConjugate
        );
    }

    #[test]
    fn reject_invalid_norm_involution() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1"]

            [norm]
            primary_involution = "invalid_involution"
            "#,
        );

        assert!(matches!(result, Err(ParseError::InvalidValue { .. })));
    }

    #[test]
    fn parse_sparse_type() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "conformal3"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3", "e4"]
            negative = ["e5"]

            [types.Line]
            grades = [3]
            description = "Line (circle through infinity)"
            field_map = [
                { name = "vx", blade = "e145" },
                { name = "vy", blade = "e245" },
                { name = "vz", blade = "e345" },
                { name = "mx", blade = "e235" },
                { name = "my", blade = "e135" },
                { name = "mz", blade = "e125" }
            ]
            "#,
        )
        .unwrap();

        let line = spec.types.iter().find(|t| t.name == "Line").unwrap();
        assert!(line.is_sparse, "Line should be marked as sparse");
        assert_eq!(line.grades, vec![3]);
        assert_eq!(line.fields.len(), 6, "Line should have 6 fields");

        // Verify blade indices were parsed correctly
        // e145 -> bits 0, 3, 4 -> 0b11001 = 25
        assert_eq!(
            line.fields[0].blade_index, 25,
            "e145 should be blade index 25"
        );
        // e235 -> bits 1, 2, 4 -> 0b10110 = 22
        assert_eq!(
            line.fields[3].blade_index, 22,
            "e235 should be blade index 22"
        );
        // e125 -> bits 0, 1, 4 -> 0b10011 = 19
        assert_eq!(
            line.fields[5].blade_index, 19,
            "e125 should be blade index 19"
        );
    }

    #[test]
    fn reject_field_map_grade_mismatch() {
        // Test that a blade grade mismatch in field_map is rejected
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Bad]
            grades = [2]
            field_map = [
                { name = "xyz", blade = "e123" }
            ]
            "#,
        );

        assert!(matches!(
            result,
            Err(ParseError::FieldMapGradeMismatch { blade_grade: 3, .. })
        ));
    }

    #[test]
    fn completeness_check_passes_for_complete_algebra() {
        // Euclidean 2D has all required types
        let result = parse_spec(
            r#"
            [algebra]
            name = "euclidean2"
            complete = true

            [signature]
            positive = ["e1", "e2"]

            [types.Scalar]
            grades = [0]
            field_map = [{ name = "s", blade = "s" }]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "y", blade = "e2" }
            ]

            [types.Bivector]
            grades = [2]
            field_map = [{ name = "b", blade = "e12" }]

            [types.Rotor]
            grades = [0, 2]
            field_map = [
                { name = "s", blade = "s" },
                { name = "xy", blade = "e12" }
            ]
            "#,
        );

        assert!(
            result.is_ok(),
            "Complete euclidean2 algebra should pass: {:?}",
            result
        );
    }

    #[test]
    fn completeness_check_fails_for_incomplete_algebra() {
        // Missing Rotor type means Vector × Vector has no output
        let result = parse_spec(
            r#"
            [algebra]
            name = "incomplete"
            complete = true

            [signature]
            positive = ["e1", "e2"]

            [types.Scalar]
            grades = [0]
            field_map = [{ name = "s", blade = "s" }]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "y", blade = "e2" }
            ]

            [types.Bivector]
            grades = [2]
            field_map = [{ name = "b", blade = "e12" }]
            "#,
        );

        assert!(matches!(result, Err(ParseError::IncompleteAlgebra { .. })));
        if let Err(ParseError::IncompleteAlgebra { count, details, .. }) = result {
            assert!(count > 0, "Should have missing products");
            assert!(details.contains("[0, 2]"), "Should mention grades [0, 2]");
        }
    }

    #[test]
    fn completeness_check_enabled_by_default() {
        // Same incomplete algebra but without complete = true
        // Since default is complete = true, this should fail
        let result = parse_spec(
            r#"
            [algebra]
            name = "incomplete"

            [signature]
            positive = ["e1", "e2"]

            [types.Scalar]
            grades = [0]
            field_map = [{ name = "s", blade = "s" }]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "y", blade = "e2" }
            ]

            [types.Bivector]
            grades = [2]
            field_map = [{ name = "b", blade = "e12" }]
            "#,
        );

        // Should fail because complete defaults to true
        assert!(matches!(result, Err(ParseError::IncompleteAlgebra { .. })));
    }

    #[test]
    fn completeness_check_can_be_disabled() {
        // Same incomplete algebra with complete = false
        let result = parse_spec(
            r#"
            [algebra]
            name = "incomplete"
            complete = false

            [signature]
            positive = ["e1", "e2"]

            [types.Scalar]
            grades = [0]
            field_map = [{ name = "s", blade = "s" }]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "y", blade = "e2" }
            ]

            [types.Bivector]
            grades = [2]
            field_map = [{ name = "b", blade = "e12" }]
            "#,
        );

        // Should succeed because complete = false explicitly disables the check
        assert!(
            result.is_ok(),
            "Incomplete algebra with complete=false should succeed: {:?}",
            result
        );
    }

    #[test]
    fn require_field_map_for_types() {
        // Test that types without field_map are rejected
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"
            complete = false

            [signature]
            positive = ["e1", "e2"]

            [types.Vector]
            grades = [1]
            "#,
        );

        assert!(matches!(result, Err(ParseError::MissingFieldMap { .. })));
    }
}
