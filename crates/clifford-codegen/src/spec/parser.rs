//! Parser for algebra specifications.
//!
//! Converts raw TOML structures to validated IR types.

use std::collections::{HashMap, HashSet};

use crate::algebra::{Algebra, binomial};
use crate::discovery::{ProductType, infer_all_products};

use super::error::ParseError;
use super::ir::{
    AlgebraSpec, BasisVector, FieldSpec, GenerationOptions, ProductEntry, ProductsSpec,
    SignatureSpec, TypeSpec, VersorSpec,
};
use super::raw::{RawAlgebraSpec, RawSignature, RawTypeSpec};

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
///
/// [signature]
/// positive = ["e1", "e2"]
///
/// [types.Vector]
/// grades = [1]
/// "#).unwrap();
///
/// assert_eq!(spec.name, "euclidean2");
/// assert_eq!(spec.signature.dim(), 2);
/// ```
pub fn parse_spec(toml_content: &str) -> Result<AlgebraSpec, ParseError> {
    let raw: RawAlgebraSpec = toml::from_str(toml_content)?;

    // Build signature
    let signature = parse_signature(&raw.signature)?;

    // Build blade name map
    let blade_names = parse_blade_names(&raw.blades, &signature)?;

    // Build types
    let types = parse_types(&raw.types, &signature, &blade_names)?;

    // Auto-infer products from types (products section in TOML is ignored)
    let products = infer_products_from_types(&types, &signature);

    // Validate the complete specification
    validate_spec(&types)?;

    Ok(AlgebraSpec {
        name: raw.algebra.name,
        module_path: raw.algebra.module_path,
        description: raw.algebra.description,
        signature,
        blade_names,
        types,
        products,
        options: GenerationOptions {
            generate_serde: raw.options.generate_serde,
            generate_arbitrary: raw.options.generate_arbitrary,
            generate_tests: raw.options.generate_tests,
        },
    })
}

/// Parses the signature section.
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

    // Build basis vectors
    let mut basis = Vec::with_capacity(dim);
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

    Ok(SignatureSpec { basis, p, q, r })
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

/// Parses a blade name like "e12" or "e123" into its index.
fn parse_blade_index(name: &str, dim: usize) -> Result<usize, ParseError> {
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

    // Each digit is a 1-based basis vector index
    let mut index = 0usize;
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

        // Convert to 0-based and set bit
        let bit = digit - 1;
        if (index >> bit) & 1 == 1 {
            // Duplicate index in blade name
            return Err(ParseError::InvalidBladeName {
                name: name.to_string(),
            });
        }
        index |= 1 << bit;
    }

    Ok(index)
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
    sig: &SignatureSpec,
    blade_names: &HashMap<usize, String>,
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

    // Compute expected field count
    let expected_fields: usize = raw.grades.iter().map(|&g| binomial(dim, g)).sum();

    // Build fields
    let fields = if raw.fields.is_empty() {
        // Generate default fields from blade names
        generate_default_fields(&raw.grades, dim, sig, blade_names)
    } else {
        // Use provided fields
        if raw.fields.len() != expected_fields {
            return Err(ParseError::FieldCountMismatch {
                type_name: name.to_string(),
                expected: expected_fields,
                got: raw.fields.len(),
            });
        }
        build_fields_from_names(&raw.fields, &raw.grades, dim, name, blade_names, sig)?
    };

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

    // Constraints are now inferred automatically during code generation (Phase 4)
    // No user-defined constraints are read from TOML
    let constraints = Vec::new();

    // Parse versor information
    let versor = if raw.versor {
        Some(VersorSpec {
            // is_unit will be determined by inferred constraints in Phase 4
            is_unit: false,
            sandwich_targets: raw
                .sandwich
                .as_ref()
                .map(|s| s.targets.clone())
                .unwrap_or_default(),
        })
    } else {
        None
    };

    Ok(TypeSpec {
        name: name.to_string(),
        grades: raw.grades.clone(),
        description: raw.description.clone(),
        fields,
        alias_of: raw.alias_of.clone(),
        constraints,
        versor,
    })
}

/// Generates default field specs from blade indices.
fn generate_default_fields(
    grades: &[usize],
    dim: usize,
    sig: &SignatureSpec,
    blade_names: &HashMap<usize, String>,
) -> Vec<FieldSpec> {
    let mut fields = Vec::new();

    for &grade in grades {
        for blade_index in 0usize..(1 << dim) {
            if blade_index.count_ones() as usize == grade {
                let name = blade_names
                    .get(&blade_index)
                    .cloned()
                    .unwrap_or_else(|| default_blade_name(blade_index, sig));

                fields.push(FieldSpec {
                    name,
                    blade_index,
                    grade,
                });
            }
        }
    }

    fields
}

/// Generates a default blade name from basis vector names.
fn default_blade_name(blade_index: usize, sig: &SignatureSpec) -> String {
    if blade_index == 0 {
        return "s".to_string();
    }

    let mut name = String::new();
    for basis in &sig.basis {
        if (blade_index >> basis.index) & 1 == 1 {
            name.push_str(&basis.name);
        }
    }
    name
}

/// Builds field specs from provided field names.
///
/// # Canonical Ordering Requirement
///
/// Fields in TOML must be listed in **canonical blade order**:
/// 1. First by grade (ascending)
/// 2. Within each grade, by blade index (ascending)
///
/// Blade indices follow the bitmask convention where bit `i` indicates
/// the presence of basis vector `eᵢ`. For example, in 3D:
///
/// | Index | Binary | Blade | Grade |
/// |-------|--------|-------|-------|
/// | 0     | 000    | 1     | 0     |
/// | 1     | 001    | e₁    | 1     |
/// | 2     | 010    | e₂    | 1     |
/// | 3     | 011    | e₁₂   | 2     |
/// | 4     | 100    | e₃    | 1     |
/// | 5     | 101    | e₁₃   | 2     |
/// | 6     | 110    | e₂₃   | 2     |
/// | 7     | 111    | e₁₂₃  | 3     |
///
/// Builds field specs from provided field names.
///
/// This function tries to compute the correct blade_index for each field:
/// 1. First, looks up the field name in the blade_names mapping (from [blades] section)
/// 2. If not found, looks up by default blade name (e.g., "s", "e12", "e123")
/// 3. If still not found, falls back to assuming canonical blade order
///
/// The fallback to canonical order maintains backward compatibility with
/// algebras that use custom field names without explicit blade mappings.
fn build_fields_from_names(
    names: &[String],
    grades: &[usize],
    dim: usize,
    type_name: &str,
    blade_names: &HashMap<usize, String>,
    sig: &SignatureSpec,
) -> Result<Vec<FieldSpec>, ParseError> {
    // Build inverted mapping: field_name -> blade_index
    let name_to_index: HashMap<String, usize> = blade_names
        .iter()
        .map(|(&idx, name)| (name.clone(), idx))
        .collect();

    // Also add default blade names (e.g., "s" for scalar, "e123" for trivector)
    let mut name_to_index_with_defaults: HashMap<String, usize> = name_to_index;
    for blade_index in 0usize..(1 << dim) {
        let default_name = default_blade_name(blade_index, sig);
        name_to_index_with_defaults
            .entry(default_name)
            .or_insert(blade_index);
    }

    // Collect expected blade indices in canonical order for fallback
    let mut blade_indices_canonical = Vec::new();
    for &grade in grades {
        for blade_index in 0usize..(1 << dim) {
            if blade_index.count_ones() as usize == grade {
                blade_indices_canonical.push((blade_index, grade));
            }
        }
    }

    if names.len() != blade_indices_canonical.len() {
        return Err(ParseError::FieldCountMismatch {
            type_name: type_name.to_string(),
            expected: blade_indices_canonical.len(),
            got: names.len(),
        });
    }

    // Build field specs - try to look up name, otherwise use canonical order
    let mut fields = Vec::with_capacity(names.len());
    for (i, name) in names.iter().enumerate() {
        let (blade_index, grade) = if let Some(&idx) = name_to_index_with_defaults.get(name) {
            // Found in mapping - use the correct blade_index
            let grade = idx.count_ones() as usize;
            (idx, grade)
        } else {
            // Not found - fall back to canonical order (maintains backward compatibility)
            blade_indices_canonical[i]
        };

        fields.push(FieldSpec {
            name: name.clone(),
            blade_index,
            grade,
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
fn infer_products_from_types(types: &[TypeSpec], signature: &SignatureSpec) -> ProductsSpec {
    // Build algebra for product computation
    let algebra = Algebra::new(signature.p, signature.q, signature.r);

    // Build entity list for inference
    let entities: Vec<(String, Vec<usize>)> = types
        .iter()
        .filter(|t| t.alias_of.is_none())
        .map(|t| (t.name.clone(), t.grades.clone()))
        .collect();

    // Infer all standard product types
    let geometric_table = infer_all_products(&entities, ProductType::Geometric, &algebra);
    let exterior_table = infer_all_products(&entities, ProductType::Exterior, &algebra);
    let inner_table = infer_all_products(&entities, ProductType::Inner, &algebra);
    let left_contraction_table =
        infer_all_products(&entities, ProductType::LeftContraction, &algebra);
    let right_contraction_table =
        infer_all_products(&entities, ProductType::RightContraction, &algebra);
    let regressive_table = infer_all_products(&entities, ProductType::Regressive, &algebra);
    let scalar_table = infer_all_products(&entities, ProductType::Scalar, &algebra);
    let antigeometric_table = infer_all_products(&entities, ProductType::Antigeometric, &algebra);
    let antiscalar_table = infer_all_products(&entities, ProductType::Antiscalar, &algebra);

    // Convert inferred products to ProductEntry format
    // Skip products that don't have matching entity types
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
                    output_constrained: false, // No wrapper types
                }
            })
            .collect()
    };

    // Note: Interior product is the symmetric inner product.
    // Left contraction is more commonly used in GA, but we provide both.

    ProductsSpec {
        geometric: convert_entries(geometric_table),
        exterior: convert_entries(exterior_table),
        interior: convert_entries(inner_table), // Symmetric inner product
        left_contraction: convert_entries(left_contraction_table),
        right_contraction: convert_entries(right_contraction_table),
        regressive: convert_entries(regressive_table),
        scalar: convert_entries(scalar_table),
        antigeometric: convert_entries(antigeometric_table),
        antiscalar: convert_entries(antiscalar_table),
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

            [types.Vector]
            grades = [1]
            fields = ["x", "y"]

            [types.Bivector]
            grades = [2]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]
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
    fn reject_field_count_mismatch() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Vector]
            grades = [1]
            fields = ["x", "y"]
            "#,
        );

        assert!(matches!(
            result,
            Err(ParseError::FieldCountMismatch {
                expected: 3,
                got: 2,
                ..
            })
        ));
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
            fields = ["x", "x"]
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
            !spec.products.exterior.is_empty(),
            "Exterior products should be inferred"
        );
        assert!(
            !spec.products.left_contraction.is_empty(),
            "Left contraction products should be inferred"
        );
        assert!(
            !spec.products.interior.is_empty(),
            "Interior products should be inferred"
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
                },
                FieldSpec {
                    name: "xz".to_string(),
                    blade_index: 5,
                    grade: 2,
                },
                FieldSpec {
                    name: "yz".to_string(),
                    blade_index: 6,
                    grade: 2,
                },
            ],
            alias_of: None,
            constraints: vec![],
            versor: None,
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
                }, // Wrong!
                FieldSpec {
                    name: "xz".to_string(),
                    blade_index: 5,
                    grade: 2,
                }, // Wrong!
                FieldSpec {
                    name: "xy".to_string(),
                    blade_index: 3,
                    grade: 2,
                }, // Wrong!
            ],
            alias_of: None,
            constraints: vec![],
            versor: None,
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
                },
                FieldSpec {
                    name: "xy".to_string(),
                    blade_index: 3,
                    grade: 2,
                },
                FieldSpec {
                    name: "xz".to_string(),
                    blade_index: 5,
                    grade: 2,
                },
                FieldSpec {
                    name: "yz".to_string(),
                    blade_index: 6,
                    grade: 2,
                },
            ],
            alias_of: None,
            constraints: vec![],
            versor: None,
        };
        assert!(super::validate_canonical_field_order(&valid_rotor));
    }
}
