//! Parser for algebra specifications.
//!
//! Converts raw TOML structures to validated IR types.

use std::collections::{HashMap, HashSet};

use crate::algebra::binomial;

use super::error::ParseError;
use super::ir::{
    AlgebraSpec, BasisVector, ConstraintKind, ConstraintSpec, FieldSpec, GenerationOptions,
    ProductEntry, ProductsSpec, SignatureSpec, TypeSpec,
};
use super::raw::{RawAlgebraSpec, RawConstraint, RawProducts, RawSignature, RawTypeSpec};

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

    // Build products
    let products = parse_products(&raw.products, &types);

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
        build_fields_from_names(&raw.fields, &raw.grades, dim, name)?
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

    // Parse constraints
    let constraints = parse_constraints(name, &raw.constraints)?;

    // Check self-alias
    if let Some(alias) = &raw.alias_of {
        if alias == name {
            return Err(ParseError::SelfAlias {
                type_name: name.to_string(),
            });
        }
    }

    Ok(TypeSpec {
        name: name.to_string(),
        grades: raw.grades.clone(),
        description: raw.description.clone(),
        fields,
        alias_of: raw.alias_of.clone(),
        constraints,
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
fn build_fields_from_names(
    names: &[String],
    grades: &[usize],
    dim: usize,
    type_name: &str,
) -> Result<Vec<FieldSpec>, ParseError> {
    // Collect blade indices for these grades in order
    let mut blade_indices = Vec::new();
    for &grade in grades {
        for blade_index in 0usize..(1 << dim) {
            if blade_index.count_ones() as usize == grade {
                blade_indices.push((blade_index, grade));
            }
        }
    }

    if names.len() != blade_indices.len() {
        return Err(ParseError::FieldCountMismatch {
            type_name: type_name.to_string(),
            expected: blade_indices.len(),
            got: names.len(),
        });
    }

    Ok(names
        .iter()
        .zip(blade_indices.iter())
        .map(|(name, &(blade_index, grade))| FieldSpec {
            name: name.clone(),
            blade_index,
            grade,
        })
        .collect())
}

/// Parses constraint specifications.
fn parse_constraints(
    type_name: &str,
    raw: &HashMap<String, RawConstraint>,
) -> Result<Vec<ConstraintSpec>, ParseError> {
    let mut constraints = Vec::new();

    for (kind_str, raw_constraint) in raw {
        let kind = parse_constraint_kind(type_name, kind_str)?;
        let wrapper_name = format!("{}{}", kind.wrapper_prefix(), type_name);

        constraints.push(ConstraintSpec {
            kind,
            wrapper_name,
            condition: raw_constraint.condition.clone(),
        });
    }

    Ok(constraints)
}

/// Parses a constraint kind string.
fn parse_constraint_kind(type_name: &str, kind: &str) -> Result<ConstraintKind, ParseError> {
    match kind {
        "unit" => Ok(ConstraintKind::Unit),
        "nonzero" => Ok(ConstraintKind::NonZero),
        "normalized" => Ok(ConstraintKind::Normalized),
        "null" => Ok(ConstraintKind::Null),
        "ideal" => Ok(ConstraintKind::Ideal),
        _ => Err(ParseError::UnknownConstraintKind {
            type_name: type_name.to_string(),
            kind: kind.to_string(),
        }),
    }
}

/// Parses the products section.
fn parse_products(raw: &RawProducts, types: &[TypeSpec]) -> ProductsSpec {
    // Collect all type names including constrained wrapper names
    let mut all_types: HashSet<&str> = types.iter().map(|t| t.name.as_str()).collect();
    for ty in types {
        for constraint in &ty.constraints {
            all_types.insert(&constraint.wrapper_name);
        }
    }

    // Helper to check if a type is constrained
    let is_constrained = |name: &str| -> bool {
        types
            .iter()
            .flat_map(|t| &t.constraints)
            .any(|c| c.wrapper_name == name)
    };

    ProductsSpec {
        geometric: parse_product_entries(&raw.geometric, is_constrained),
        outer: parse_product_entries(&raw.outer, is_constrained),
        left_contraction: parse_product_entries(&raw.left_contraction, is_constrained),
        right_contraction: parse_product_entries(&raw.right_contraction, is_constrained),
        regressive: parse_product_entries(&raw.regressive, is_constrained),
        scalar: parse_product_entries(&raw.scalar, is_constrained),
    }
}

/// Parses product entries from a map of "Lhs_Rhs" -> "Output".
fn parse_product_entries(
    entries: &HashMap<String, String>,
    is_constrained: impl Fn(&str) -> bool,
) -> Vec<ProductEntry> {
    entries
        .iter()
        .filter_map(|(key, output)| {
            // Parse "Lhs_Rhs" format
            let parts: Vec<&str> = key.split('_').collect();
            if parts.len() != 2 {
                return None;
            }
            Some(ProductEntry {
                lhs: parts[0].to_string(),
                rhs: parts[1].to_string(),
                output: output.clone(),
                output_constrained: is_constrained(output),
            })
        })
        .collect()
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

            [types.Vector]
            grades = [1]
            fields = ["x", "y"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]
            "#,
        )
        .unwrap();

        assert_eq!(spec.types.len(), 2);

        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();
        assert_eq!(vector.grades, vec![1]);
        assert_eq!(vector.fields.len(), 2);

        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.grades, vec![0, 2]);
        assert_eq!(rotor.fields.len(), 2);
    }

    #[test]
    fn parse_with_constraints() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]

            [types.Rotor.constraints.unit]
            "#,
        )
        .unwrap();

        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.constraints.len(), 1);

        let unit = &rotor.constraints[0];
        assert_eq!(unit.kind, ConstraintKind::Unit);
        assert_eq!(unit.wrapper_name, "UnitRotor");
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
}
