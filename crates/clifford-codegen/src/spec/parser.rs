//! Parser for algebra specifications.
//!
//! Converts raw TOML structures to validated IR types.

use std::collections::{HashMap, HashSet};

use crate::algebra::{Algebra, binomial};
use crate::discovery::{ProductType, infer_all_products};

use super::error::ParseError;
use super::ir::{
    AlgebraSpec, BasisVector, FieldSpec, GenerationOptions, ProductEntry, ProductsSpec,
    SignConvention, SignatureSpec, TypeSpec, UserConstraint, normalize_constraint_expr,
};
use super::raw::{RawAlgebraSpec, RawSignature, RawTypeSpec, RawUserConstraint};

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

    // Check self-alias
    if let Some(alias) = &raw.alias_of {
        if alias == name {
            return Err(ParseError::SelfAlias {
                type_name: name.to_string(),
            });
        }
    }

    // Validate geometric_solve_for field if present
    if let Some(ref solve_for) = raw.geometric_solve_for {
        let valid_field = fields.iter().any(|f| &f.name == solve_for);
        if !valid_field {
            return Err(ParseError::InvalidSolveFor {
                type_name: name.to_string(),
                field: solve_for.clone(),
            });
        }
        // geometric_solve_for requires geometric_constraint
        if raw.geometric_constraint.is_none() {
            return Err(ParseError::SolveForWithoutConstraint {
                type_name: name.to_string(),
            });
        }
    }

    // Validate antiproduct_solve_for field if present
    if let Some(ref solve_for) = raw.antiproduct_solve_for {
        let valid_field = fields.iter().any(|f| &f.name == solve_for);
        if !valid_field {
            return Err(ParseError::InvalidSolveFor {
                type_name: name.to_string(),
                field: solve_for.clone(),
            });
        }
        // antiproduct_solve_for requires antiproduct_constraint
        if raw.antiproduct_constraint.is_none() {
            return Err(ParseError::SolveForWithoutConstraint {
                type_name: name.to_string(),
            });
        }
    }

    // Check: if constraints are independent, both solve_for fields must differ
    let constraints_independent = match (&raw.geometric_constraint, &raw.antiproduct_constraint) {
        (Some(gc), Some(ac)) => normalize_constraint_expr(gc) != normalize_constraint_expr(ac),
        _ => false,
    };

    if constraints_independent {
        // Independent constraints require two different solve_for fields
        match (&raw.geometric_solve_for, &raw.antiproduct_solve_for) {
            (Some(gsf), Some(asf)) if gsf == asf => {
                return Err(ParseError::IndependentConstraintsSameSolveFor {
                    type_name: name.to_string(),
                });
            }
            (None, _) | (_, None) => {
                return Err(ParseError::IndependentConstraintsMissingSolveFor {
                    type_name: name.to_string(),
                });
            }
            _ => {} // Valid: different solve_for fields
        }
    }

    // Parse and validate user constraints
    let constraints = parse_user_constraints(&raw.constraints, &fields, name)?;

    Ok(TypeSpec {
        name: name.to_string(),
        grades: raw.grades.clone(),
        description: raw.description.clone(),
        fields,
        alias_of: raw.alias_of.clone(),
        geometric_constraint: raw.geometric_constraint.clone(),
        antiproduct_constraint: raw.antiproduct_constraint.clone(),
        geometric_solve_for: raw.geometric_solve_for.clone(),
        antiproduct_solve_for: raw.antiproduct_solve_for.clone(),
        constraints,
    })
}

/// Parses and validates user-defined constraints.
fn parse_user_constraints(
    raw_constraints: &[RawUserConstraint],
    fields: &[FieldSpec],
    type_name: &str,
) -> Result<Vec<UserConstraint>, ParseError> {
    let mut constraints = Vec::with_capacity(raw_constraints.len());

    for raw in raw_constraints {
        // Validate solve_for field if present
        if let Some(ref solve_for) = raw.solve_for {
            let valid_field = fields.iter().any(|f| &f.name == solve_for);
            if !valid_field {
                return Err(ParseError::InvalidSolveFor {
                    type_name: type_name.to_string(),
                    field: solve_for.clone(),
                });
            }
        }

        // Parse sign convention
        let sign = match raw.sign.to_lowercase().as_str() {
            "positive" => SignConvention::Positive,
            "negative" => SignConvention::Negative,
            _ => {
                return Err(ParseError::InvalidSignConvention {
                    type_name: type_name.to_string(),
                    constraint_name: raw.name.clone(),
                    sign: raw.sign.clone(),
                });
            }
        };

        // Detect if constraint has domain restrictions
        // A constraint has domain restrictions if solving requires sqrt
        // (detected by quadratic terms in the expression)
        let has_domain_restriction = detect_domain_restriction(&raw.expression, &raw.solve_for);

        constraints.push(UserConstraint {
            name: raw.name.clone(),
            description: raw.description.clone(),
            expression: raw.expression.clone(),
            solve_for: raw.solve_for.clone(),
            sign,
            enforce: raw.enforce.clone(),
            has_domain_restriction,
        });
    }

    Ok(constraints)
}

/// Detects if a constraint expression has domain restrictions when solved for a variable.
///
/// Returns true if:
/// - The constraint involves squared terms of the solve_for variable
/// - Solving would require taking a square root (potentially undefined for some inputs)
fn detect_domain_restriction(expression: &str, solve_for: &Option<String>) -> bool {
    let Some(var) = solve_for else {
        return false;
    };

    // Check if the variable appears squared in the expression
    // Patterns: "var*var", "var^2", "var**2"
    let squared_patterns = [
        format!("{}*{}", var, var),
        format!("{}^2", var),
        format!("{}**2", var),
    ];

    for pattern in &squared_patterns {
        if expression.contains(pattern) {
            return true;
        }
    }

    false
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

/// Infers products automatically from types.
///
/// Products are always auto-inferred from the defined types.
/// The TOML no longer needs (or supports) a products section.
fn infer_products_from_types(types: &[TypeSpec], signature: &SignatureSpec) -> ProductsSpec {
    // Build algebra for product computation
    let algebra = Algebra::new(signature.p, signature.q, signature.r);

    // Build entity list for inference
    let entities: Vec<(String, Vec<usize>)> = types
        .iter()
        .filter(|t| t.alias_of.is_none())
        .map(|t| (t.name.clone(), t.grades.clone()))
        .collect();

    // Infer products for geometric and outer only
    // Left contraction and other products are only generated if explicitly specified
    let geometric_table = infer_all_products(&entities, ProductType::Geometric, &algebra);
    let outer_table = infer_all_products(&entities, ProductType::Outer, &algebra);

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

    ProductsSpec {
        geometric: convert_entries(geometric_table),
        outer: convert_entries(outer_table),
        left_contraction: vec![],  // Only generated if explicitly specified
        right_contraction: vec![], // Not commonly used
        regressive: vec![],        // Not commonly used
        scalar: vec![],            // Can be derived from geometric
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
    fn parse_with_constraint() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Scalar]
            grades = [0]

            [types.Vector]
            grades = [1]

            [types.Bivector]
            grades = [2]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]
            geometric_constraint = "s * s + xy * xy = 1"
            "#,
        )
        .unwrap();

        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(
            rotor.geometric_constraint,
            Some("s * s + xy * xy = 1".to_string())
        );
    }

    #[test]
    fn parse_with_geometric_solve_for() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Motor]
            grades = [0, 2, 4]
            fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
            geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
            geometric_solve_for = "e0123"
            "#,
        )
        .unwrap();

        let motor = spec.types.iter().find(|t| t.name == "Motor").unwrap();
        assert_eq!(motor.geometric_solve_for, Some("e0123".to_string()));
    }

    #[test]
    fn reject_invalid_geometric_solve_for() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]
            geometric_constraint = "s * s + xy * xy = 1"
            geometric_solve_for = "nonexistent"
            "#,
        );

        assert!(matches!(result, Err(ParseError::InvalidSolveFor { .. })));
    }

    #[test]
    fn reject_geometric_solve_for_without_constraint() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]
            geometric_solve_for = "xy"
            "#,
        );

        assert!(matches!(
            result,
            Err(ParseError::SolveForWithoutConstraint { .. })
        ));
    }

    #[test]
    fn parse_independent_constraints() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.TestType]
            grades = [0, 2, 4]
            fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
            geometric_constraint = "s + e0123 = 0"
            antiproduct_constraint = "e12 + e03 = 0"
            geometric_solve_for = "e0123"
            antiproduct_solve_for = "e03"
            "#,
        )
        .unwrap();

        let tt = spec.types.iter().find(|t| t.name == "TestType").unwrap();
        assert_eq!(tt.geometric_solve_for, Some("e0123".to_string()));
        assert_eq!(tt.antiproduct_solve_for, Some("e03".to_string()));
    }

    #[test]
    fn reject_independent_constraints_same_solve_for() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.TestType]
            grades = [0, 2, 4]
            fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
            geometric_constraint = "s + e0123 = 0"
            antiproduct_constraint = "e12 + e03 = 0"
            geometric_solve_for = "e0123"
            antiproduct_solve_for = "e0123"
            "#,
        );

        assert!(matches!(
            result,
            Err(ParseError::IndependentConstraintsSameSolveFor { .. })
        ));
    }

    #[test]
    fn reject_independent_constraints_missing_solve_for() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.TestType]
            grades = [0, 2, 4]
            fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
            geometric_constraint = "s + e0123 = 0"
            antiproduct_constraint = "e12 + e03 = 0"
            geometric_solve_for = "e0123"
            "#,
        );

        assert!(matches!(
            result,
            Err(ParseError::IndependentConstraintsMissingSolveFor { .. })
        ));
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
    fn parse_user_constraints() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Motor]
            grades = [0, 2, 4]
            fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
            geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
            geometric_solve_for = "e0123"

            [[types.Motor.constraints]]
            name = "unit"
            expression = "s*s + e12*e12 + e13*e13 + e23*e23 = 1"
            solve_for = "s"
            sign = "positive"
            "#,
        )
        .unwrap();

        let motor = spec.types.iter().find(|t| t.name == "Motor").unwrap();
        assert_eq!(motor.constraints.len(), 1);
        assert_eq!(motor.constraints[0].name, "unit");
        assert_eq!(motor.constraints[0].solve_for, Some("s".to_string()));
        assert_eq!(motor.constraints[0].sign, SignConvention::Positive);
        assert!(motor.constraints[0].has_domain_restriction);
    }

    #[test]
    fn parse_user_constraint_with_enforce() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]

            [[types.Rotor.constraints]]
            name = "unit"
            description = "Unit rotor constraint"
            expression = "s*s + xy*xy = 1"
            enforce = "normalize"
            "#,
        )
        .unwrap();

        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.constraints.len(), 1);
        assert_eq!(rotor.constraints[0].name, "unit");
        assert_eq!(rotor.constraints[0].solve_for, None);
        assert_eq!(rotor.constraints[0].enforce, Some("normalize".to_string()));
    }

    #[test]
    fn reject_invalid_sign_convention() {
        let result = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy"]

            [[types.Rotor.constraints]]
            name = "unit"
            expression = "s*s + xy*xy = 1"
            solve_for = "s"
            sign = "invalid"
            "#,
        );

        assert!(matches!(
            result,
            Err(ParseError::InvalidSignConvention { .. })
        ));
    }

    #[test]
    fn detect_domain_restriction() {
        // Test the detection function directly
        assert!(super::detect_domain_restriction(
            "s*s + b*b = 1",
            &Some("s".to_string())
        ));
        assert!(super::detect_domain_restriction(
            "x^2 + y = 1",
            &Some("x".to_string())
        ));
        assert!(!super::detect_domain_restriction(
            "2*s*b = 0",
            &Some("s".to_string())
        ));
        assert!(!super::detect_domain_restriction("s*s = 1", &None));
    }
}
