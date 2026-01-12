# PRD-18.4: Infer Constraints from Algebra Structure

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Automatically derive constraint expressions from the algebra's signature and grade structure

## Overview

Instead of user-specified constraints, the codegen tool should infer them:

| Type | Grades | Algebra | Inferred Constraint |
|------|--------|---------|---------------------|
| Rotor | [0, 2] | Euclidean 3D | Unit norm: `s² + xy² + xz² + yz² = 1` |
| Motor | [0, 2, 4] | PGA 3D | Study: `s·e0123 + e23·e01 + e31·e02 + e12·e03 = 0` |
| Line | [2] | PGA 3D | Plücker: `e01·e23 + e02·e31 + e03·e12 = 0` |
| Flector | [1, 3] | PGA 3D | Geometric: `e1·e023 + e2·e031 + e3·e012 + e0·e123 = 0` |

## Deliverables

### Modified Files (clifford-codegen)

#### `src/algebra/constraints.rs`

Add constraint inference functions:

```rust
/// Infers the constraint for a type based on its grade structure.
///
/// Returns None if the type has no constraint (e.g., single-grade types
/// in Euclidean algebras).
pub fn infer_constraint(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
    fields: &[FieldSpec],
) -> Option<InferredConstraint> {
    // Check if type satisfies geometric constraint
    if !satisfies_geometric_constraint(grades, algebra, table) {
        return None;  // No valid constraint exists
    }

    // Derive constraint expression symbolically
    let constraint_expr = derive_constraint_expression(grades, algebra, table, fields)?;

    Some(InferredConstraint {
        name: infer_constraint_name(grades, algebra),
        expression: constraint_expr,
        constraint_type: infer_constraint_type(grades, algebra),
    })
}

/// Determines the constraint name based on algebra type and grades.
fn infer_constraint_name(grades: &[usize], algebra: &Algebra) -> &'static str {
    let dim = algebra.dim();
    let has_degenerate = algebra.signature().zero > 0;

    match (grades, has_degenerate) {
        // PGA Motor: [0, 2, 4] with degenerate basis
        (&[0, 2, 4], true) if dim == 4 => "study",

        // PGA Line: [2] with degenerate basis
        (&[2], true) if dim == 4 => "plucker",

        // PGA Flector: [1, 3] with degenerate basis
        (&[1, 3], true) if dim == 4 => "geometric",

        // Euclidean Rotor: [0, 2] without degenerate basis
        (&[0, 2], false) => "unit",

        _ => "geometric",
    }
}

/// Derives the constraint expression using Symbolica.
fn derive_constraint_expression(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
    fields: &[FieldSpec],
) -> Option<String> {
    // Use existing symbolic computation infrastructure
    // to compute u * rev(u) and extract non-scalar terms

    // For geometric constraint: scalar_part(u * rev(u)) - norm² = 0
    // For antiproduct constraint: antiscalar_part(u ⊟ antirev(u)) - antinorm² = 0

    todo!("Implement using Symbolica")
}

#[derive(Debug, Clone)]
pub struct InferredConstraint {
    pub name: &'static str,
    pub expression: String,
    pub constraint_type: ConstraintType,
}

#[derive(Debug, Clone, Copy)]
pub enum ConstraintType {
    /// Linear constraint (e.g., Study condition)
    Linear,
    /// Quadratic constraint (e.g., unit norm)
    Quadratic,
}
```

#### `src/spec/ir.rs`

Add inferred constraint to TypeSpec:

```rust
pub struct TypeSpec {
    pub name: String,
    pub grades: Vec<usize>,
    pub fields: Vec<FieldSpec>,
    pub versor: Option<VersorSpec>,
    // NEW: Inferred constraint (not user-specified)
    pub inferred_constraint: Option<InferredConstraint>,
}
```

#### `src/spec/parser.rs`

Infer constraints during parsing:

```rust
impl AlgebraSpec {
    pub fn from_raw(raw: RawAlgebraSpec) -> Result<Self, ParseError> {
        // ... existing parsing ...

        // Infer constraints for each type
        for type_spec in &mut spec.types {
            type_spec.inferred_constraint = infer_constraint(
                &type_spec.grades,
                &spec.algebra,
                &spec.product_table,
                &type_spec.fields,
            );
        }

        Ok(spec)
    }
}
```

#### `src/codegen/types.rs`

Generate `new_checked()` using inferred constraints:

```rust
fn generate_new_checked(type_spec: &TypeSpec) -> TokenStream {
    let Some(constraint) = &type_spec.inferred_constraint else {
        // No constraint - new_checked just returns Ok(Self)
        return quote! {
            pub fn new_checked(/* fields */, _tolerance: T) -> Result<Self, crate::ConstraintError> {
                Ok(Self::new_unchecked(/* fields */))
            }
        };
    };

    // Generate validation code based on constraint expression
    let residual_expr = parse_constraint_to_residual(&constraint.expression);
    let constraint_name = constraint.name;

    quote! {
        pub fn new_checked(/* fields */, tolerance: T) -> Result<Self, crate::ConstraintError> {
            let residual = #residual_expr;
            if residual.abs() > tolerance {
                return Err(crate::ConstraintError::new(
                    stringify!(#type_name),
                    #constraint_name,
                    residual.to_f64().unwrap_or(0.0),
                ));
            }
            Ok(Self::new_unchecked(/* fields */))
        }
    }
}
```

## Testing

```rust
#[test]
fn infers_study_condition_for_motor() {
    let algebra = Algebra::pga(3);
    let table = ProductTable::new(&algebra);
    let fields = motor_fields();

    let constraint = infer_constraint(&[0, 2, 4], &algebra, &table, &fields);

    assert!(constraint.is_some());
    let c = constraint.unwrap();
    assert_eq!(c.name, "study");
    assert!(c.expression.contains("e0123"));
}

#[test]
fn infers_plucker_condition_for_line() {
    let algebra = Algebra::pga(3);
    let table = ProductTable::new(&algebra);
    let fields = line_fields();

    let constraint = infer_constraint(&[2], &algebra, &table, &fields);

    assert!(constraint.is_some());
    let c = constraint.unwrap();
    assert_eq!(c.name, "plucker");
}

#[test]
fn no_constraint_for_vector() {
    let algebra = Algebra::euclidean(3);
    let table = ProductTable::new(&algebra);
    let fields = vector_fields();

    // Single-grade types in Euclidean have no algebraic constraint
    let constraint = infer_constraint(&[1], &algebra, &table, &fields);
    assert!(constraint.is_none());
}
```

## Success Criteria

1. `infer_constraint()` function implemented
2. Correctly identifies Study, Plücker, unit norm constraints
3. Derives constraint expressions symbolically
4. TypeSpec includes `inferred_constraint` field
5. `new_checked()` generated using inferred constraints
6. All existing tests pass with inferred constraints
