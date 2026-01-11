# PRD-14.2: Specification Format & Parser

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Define and implement the TOML specification format for algebras

## Overview

The specification format describes:
1. The algebra's metric signature
2. Which grade combinations become named types
3. Field names and orderings for each type
4. Constraints (unit, normalized, nonzero, null)
5. Constraint-preserving operations

The parser validates specifications and produces an intermediate representation consumed by code generators.

## Specification Schema

### Top-Level Structure

```toml
[algebra]
name = "euclidean3"                    # Required: identifier for the algebra
module_path = "euclidean::dim3"        # Optional: Rust module path
description = "3D Euclidean GA"        # Optional: doc comment

[signature]
positive = ["e1", "e2", "e3"]          # Basis vectors squaring to +1
negative = []                           # Basis vectors squaring to -1
zero = []                               # Basis vectors squaring to 0

[blades]
# Custom names for blades (optional, defaults to concatenation)
e12 = "xy"
e13 = "xz"
e23 = "yz"
e123 = "xyz"

[types.TypeName]
grades = [0, 2]                         # Which grades this type contains
description = "..."                     # Doc comment
fields = ["s", "xy", "xz", "yz"]       # Optional: custom field names

[types.TypeName.constraints.unit]
# Defines a UnitTypeName wrapper
norm = "euclidean"                      # Norm type
constructors = [...]                    # Safe constructors
preserving_ops = [...]                  # Operations preserving constraint

[options]
# Generation options
generate_nalgebra = true
generate_serde = true
generate_arbitrary = true
```

### Signature Section

```toml
[signature]
# Names are used for documentation; order determines index
positive = ["e1", "e2", "e3"]    # Indices 0, 1, 2; metric +1
negative = ["e_minus"]            # Index 3; metric -1
zero = ["e0"]                     # Index 4; metric 0
```

**Ordering Rules**:
1. Positive basis vectors come first (indices 0 to p-1)
2. Negative basis vectors next (indices p to p+q-1)
3. Zero/degenerate last (indices p+q to p+q+r-1)

This matches the standard GA convention: signature (p, q, r).

### Blades Section

```toml
[blades]
# Map canonical blade names to custom field names
# Canonical names use basis indices: e1, e2, e12, e123, etc.
e1 = "x"
e2 = "y"
e3 = "z"
e12 = "xy"
e13 = "xz"
e23 = "yz"
e123 = "xyz"
```

If not specified, blades use concatenated basis names.

### Types Section

```toml
[types.Vector]
grades = [1]                           # Contains only grade-1 blades
description = "Grade-1 vector"
fields = ["x", "y", "z"]              # Must match blade count for grade 1

[types.Rotor]
grades = [0, 2]                        # Contains grades 0 and 2
description = "Rotation element"
fields = ["s", "xy", "xz", "yz"]      # scalar + 3 bivector components
alias_of = "Even"                      # Optional: same storage as Even
```

**Field Ordering Convention**:
1. Fields are ordered by grade (ascending)
2. Within a grade, fields are ordered by blade index (ascending)

If `fields` is omitted, names come from the `[blades]` section.

### Constraints Section

```toml
[types.Rotor.constraints.unit]
# Generates: UnitRotor<T>
norm = "euclidean"                     # ||R||² = s² + xy² + xz² + yz²
constructors = [
    "identity()",
    "from_angle_plane(angle: T, plane: Bivector<T>)",
    "from_vectors(a: Vector<T>, b: Vector<T>)",
]
preserving_ops = ["compose", "inverse", "slerp"]

[types.Vector.constraints.nonzero]
# Generates: NonZeroVector<T>
condition = "norm_squared > epsilon"
constructors = []                       # Only try_from

[types.Point.constraints.null]
# For CGA: generates NullPoint<T> (or just use Point with constraint)
condition = "inner_self == 0"
constructors = [
    "from_euclidean(x: T, y: T, z: T)",
    "origin()",
]
preserving_ops = ["sandwich_by_motor"]
```

**Constraint Types**:

| Type | Meaning | Generated Name |
|------|---------|----------------|
| `unit` | ‖x‖ = 1 | `Unit{Type}` |
| `nonzero` | ‖x‖ ≠ 0 | `NonZero{Type}` |
| `normalized` | Context-dependent canonical form | `Normalized{Type}` |
| `null` | x · x = 0 | Uses base name or custom |
| `ideal` | In ideal subspace | `Ideal{Type}` |

**Norm Types**:

| Name | Formula |
|------|---------|
| `euclidean` | s² + x² + y² + ... |
| `motor` | PGA motor norm (scalar + Euclidean bivector part) |
| `cga` | CGA-specific norm |

## Intermediate Representation

The parser produces this IR:

```rust
/// Parsed algebra specification.
#[derive(Debug, Clone)]
pub struct AlgebraSpec {
    /// Algebra identifier
    pub name: String,
    /// Rust module path
    pub module_path: Option<String>,
    /// Documentation
    pub description: Option<String>,
    /// Metric signature
    pub signature: SignatureSpec,
    /// Blade name mappings
    pub blade_names: HashMap<usize, String>,
    /// Type definitions
    pub types: Vec<TypeSpec>,
    /// Generation options
    pub options: GenerationOptions,
}

#[derive(Debug, Clone)]
pub struct SignatureSpec {
    /// Basis vector names and metrics
    pub basis: Vec<BasisVector>,
    /// p (positive), q (negative), r (zero) counts
    pub p: usize,
    pub q: usize,
    pub r: usize,
}

#[derive(Debug, Clone)]
pub struct BasisVector {
    pub name: String,
    pub index: usize,
    pub metric: i8,  // +1, -1, or 0
}

#[derive(Debug, Clone)]
pub struct TypeSpec {
    /// Type name (e.g., "Rotor")
    pub name: String,
    /// Grades contained
    pub grades: Vec<usize>,
    /// Documentation
    pub description: Option<String>,
    /// Field names in order
    pub fields: Vec<FieldSpec>,
    /// Type this aliases (same storage)
    pub alias_of: Option<String>,
    /// Constraints
    pub constraints: Vec<ConstraintSpec>,
}

#[derive(Debug, Clone)]
pub struct FieldSpec {
    /// Field name
    pub name: String,
    /// Blade index this field holds
    pub blade_index: usize,
    /// Grade of this blade
    pub grade: usize,
}

#[derive(Debug, Clone)]
pub struct ConstraintSpec {
    /// Constraint kind
    pub kind: ConstraintKind,
    /// Generated wrapper name (e.g., "UnitRotor")
    pub wrapper_name: String,
    /// Norm type for unit constraints
    pub norm_type: Option<NormType>,
    /// Condition expression for other constraints
    pub condition: Option<String>,
    /// Constructor signatures
    pub constructors: Vec<ConstructorSpec>,
    /// Operations that preserve this constraint
    pub preserving_ops: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintKind {
    Unit,
    NonZero,
    Normalized,
    Null,
    Ideal,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormType {
    Euclidean,
    Motor,
    Cga,
    Custom,
}

#[derive(Debug, Clone)]
pub struct ConstructorSpec {
    pub name: String,
    pub params: Vec<ParamSpec>,
}

#[derive(Debug, Clone)]
pub struct ParamSpec {
    pub name: String,
    pub ty: String,
}

#[derive(Debug, Clone, Default)]
pub struct GenerationOptions {
    pub generate_nalgebra: bool,
    pub generate_serde: bool,
    pub generate_arbitrary: bool,
    pub generate_tests: bool,
}
```

## Parser Implementation

```rust
use serde::Deserialize;
use std::collections::HashMap;

/// Raw TOML structure (matches file format exactly)
#[derive(Debug, Deserialize)]
struct RawAlgebraSpec {
    algebra: RawAlgebraInfo,
    signature: RawSignature,
    #[serde(default)]
    blades: HashMap<String, String>,
    #[serde(default)]
    types: HashMap<String, RawTypeSpec>,
    #[serde(default)]
    options: RawOptions,
}

#[derive(Debug, Deserialize)]
struct RawAlgebraInfo {
    name: String,
    module_path: Option<String>,
    description: Option<String>,
}

#[derive(Debug, Deserialize)]
struct RawSignature {
    #[serde(default)]
    positive: Vec<String>,
    #[serde(default)]
    negative: Vec<String>,
    #[serde(default)]
    zero: Vec<String>,
}

#[derive(Debug, Deserialize)]
struct RawTypeSpec {
    grades: Vec<usize>,
    description: Option<String>,
    #[serde(default)]
    fields: Vec<String>,
    alias_of: Option<String>,
    #[serde(default)]
    constraints: HashMap<String, RawConstraint>,
}

#[derive(Debug, Deserialize)]
struct RawConstraint {
    norm: Option<String>,
    condition: Option<String>,
    #[serde(default)]
    constructors: Vec<String>,
    #[serde(default)]
    preserving_ops: Vec<String>,
}

#[derive(Debug, Deserialize, Default)]
struct RawOptions {
    #[serde(default)]
    generate_nalgebra: bool,
    #[serde(default)]
    generate_serde: bool,
    #[serde(default = "default_true")]
    generate_arbitrary: bool,
    #[serde(default = "default_true")]
    generate_tests: bool,
}

fn default_true() -> bool { true }

/// Parses a TOML specification into the IR.
pub fn parse_spec(toml_content: &str) -> Result<AlgebraSpec, ParseError> {
    let raw: RawAlgebraSpec = toml::from_str(toml_content)?;

    // Build signature
    let signature = parse_signature(&raw.signature)?;

    // Build blade name map
    let blade_names = parse_blade_names(&raw.blades, &signature)?;

    // Build types
    let types = parse_types(&raw.types, &signature, &blade_names)?;

    // Validate
    validate_spec(&signature, &types)?;

    Ok(AlgebraSpec {
        name: raw.algebra.name,
        module_path: raw.algebra.module_path,
        description: raw.algebra.description,
        signature,
        blade_names,
        types,
        options: GenerationOptions {
            generate_nalgebra: raw.options.generate_nalgebra,
            generate_serde: raw.options.generate_serde,
            generate_arbitrary: raw.options.generate_arbitrary,
            generate_tests: raw.options.generate_tests,
        },
    })
}
```

## Validation

The parser validates:

### 1. Signature Validation

```rust
fn validate_signature(sig: &SignatureSpec) -> Result<(), ParseError> {
    // Check dimension is reasonable
    let dim = sig.p + sig.q + sig.r;
    if dim > 6 {
        return Err(ParseError::DimensionTooLarge(dim));
    }
    if dim == 0 {
        return Err(ParseError::EmptySignature);
    }

    // Check unique basis names
    let mut names = HashSet::new();
    for basis in &sig.basis {
        if !names.insert(&basis.name) {
            return Err(ParseError::DuplicateBasisName(basis.name.clone()));
        }
    }

    Ok(())
}
```

### 2. Type Validation

```rust
fn validate_type(ty: &TypeSpec, sig: &SignatureSpec) -> Result<(), ParseError> {
    let dim = sig.p + sig.q + sig.r;

    // Check grades are valid
    for &grade in &ty.grades {
        if grade > dim {
            return Err(ParseError::InvalidGrade {
                type_name: ty.name.clone(),
                grade,
                max: dim,
            });
        }
    }

    // Check field count matches
    let expected_fields: usize = ty.grades.iter()
        .map(|&g| binomial(dim, g))
        .sum();

    if ty.fields.len() != expected_fields {
        return Err(ParseError::FieldCountMismatch {
            type_name: ty.name.clone(),
            expected: expected_fields,
            got: ty.fields.len(),
        });
    }

    // Check field names are unique
    let mut field_names = HashSet::new();
    for field in &ty.fields {
        if !field_names.insert(&field.name) {
            return Err(ParseError::DuplicateFieldName {
                type_name: ty.name.clone(),
                field: field.name.clone(),
            });
        }
    }

    Ok(())
}
```

### 3. Constraint Validation

```rust
fn validate_constraints(ty: &TypeSpec) -> Result<(), ParseError> {
    for constraint in &ty.constraints {
        // Unit constraint requires norm type
        if constraint.kind == ConstraintKind::Unit && constraint.norm_type.is_none() {
            return Err(ParseError::MissingNormType {
                type_name: ty.name.clone(),
            });
        }

        // Parse constructor signatures
        for ctor in &constraint.constructors {
            parse_constructor_signature(&ctor.name)?;
        }
    }

    Ok(())
}
```

## Error Handling

```rust
#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("TOML parse error: {0}")]
    Toml(#[from] toml::de::Error),

    #[error("Dimension {0} exceeds maximum of 6")]
    DimensionTooLarge(usize),

    #[error("Signature must have at least one basis vector")]
    EmptySignature,

    #[error("Duplicate basis vector name: {0}")]
    DuplicateBasisName(String),

    #[error("Type '{type_name}' has invalid grade {grade} (max: {max})")]
    InvalidGrade {
        type_name: String,
        grade: usize,
        max: usize,
    },

    #[error("Type '{type_name}' has {got} fields but needs {expected}")]
    FieldCountMismatch {
        type_name: String,
        expected: usize,
        got: usize,
    },

    #[error("Type '{type_name}' has duplicate field name: {field}")]
    DuplicateFieldName {
        type_name: String,
        field: String,
    },

    #[error("Unit constraint for '{type_name}' requires norm type")]
    MissingNormType {
        type_name: String,
    },

    #[error("Invalid blade name: {0}")]
    InvalidBladeName(String),

    #[error("Unknown type reference: {0}")]
    UnknownType(String),
}
```

## Bundled Specifications

The generator includes specifications for common algebras:

### euclidean2.toml

```toml
[algebra]
name = "euclidean2"
module_path = "euclidean::dim2"
description = "2D Euclidean Geometric Algebra"

[signature]
positive = ["e1", "e2"]

[blades]
e1 = "x"
e2 = "y"
e12 = "xy"

[types.Scalar]
grades = [0]
description = "Grade-0 scalar"
fields = ["s"]

[types.Vector]
grades = [1]
description = "Grade-1 vector"
fields = ["x", "y"]

[types.Vector.constraints.unit]
norm = "euclidean"
constructors = ["unit_x()", "unit_y()"]

[types.Vector.constraints.nonzero]
condition = "norm_squared > epsilon"

[types.Bivector]
grades = [2]
description = "Grade-2 pseudoscalar"
fields = ["xy"]

[types.Rotor]
grades = [0, 2]
description = "2D rotation element"
fields = ["s", "xy"]

[types.Rotor.constraints.unit]
norm = "euclidean"
constructors = ["identity()", "from_angle(angle: T)"]
preserving_ops = ["compose", "inverse"]

[types.Even]
grades = [0, 2]
description = "Even subalgebra"
fields = ["s", "xy"]
alias_of = "Rotor"

[types.Full]
grades = [0, 1, 2]
description = "Full multivector"
fields = ["s", "x", "y", "xy"]
```

### euclidean3.toml

```toml
[algebra]
name = "euclidean3"
module_path = "euclidean::dim3"
description = "3D Euclidean Geometric Algebra"

[signature]
positive = ["e1", "e2", "e3"]

[blades]
e1 = "x"
e2 = "y"
e3 = "z"
e12 = "xy"
e13 = "xz"
e23 = "yz"
e123 = "xyz"

[types.Scalar]
grades = [0]
fields = ["s"]

[types.Vector]
grades = [1]
description = "Grade-1 vector: e₁, e₂, e₃"
fields = ["x", "y", "z"]

[types.Vector.constraints.unit]
norm = "euclidean"
constructors = ["unit_x()", "unit_y()", "unit_z()"]

[types.Vector.constraints.nonzero]
condition = "norm_squared > epsilon"

[types.Bivector]
grades = [2]
description = "Grade-2 bivector: e₁₂, e₁₃, e₂₃"
fields = ["xy", "xz", "yz"]

[types.Bivector.constraints.unit]
norm = "euclidean"
constructors = ["unit_xy()", "unit_xz()", "unit_yz()"]

[types.Trivector]
grades = [3]
description = "Grade-3 pseudoscalar: e₁₂₃"
fields = ["xyz"]

[types.Rotor]
grades = [0, 2]
description = "3D rotation element"
fields = ["s", "xy", "xz", "yz"]

[types.Rotor.constraints.unit]
norm = "euclidean"
constructors = [
    "identity()",
    "from_angle_plane(angle: T, plane: Bivector<T>)",
    "from_angle_axis(angle: T, axis: Vector<T>)",
    "from_vectors(a: Vector<T>, b: Vector<T>)",
]
preserving_ops = ["compose", "inverse", "slerp"]

[types.Even]
grades = [0, 2]
description = "Even subalgebra"
fields = ["s", "xy", "xz", "yz"]
alias_of = "Rotor"

[types.Odd]
grades = [1, 3]
description = "Odd subalgebra"
fields = ["x", "y", "z", "xyz"]

[types.Full]
grades = [0, 1, 2, 3]
description = "Full multivector"
fields = ["s", "x", "y", "z", "xy", "xz", "yz", "xyz"]
```

## Testing

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_euclidean3() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();

        assert_eq!(spec.name, "euclidean3");
        assert_eq!(spec.signature.p, 3);
        assert_eq!(spec.signature.q, 0);
        assert_eq!(spec.signature.r, 0);

        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.grades, vec![0, 2]);
        assert_eq!(rotor.fields.len(), 4); // s, xy, xz, yz
    }

    #[test]
    fn reject_invalid_grade() {
        let toml = r#"
        [algebra]
        name = "test"

        [signature]
        positive = ["e1", "e2"]

        [types.Bad]
        grades = [5]  # Invalid: max grade is 2
        "#;

        let err = parse_spec(toml).unwrap_err();
        assert!(matches!(err, ParseError::InvalidGrade { grade: 5, .. }));
    }

    #[test]
    fn reject_field_count_mismatch() {
        let toml = r#"
        [algebra]
        name = "test"

        [signature]
        positive = ["e1", "e2", "e3"]

        [types.Vector]
        grades = [1]
        fields = ["x", "y"]  # Missing z
        "#;

        let err = parse_spec(toml).unwrap_err();
        assert!(matches!(err, ParseError::FieldCountMismatch { expected: 3, got: 2, .. }));
    }
}
```

## Deliverables

- [ ] TOML schema definition
- [ ] Raw TOML deserializer (serde)
- [ ] IR types (`AlgebraSpec`, `TypeSpec`, etc.)
- [ ] Parser with validation
- [ ] Error types with helpful messages
- [ ] Bundled specifications (euclidean2, euclidean3, pga2, pga3, cga2, cga3)
- [ ] Parser tests

## Dependencies

- `toml` - TOML parsing
- `serde` - Deserialization
- `thiserror` - Error handling

## Success Criteria

1. All bundled specifications parse without error
2. Invalid specifications produce clear error messages
3. IR contains all information needed for code generation
