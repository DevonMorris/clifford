# PRD-16.1: Blade Ordering Audit and Validation

**Status**: Draft
**Parent**: PRD-16
**Goal**: Ensure consistent blade ordering between type specifications, generated code, and algebraic computations

## Problem Statement

The code generator uses a bitmask-based blade representation where each blade is identified by an integer index. For an n-dimensional algebra:

- Blade index 0 = scalar (1)
- Blade index 1 = e₁
- Blade index 2 = e₂
- Blade index 3 = e₁e₂ = e₁₂
- Blade index 4 = e₃
- Blade index 5 = e₁e₃ = e₁₃
- ...and so on

This bitmask representation is **canonical** - bit `i` indicates presence of basis vector `eᵢ`. However, there are potential ordering inconsistencies between:

1. **TOML field order**: The order of field names in the `fields = [...]` array
2. **Canonical blade order**: The order blades appear when iterating by index
3. **Product computation order**: The order terms are accumulated during product generation

### Current Implementation

In `parser.rs`, `build_fields_from_names()` assigns blade indices canonically:

```rust
fn build_fields_from_names(names: &[String], grades: &[usize], dim: usize, ...) {
    let mut blade_indices = Vec::new();
    for &grade in grades {
        for blade_index in 0usize..(1 << dim) {
            if blade_index.count_ones() as usize == grade {
                blade_indices.push((blade_index, grade));
            }
        }
    }
    // Zip with names in order provided
    names.iter().zip(blade_indices.iter())...
}
```

This assumes the TOML `fields` array is in canonical order. If a user writes:

```toml
[types.Bivector]
fields = ["yz", "xz", "xy"]  # Wrong order!
```

The parser will silently assign:
- `yz` → blade index 3 (e₁₂) - WRONG, should be e₂₃
- `xz` → blade index 5 (e₁₃) - WRONG, should be e₁₃
- `xy` → blade index 6 (e₂₃) - WRONG, should be e₁₂

### Consequences

1. **Sign errors in products**: If a field is assigned the wrong blade index, the computed product sign will be incorrect
2. **Incorrect conversions**: Multivector conversions will place coefficients in wrong positions
3. **Silent failures**: No validation catches this - tests may pass by coincidence

### Investigation Needed

The current implementation in `products.rs` uses `field.blade_index` correctly in `compute_terms()`. However, we need to verify:

1. **Parser correctness**: Does `build_fields_from_names()` always produce canonical ordering?
2. **Grade ordering**: For multi-grade types, are grades processed in ascending order?
3. **Degenerate basis handling**: In PGA (e₀ squares to 0), is the blade ordering correct?

## Solution

### Phase 1: Add Validation

Add validation in the parser to detect non-canonical field ordering:

```rust
/// Validates that provided field names match canonical blade ordering.
fn validate_field_ordering(
    names: &[String],
    grades: &[usize],
    dim: usize,
    type_name: &str,
) -> Result<(), ParseError> {
    // Build expected blade indices in canonical order
    let mut expected_indices = Vec::new();
    for &grade in grades {
        for blade_index in 0usize..(1 << dim) {
            if blade_index.count_ones() as usize == grade {
                expected_indices.push(blade_index);
            }
        }
    }

    if names.len() != expected_indices.len() {
        return Err(ParseError::FieldCountMismatch { ... });
    }

    // Names must be provided in canonical index order
    // (we can't validate the actual names without blade name mapping,
    // but we validate the count and let build_fields_from_names handle assignment)

    Ok(())
}
```

### Phase 2: Add Documentation

Document the canonical ordering requirement clearly in:

1. **TOML spec format docs**: Explain that fields must be in canonical order
2. **Error messages**: Provide helpful errors when ordering is wrong
3. **Generated code comments**: Include blade index in struct field comments

### Phase 3: Add Verification Tool

Add a `verify` subcommand that checks existing generated code:

```bash
cargo run --package clifford-codegen -- verify algebras/projective3.toml
```

This command will:

1. Parse the algebra specification
2. Compute products symbolically
3. Parse the generated Rust code (or regenerate fresh)
4. Compare expected vs actual formulas
5. Report any discrepancies

```rust
// In main.rs
Commands::Verify { spec_file } => {
    let spec = parse_spec(&std::fs::read_to_string(&spec_file)?)?;
    let algebra = build_algebra(&spec);

    // Verify each product
    for entry in &spec.products.geometric {
        let expected = compute_product_symbolically(&spec, &algebra, entry);
        let generated = parse_generated_product(&spec, entry)?;

        if expected != generated {
            eprintln!("Mismatch in {} * {} -> {}:", entry.lhs, entry.rhs, entry.output);
            eprintln!("  Expected: {:?}", expected);
            eprintln!("  Got:      {:?}", generated);
        }
    }
}
```

### Phase 4: Test Coverage

Add property-based tests that verify blade ordering:

```rust
#[test]
fn blade_indices_are_canonical() {
    for toml_file in glob("algebras/*.toml") {
        let spec = parse_spec(&std::fs::read_to_string(toml_file)?)?;

        for ty in &spec.types {
            let mut prev_grade = 0;
            let mut prev_index_in_grade = 0;

            for field in &ty.fields {
                // Fields must be ordered by grade first
                assert!(field.grade >= prev_grade,
                    "Type {} field {} has grade {} after grade {}",
                    ty.name, field.name, field.grade, prev_grade);

                // Within a grade, blade indices must be ascending
                if field.grade == prev_grade {
                    assert!(field.blade_index > prev_index_in_grade,
                        "Type {} field {} has non-ascending blade index",
                        ty.name, field.name);
                }

                prev_grade = field.grade;
                prev_index_in_grade = field.blade_index;
            }
        }
    }
}
```

## Canonical Ordering Specification

For reference, the canonical blade ordering for common dimensions:

### 2D Euclidean (Cl(2,0,0))

| Index | Binary | Blade | Grade |
|-------|--------|-------|-------|
| 0 | 00 | 1 | 0 |
| 1 | 01 | e₁ | 1 |
| 2 | 10 | e₂ | 1 |
| 3 | 11 | e₁₂ | 2 |

### 3D Euclidean (Cl(3,0,0))

| Index | Binary | Blade | Grade |
|-------|--------|-------|-------|
| 0 | 000 | 1 | 0 |
| 1 | 001 | e₁ | 1 |
| 2 | 010 | e₂ | 1 |
| 3 | 011 | e₁₂ | 2 |
| 4 | 100 | e₃ | 1 |
| 5 | 101 | e₁₃ | 2 |
| 6 | 110 | e₂₃ | 2 |
| 7 | 111 | e₁₂₃ | 3 |

### 3D PGA (Cl(3,0,1))

| Index | Binary | Blade | Grade |
|-------|--------|-------|-------|
| 0 | 0000 | 1 | 0 |
| 1 | 0001 | e₁ | 1 |
| 2 | 0010 | e₂ | 1 |
| 3 | 0011 | e₁₂ | 2 |
| 4 | 0100 | e₃ | 1 |
| 5 | 0101 | e₁₃ | 2 |
| 6 | 0110 | e₂₃ | 2 |
| 7 | 0111 | e₁₂₃ | 3 |
| 8 | 1000 | e₀ | 1 |
| 9 | 1001 | e₀₁ | 2 |
| 10 | 1010 | e₀₂ | 2 |
| 11 | 1011 | e₀₁₂ | 3 |
| 12 | 1100 | e₀₃ | 2 |
| 13 | 1101 | e₀₁₃ | 3 |
| 14 | 1110 | e₀₂₃ | 3 |
| 15 | 1111 | e₀₁₂₃ | 4 |

**Note**: The degenerate basis vector e₀ is assigned the highest bit (index 3 in 4D). This matches the TOML signature order: `positive = ["e1", "e2", "e3"]`, `zero = ["e0"]`.

## Deliverables

- [ ] Add `validate_field_ordering()` function in `parser.rs`
- [ ] Add `ParseError::NonCanonicalFieldOrdering` variant
- [ ] Update documentation with canonical ordering rules
- [ ] Add `verify` subcommand to CLI
- [ ] Add property-based tests for blade ordering
- [ ] Audit all existing algebra specs in `algebras/`

## Success Criteria

1. Parser rejects TOML specs with non-canonical field ordering
2. `verify` command passes for all existing algebras
3. Documentation clearly explains the ordering requirement
4. All tests pass

## Files Changed

| File | Action | Description |
|------|--------|-------------|
| `crates/clifford-codegen/src/spec/parser.rs` | Update | Add validation |
| `crates/clifford-codegen/src/spec/error.rs` | Update | Add error variant |
| `crates/clifford-codegen/src/main.rs` | Update | Add verify subcommand |
| `algebras/*.toml` | Audit | Verify all specs have correct ordering |

## Testing

```rust
#[test]
fn reject_non_canonical_ordering() {
    let result = parse_spec(r#"
        [algebra]
        name = "bad"

        [signature]
        positive = ["e1", "e2", "e3"]

        [types.Bivector]
        grades = [2]
        fields = ["yz", "xz", "xy"]  # Wrong order: should be xy, xz, yz
    "#);

    assert!(matches!(result, Err(ParseError::NonCanonicalFieldOrdering { .. })));
}
```
