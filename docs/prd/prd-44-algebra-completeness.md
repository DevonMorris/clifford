# PRD-44: Algebra Completeness Checking

**Status**: Draft
**Goal**: Ensure algebras define output types for all possible products

## Problem Statement

Currently, when a product between two types produces grades that don't match any defined type, the product is **silently skipped**. This leads to surprising behavior where some products simply don't exist.

Running `report_missing_product_matches` shows **132 missing products** across all algebras:

| Algebra | Missing | Common Missing Grades |
|---------|---------|----------------------|
| euclidean3 | 9 | [1, 3] (odd elements) |
| projective3 | 22 | [0, 2], [2, 4] |
| conformal3 | 36 | [1, 3], [0, 2], [2, 4], [1, 3, 5] |
| quaternion | 7 | [0, 2], [0, 1], [1, 2] |
| dualquat | 19 | [0, 2], [1, 3], [1, 2, 3] |
| minkowski3 | 24 | [1, 3], [0, 2], [2, 4] |
| others | 15 | various |

Example: In euclidean3, `Bivector × Vector → grades [1, 3]` but there's no `Odd` type defined, so this product doesn't exist.

## Proposed Solution

### Option A: Error on Missing Products (Recommended)

Add a `complete = true` flag to algebra specs that enforces completeness:

```toml
[algebra]
name = "euclidean3"
complete = true  # Error if any products lack output types
```

During codegen:
1. Compute all products between defined types
2. Check if each non-zero product has a matching output type
3. **Error** with a list of missing types if `complete = true`

Error message example:
```
error: Algebra 'euclidean3' is incomplete. Missing types for product outputs:
  - grades [1, 3]: Bivector × Vector, Rotor × Vector, Vector × Bivector, ...

Add type definitions or set `complete = false` to allow partial algebras.
```

### Option B: Warn on Missing Products

Add a `--warn-incomplete` flag to codegen CLI:

```bash
clifford-codegen --warn-incomplete algebras/euclidean3.toml
```

Emits warnings but doesn't fail:
```
warning: euclidean3: 9 products have no output type
  - Bivector × Vector → [1, 3] (no type)
  - Rotor × Vector → [1, 3] (no type)
  ...
```

### Option C: Auto-Generate Missing Types

Automatically create types for missing grade combinations:

```rust
// If Bivector × Vector → [1, 3] has no type, auto-generate:
[types.Grades_1_3]
grades = [1, 3]
fields = ["e1", "e2", "e3", "e123"]  // Auto-named
```

**Not recommended** - type names should be meaningful (e.g., `Odd`, `Even`, `Flector`).

## Implementation Plan

### Step 1: Add completeness check function

```rust
/// Checks if all products have matching output types.
/// Returns list of missing grade combinations.
pub fn check_algebra_completeness(
    spec: &AlgebraSpec,
    algebra: &Algebra,
) -> Vec<MissingProduct> {
    let mut missing = Vec::new();

    for product_type in ProductType::all() {
        let table = infer_all_products(&entities, *product_type, algebra);

        for (lhs, rhs, result) in &table.entries {
            if !result.is_zero && result.matching_entity.is_none() {
                missing.push(MissingProduct {
                    lhs: lhs.clone(),
                    rhs: rhs.clone(),
                    product_type: *product_type,
                    output_grades: result.output_grades.clone(),
                });
            }
        }
    }

    missing
}

pub struct MissingProduct {
    pub lhs: String,
    pub rhs: String,
    pub product_type: ProductType,
    pub output_grades: Vec<usize>,
}
```

### Step 2: Add `complete` flag to spec

```rust
// In raw.rs
pub struct RawAlgebraInfo {
    pub name: String,
    pub module_path: Option<String>,
    pub description: Option<String>,
    #[serde(default)]
    pub complete: bool,  // New field
}

// In ir.rs
pub struct AlgebraSpec {
    // ... existing fields ...
    pub complete: bool,
}
```

### Step 3: Check completeness during parsing

```rust
// In parser.rs
pub fn parse_spec(toml_content: &str) -> Result<AlgebraSpec, ParseError> {
    let raw: RawAlgebraSpec = toml::from_str(toml_content)?;

    // ... existing parsing ...

    // Check completeness if requested
    if raw.algebra.complete {
        let missing = check_algebra_completeness(&spec, &algebra);
        if !missing.is_empty() {
            return Err(ParseError::IncompleteAlgebra { missing });
        }
    }

    Ok(spec)
}
```

### Step 4: Add error type

```rust
// In error.rs
#[error("algebra '{name}' is incomplete: {count} products have no output type")]
IncompleteAlgebra {
    name: String,
    count: usize,
    missing: Vec<MissingProduct>,
},
```

### Step 5: Update existing algebras

For each algebra, either:
1. Add `complete = true` and define missing types, OR
2. Keep `complete = false` (default) for partial algebras

## Recommended Type Additions

To achieve completeness, algebras typically need:

| Algebra | Missing Type | Grades | Description |
|---------|--------------|--------|-------------|
| euclidean3 | Odd | [1, 3] | Vector + Trivector |
| projective3 | Bivector_Quadvector | [2, 4] | Mixed even grades |
| projective3 | Scalar_Bivector | [0, 2] | Even part of flector |
| quaternion | Scalar_Imaginary | [0, 1] | Full quaternion |
| quaternion | Scalar_Bivector | [0, 2] | Even quaternion |

## Verification

1. `cargo build` with `complete = true` fails for incomplete algebras
2. `cargo build` with `complete = false` succeeds (current behavior)
3. Adding missing types makes `complete = true` pass
4. Error message lists all missing grade combinations

## Future Work

- PRD-45: Blade-Level Product Inference (for sparse types)
- Suggested type generator based on missing grades
- Algebraically closed subalgebra detection

## References

- `crates/clifford-codegen/src/discovery/products.rs`
- `crates/clifford-codegen/src/spec/parser.rs`
- `report_missing_product_matches` test
