# PRD-45: Blade-Level Product Inference

**Status**: Draft
**Depends on**: PRD-44 (Algebra Completeness Checking)
**Goal**: Enable product generation for sparse types by inferring products at the blade level rather than grade level

## Problem Statement

The current product inference system operates at the **grade level**:

```
Vector (grade 1) × Vector (grade 1) → Scalar + Bivector (grades 0, 2)
```

This works for full types that span all blades of their grades. However, it fails for **sparse types** (PRD-43) that use only a subset of blades within a grade:

- `Line` uses 6 of 10 grade-3 blades in CGA
- `Plane` uses 4 of 5 grade-4 blades in CGA

Grade-level inference would incorrectly assume all 10 (or 5) blades are present, producing wrong results.

### Current Limitation

Sparse types are excluded from automatic product inference in `parser.rs`:

```rust
let entities: Vec<(String, Vec<usize>)> = types
    .iter()
    .filter(|t| t.alias_of.is_none() && !t.is_sparse)  // Sparse types excluded
    .map(|t| (t.name.clone(), t.grades.clone()))
    .collect();
```

This means:
- `Line × Scalar` → ❌ No product generated (should produce `Line`)
- `Plane × RoundPoint` → ❌ No product generated

## Current Product Coverage

Running `report_missing_product_matches` test shows **132 missing products** across all algebras:

| Algebra | Missing Products | Common Patterns |
|---------|-----------------|-----------------|
| euclidean3 | 9 | grades [1, 3] (odd elements) |
| projective3 | 22 | grades [0, 2], [2, 4] |
| conformal3 | 36 | grades [1, 3], [0, 2], [2, 4], [1, 3, 5] |
| quaternion | 7 | grades [0, 2], [0, 1], [1, 2] |
| dualquat | 19 | grades [0, 2], [1, 3], [1, 2, 3] |
| minkowski3 | 24 | grades [1, 3], [0, 2], [2, 4] |
| others | 15 | various |

Most missing products are for multi-grade outputs that don't have corresponding types defined.

## Proposed Solution

### Phase 1: Blade-Level Inference Infrastructure

Add blade-level inference alongside the existing grade-level inference:

```rust
/// Entity representation with exact blade set
pub struct EntityBladeSet {
    pub name: String,
    pub blades: BTreeSet<usize>,  // Exact blade indices
    pub grades: Vec<usize>,       // For compatibility
    pub is_sparse: bool,
}

/// Infer product output at blade level
pub fn infer_product_blades(
    lhs_blades: &[usize],
    rhs_blades: &[usize],
    product_type: ProductType,
    known_entities: &[EntityBladeSet],
    table: &ProductTable,
) -> ProductResult {
    // Compute exact output blades
    let mut output_blades = BTreeSet::new();
    for &a in lhs_blades {
        for &b in rhs_blades {
            let (sign, result) = compute_blade_product(a, b, product_type, table);
            if sign != 0 {
                output_blades.insert(result);
            }
        }
    }

    // Find matching entity by exact blade set
    if let Some(entity) = known_entities
        .iter()
        .find(|e| e.blades == output_blades)
    {
        return ProductResult {
            matching_entity: Some(entity.name.clone()),
            ..
        };
    }

    ProductResult { matching_entity: None, .. }
}
```

### Phase 2: Hybrid Inference Strategy

Use blade-level inference when any operand is sparse, grade-level otherwise:

```rust
fn infer_products_from_types(types: &[TypeSpec], signature: &SignatureSpec) -> ProductsSpec {
    // Build blade-level entity list
    let blade_entities: Vec<EntityBladeSet> = types
        .iter()
        .filter(|t| t.alias_of.is_none())
        .map(|t| EntityBladeSet {
            name: t.name.clone(),
            blades: t.fields.iter().map(|f| f.blade_index).collect(),
            grades: t.grades.clone(),
            is_sparse: t.is_sparse,
        })
        .collect();

    // For each type pair, use appropriate inference
    for lhs in &blade_entities {
        for rhs in &blade_entities {
            if lhs.is_sparse || rhs.is_sparse {
                // Use blade-level inference
                infer_product_blades(&lhs.blades, &rhs.blades, ...)
            } else {
                // Use existing grade-level inference (faster)
                infer_product(&lhs.grades, &rhs.grades, ...)
            }
        }
    }
}
```

### Phase 3: Exact Blade Matching

Entity matching must work at the blade level:

```rust
// Grade-level matching (current)
entities.iter().find(|(grades, _)| grades == &output_grades)

// Blade-level matching (new)
entities.iter().find(|e| e.blades == output_blades)
```

This enables:
- `Line × Scalar → Line` (output blades match Line's blade set exactly)
- `Plane × Scalar → Plane`
- `Line × Line → ???` (may or may not match an entity)

## Implementation Plan

### Step 1: Add EntityBladeSet struct

Add to `discovery/products.rs`:
- `EntityBladeSet` struct with name, blades, grades, is_sparse
- Helper to build from `TypeSpec`

### Step 2: Implement blade-level product computation

Add to `discovery/products.rs`:
- `compute_blade_product()` function that handles all product types
- `infer_product_blades()` main entry point

### Step 3: Update entity matching

Modify `infer_product_blades()` to:
- Match by exact blade set, not grades
- Return `None` if no exact match (product won't be generated)

### Step 4: Hybrid inference in parser

Update `infer_products_from_types()` in `spec/parser.rs`:
- Build blade-level entity list
- Use blade-level inference when any operand is sparse
- Keep grade-level inference for non-sparse (performance)

### Step 5: Re-enable sparse types in product inference

Remove the `!t.is_sparse` filter from entity building.

### Step 6: Add tests

- Test sparse × non-sparse products
- Test sparse × sparse products
- Test products with no matching entity
- Verify existing grade-level inference unchanged

## Expected Outcomes

After implementation:

| Product | Before | After |
|---------|--------|-------|
| Line × Scalar | ❌ | ✓ → Line |
| Scalar × Line | ❌ | ✓ → Line |
| Plane × Scalar | ❌ | ✓ → Plane |
| FlatPoint × Scalar | ❌ | ✓ → FlatPoint |
| Line × Line | ❌ | ✓ or ❌ (depends on output matching) |
| Line × Circle | ❌ | ✓ or ❌ (depends on output matching) |

## Error Handling

When a product is computed but no matching entity exists:

**Current behavior**: Product is silently skipped (no code generated)

**Proposed enhancement** (optional, separate PRD):
- Add `--warn-missing-products` flag to codegen
- Emit warning during build showing which products have no output type
- Add `#[cfg(feature = "strict")]` mode that fails build if products are missing

## Performance Considerations

Blade-level inference is O(|blades_a| × |blades_b|) per product vs O(|grades_a| × |grades_b|) for grade-level. In practice:
- Grade-level: ~2-5 iterations (few grades per type)
- Blade-level: ~4-100 iterations (varies by type)

Mitigation: Only use blade-level when sparse types are involved.

## Testing

1. Unit tests for `infer_product_blades()`
2. Integration tests with conformal3 sparse types
3. Verify no regression for existing algebras
4. Run `report_missing_product_matches` before/after

## Future Work

- Optimized sparse × sparse product generation
- Automatic detection of sparse output types
- Sparse type constraint propagation

## References

- PRD-43: Sparse Blade Types for Subspace Constraints
- PRD-44: Algebra Completeness Checking
- `crates/clifford-codegen/src/discovery/products.rs`
- `crates/clifford-codegen/src/spec/parser.rs`
