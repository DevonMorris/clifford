# PRD-17.8: Use Constrained Constructors in Generated Products

**Status**: Draft
**Parent**: [PRD-17](prd-17-codegen-products.md)
**Goal**: Replace `new_unchecked` with `new` in generated products to maintain type invariants via existing Symbolica-based constraint solving

## Problem Statement

Generated product functions currently use `new_unchecked` for constrained output types:

```rust
// Current generated code (problematic)
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...
    Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
}
```

**The problem:** Due to floating-point arithmetic errors, the output may not satisfy the type's constraints, and `new_unchecked` bypasses all constraint solving.

## Solution

Use the existing `new()` constructor which applies Symbolica-based constraint solving:

```rust
// Correct generated code
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...
    Motor::new(s, e23, e31, e12, e01, e02, e03)  // Symbolica solves for e0123
}
```

The constraint system already:
1. Parses constraint expressions from TOML
2. Uses Symbolica to solve for constrained fields
3. Generates `new()` constructors that only take free fields
4. Automatically computes constrained field values

**No changes to the constraint system are needed** - just use `new()` instead of `new_unchecked()` in generated products.

## Implementation

### Update Product Generation

In `products.rs`, change output construction:

```rust
fn generate_product_function(&self, ..., output: &TypeSpec) -> TokenStream {
    let terms = self.compute_terms(...);

    if output.has_constraints() {
        // Use new() which applies constraint solving
        let free_fields = output.free_fields();
        quote! {
            pub fn #fn_name<T: Float>(a: &#lhs<T>, b: &#rhs<T>) -> #output<T> {
                #terms
                #output::new(#(#free_fields),*)
            }
        }
    } else {
        // Unconstrained types can use direct construction
        let all_fields = output.all_fields();
        quote! {
            pub fn #fn_name<T: Float>(a: &#lhs<T>, b: &#rhs<T>) -> #output<T> {
                #terms
                #output::new(#(#all_fields),*)
            }
        }
    }
}
```

### Affected Products

Any product that outputs a constrained type:

| Output Type | Constraints | Current | Fixed |
|-------------|-------------|---------|-------|
| Motor | unit, study | `new_unchecked(8 fields)` | `new(7 fields)` |
| Line | plücker | `new_unchecked(6 fields)` | `new(5 fields)` |
| Flector | geometric | `new_unchecked(8 fields)` | `new(7 fields)` |

## Deliverables

- [ ] Update `generate_product_function()` to use `new()` for constrained outputs
- [ ] Update `generate_sandwich_function()` similarly
- [ ] Update `generate_antisandwich_function()` similarly
- [ ] Regenerate all algebras
- [ ] Add stress tests verifying constraints maintained after many operations

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/codegen/products.rs` | Use `new()` instead of `new_unchecked()` |
| `src/generated/*/products.rs` | Regenerate |

## Testing Strategy

```rust
proptest! {
    #[test]
    fn motor_composition_maintains_constraints(
        motors in prop::collection::vec(any::<Motor<f64>>(), 100..200)
    ) {
        let mut result = Motor::identity();
        for m in motors {
            result = geometric_motor_motor(&result, &m);
        }

        // Constraints should be satisfied
        prop_assert!(result.study_residual().abs() < 1e-10);
        prop_assert!((result.weight_norm() - 1.0).abs() < 1e-10);
    }
}
```

## Success Criteria

1. No `new_unchecked` in generated products for constrained types
2. All constraint residuals < 1e-10 after 1000+ operations
3. Generated products use `new()` which applies Symbolica constraint solving
