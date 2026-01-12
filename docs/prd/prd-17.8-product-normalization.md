# PRD-17.8: Use Constrained Constructors in Generated Products

**Status**: Draft
**Parent**: [PRD-17](prd-17-codegen-products.md)
**Goal**: Replace `new_unchecked` with `new` in generated products to maintain type invariants via existing constraint solving

## Problem Statement

Generated product functions currently use `new_unchecked` for constrained output types:

```rust
// Current generated code (problematic)
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...
    Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
}
```

**The problem:** Due to floating-point arithmetic errors, the output may not satisfy the type's constraints, and `new_unchecked` bypasses all validation and constraint solving.

## Solution

Use the existing `new()` constructor which already applies constraint solving:

```rust
// Correct generated code
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...
    Motor::new(s, e23, e31, e12, e01, e02, e03)  // e0123 solved from Study constraint
}
```

The `new()` constructor:
1. Takes only the free (unconstrained) fields
2. Solves for constrained fields using the constraint expressions
3. Guarantees the output satisfies all type invariants

### How Constraint Solving Already Works

For `Motor`, the TOML specifies:

```toml
[[types.Motor.constraints]]
name = "study"
expression = "s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0"
solve_for = "e0123"
```

The generated `new()` constructor solves for `e0123`:

```rust
impl<T: Float> Motor<T> {
    pub fn new(s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T) -> Self {
        // Solve: e0123 = -(e23*e01 + e31*e02 + e12*e03) / s
        let e0123 = -(e23 * e01 + e31 * e02 + e12 * e03) / s;
        Self { s, e23, e31, e12, e01, e02, e03, e0123 }
    }
}
```

### What About Unit Norm?

The unit norm constraint is different - it's a normalization constraint, not a projection constraint:

```toml
[[types.Motor.constraints]]
name = "unit"
expression = "s*s + e23*e23 + e31*e31 + e12*e12 = 1"
solve_for = "s"
sign = "positive"
```

For this, we need to normalize *before* calling `new()`:

```rust
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...

    // Normalize weight (unit constraint)
    let wn = (s * s + e23 * e23 + e31 * e31 + e12 * e12).sqrt();
    let s = s / wn;
    let e23 = e23 / wn;
    let e31 = e31 / wn;
    let e12 = e12 / wn;
    let e01 = e01 / wn;
    let e02 = e02 / wn;
    let e03 = e03 / wn;

    // new() will solve for e0123 from Study constraint
    Motor::new(s, e23, e31, e12, e01, e02, e03)
}
```

## Implementation Plan

### Phase 1: Categorize Constraint Types

Constraints fall into two categories:

| Type | Description | Example | Handling |
|------|-------------|---------|----------|
| **Projection** | Solve for one field from others | Study: `e0123 = f(s, e23, ...)` | Use `new()` constructor |
| **Normalization** | Scale all fields to satisfy | Unit: `s² + e23² + ... = 1` | Normalize before `new()` |

Update constraint spec to indicate type:

```toml
[[types.Motor.constraints]]
name = "unit"
expression = "s*s + e23*e23 + e31*e31 + e12*e12 = 1"
constraint_type = "normalization"  # NEW
normalize_fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03"]

[[types.Motor.constraints]]
name = "study"
expression = "s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0"
constraint_type = "projection"  # NEW (default)
solve_for = "e0123"
```

### Phase 2: Update Product Generation

Modify product generation to:

1. Identify normalization constraints on output type
2. Generate normalization code for those constraints
3. Call `new()` instead of `new_unchecked()`

```rust
fn generate_product_function(&self, ..., output: &TypeSpec) -> TokenStream {
    let terms = self.compute_terms(...);

    // Get normalization constraints
    let norm_constraints = output.constraints
        .iter()
        .filter(|c| c.constraint_type == ConstraintType::Normalization);

    // Generate normalization code
    let normalization = self.generate_normalization(norm_constraints);

    // Get free fields (not solved by projection constraints)
    let free_fields = output.free_fields();

    quote! {
        pub fn #fn_name<T: Float>(a: &#lhs<T>, b: &#rhs<T>) -> #output<T> {
            #terms
            #normalization
            #output::new(#(#free_fields),*)
        }
    }
}
```

### Phase 3: Handle Edge Cases

**Types with only projection constraints (e.g., Line with Plücker):**
```rust
// Just use new(), no normalization needed
pub fn exterior_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    // ... compute terms ...
    Line::new(e01, e02, e23, e31, e12)  // e03 solved from Plücker
}
```

**Types with no constraints (e.g., Point, Plane):**
```rust
// Can use new() directly (same as new_unchecked for unconstrained types)
pub fn exterior_line_point<T: Float>(a: &Line<T>, b: &Point<T>) -> Plane<T> {
    // ... compute terms ...
    Plane::new(e023, e031, e012, e123)
}
```

**Types with multiple normalization constraints:**
```rust
// Apply normalizations in order specified by constraint ordering
```

## Deliverables

- [ ] Add `constraint_type` field to constraint spec (default: "projection")
- [ ] Add `normalize_fields` field for normalization constraints
- [ ] Update product generation to use `new()` instead of `new_unchecked()`
- [ ] Generate normalization code for normalization constraints
- [ ] Update existing TOML files with constraint types
- [ ] Regenerate all algebras
- [ ] Add stress tests verifying constraints maintained after many operations

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/spec/ir.rs` | Add ConstraintType enum, normalize_fields |
| `crates/clifford-codegen/src/spec/raw.rs` | Add raw constraint_type field |
| `crates/clifford-codegen/src/spec/parser.rs` | Parse constraint_type |
| `crates/clifford-codegen/src/codegen/products.rs` | Use new(), add normalization |
| `algebras/*.toml` | Add constraint_type to constraints |
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
3. Generated products use `new()` which applies constraint solving
4. Normalization constraints are applied before `new()` call
