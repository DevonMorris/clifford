# PRD-17.8: Product Normalization for Constrained Types

**Status**: Draft
**Parent**: [PRD-17](prd-17-codegen-products.md)
**Goal**: Remove `new_unchecked` from generated products and add automatic normalization to maintain type invariants

## Problem Statement

Generated product functions currently use `new_unchecked` for constrained output types:

```rust
// Current generated code (problematic)
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...
    Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
}
```

**The problem:** Due to floating-point arithmetic errors, the output may not satisfy the type's constraints:

1. **Study condition drift**: `s * e0123 + e23 * e01 + e31 * e02 + e12 * e03 = 0` may become non-zero
2. **Unit norm drift**: `s² + e23² + e31² + e12² = 1` may drift from unity
3. **Plücker constraint drift**: `e01 * e23 + e02 * e31 + e03 * e12 = 0` for lines
4. **Accumulation**: Errors compound with repeated operations (e.g., animation loops)

### Concrete Example

```rust
// After many compositions, errors accumulate
let mut motor = Motor::identity();
for _ in 0..1000 {
    motor = motor.compose(&small_rotation);
}

// motor may no longer satisfy Study condition
// motor.study_residual() could be 1e-6 or worse
// motor.weight_norm() could be 1.0001 or 0.9999
```

### Current Workaround

Users must manually call `normalized()` after operations:

```rust
// Current user code (verbose, easy to forget)
let result = motor1.compose(&motor2).normalized();
```

This is:
- Easy to forget
- Inconsistent (some code paths normalize, others don't)
- Not enforced by the type system

## Solution

### Option A: Normalize in Generated Products (Recommended)

Generated products for constrained types should automatically normalize:

```rust
// Generated code with normalization
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...

    // Normalize to maintain invariants
    let wn = (s * s + e23 * e23 + e31 * e31 + e12 * e12).sqrt();
    let s = s / wn;
    let e23 = e23 / wn;
    let e31 = e31 / wn;
    let e12 = e12 / wn;
    let e01 = e01 / wn;
    let e02 = e02 / wn;
    let e03 = e03 / wn;

    // Project e0123 to satisfy Study condition
    let e0123 = -(e23 * e01 + e31 * e02 + e12 * e03) / s;

    Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
}
```

**Pros:**
- Automatic, no user intervention needed
- Invariants always maintained
- Small overhead per operation (sqrt + division)

**Cons:**
- Slightly slower (but typically negligible)
- May not be desired for intermediate calculations

### Option B: Separate Normalized and Raw Products

Generate both variants:

```rust
// Raw product (for advanced users, intermediate calculations)
pub fn geometric_motor_motor_raw<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute without normalization ...
    Motor::new_unchecked(...)
}

// Normalized product (default, maintains invariants)
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    geometric_motor_motor_raw(a, b).normalized()
}
```

**Pros:**
- Flexibility for advanced users
- Clear intent in function names

**Cons:**
- API surface doubles
- Users might accidentally use `_raw` variants

### Option C: Configurable via TOML

Allow per-type configuration:

```toml
[types.Motor]
grades = [0, 2, 4]
versor = true
constraints = [...]

# New: normalization policy for products
normalization = "always"  # or "never" or "explicit"
```

## Implementation Plan

### Phase 1: Add Normalization Functions to Generated Types

For each constrained type, generate a `normalize()` method that enforces all constraints:

```rust
impl<T: Float> Motor<T> {
    /// Normalize to satisfy all constraints (unit norm + Study condition).
    ///
    /// This should be called after operations that may cause drift.
    #[inline]
    pub fn normalize(&self) -> Self {
        // 1. Normalize weight norm
        let wn = self.weight_norm();
        if wn < T::epsilon() {
            return *self;
        }

        let s = self.s() / wn;
        let e23 = self.e23() / wn;
        let e31 = self.e31() / wn;
        let e12 = self.e12() / wn;
        let e01 = self.e01() / wn;
        let e02 = self.e02() / wn;
        let e03 = self.e03() / wn;

        // 2. Project e0123 to satisfy Study condition
        let e0123 = if s.abs() > T::epsilon() {
            -(e23 * e01 + e31 * e02 + e12 * e03) / s
        } else {
            T::zero()
        };

        Self::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
    }
}
```

### Phase 2: Update Product Generation

Modify `generate_product_function()` to optionally include normalization:

```rust
fn generate_product_function(
    &self,
    lhs: &TypeSpec,
    rhs: &TypeSpec,
    output: &TypeSpec,
    kind: ProductKind,
) -> TokenStream {
    let terms = self.compute_terms(lhs, rhs, output, kind);
    let construction = self.generate_construction(&terms, output);

    // If output type has constraints, wrap in normalization
    if output.has_constraints() && self.should_normalize(kind) {
        quote! {
            pub fn #fn_name<T: Float>(a: &#lhs_ty<T>, b: &#rhs_ty<T>) -> #output_ty<T> {
                #construction.normalize()
            }
        }
    } else {
        quote! {
            pub fn #fn_name<T: Float>(a: &#lhs_ty<T>, b: &#rhs_ty<T>) -> #output_ty<T> {
                #construction
            }
        }
    }
}
```

### Phase 3: TOML Configuration

Add normalization policy to type specs:

```toml
[types.Motor]
grades = [0, 2, 4]
versor = true

# Normalization policy for this type
[types.Motor.normalization]
# When to normalize products that output this type
products = "always"  # "always", "never", "configurable"

# Constraints to enforce during normalization
constraints = ["unit", "study"]

[types.Line]
grades = [2]

[types.Line.normalization]
products = "never"  # Lines don't need normalization typically
```

### Phase 4: Update Constraint Specification

Extend constraint spec to include normalization formula:

```toml
[[types.Motor.constraints]]
name = "unit"
expression = "s*s + e23*e23 + e31*e31 + e12*e12 = 1"
solve_for = "s"
sign = "positive"
# New: normalization formula
normalize = "divide_by_sqrt"  # Divide all components by sqrt of LHS

[[types.Motor.constraints]]
name = "study"
expression = "s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0"
solve_for = "e0123"
# New: normalization happens after unit constraint
normalize_order = 2
```

## Affected Products

Products that output constrained types and need normalization:

| Product | Output Type | Constraints |
|---------|-------------|-------------|
| `geometric_motor_motor` | Motor | unit + study |
| `geometric_flector_flector` | Motor | unit + study |
| `sandwich_motor_motor` | Motor | unit + study |
| `antisandwich_motor_motor` | Motor | unit + study |
| `geometric_rotor_rotor` | Rotor | unit |
| `sandwich_rotor_*` | Rotor | unit |
| `exterior_point_point` | Line | plücker |

## Performance Considerations

### Cost of Normalization

Per motor normalization:
- 4 multiplications (for weight norm squared)
- 1 sqrt
- 8 divisions (normalize components)
- 4 multiplications + 2 additions (for Study projection)

**Total:** ~20 FLOPs per normalization

For typical use cases (transforms at 60fps), this is negligible.

### When NOT to Normalize

Some operations don't need normalization:
- Intermediate calculations that will be normalized at the end
- Read-only inspections (no modification)
- Known-clean inputs (e.g., identity motor)

The `_raw` variant (Option B) allows bypassing for these cases.

## Deliverables

- [ ] Add normalization policy to TypeSpec
- [ ] Generate `normalize()` method for constrained types
- [ ] Update product generation to call normalize when needed
- [ ] Add TOML configuration for normalization policy
- [ ] Generate `_raw` variants for advanced users (optional)
- [ ] Update existing extension methods to use generated products
- [ ] Add tests verifying constraints are maintained after many operations
- [ ] Document normalization behavior in generated code

## Testing Strategy

### Stress Tests

```rust
proptest! {
    #[test]
    fn motor_composition_maintains_study_condition(
        motors in prop::collection::vec(any::<Motor<f64>>(), 100..200)
    ) {
        let mut result = Motor::identity();
        for m in motors {
            result = result.compose(&m);
        }

        // Should still satisfy constraints
        prop_assert!(result.study_residual().abs() < 1e-10);
        prop_assert!((result.weight_norm() - 1.0).abs() < 1e-10);
    }
}
```

### Comparison Tests

```rust
#[test]
fn normalized_vs_raw_products() {
    let a = Motor::from_rotation_x(0.1);
    let b = Motor::from_translation(1.0, 2.0, 3.0);

    // Both should give same result for clean inputs
    let normalized = geometric_motor_motor(&a, &b);
    let raw = geometric_motor_motor_raw(&a, &b);

    assert!(abs_diff_eq!(normalized, raw, epsilon = 1e-14));
}
```

## Success Criteria

1. Generated products maintain type constraints within floating-point tolerance
2. No `new_unchecked` in products that output constrained types (except in `_raw` variants)
3. Stress tests pass with 1000+ sequential operations
4. Performance regression < 5% for typical workloads

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/codegen/products.rs` | Add normalization to product generation |
| `crates/clifford-codegen/src/codegen/types.rs` | Generate normalize() methods |
| `crates/clifford-codegen/src/spec/ir.rs` | Add normalization policy to TypeSpec |
| `crates/clifford-codegen/src/spec/raw.rs` | Add raw normalization fields |
| `crates/clifford-codegen/src/spec/parser.rs` | Parse normalization config |
| `algebras/*.toml` | Add normalization configuration |
| `src/generated/*/products.rs` | Regenerate with normalization |
| `src/generated/*/types.rs` | Regenerate with normalize() methods |

## References

- [Geometric Algebra for Computer Science](https://geometricalgebra.org/) - Numerical stability
- [Quaternion normalization](https://www.3dgep.com/understanding-quaternions/) - Similar problem in quaternion math
- Current `Motor::normalized()` implementation in extensions.rs
