# PRD-17.8: Use Constrained Constructors in Generated Products

**Status**: Won't Do (Investigated)
**Parent**: [PRD-17](prd-17-codegen-products.md)
**Goal**: ~~Replace `new_unchecked` with `new` in generated products to maintain type invariants via existing Symbolica-based constraint solving~~

## Summary

This PRD proposed using `new()` instead of `new_unchecked()` for constrained output types in generated products. After investigation and implementation, **this approach does not work** because it fundamentally misunderstands the relationship between algebraic products and type constraints.

## Problem Statement (Original)

Generated product functions use `new_unchecked` for constrained output types. Due to floating-point arithmetic errors, the output may not satisfy the type's constraints.

```rust
// Current generated code
pub fn geometric_motor_motor<T: Float>(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    // ... compute product terms ...
    Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
}
```

## Why This Approach Fails

### Root Cause: Constraints vs. Product Outputs

The constraint system (e.g., Motor's Study condition: `s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0`) applies to **normalized, valid instances** of a type, not to **arbitrary product results**.

When computing products like `Motor ∧ Motor` (exterior product):
1. The algebraic result naturally produces `e0123 = 0` (wedge products don't generate grade-4 from these inputs)
2. Using `new()` would apply constraint solving, computing `e0123 = -(e23*e01 + e31*e02 + e12*e03) / s`
3. This adds a **non-zero** `e0123` term that shouldn't exist

### Concrete Example

```
Input: Two Motors a, b
Exterior product (correct): 0.91*1 + 0.40*e23 + 0.40*e34
With constraint solving (wrong): 0.91*1 + 0.40*e23 + 0.40*e34 + -0.18*e1234
```

The constraint solver incorrectly "fixes" the result by adding a term that the algebraic operation never produced.

### Key Insight

**Product outputs are mathematically correct as computed.** The constraint system is for:
- Validating user-constructed instances
- Ensuring types created via factory methods are valid
- Computing derived fields during explicit normalization

**NOT for:**
- Modifying algebraic product results
- "Correcting" intermediate computations

## Resolution

Keep using `new_unchecked()` for constrained types in product outputs because:

1. **Product outputs are algebraically correct** - The formulas are derived from the multiplication table
2. **Constraint solving modifies results** - It changes values that should remain as computed
3. **Constraints apply to normalized instances** - Not to every multivector with the same grade structure

The correct place for constraint enforcement is:
- `Motor::normalize()` - explicit re-normalization when drift accumulates
- Factory methods like `Motor::from_translation()` - guarantee validity by construction
- `Motor::try_from_components()` - validate user inputs

## Lessons Learned

1. **Type constraints are not always enforced** - A `Motor` returned from a product is a "motor-shaped" multivector, not necessarily a valid motor satisfying all constraints
2. **new_unchecked is appropriate for products** - It represents "construct without validation because the computation is correct"
3. **Floating-point drift is handled by explicit normalization** - Users should call `normalize()` periodically, not rely on automatic constraint solving

## Alternative Approaches (For Future Consideration)

If constraint drift becomes a real problem:

1. **Periodic renormalization** - Have users call `motor.normalize()` after many operations
2. **Numerical stability improvements** - Use compensated summation or higher precision for accumulation
3. **Domain-specific products** - Create `geometric_motor_motor_normalized()` variants that explicitly normalize output
