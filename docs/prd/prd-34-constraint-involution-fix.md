# PRD-34: Fix Constraint System to Use Correct Involution

**Status**: Draft
**Goal**: Fix the constraint derivation system to only apply versor constraints to versor types, and use the algebra's configured involution where appropriate

## Problem Statement

The constraint derivation system currently has a fundamental issue: it always uses `Reverse` to compute constraints, regardless of the algebra's configured involution. This causes incorrect constraints to be inferred for types that don't use reverse.

### Example: Hyperbolic Numbers

Hyperbolic numbers (Cl(1,0,0)) use **grade involution** for their norm:
```toml
[norm]
primary_involution = "grade_involution"
```

However, the constraint system computes `u * reverse(u)` and infers a constraint:
```rust
// Generated for Hyperbolic type - WRONG!
pub fn new_checked(real: T, hyp: T, tolerance: T) -> Result<Self, &'static str> {
    let expected = (T::zero()) / (real);
    let actual = hyp;
    if (actual - expected).abs() > tolerance {
        return Err("Hyperbolic constraint");
    }
    Ok(Self::new_unchecked(real, hyp))
}
```

This constraint forces `hyp == 0`, which is completely wrong. All hyperbolic numbers `a + bj` should be valid.

### Root Cause

The constraint system was designed for **versors** - elements that satisfy `v * reverse(v) = scalar`. This is the defining property of versors (rotors, motors, etc.).

The issue is that the constraint derivation code (`derive_geometric_constraint`) unconditionally uses `reverse`:

```rust
// From constraint_derive.rs line 51-53:
/// Derives the geometric constraint for a type.
///
/// Computes `u * reverse(u)` symbolically and returns the non-scalar
/// terms that must equal zero.
```

This doesn't make sense for:
1. Non-versor types (like general multivectors)
2. Types in algebras that use a different involution for their norm
3. Simple algebras like 1D number systems (complex, hyperbolic, dual)

## Proposed Solution

### Understanding the Problem Better

The constraint system computes `u * reverse(u)` and checks if non-scalar terms must be zero. But this hardcodes **reverse** when the algebra may use a different involution.

For hyperbolic numbers using **grade involution**:
- `z * involute(z) = (a + bj)(a - bj) = a² - b²` → **scalar only** → no constraint needed!

But the current code uses reverse:
- `z * reverse(z) = (a + bj)(a + bj) = a² + 2abj + b²` → has grade-1 term → incorrectly derives `ab = 0`

### Solution: Use Algebra's Configured Involution (Recommended)

The constraint derivation should use the algebra's `primary_involution` instead of hardcoded reverse.

| Algebra | Involution | `u * involution(u)` | Constraint |
|---------|------------|---------------------|------------|
| Hyperbolic | grade_involution | `a² - b²` (scalar) | None |
| Complex | clifford_conjugate | `a² + b²` (scalar) | None |
| Dual | grade_involution | `a²` (scalar) | None |
| PGA | reverse (default) | Plücker terms | Preserved |
| Euclidean | reverse (default) | Versor terms | Preserved |

This is more principled than checking `dim <= 1` because:
1. It respects the algebra's mathematical structure
2. It would work correctly for any future algebra with non-reverse involution
3. It aligns constraint derivation with norm computation

## Implementation Plan

### Step 1: Update `ConstraintDeriver` to Accept Involution Kind

The `ConstraintDeriver` in `constraint_derive.rs` needs access to the algebra's involution:

```rust
pub struct ConstraintDeriver<'a> {
    algebra: &'a Algebra,
    table: ProductTable<'a>,
    involution: InvolutionKind,  // NEW: from NormSpec
}

impl<'a> ConstraintDeriver<'a> {
    pub fn new(algebra: &'a Algebra, involution: InvolutionKind) -> Self {
        Self {
            algebra,
            table: ProductTable::new(algebra),
            involution,
        }
    }
}
```

### Step 2: Update `compute_product_reverse_at_grade` to Use Involution

Rename to `compute_product_involution_at_grade` and use the configured involution:

```rust
fn compute_product_involution_at_grade(
    &self,
    ty: &TypeSpec,
    output_grade: usize,
    symbols: &HashMap<usize, Atom>,
) -> Atom {
    // Use self.involution instead of hardcoded reverse_sign
    let inv_sign = match self.involution {
        InvolutionKind::Reverse => reverse_sign(grade),
        InvolutionKind::GradeInvolution => grade_involution_sign(grade),
        InvolutionKind::CliffordConjugate => clifford_conjugate_sign(grade),
    };
    // ... rest of computation
}
```

### Step 3: Update Call Sites to Pass Involution

Update `TypeGenerator` and `TraitsGenerator` to pass the algebra's involution:

```rust
// In types.rs
let deriver = ConstraintDeriver::new(&self.algebra, self.spec.norm.primary_involution);
```

### Step 4: Add Involution Sign Helper Functions

Add helper functions for each involution type in `algebra/mod.rs`:

```rust
/// Grade involution sign: (-1)^k
pub fn grade_involution_sign(grade: usize) -> i8 {
    if grade % 2 == 0 { 1 } else { -1 }
}

/// Clifford conjugate sign: (-1)^(k(k+1)/2)
pub fn clifford_conjugate_sign(grade: usize) -> i8 {
    if (grade * (grade + 1) / 2) % 2 == 0 { 1 } else { -1 }
}
```

### Step 5: Regenerate All Algebras

Run `cargo build` to regenerate:
- Hyperbolic numbers: no constraint (grade_involution → scalar result)
- Complex numbers: no constraint (clifford_conjugate → scalar result)
- Dual numbers: no constraint (grade_involution → scalar result)
- Rotors/Motors: constraint preserved (reverse → versor condition)
- Lines in PGA: constraint preserved (reverse → Plücker constraint)

### Step 6: Add Tests

Verify that:
- Hyperbolic `Hyperbolic` type has no `new_checked` method
- Complex `Complex` type has no `new_checked` method
- Dual `Dual` type has no `new_checked` method
- `Motor::new_checked()` still validates the versor condition
- `Line::new_checked()` still validates the Plücker constraint

## Files to Modify

1. `crates/clifford-codegen/src/symbolic/constraint_derive.rs`
   - Add `involution: InvolutionKind` field to `ConstraintDeriver`
   - Update `new()` to accept involution parameter
   - Rename `compute_product_reverse_at_grade` → `compute_product_involution_at_grade`
   - Use involution-specific sign instead of hardcoded `reverse_sign`

2. `crates/clifford-codegen/src/algebra/mod.rs`
   - Add `grade_involution_sign()` function
   - Add `clifford_conjugate_sign()` function

3. `crates/clifford-codegen/src/codegen/types.rs`
   - Pass `spec.norm.primary_involution` when creating `ConstraintDeriver`

4. `crates/clifford-codegen/src/codegen/traits.rs`
   - Pass `spec.norm.primary_involution` when creating `ConstraintDeriver`

## Verification

1. `cargo build` - regenerate all algebras
2. `cargo nextest run` - all tests pass
3. Verify algebras with non-reverse involution have no constraints
4. Verify PGA types (Line, Motor) still have correct constraints

## Success Criteria

1. Hyperbolic numbers have no geometric constraint (grade_involution)
2. Complex numbers have no geometric constraint (clifford_conjugate)
3. Dual numbers have no geometric constraint (grade_involution)
4. Versors (Rotor, Motor) retain their constraint (reverse)
5. Lines retain Plücker constraint (reverse)
6. All tests pass
7. No clippy warnings

## Future Considerations

1. **Normed constraint validation**: In the future, we might want a separate constraint system for `Normed` types that validates `norm_squared() > 0` for unit normalization. This would be separate from the versor constraint.

2. **Custom constraints**: Some algebras might need custom constraints specified in TOML. This could be added later.

## References

- [Geometric Algebra - Versors](https://en.wikipedia.org/wiki/Geometric_algebra#Versors)
- PRD-31: Hyperbolic Numbers (introduced involution configuration)
- `crates/clifford-codegen/src/symbolic/constraint_derive.rs` - constraint derivation code
