# PRD-33: Dual Numbers Algebra

**Status**: Draft
**Goal**: Generate dual numbers algebra Cl(0,0,1) via code generation

## Background

Dual numbers are a 2-dimensional algebra over the reals with basis `{1, ε}` where `ε² = 0`. They arise naturally as the Clifford algebra Cl(0,0,1) with one degenerate (null) basis vector.

### Primary Application: Automatic Differentiation

Dual numbers are the foundation of **forward-mode automatic differentiation**. When computing `f(a + bε)`, the result is `f(a) + b·f'(a)·ε`. This means:

```
f(x + ε) = f(x) + f'(x)·ε
```

The dual part automatically carries the derivative!

### Comparison with Other 1D Clifford Algebras

| Property | Complex Cl(0,1,0) | Hyperbolic Cl(1,0,0) | Dual Cl(0,0,1) |
|----------|-------------------|----------------------|----------------|
| Signature | (0,1,0) | (1,0,0) | (0,0,1) |
| Unit squares to | `i² = -1` | `j² = +1` | `ε² = 0` |
| Norm formula | `a² + b²` | `a² - b²` | `a²` |
| Zero divisors | None | Yes: `1±j` | Yes: all `bε` |
| Division algebra | Yes | No | No |

## Mathematical Definitions

### Algebra Structure

The Clifford algebra Cl(0,0,1) has:
- **Dimension**: 1 (one basis vector)
- **Total blades**: 2 (scalar and the unit `ε`)
- **Metric**: `ε² = 0` (degenerate/null)

### Elements

| Grade | Blade | Description |
|-------|-------|-------------|
| 0 | `1` | Scalar (real part) |
| 1 | `ε` | Dual unit (nilpotent) |

A general dual number is: `z = a + bε` where `a, b ∈ ℝ`.

### Operations

**Addition**: `(a + bε) + (c + dε) = (a+c) + (b+d)ε`

**Multiplication** (geometric product):
```
(a + bε)(c + dε) = ac + adε + bcε + bdε²
                 = ac + adε + bcε + 0
                 = ac + (ad + bc)ε
```

Note: The product rule for dual numbers mirrors the product rule for derivatives!

**Involutions**:

For dual numbers, all three involutions (reverse, grade involution, Clifford conjugate) give the same result for grades 0 and 1:
- Grade 0: unchanged
- Grade 1: negated (for grade involution and Clifford conjugate) or unchanged (for reverse)

Actually, let's check:
| Involution | Grade 0 | Grade 1 | Formula |
|------------|---------|---------|---------|
| Reverse | +1 | +1 | `(-1)^(k(k-1)/2)` |
| Grade involution | +1 | -1 | `(-1)^k` |
| Clifford conjugate | +1 | -1 | `(-1)^(k(k+1)/2)` |

For the "conjugate" `a - bε`, we need grade involution or Clifford conjugate.

**Norm** (using any conjugating involution):
```
z * conjugate(z) = (a + bε)(a - bε)
                 = a² - abε + abε - b²ε²
                 = a² - b²·0
                 = a²
```

The norm of a dual number is simply the square of its real part! This is a **degenerate norm** - all pure dual numbers `bε` have norm zero.

## API Design

### Types

```rust
/// Scalar (grade-0 element)
pub struct Scalar<T> {
    pub s: T,
}

/// Pure dual unit (ε where ε² = 0)
pub struct DualUnit<T> {
    pub eps: T,  // or "d" for dual
}

/// Dual number: a + bε where ε² = 0
pub struct Dual<T> {
    pub real: T,
    pub dual: T,
}
```

### Traits

The dual algebra will implement:
- `Add`, `Sub`, `Mul`, `Neg` - basic arithmetic
- `GeometricProduct` - same as `Mul` for this algebra
- `Reverse` - identity (doesn't change grade-1 in this case)
- `Involute` - dual conjugate (negates grade-1)
- `Normed` - degenerate norm: `norm_squared() = real²`

### Norm Configuration

```toml
[norm]
primary_involution = "grade_involution"
```

Note: For Cl(0,0,1), the metric contribution is 0, so `z * involute(z) = a² - b²·0 = a²`.

## Implementation Plan

### Step 1: Create `algebras/dual.toml`

```toml
# Dual Numbers
# Cl(0,0,1) - foundation for automatic differentiation

[algebra]
name = "dual"
module_path = "dual"
description = "Dual numbers Cl(0,0,1)"

[signature]
positive = []
negative = []
zero = ["e1"]

[norm]
primary_involution = "grade_involution"

[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"

[types.DualUnit]
grades = [1]
description = "Pure dual unit (ε where ε² = 0)"
fields = ["eps"]

[types.Dual]
grades = [0, 1]
description = "Dual number: a + bε where ε² = 0"
fields = ["real", "dual"]
```

### Step 2: Add to `build.rs`

Add `"dual"` to the `ALGEBRAS` array.

### Step 3: Create Module Structure

Create `src/specialized/dual/mod.rs` with documentation about automatic differentiation.

### Step 4: Add Signature Type

Add `Cl0_0_1` signature to `src/signature/` with metric returning 0.

### Step 5: Regenerate

Run `cargo build` to generate the dual algebra.

### Step 6: Add Tests

Property-based tests for:
- `ε * ε = 0` (nilpotent)
- Multiplication follows product rule pattern
- `norm_squared()` equals `real²`
- Pure dual numbers have zero norm

## Testing Strategy

### Unit Tests

```rust
#[test]
fn dual_unit_is_nilpotent() {
    // ε² = 0
    let eps = DualUnit::new(1.0);
    let product = eps.geometric_product(&eps);
    assert!(abs_diff_eq!(product.s(), 0.0, epsilon = 1e-10));
}

#[test]
fn dual_multiplication_product_rule() {
    // (a + bε)(c + dε) = ac + (ad + bc)ε
    let z1 = Dual::new(2.0, 3.0);  // 2 + 3ε
    let z2 = Dual::new(5.0, 7.0);  // 5 + 7ε
    let product = z1.geometric_product(&z2);

    // real part: 2 * 5 = 10
    assert!(abs_diff_eq!(product.real(), 10.0, epsilon = 1e-10));
    // dual part: 2*7 + 3*5 = 14 + 15 = 29
    assert!(abs_diff_eq!(product.dual(), 29.0, epsilon = 1e-10));
}

#[test]
fn dual_norm_is_degenerate() {
    // |a + bε|² = a²
    let z = Dual::new(3.0, 100.0);  // dual part doesn't matter
    assert!(abs_diff_eq!(z.norm_squared(), 9.0, epsilon = 1e-10));
}
```

### Property-Based Tests

```rust
proptest! {
    #[test]
    fn dual_unit_squared_is_zero(a in any::<f64>()) {
        let eps = DualUnit::new(a);
        let product = eps.geometric_product(&eps);
        prop_assert!(abs_diff_eq!(product.s(), 0.0, epsilon = 1e-10));
    }

    #[test]
    fn norm_ignores_dual_part(
        real in any::<f64>(),
        dual in any::<f64>(),
    ) {
        let z = Dual::new(real, dual);
        prop_assert!(relative_eq!(z.norm_squared(), real * real, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn multiplication_is_commutative(
        z1 in any::<Dual<f64>>(),
        z2 in any::<Dual<f64>>(),
    ) {
        let p1 = z1.geometric_product(&z2);
        let p2 = z2.geometric_product(&z1);
        prop_assert!(relative_eq!(p1, p2, epsilon = RELATIVE_EQ_EPS));
    }
}
```

### Automatic Differentiation Test

```rust
#[test]
fn autodiff_polynomial() {
    // f(x) = x² at x = 3
    // f(3 + ε) = (3 + ε)² = 9 + 6ε
    // So f(3) = 9, f'(3) = 6
    let x = Dual::new(3.0, 1.0);  // x + ε represents "x with derivative 1"
    let f_x = x.geometric_product(&x);  // x²

    assert!(abs_diff_eq!(f_x.real(), 9.0, epsilon = 1e-10));  // f(3) = 9
    assert!(abs_diff_eq!(f_x.dual(), 6.0, epsilon = 1e-10));  // f'(3) = 6
}
```

## Dependencies

- PRD-31 (Hyperbolic Numbers): The involution system with configurable norm

## Success Criteria

1. `algebras/dual.toml` created with correct specification
2. Code generation produces correct `Dual`, `Scalar`, `DualUnit` types
3. `ε * ε = 0` (nilpotent property)
4. `Normed` trait computes `norm_squared() = real²`
5. Multiplication matches dual number arithmetic
6. All verification passes: `cargo fmt && cargo clippy && cargo nextest run && cargo doc --no-deps`

## Future Enhancements

1. **Autodiff utilities**: Helper functions for computing derivatives
2. **Higher-order duals**: Nested dual numbers for higher derivatives
3. **Interop with autodiff crates**: Conversion traits to/from existing autodiff libraries
4. **Jet types**: Dual numbers with vector-valued dual parts for gradients

## References

- [Wikipedia - Dual Numbers](https://en.wikipedia.org/wiki/Dual_number)
- [Wikipedia - Automatic Differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)
- [Dual Numbers & Automatic Differentiation](https://blog.demofox.org/2014/12/30/dual-numbers-automatic-differentiation/)
