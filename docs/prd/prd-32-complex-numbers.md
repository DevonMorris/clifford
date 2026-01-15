# PRD-32: Complex Numbers Algebra

**Status**: Draft
**Goal**: Generate complex numbers algebra Cl(0,1,0) via code generation

## Background

Complex numbers are one of the most fundamental mathematical structures, and they naturally arise as the Clifford algebra Cl(0,1,0). This algebra has:

- One basis vector `e₁` with `e₁² = -1`
- Elements of the form `a + bi` where `i² = -1`
- A positive-definite norm: `|z|² = a² + b²`

### Comparison with Hyperbolic Numbers

| Property | Complex Cl(0,1,0) | Hyperbolic Cl(1,0,0) |
|----------|-------------------|----------------------|
| Signature | (0,1,0) | (1,0,0) |
| Unit squares to | `i² = -1` | `j² = +1` |
| Norm formula | `z * reverse(z) = a² + b²` | `z * involute(z) = a² - b²` |
| Norm type | Positive-definite | Indefinite |
| Zero divisors | None | Yes: `(1+j)(1-j) = 0` |
| Involution for norm | Reverse | Grade involution |

### Why Include Complex Numbers?

1. **Completeness**: Demonstrates the codegen system handles both positive and negative metric signatures
2. **Educational**: Shows how Cl(0,1,0) naturally produces complex numbers
3. **Verification**: Complex number properties are well-known, providing excellent test cases
4. **Foundation**: Required for understanding higher algebras like quaternions Cl(0,2,0)

## Mathematical Definitions

### Algebra Structure

The Clifford algebra Cl(0,1,0) has:
- **Dimension**: 1 (one basis vector)
- **Total blades**: 2 (scalar and the unit `e₁`)
- **Metric**: `e₁² = -1`

### Elements

| Grade | Blade | Description |
|-------|-------|-------------|
| 0 | `1` | Scalar (real part) |
| 1 | `e₁` | Imaginary unit `i` |

A general complex number is: `z = a + bi` where `a, b ∈ ℝ`.

### Operations

**Addition**: `(a + bi) + (c + di) = (a+c) + (b+d)i`

**Multiplication** (geometric product):
```
(a + bi)(c + di) = ac + adi + bci + bdi²
                 = ac + adi + bci - bd
                 = (ac - bd) + (ad + bc)i
```

**Involutions**:

For complex conjugation, we need an involution that:
- Leaves grade-0 (real part) unchanged
- Negates grade-1 (imaginary part)

| Involution | Grade 0 | Grade 1 | Formula |
|------------|---------|---------|---------|
| Reverse | +1 | +1 | `(-1)^(k(k-1)/2)` |
| Grade involution | +1 | -1 | `(-1)^k` |
| Clifford conjugate | +1 | -1 | `(-1)^(k(k+1)/2)` |

Both **grade involution** and **Clifford conjugate** give complex conjugation for this algebra.
We use Clifford conjugate as it's the conventional choice for "conjugation" operations.

**Conjugate**: `conjugate(a + bi) = a - bi`

**Norm** (using Clifford conjugate):
```
z * conjugate(z) = (a + bi)(a - bi)
                 = a² - abi + abi - b²i²
                 = a² + b²
```

## API Design

### Types

```rust
/// Scalar (grade-0 element)
pub struct Scalar<T> {
    pub s: T,
}

/// Pure imaginary unit (i where i² = -1)
pub struct ImagUnit<T> {
    pub i: T,
}

/// Complex number: a + bi where i² = -1
pub struct Complex<T> {
    pub real: T,
    pub imag: T,
}
```

### Traits

The complex algebra will implement:
- `Add`, `Sub`, `Mul`, `Neg` - basic arithmetic
- `GeometricProduct` - same as `Mul` for this algebra
- `Reverse` - identity for this algebra (doesn't flip grade-1)
- `Involute` - complex conjugate (flips grade-1)
- `Normed` - using Clifford conjugate: `norm_squared() = real² + imag²`

### Norm Configuration

The TOML will specify Clifford conjugate as the involution for norm:

```toml
[norm]
primary_involution = "clifford_conjugate"
```

This makes `norm_squared()` compute `z * clifford_conjugate(z) = a² + b²`.

## Implementation Plan

### Step 1: Create `algebras/complex.toml`

```toml
# Complex Numbers
# Cl(0,1,0) with Clifford conjugate for norm

[algebra]
name = "complex"
module_path = "complex"
description = "Complex numbers Cl(0,1,0)"

[signature]
positive = []
negative = ["e1"]
zero = []

[norm]
primary_involution = "clifford_conjugate"

[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"

[types.ImagUnit]
grades = [1]
description = "Pure imaginary unit (i where i² = -1)"
fields = ["i"]

[types.Complex]
grades = [0, 1]
description = "Complex number: a + bi where i² = -1"
fields = ["real", "imag"]
```

### Step 2: Add to `build.rs`

Add `"complex"` to the `ALGEBRAS` array.

### Step 3: Create Module Structure

Create `src/specialized/complex/mod.rs` with documentation.

### Step 4: Add Signature Type

Add `Cl0_1_0` signature to `src/signature/` with appropriate metric.

### Step 5: Regenerate

Run `cargo build` to generate the complex algebra.

### Step 6: Add Tests

Property-based tests for:
- Multiplication matches expected complex multiplication
- `norm_squared()` equals `real² + imag²`
- `norm_squared()` is always non-negative
- Multiplication is commutative (complex numbers are commutative)

## Testing Strategy

### Unit Tests

```rust
#[test]
fn complex_multiplication() {
    // (1 + 2i) * (3 + 4i) = 3 + 4i + 6i + 8i² = 3 + 10i - 8 = -5 + 10i
    let z1 = Complex::new(1.0, 2.0);
    let z2 = Complex::new(3.0, 4.0);
    let product = z1.geometric_product(&z2);

    assert!(abs_diff_eq!(product.real, -5.0, epsilon = 1e-10));
    assert!(abs_diff_eq!(product.imag, 10.0, epsilon = 1e-10));
}

#[test]
fn complex_norm() {
    // |3 + 4i|² = 9 + 16 = 25
    let z = Complex::new(3.0, 4.0);
    assert!(abs_diff_eq!(z.norm_squared(), 25.0, epsilon = 1e-10));
}
```

### Property-Based Tests

```rust
proptest! {
    #[test]
    fn norm_is_non_negative(z in any::<Complex<f64>>()) {
        prop_assert!(z.norm_squared() >= 0.0);
    }

    #[test]
    fn multiplication_is_commutative(
        z1 in any::<Complex<f64>>(),
        z2 in any::<Complex<f64>>(),
    ) {
        let p1 = z1.geometric_product(&z2);
        let p2 = z2.geometric_product(&z1);
        prop_assert!(relative_eq!(p1, p2, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn norm_of_product(
        z1 in any::<Complex<f64>>(),
        z2 in any::<Complex<f64>>(),
    ) {
        // |z1 * z2|² = |z1|² * |z2|²
        let product = z1.geometric_product(&z2);
        let expected = z1.norm_squared() * z2.norm_squared();
        prop_assert!(relative_eq!(product.norm_squared(), expected, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }
}
```

## Dependencies

- PRD-31 (Hyperbolic Numbers): The involution system with `InvolutionKind::CliffordConjugate` support

## Success Criteria

1. `algebras/complex.toml` created with correct specification
2. Code generation produces correct `Complex`, `Scalar`, `ImagUnit` types
3. `Normed` trait uses Clifford conjugate: `norm_squared() = real² + imag²`
4. All property tests pass, especially `norm_of_product`
5. Multiplication matches standard complex multiplication
6. All verification passes: `cargo fmt && cargo clippy && cargo nextest run && cargo doc --no-deps`

## Future Enhancements

1. **Quaternions**: Cl(0,2,0) - four-dimensional division algebra
2. **Dual numbers**: Cl(0,0,1) - automatic differentiation
3. **Interop with `num-complex`**: Conversion traits to/from standard complex crate

## References

- [Wikipedia - Complex Numbers](https://en.wikipedia.org/wiki/Complex_number)
- [Wikipedia - Clifford Algebra](https://en.wikipedia.org/wiki/Clifford_algebra)
- [Geometric Algebra Primer](http://www.jaapsuter.com/geometric-algebra.pdf)
