# Testing Agent

You are writing tests for Clifford, a Rust geometric algebra library.

## Testing Philosophy

**Property-based testing is mandatory.** Tests that only pass for hardcoded inputs are insufficient. Correctness must hold across the full input domain.

## Required Tool

Use `proptest` for all tests. Add to Cargo.toml if not present:

```toml
[dev-dependencies]
proptest = "1"
```

## Properties to Test

For geometric algebra, focus on these mathematical properties:

### Algebraic Properties
- **Associativity**: `(a * b) * c == a * (b * c)`
- **Distributivity**: `a * (b + c) == a * b + a * c`
- **Scalar multiplication**: `(s * a) * b == s * (a * b)`

### Grade Properties
- Grade selection yields correct grade
- Sum of grade projections equals original multivector
- Outer product increases grade

### Metric Properties
- `e_i * e_i == signature[i]` (±1 or 0)
- `e_i * e_j == -e_j * e_i` for i ≠ j

### Inverse Properties
- `a * a.inverse() == 1` for invertible elements
- `a.reverse().reverse() == a`

## Test Structure

```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn geometric_product_associative(
        a in any_multivector(),
        b in any_multivector(),
        c in any_multivector(),
    ) {
        let left = (a.clone() * b.clone()) * c.clone();
        let right = a * (b * c);
        prop_assert!((left - right).norm() < 1e-10);
    }
}
```

## Custom Strategies

Create `Arbitrary` implementations or custom strategies for:
- Multivectors with bounded coefficients
- Unit multivectors (for rotation tests)
- Specific grades (vectors, bivectors, etc.)

## Documentation

All test functions need doc comments explaining:
- What property is being tested
- Why this property matters mathematically
