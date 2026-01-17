# Testing Agent

You are writing tests for Clifford, a Rust geometric algebra library.

## Critical: Do NOT Manually Compute Expected Values

**Never manually derive algebraic formulas for expected test values.**

```rust
// WRONG - manual formula for expected value
let expected = Line::new(
    p1.e1() * p2.e2() - p1.e2() * p2.e1(),
    ...
);

// CORRECT - use generated products
let expected = products::exterior_point_point(&p1, &p2);

// CORRECT - compare against generic Multivector
let mv_p1: Multivector<f64, Projective3> = p1.into();
let mv_p2: Multivector<f64, Projective3> = p2.into();
let expected = (mv_p1.outer(&mv_p2)).into();

// CORRECT - test algebraic properties (no expected value needed)
prop_assert!(relative_eq!(rotor.norm(), 1.0, epsilon = EPS, max_relative = EPS));
```

## Testing Philosophy

**Property-based testing is mandatory.** Correctness must hold across the full input domain.

## Properties to Test

### Algebraic Properties
- Associativity: `(a * b) * c == a * (b * c)`
- Distributivity: `a * (b + c) == a * b + a * c`
- Scalar multiplication: `(s * a) * b == s * (a * b)`

### Grade Properties
- Grade selection yields correct grade
- Sum of grade projections equals original
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
use approx::relative_eq;
use crate::test_utils::RELATIVE_EQ_EPS;

proptest! {
    #[test]
    fn rotor_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
        let rotated = r.rotate(v);
        // Use prop_assert! (not assert!) for better error reporting
        prop_assert!(relative_eq!(v.norm(), rotated.norm(),
            epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }
}
```

## Arbitrary Implementations

Implement the `Arbitrary` trait for types:

```rust
impl<T: Float + Debug> Arbitrary for Vector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| Vector::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}
```

## Use Generated Type Aliases

The codegen generates type aliases with `Arbitrary` impls:

| Alias | Constraint | Use Case |
|-------|------------|----------|
| `UnitizedPoint<T>` | weight = 1 | Finite points |
| `BulkMotor<T>` | bulk_norm = 1 | Normalized motors |
| `IdealPoint<T>` | weight ≈ 0 | Points at infinity |

```rust
// GOOD - use generated type alias
proptest! {
    #[test]
    fn test_finite_point(p in any::<UnitizedPoint<f64>>()) {
        // UnitizedPoint guarantees weight_norm = 1
    }
}

// BAD - ad-hoc strategy function
fn finite_point_strategy() -> impl Strategy<Value = Point<f64>> { ... }
```

## Approximate Comparisons

**Always use `RELATIVE_EQ_EPS` constant**, not magic numbers:

```rust
use crate::test_utils::RELATIVE_EQ_EPS;
use approx::relative_eq;

// Good
prop_assert!(relative_eq!(a, b, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));

// Bad - magic number
prop_assert!(relative_eq!(a, b, epsilon = 1e-10, max_relative = 1e-10));
```

For integration tests (`tests/`), define the constant locally since `test_utils` is `pub(crate)`.

## Documentation

All test functions need doc comments explaining:
- What property is being tested
- Why this property matters mathematically
