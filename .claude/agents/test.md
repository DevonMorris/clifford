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

## Implementing Arbitrary

**Always implement the `Arbitrary` trait** for types instead of free functions. This enables using `any::<Type>()` which is cleaner and more idiomatic.

```rust
use proptest::prelude::*;
use proptest::arbitrary::{Arbitrary, StrategyFor};
use proptest::strategy::{BoxedStrategy, Strategy};

impl Arbitrary for Vec3<f64> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-1e10..1e10, -1e10..1e10, -1e10..1e10)
            .prop_map(|(x, y, z)| Vec3::new(x, y, z))
            .boxed()
    }
}

// For constrained variants, create wrapper types:
#[derive(Debug, Clone)]
struct UnitVec3(Vec3<f64>);

impl Arbitrary for UnitVec3 {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Vec3<f64>>()
            .prop_filter("non-zero", |v| v.norm_squared() > 1e-10)
            .prop_map(|v| UnitVec3(v.normalized()))
            .boxed()
    }
}
```

## Test Structure

```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn geometric_product_associative(
        a in any::<Multivector<f64, Euclidean3>>(),
        b in any::<Multivector<f64, Euclidean3>>(),
        c in any::<Multivector<f64, Euclidean3>>(),
    ) {
        let left = (a.clone() * b.clone()) * c.clone();
        let right = a * (b * c);
        prop_assert!((left - right).norm() < 1e-10);
    }
}
```

## Strategy Guidelines

- Implement `Arbitrary` for base types (`Vec3<f64>`, `Bivec3<f64>`, etc.)
- Use wrapper types for constrained values (`UnitVec3`, `NonZeroVec3`)
- Chain with `.prop_filter()` and `.prop_map()` for derived strategies
- Use `any::<Type>()` in tests, not free functions

## Arbitrary Module Structure

Each module has an `arbitrary` submodule with Arbitrary implementations and wrapper types:

```
src/specialized/ga3d/
├── mod.rs
├── types.rs
├── ops.rs
└── arbitrary.rs  # Contains Arbitrary impls and wrapper types
```

### Key Points

1. **Feature gating**: Arbitrary modules use `#[cfg(any(test, feature = "proptest-support"))]`
2. **Import path**: Use `crate::specialized::ga3d::arbitrary::{NonZeroVec3, UnitVec3, UnitRotor3}`
3. **Wrapper types are public**: All wrapper types (NonZeroVec3, UnitVec3, etc.) have `pub` visibility
4. **New types need Arbitrary**: When adding new types, add `impl Arbitrary` in the module's `arbitrary.rs`

### Available Wrapper Types

| Module | Types |
|--------|-------|
| `ga2d::arbitrary` | `NonZeroVec2`, `UnitVec2`, `UnitRotor2` |
| `ga3d::arbitrary` | `NonZeroVec3`, `UnitVec3`, `UnitBivec3`, `UnitRotor3` |
| `algebra::arbitrary` | `VectorE3`, `NonZeroVectorE3`, `UnitVectorE3` |

## Approximate Comparisons

**Always use the `approx` crate** for floating-point comparisons. Never hand-roll comparisons.

```rust
use approx::abs_diff_eq;

// Good: use approx macros
prop_assert!(abs_diff_eq!(rotated.norm(), v.norm(), epsilon = 1e-9));

// Avoid: hand-rolled comparisons
prop_assert!((rotated.norm() - v.norm()).abs() < 1e-9);
```

Available comparison methods:
- `abs_diff_eq!` - absolute difference comparison
- `relative_eq!` - relative difference comparison
- `ulps_eq!` - units in last place comparison

All types (`Vec3`, `Bivec3`, `Rotor3`, etc.) implement `AbsDiffEq`, `RelativeEq`, and `UlpsEq` for both f32 and f64.

## Wrapper Type Ergonomics

Wrapper types (`UnitVec3`, `NonZeroVec3`, etc.) implement:
- `Deref` - auto-dereference to inner type: `wrapper.method()` works directly
- `AsRef` - borrow as inner type: `wrapper.as_ref()`
- `From` - convert to inner type: `Vec3::from(wrapper)` or `wrapper.into()`
- `into_inner()` - consume and return inner: `wrapper.into_inner()`

## Documentation

All test functions need doc comments explaining:
- What property is being tested
- Why this property matters mathematically
