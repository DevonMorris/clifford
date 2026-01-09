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

// All types use generic Arbitrary impls with Float::from_f64() for conversion
impl<T: Float + Debug> Arbitrary for Vector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate f64 values and convert to T
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| Vector::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
            .boxed()
    }
}

// Wrapper types use where clauses requiring the inner type to be Arbitrary
#[derive(Debug, Clone)]
struct UnitVector<T: Float>(Vector<T>);

impl<T> Arbitrary for UnitVector<T>
where
    T: Float + Debug,
    Vector<T>: Arbitrary + Debug,
    <Vector<T> as Arbitrary>::Strategy: 'static,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        let threshold = T::from_f64(1e-6);  // Use Float::from_f64 for conversion
        any::<Vector<T>>()
            .prop_filter("non-zero", move |v| v.norm_squared() > threshold)
            .prop_map(|v| UnitVector(v.normalized()))
            .boxed()
    }
}

// Usage: specify the float type explicitly
// any::<Vector<f64>>() or any::<Vector<f32>>()
// any::<UnitVector<f64>>() or any::<UnitVector<f32>>()
```

## Test Structure

```rust
use proptest::prelude::*;
use approx::abs_diff_eq;

proptest! {
    #[test]
    fn geometric_product_associative(
        a in any::<Multivector<f64, Euclidean3>>(),
        b in any::<Multivector<f64, Euclidean3>>(),
        c in any::<Multivector<f64, Euclidean3>>(),
    ) {
        let left = (a.clone() * b.clone()) * c.clone();
        let right = a * (b * c);
        prop_assert!(abs_diff_eq!(left, right, epsilon = 1e-10));
    }
}
```

## Use `prop_assert!` in Proptest Blocks

**Always use `prop_assert!` instead of `assert!`** inside `proptest!` blocks. This provides better error reporting with counterexamples and enables proptest's shrinking behavior.

```rust
proptest! {
    #[test]
    fn rotor_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
        let rotated = r.rotate(v);
        // Good: prop_assert! for better proptest integration
        prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));

        // Avoid: assert! loses proptest's shrinking and reporting benefits
        // assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Strategy Guidelines

- All `Arbitrary` impls are generic over `T: Float + Debug`
- Use `Float::from_f64()` to convert f64 range values to the target type
- Always specify the float type explicitly: `any::<Vector<f64>>()` not `any::<Vector>()`
- Internal tests should use f64 for consistency
- Wrapper types use where clauses requiring the inner type to be Arbitrary
- Chain with `.prop_filter()` and `.prop_map()` for derived strategies
- Use `any::<Type>()` in tests, not free functions

## Arbitrary Module Structure

Each module has an `arbitrary` submodule with Arbitrary implementations and wrapper types:

```
src/specialized/euclidean::dim3/
├── mod.rs
├── types.rs
├── ops.rs
└── arbitrary.rs  # Contains Arbitrary impls and wrapper types
```

### Key Points

1. **Feature gating**: Arbitrary modules use `#[cfg(any(test, feature = "proptest-support"))]`
2. **Import path**: Use `crate::specialized::euclidean::dim3::arbitrary::{NonZeroVector, UnitVector, UnitRotor}`
3. **All types are generic**: All Arbitrary impls use `impl<T: Float + Debug> Arbitrary for Type<T>`
4. **New types need Arbitrary**: When adding new types, add generic `impl Arbitrary` in the module's `arbitrary.rs`
5. **Use `Float::from_f64()`**: For converting f64 range values and threshold constants

### Available Types (all generic over Float)

| Module | Base Types | Wrapper Types |
|--------|-----------|---------------|
| `euclidean::dim2::arbitrary` | `Vector<T>` | `NonZeroVec2<T>`, `UnitVec2<T>`, `UnitRotor2<T>` |
| `euclidean::dim3::arbitrary` | `Vector<T>`, `Bivector<T>` | `NonZeroVector<T>`, `UnitVector<T>`, `UnitBivec3<T>`, `UnitRotor<T>` |
| `algebra::arbitrary` | `Multivector<T, S>` | `VectorE3`, `NonZeroVectorE3`, `UnitVectorE3` (f64 only) |

## Approximate Comparisons

**Always use the `approx` crate** for floating-point comparisons. Never hand-roll comparisons.

**Always use `ABS_DIFF_EQ_EPS` constant** instead of magic numbers. The constant is defined in `src/lib.rs::test_utils`:

```rust
use crate::test_utils::ABS_DIFF_EQ_EPS;
use approx::abs_diff_eq;

// Good: use approx macros with standard constant
prop_assert!(abs_diff_eq!(rotated.norm(), v.norm(), epsilon = ABS_DIFF_EQ_EPS));

// Avoid: magic numbers
prop_assert!(abs_diff_eq!(rotated.norm(), v.norm(), epsilon = 1e-9));

// Avoid: hand-rolled comparisons
prop_assert!((rotated.norm() - v.norm()).abs() < 1e-9);
```

For integration tests (`tests/` directory), define the constant locally since `test_utils` is `pub(crate)`.

Available comparison methods:
- `abs_diff_eq!` - absolute difference comparison
- `relative_eq!` - relative difference comparison
- `ulps_eq!` - units in last place comparison

All types (`Vector`, `Bivector`, `Rotor`, etc.) implement `AbsDiffEq`, `RelativeEq`, and `UlpsEq` for both f32 and f64.

## Wrapper Type Ergonomics

Wrapper types (`UnitVector<T>`, `NonZeroVector<T>`, etc.) implement:
- `Deref` - auto-dereference to inner type: `wrapper.method()` works directly
- `AsRef` - borrow as inner type: `wrapper.as_ref()`
- `From` - convert to inner type: `Vector::from(wrapper)` or `wrapper.into()`
- `into_inner()` - consume and return inner: `wrapper.into_inner()`

**Note**: Always specify the float type when using wrapper types: `any::<UnitVector<f64>>()`, not `any::<UnitVector>()`.

## Documentation

All test functions need doc comments explaining:
- What property is being tested
- Why this property matters mathematically
