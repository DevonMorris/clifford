# PRD-2: Core Multivector

**Status**: Pending
**Goal**: Generic multivector type with geometric product

## Deliverables

### 1. Multivector Type (`src/algebra/multivector.rs`)

```rust
use crate::{Float, Signature};
use std::marker::PhantomData;

pub const MAX_BLADES: usize = 64;

/// A multivector in a Clifford algebra with signature S.
///
/// # Type Parameters
/// - `T`: Scalar type (f32, f64)
/// - `S`: Metric signature
///
/// # Storage
/// Dense array of 2^N coefficients in canonical order.
#[derive(Clone, Debug)]
pub struct Multivector<T: Float, S: Signature> {
    coeffs: [T; MAX_BLADES],
    _signature: PhantomData<S>,
}
```

### 2. Constructors

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Zero multivector
    pub fn zero() -> Self;

    /// Scalar (grade 0)
    pub fn one() -> Self;
    pub fn scalar(value: T) -> Self;

    /// Basis blade by index
    pub fn basis_blade(index: usize) -> Self;

    /// Basis vector e_i (grade 1)
    pub fn basis_vector(i: usize) -> Self;

    /// Vector from components
    pub fn vector(components: &[T]) -> Self;
}
```

### 3. Accessors

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Get coefficient at index
    pub fn get(&self, index: usize) -> T;

    /// Set coefficient at index
    pub fn set(&mut self, index: usize, value: T);

    /// Scalar part (grade 0)
    pub fn scalar_part(&self) -> T;

    /// Check if approximately zero
    pub fn is_zero(&self, epsilon: T) -> bool;
}
```

### 4. Geometric Product (`src/algebra/ops/geometric.rs`)

```rust
impl<T: Float, S: Signature> Mul for Multivector<T, S> {
    type Output = Self;

    /// Geometric product of two multivectors.
    ///
    /// # Complexity
    /// O(4^n) where n is dimension
    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = Self::zero();
        for i in 0..S::NUM_BLADES {
            for j in 0..S::NUM_BLADES {
                let (sign, k) = basis_product(i, j, S::metric);
                if sign != 0 {
                    result.coeffs[k] += T::from_i8(sign) * self.coeffs[i] * rhs.coeffs[j];
                }
            }
        }
        result
    }
}
```

### 5. Unary Operations (`src/algebra/unary.rs`)

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Reverse: reverses order of basis vectors in each blade.
    /// For grade k: sign = (-1)^(k(k-1)/2)
    pub fn reverse(&self) -> Self;

    /// Grade involution: negates odd-grade parts.
    /// For grade k: sign = (-1)^k
    pub fn involute(&self) -> Self;

    /// Clifford conjugate: reverse + involute
    pub fn conjugate(&self) -> Self;
}
```

### 6. Norms (`src/algebra/norms.rs`)

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Squared norm: scalar part of self * self.reverse()
    pub fn norm_squared(&self) -> T;

    /// Norm: sqrt(|norm_squared|)
    pub fn norm(&self) -> T;

    /// Normalize to unit norm
    pub fn normalize(&self) -> Option<Self>;

    /// Multiplicative inverse: reverse / norm_squared
    pub fn inverse(&self) -> Option<Self>;
}
```

## Files to Create

- `src/algebra/mod.rs`
- `src/algebra/multivector.rs`
- `src/algebra/ops/mod.rs`
- `src/algebra/ops/geometric.rs`
- `src/algebra/unary.rs`
- `src/algebra/norms.rs`

## Testing (proptest)

```rust
proptest! {
    #[test]
    fn geometric_product_associative(
        a in arb_multivector::<f64, Euclidean3>(),
        b in arb_multivector::<f64, Euclidean3>(),
        c in arb_multivector::<f64, Euclidean3>(),
    ) {
        let lhs = (a.clone() * b.clone()) * c.clone();
        let rhs = a * (b * c);
        prop_assert!(lhs.approx_eq(&rhs, 1e-10));
    }

    #[test]
    fn geometric_product_distributive(
        a in arb_multivector::<f64, Euclidean3>(),
        b in arb_multivector::<f64, Euclidean3>(),
        c in arb_multivector::<f64, Euclidean3>(),
    ) {
        let lhs = a.clone() * (b.clone() + c.clone());
        let rhs = a.clone() * b + a * c;
        prop_assert!(lhs.approx_eq(&rhs, 1e-10));
    }

    #[test]
    fn reverse_involutory(a in arb_multivector::<f64, Euclidean3>()) {
        prop_assert!(a.reverse().reverse().approx_eq(&a, 1e-10));
    }

    #[test]
    fn reverse_antimorphism(
        a in arb_multivector::<f64, Euclidean3>(),
        b in arb_multivector::<f64, Euclidean3>(),
    ) {
        let lhs = (a.clone() * b.clone()).reverse();
        let rhs = b.reverse() * a.reverse();
        prop_assert!(lhs.approx_eq(&rhs, 1e-10));
    }

    #[test]
    fn inverse_property(a in arb_invertible::<f64, Euclidean3>()) {
        if let Some(inv) = a.inverse() {
            let product = a.clone() * inv;
            prop_assert!(product.approx_eq(&Multivector::one(), 1e-10));
        }
    }
}
```

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - all proptest properties pass
- [ ] `cargo clippy` - no warnings
- [ ] All items have rustdoc with mathematical explanations
