# PRD-1: Foundation Layer

**Status**: Pending
**Goal**: Core traits and infrastructure

## Deliverables

### 1. Float Trait (`src/scalar/float.rs`)

```rust
pub trait Float:
    Copy + Clone + Default + PartialEq + PartialOrd + Debug
    + Add<Output = Self> + Sub<Output = Self>
    + Mul<Output = Self> + Div<Output = Self> + Neg<Output = Self>
    + AddAssign + SubAssign + MulAssign + DivAssign
    + 'static
{
    const ZERO: Self;
    const ONE: Self;
    const EPSILON: Self;

    fn abs(self) -> Self;
    fn sqrt(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn atan2(y: Self, x: Self) -> Self;
    fn approx_eq(self, other: Self, epsilon: Self) -> bool;
}
```

Implement for `f32` and `f64`.

### 2. Signature Trait (`src/signature/metric.rs`)

```rust
pub trait Signature: Copy + Clone + Default + 'static {
    /// Basis vectors squaring to +1
    const P: usize;
    /// Basis vectors squaring to -1
    const Q: usize;
    /// Basis vectors squaring to 0 (degenerate)
    const R: usize;
    /// Total dimension: P + Q + R
    const DIM: usize;
    /// Total basis blades: 2^DIM
    const NUM_BLADES: usize;

    /// Returns +1, -1, or 0 for basis vector i
    fn metric(i: usize) -> i8;
}
```

### 3. Euclidean Signatures (`src/signature/euclidean.rs`)

```rust
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Euclidean2;

impl Signature for Euclidean2 {
    const P: usize = 2;
    const Q: usize = 0;
    const R: usize = 0;
    const DIM: usize = 2;
    const NUM_BLADES: usize = 4;

    fn metric(_i: usize) -> i8 { 1 }
}

// Similarly: Euclidean3, Euclidean4
```

### 4. Grade Utilities (`src/basis/index.rs`)

```rust
/// Returns the grade (number of basis vectors) of a blade index.
/// Grade = number of 1-bits in binary representation.
#[inline]
pub const fn grade_of_blade(index: usize) -> usize {
    index.count_ones() as usize
}

/// Returns the number of blades of grade k in dimension n.
/// This is the binomial coefficient C(n, k).
pub const fn blades_of_grade(dim: usize, grade: usize) -> usize {
    // binomial(dim, grade)
}
```

### 5. Basis Product (`src/basis/blade.rs`)

```rust
/// Computes the sign and result index when multiplying basis blades.
///
/// # Arguments
/// * `a` - Index of first basis blade
/// * `b` - Index of second basis blade
/// * `metric` - Function returning +1, -1, or 0 for each basis vector
///
/// # Returns
/// `(sign, result_index)` where sign is -1, 0, or 1
pub fn basis_product<F>(a: usize, b: usize, metric: F) -> (i8, usize)
where
    F: Fn(usize) -> i8,
{
    let result_index = a ^ b;
    let sign = compute_sign(a, b, metric);
    (sign, result_index)
}

/// Count swaps needed to reorder and apply metric.
fn compute_sign<F>(a: usize, b: usize, metric: F) -> i8
where
    F: Fn(usize) -> i8,
{
    // Implementation using bit manipulation
}
```

## Files to Create

- `src/scalar/mod.rs`
- `src/scalar/float.rs`
- `src/signature/mod.rs`
- `src/signature/metric.rs`
- `src/signature/euclidean.rs`
- `src/basis/mod.rs`
- `src/basis/index.rs`
- `src/basis/blade.rs`

## Testing (proptest)

```rust
proptest! {
    #[test]
    fn grade_is_popcount(index in 0usize..64) {
        prop_assert_eq!(grade_of_blade(index), index.count_ones() as usize);
    }

    #[test]
    fn basis_product_result_is_xor(a in 0usize..16, b in 0usize..16) {
        let (_, result) = basis_product(a, b, |_| 1);
        prop_assert_eq!(result, a ^ b);
    }

    #[test]
    fn euclidean_metric_is_positive(i in 0usize..4) {
        prop_assert_eq!(Euclidean3::metric(i.min(2)), 1);
    }
}
```

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - proptest properties pass
- [ ] `cargo clippy` - no warnings
- [ ] All items have rustdoc
