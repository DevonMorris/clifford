# PRD-3: Derived Products & Grade Operations

**Status**: Pending
**Goal**: Complete product suite and grade operations

## Deliverables

### 1. Outer Product (`src/algebra/ops/outer.rs`)

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Outer (wedge) product.
    ///
    /// # Mathematical Definition
    /// a ∧ b extracts the grade-raising part of the geometric product.
    /// For a k-blade and j-blade: grade(a ∧ b) = k + j
    ///
    /// # Properties
    /// - Anticommutative for vectors: a ∧ b = -(b ∧ a)
    /// - Associative: (a ∧ b) ∧ c = a ∧ (b ∧ c)
    pub fn outer(&self, other: &Self) -> Self;
}
```

### 2. Inner Product (`src/algebra/ops/inner.rs`)

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Inner (dot) product - left contraction.
    ///
    /// # Mathematical Definition
    /// a · b extracts the grade-lowering part.
    /// For k-blade A and j-blade B (k ≤ j): grade(A · B) = j - k
    pub fn inner(&self, other: &Self) -> Self;

    /// Left contraction (⌋)
    pub fn left_contract(&self, other: &Self) -> Self;

    /// Right contraction (⌊)
    pub fn right_contract(&self, other: &Self) -> Self;
}
```

### 3. Regressive Product (`src/algebra/ops/regressive.rs`)

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Regressive (vee) product.
    ///
    /// # Mathematical Definition
    /// a ∨ b = (a* ∧ b*)* where * is the dual
    ///
    /// The regressive product is grade-lowering and represents
    /// the "meet" of geometric objects in PGA.
    pub fn regressive(&self, other: &Self) -> Self;
}
```

### 4. Sandwich Product (`src/algebra/ops/sandwich.rs`)

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Sandwich product: self * x * self.reverse()
    ///
    /// Used for transformations:
    /// - Rotation: R * v * R̃
    /// - Reflection: n * v * n
    pub fn sandwich(&self, x: &Self) -> Self;
}
```

### 5. Grade Operations (`src/basis/grade.rs`)

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Extract grade-k part: ⟨M⟩_k
    pub fn grade_select(&self, k: usize) -> Self;

    /// Extract multiple grades
    pub fn grade_select_many(&self, grades: &[usize]) -> Self;

    /// Even part (grades 0, 2, 4, ...)
    pub fn even(&self) -> Self;

    /// Odd part (grades 1, 3, 5, ...)
    pub fn odd(&self) -> Self;

    /// Decompose into list of homogeneous parts
    pub fn decompose(&self) -> Vec<Self>;

    /// Get the grade of this multivector (if homogeneous)
    pub fn grade(&self) -> Option<usize>;
}
```

### 6. Dual Operations

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Hodge dual: A* = A · I⁻¹
    pub fn dual(&self) -> Self;

    /// Undual: inverse of dual
    pub fn undual(&self) -> Self;

    /// Unit pseudoscalar I = e₁e₂...eₙ
    pub fn pseudoscalar() -> Self;
}
```

## Files to Create

- `src/algebra/ops/outer.rs`
- `src/algebra/ops/inner.rs`
- `src/algebra/ops/regressive.rs`
- `src/algebra/ops/sandwich.rs`
- `src/basis/grade.rs`

## Testing (proptest)

```rust
proptest! {
    #[test]
    fn wedge_anticommutative_vectors(
        a in arb_vector::<f64, Euclidean3>(),
        b in arb_vector::<f64, Euclidean3>(),
    ) {
        let lhs = a.outer(&b);
        let rhs = b.outer(&a).neg();
        prop_assert!(lhs.approx_eq(&rhs, 1e-10));
    }

    #[test]
    fn wedge_associative(
        a in arb_multivector::<f64, Euclidean3>(),
        b in arb_multivector::<f64, Euclidean3>(),
        c in arb_multivector::<f64, Euclidean3>(),
    ) {
        let lhs = a.outer(&b).outer(&c);
        let rhs = a.outer(&b.outer(&c));
        prop_assert!(lhs.approx_eq(&rhs, 1e-10));
    }

    #[test]
    fn grade_decomposition_complete(a in arb_multivector::<f64, Euclidean3>()) {
        let sum: Multivector<f64, Euclidean3> = (0..=3)
            .map(|k| a.grade_select(k))
            .fold(Multivector::zero(), |acc, x| acc + x);
        prop_assert!(sum.approx_eq(&a, 1e-10));
    }

    #[test]
    fn dual_undual_inverse(a in arb_multivector::<f64, Euclidean3>()) {
        let roundtrip = a.dual().undual();
        prop_assert!(roundtrip.approx_eq(&a, 1e-10));
    }

    #[test]
    fn outer_increases_grade(
        a in arb_vector::<f64, Euclidean3>(),
        b in arb_vector::<f64, Euclidean3>(),
    ) {
        let wedge = a.outer(&b);
        // Result should be grade 2 (bivector)
        prop_assert!(wedge.grade_select(0).is_zero(1e-10));
        prop_assert!(wedge.grade_select(1).is_zero(1e-10));
        // grade 2 may be non-zero
    }
}
```

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - all proptest properties pass
- [ ] `cargo clippy` - no warnings
- [ ] All items have rustdoc with mathematical notation
