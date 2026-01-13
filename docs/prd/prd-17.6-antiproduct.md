# PRD-17.6: Complete Anti-Product Suite for PGA

**Status**: Draft
**Parent**: [PRD-17](prd-17-codegen-products.md)
**Goal**: Complete the anti-product suite including antiproduct, antiscalar, antireverse, dual/antidual

## Current State

The codegen already implements core anti-product functionality:

| Operation | Status | Location |
|-----------|--------|----------|
| `antiproduct()` | ✅ Implemented | `table.rs:318-336` |
| `antisandwich_*_*()` | ✅ Generated | `products.rs:603-692` |
| Antireverse sign | ✅ Used in antisandwich | `products.rs:926-980` |
| `complement()` | ✅ Implemented | `table.rs:261-273` |

**What's missing:**

| Operation | Status | Notes |
|-----------|--------|-------|
| `antigeometric_*_*()` | ❌ Missing | Standalone antiproduct functions |
| `antiscalar_*_*()` | ❌ Missing | Grade-n part of antigeometric |
| `antiwedge_*_*()` | ❌ Missing | Exterior antiproduct (regressive) |
| `reverse_*()` | ❌ Missing | Standalone reverse operation |
| `antireverse_*()` | ❌ Missing | Standalone antireverse operation |
| `dual_*()` | ❌ Missing | Bulk dual operation |
| `antidual_*()` | ❌ Missing | Weight dual operation |
| `Multivector::antiproduct()` | ❌ Missing | Generic method |
| `Multivector::antiwedge()` | ❌ Missing | Generic method |
| `Multivector::dual()` | ❌ Missing | Generic method |
| `Multivector::antidual()` | ❌ Missing | Generic method |

## Background: Anti-Operations from RGA Wiki

Reference: [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/)

### Geometric Antiproduct (⊛)

The geometric antiproduct is dual to the geometric product:

```
a ⊛ b = ∁(∁a * ∁b)
```

Where `∁` is the complement (right complement: `u ∧ ū = I`).

**De Morgan Laws:**
```
∁(a * b) = ∁a ⊛ ∁b
∁(a ⊛ b) = ∁a * ∁b
```

**Grade rules:** For antigrade = n - grade:
```
antigrade(a ⊛ b) follows geometric product rules for antigrade(a) * antigrade(b)
```

### Antiscalar Product

The antiscalar is the grade-n (pseudoscalar grade) part of the antigeometric product:

```
antiscalar(a, b) = ⟨a ⊛ b⟩_n
```

This is dual to the scalar product `⟨a * b⟩_0`.

### Antiwedge Product (∨)

The antiwedge (regressive product) is dual to the exterior product:

```
a ∨ b = ∁(∁a ∧ ∁b)
```

**De Morgan Laws:**
```
∁(a ∧ b) = ∁a ∨ ∁b
∁(a ∨ b) = ∁a ∧ ∁b
```

**Grade rule:**
```
antigrade(a ∨ b) = antigrade(a) + antigrade(b)
```

### Reverse and Antireverse

**Reverse** (for grade-k blade):
```
ũ = (-1)^(k(k-1)/2) * u
```

**Antireverse** (for grade-k blade in n-dimensional algebra):
```
u̲ = (-1)^((n-k)(n-k-1)/2) * u
```

**Relationship:**
```
u̲ = (-1)^(k*(n-k)) * (-1)^(n(n-1)/2) * ũ
```

### Dual and Antidual

**Bulk Dual (★):**
```
u★ = ũ ⊡ 𝟙
```
Where `⊡` is the "wedge dot" (scalar part of product with pseudoscalar).

**Weight Dual (☆):**
```
u☆ = ũ ⊗ 1
```
Where `⊗` uses the antiproduct.

**Double dual:**
```
u★★ = u☆☆ = (-1)^(gr(u)·ag(u)) · det(𝔤) · u
```

In PGA (det = 0): `u★★ = u☆☆ = 0`

## Implementation Plan

### Phase 1: Unary Operations in Codegen

Add generation for unary operations:

```rust
// In codegen/unary.rs (new file)

pub enum UnaryKind {
    Reverse,      // ũ = (-1)^(k(k-1)/2) * u
    Antireverse,  // u̲ = (-1)^((n-k)(n-k-1)/2) * u
    Dual,         // u★ - bulk dual
    Antidual,     // u☆ - weight dual
    Complement,   // ∁u - right complement
}
```

### Phase 2: Binary Antiproducts

Add standalone antiproduct generation:

```rust
pub enum ProductKind {
    // Existing
    Geometric,
    Exterior,
    LeftContraction,
    Regressive,
    Scalar,
    // New
    Antigeometric,  // a ⊛ b = ∁(∁a * ∁b)
    Antiwedge,      // a ∨ b = ∁(∁a ∧ ∁b) [same as Regressive]
    Antiscalar,     // ⟨a ⊛ b⟩_n
}
```

Note: `Antiwedge` is mathematically equivalent to `Regressive`, but we may want both names for clarity.

### Phase 3: Multivector Methods

Add methods to generic Multivector:

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Right complement: u ∧ ū = I
    pub fn complement(&self) -> Self { ... }

    /// Reverse: ũ = (-1)^(k(k-1)/2) * u per grade
    pub fn reverse(&self) -> Self { ... }

    /// Antireverse: u̲ = (-1)^((n-k)(n-k-1)/2) * u per grade
    pub fn antireverse(&self) -> Self { ... }

    /// Geometric antiproduct: a ⊛ b = ∁(∁a * ∁b)
    pub fn antiproduct(&self, other: &Self) -> Self {
        self.complement().geometric(&other.complement()).complement()
    }

    /// Antiwedge (regressive): a ∨ b = ∁(∁a ∧ ∁b)
    pub fn antiwedge(&self, other: &Self) -> Self {
        self.complement().exterior(&other.complement()).complement()
    }

    /// Bulk dual: u★ = reverse(u) contracted with pseudoscalar
    pub fn dual(&self) -> Self { ... }

    /// Weight antidual: u☆ = reverse(u) anticontracted with scalar
    pub fn antidual(&self) -> Self { ... }
}
```

### Phase 4: TOML Configuration

Add product sections:

```toml
[products]
geometric = true
exterior = true
left_contraction = true
scalar = true
# New
antigeometric = true
antiwedge = true
antiscalar = true

[unary]
reverse = true
antireverse = true
dual = true
antidual = true
complement = true
```

### Phase 5: Antisandwich Already Works

The antisandwich implementation already exists and uses the antiproduct correctly. Verify it works for all cases:

```rust
// Already generated: antisandwich_motor_point, antisandwich_flector_point, etc.
// These use table.antiproduct() internally
```

## Computing Antiproduct Tables

The antiproduct table can be computed directly from the geometric product table:

```rust
impl ProductTable {
    /// Precompute antiproduct table: a ⊛ b = ∁(∁a * ∁b)
    fn compute_antiproduct_table(&self) -> Vec<(i8, usize)> {
        let n = self.num_blades();
        let mut table = vec![(0i8, 0usize); n * n];

        for a in 0..n {
            for b in 0..n {
                // Get complements
                let (sign_ca, comp_a) = self.complement(a);
                let (sign_cb, comp_b) = self.complement(b);

                // Geometric product of complements
                let (sign_prod, prod) = self.geometric(comp_a, comp_b);
                if sign_prod == 0 {
                    table[a * n + b] = (0, 0);
                    continue;
                }

                // Complement of result
                let (sign_result, result) = self.complement(prod);

                // Total sign
                let total_sign = sign_ca * sign_cb * sign_prod * sign_result;
                table[a * n + b] = (total_sign, result);
            }
        }

        table
    }
}
```

## Deliverables

### Phase 1: Unary Operations
- [ ] Add `UnaryKind` enum in codegen
- [ ] Implement `generate_reverse()` function
- [ ] Implement `generate_antireverse()` function
- [ ] Implement `generate_complement()` function
- [ ] Implement `generate_dual()` function
- [ ] Implement `generate_antidual()` function
- [ ] Update TOML parser for `[unary]` section

### Phase 2: Binary Products
- [ ] Add `ProductKind::Antigeometric` variant
- [ ] Add `ProductKind::Antiscalar` variant
- [ ] Implement antigeometric product generation
- [ ] Implement antiscalar product generation
- [ ] Update TOML parser for antiproduct sections

### Phase 3: Multivector Methods
- [ ] Add `Multivector::complement()` method
- [ ] Add `Multivector::antireverse()` method
- [ ] Add `Multivector::antiproduct()` method
- [ ] Add `Multivector::antiwedge()` method
- [ ] Add `Multivector::dual()` method
- [ ] Add `Multivector::antidual()` method
- [ ] Add property-based tests for De Morgan laws
- [ ] Add property-based tests for double-dual identities

### Phase 4: PGA Verification
- [ ] Verify antisandwich produces correct transformations
- [ ] Verify generated antigeometric matches Multivector::antiproduct
- [ ] Re-enable nalgebra comparison tests (currently ignored)

## Testing Strategy

### Algebraic Property Tests

```rust
proptest! {
    // De Morgan law: ∁(a * b) = ∁a ⊛ ∁b
    #[test]
    fn de_morgan_geometric_antiproduct(a in any::<Mv>(), b in any::<Mv>()) {
        let lhs = (a.geometric(&b)).complement();
        let rhs = a.complement().antiproduct(&b.complement());
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = EPS));
    }

    // Antireverse property: (ab)̲ = b̲ a̲
    #[test]
    fn antireverse_reverses_product(a in any::<Mv>(), b in any::<Mv>()) {
        let lhs = (a.antiproduct(&b)).antireverse();
        let rhs = b.antireverse().antiproduct(&a.antireverse());
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = EPS));
    }
}
```

### Transformation Tests

```rust
proptest! {
    // Motor translation via antisandwich
    #[test]
    fn motor_antisandwich_translates(dx: f64, dy: f64, dz: f64, p in any::<Point<f64>>()) {
        let motor = Motor::from_translation(dx, dy, dz);
        let result = antisandwich_motor_point(&motor, &p);

        prop_assert!(abs_diff_eq!(result.x(), p.x() + dx, epsilon = EPS));
        prop_assert!(abs_diff_eq!(result.y(), p.y() + dy, epsilon = EPS));
        prop_assert!(abs_diff_eq!(result.z(), p.z() + dz, epsilon = EPS));
    }
}
```

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/codegen/unary.rs` | New - unary operation generation |
| `crates/clifford-codegen/src/codegen/products.rs` | Update - add antiproduct kinds |
| `crates/clifford-codegen/src/codegen/mod.rs` | Update - export unary module |
| `crates/clifford-codegen/src/spec/ir.rs` | Update - add unary spec |
| `crates/clifford-codegen/src/spec/parser.rs` | Update - parse unary section |
| `crates/clifford-codegen/src/spec/raw.rs` | Update - raw unary fields |
| `src/algebra/multivector.rs` | Update - add anti-methods |
| `algebras/*.toml` | Update - add unary/antiproduct config |
| `src/generated/*/products.rs` | Regenerate |

## Success Criteria

1. All De Morgan identities pass property tests
2. All antisandwich transformations match nalgebra Isometry3
3. Generated antiproducts match Multivector::antiproduct()
4. Dual/antidual satisfy double-dual identity (within floating point for non-PGA)
5. PGA nalgebra tests (currently ignored) pass

## References

- [RGA - Geometric Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_products)
- [RGA - Complements](https://rigidgeometricalgebra.org/wiki/index.php?title=Complements)
- [RGA - Duals](https://rigidgeometricalgebra.org/wiki/index.php?title=Duals)
- [RGA - Reverses](https://rigidgeometricalgebra.org/wiki/index.php?title=Reverses)
- [RGA - Motor](https://rigidgeometricalgebra.org/wiki/index.php?title=Motor)
- [Look, Ma, No Matrices!](https://enkimute.github.io/LookMaNoMatrices/)
