# PRD-24: Remove Geometric Products, Add Dot Products

**Status**: Draft
**Goal**: Remove `GeometricProduct`/`Antigeometric` traits and free functions (not type-safe), add RGA dot products

## Background

The geometric product of two blades with definite grades can produce a multivector spanning multiple grades. For example:

```rust
// Vector × Vector = Scalar + Bivector (NOT a single type)
let v1: Vector<f64> = Vector::new(1.0, 0.0, 0.0);
let v2: Vector<f64> = Vector::new(0.0, 1.0, 0.0);
let result = v1.geometric(&v2);  // What type is this?
// result has grade-0 (scalar) AND grade-2 (bivector) components
```

The current `GeometricProduct` trait forces a single `Output` type, which is mathematically incorrect for most blade combinations.

## Problem Statement

### Issue 1: Geometric Product Output is Multi-Grade

The geometric product `a × b` produces components at grades `|gr(a) - gr(b)|, |gr(a) - gr(b)| + 2, ..., gr(a) + gr(b)`.

| LHS Grade | RHS Grade | Output Grades |
|-----------|-----------|---------------|
| 1 (Vector) | 1 (Vector) | 0, 2 (Scalar + Bivector) |
| 1 (Vector) | 2 (Bivector) | 1, 3 (Vector + Trivector) |
| 2 (Bivector) | 2 (Bivector) | 0, 2, 4 (Scalar + Bivector + Quadvector) |

There's no single specialized type that captures these mixed-grade outputs.

### Issue 2: Current Implementations are Incomplete/Wrong

Looking at the generated code, `GeometricProduct` implementations exist but they **project** the result to a single grade, losing information:

```rust
// Current generated code (projective3)
impl<T: Float> GeometricProduct<Vector<T>> for Vector<T> {
    type Output = ???;  // What can this be?
    // The true result is Scalar + Bivector, not representable as one type
}
```

### Issue 3: Versors Are Different

Versors (Rotor, Motor, Flector) have sandwich products that **do** preserve grade:
- `Rotor × Vector × Rotor⁻¹` → `Vector` (grade 1 preserved)
- `Motor × Point × Motor⁻¹` → `Point` (same type out)

The `Sandwich` and `Antisandwich` traits are type-safe and should be kept.

### Issue 4: Missing RGA Dot Products

[RGA defines dot products](https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products) that **are** type-safe:
- **Dot product** `a • b`: Non-zero only when `grade(a) = grade(b)`, returns scalar
- **Antidot product** `a ⊚ b`: The antiproduct complement, also returns scalar

These are exactly what we need for measuring angles and magnitudes between same-grade elements.

## Solution

### Phase 1: Remove GeometricProduct Trait and Free Functions

Remove from `src/ops.rs`:
```rust
// REMOVE this trait definition
pub trait GeometricProduct<Rhs = Self> {
    type Output;
    fn geometric(&self, rhs: &Rhs) -> Self::Output;
}
```

Remove from codegen - stop generating:
- `geometric_*` free functions in `products.rs`
- `impl GeometricProduct<...> for ...` trait impls

### Phase 2: Remove Antigeometric Trait and Free Functions

Remove from `src/ops.rs`:
```rust
// REMOVE this trait definition
pub trait Antigeometric<Rhs = Self> {
    type Output;
    fn antigeometric(&self, rhs: &Rhs) -> Self::Output;
}
```

Remove from codegen - stop generating:
- `antigeometric_*` free functions in `products.rs`
- `impl Antigeometric<...> for ...` trait impls

### Phase 3: Add RGA Dot Product Trait

Add to `src/ops.rs`:
```rust
/// Dot product (metric inner product).
///
/// The dot product `a • b` is non-zero only when `grade(a) = grade(b)`.
/// It measures the "alignment" between same-grade elements using the metric.
///
/// Defined as: `a • b = aᵀ G b` where G is the metric exomorphism matrix.
///
/// # Properties
///
/// - Only same-grade elements produce non-zero results
/// - Symmetric: `a • b = b • a`
/// - Returns a scalar
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Dot;
///
/// let cos_angle: f64 = vector1.dot(&vector2);
/// let bivector_magnitude_sq: f64 = bivector.dot(&bivector);
/// ```
///
/// # Reference
///
/// [RGA Dot Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products)
pub trait Dot<Rhs = Self> {
    /// The scalar type.
    type Scalar;

    /// Computes the dot product `self • rhs`.
    fn dot(&self, rhs: &Rhs) -> Self::Scalar;
}
```

### Phase 4: Add RGA Antidot Product Trait

Add to `src/ops.rs`:
```rust
/// Antidot product (metric antiproduct inner product).
///
/// The antidot product `a ⊚ b` is the De Morgan dual of the dot product:
/// `a ⊚ b = ā • b̄` (complement dot complement).
///
/// Defined as: `a ⊚ b = aᵀ 𝔾 b` where 𝔾 is the metric antiexomorphism matrix.
///
/// # Properties
///
/// - Only same-antigrade elements produce non-zero results
/// - Symmetric
/// - Returns a scalar
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Antidot;
///
/// let result: f64 = plane1.antidot(&plane2);
/// ```
///
/// # Reference
///
/// [RGA Dot Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products)
pub trait Antidot<Rhs = Self> {
    /// The scalar type.
    type Scalar;

    /// Computes the antidot product `self ⊚ rhs`.
    fn antidot(&self, rhs: &Rhs) -> Self::Scalar;
}
```

### Phase 5: Update Codegen

In `crates/clifford-codegen/src/codegen/`:

1. **Remove** from `products.rs`:
   - `generate_geometric_products()` function
   - `generate_antigeometric_products()` function
   - `ProductKind::Geometric` and `ProductKind::Antigeometric` variants

2. **Remove** from `traits.rs`:
   - `generate_geometric_product_trait()` function
   - `generate_antigeometric_trait()` function

3. **Add** to `products.rs`:
   - `generate_dot_products()` function
   - `generate_antidot_products()` function
   - `ProductKind::Dot` and `ProductKind::Antidot` variants

4. **Add** to `traits.rs`:
   - `generate_dot_trait()` function
   - `generate_antidot_trait()` function

### Phase 6: Recompute Sandwich Products

Without geometric products, sandwich products must be computed using the type-safe products. The sandwich `V × X × rev(V)` can be decomposed into grade-preserving operations:

For a versor V and operand X, the sandwich product decomposes based on which grade components survive:

```rust
// Sandwich computation using type-safe products
// V × X × rev(V) where V is a versor and X is an element

// For Rotor (grades 0, 2) × Vector (grade 1) × rev(Rotor):
// The result is grade 1 (Vector), computed via:
// - dot products for same-grade contractions
// - wedge/antiwedge for grade-changing parts that cancel to target grade
```

**Implementation approach:**

The codegen will compute sandwich products by:
1. Expanding `V × X × rev(V)` symbolically using Symbolica
2. Collecting only the terms that contribute to the output grade
3. These terms naturally involve compositions of:
   - `dot(a, b)` - same-grade contractions (scalar result)
   - `antidot(a, b)` - same-antigrade contractions
   - `wedge(a, b)` - grade-raising
   - `antiwedge(a, b)` - antigrade-raising (meet)

The sandwich output is grade-preserving, so the symbolic expansion will yield a result expressible purely in terms of dot, antidot, wedge, and antiwedge operations.

**Example - Rotor sandwich on Vector (Euclidean 3D):**
```rust
// Rotor R = s + B (scalar + bivector)
// Vector v
// R × v × rev(R) = (s + B) × v × (s - B)
//                = s²v + s(Bv - vB) + B×v×(-B)
//                = s²v + 2s(B·v) + B×v×(-B)
//
// The B·v is a contraction (dot-like), the B×v×(-B) uses wedge structure
```

### Phase 7: Keep Mul Operators for Versors

The `std::ops::Mul` implementations for versor composition will be recomputed using the new product decomposition:
```rust
// KEEP - but recompute using dot/wedge decomposition
impl<T: Float> Mul<Rotor<T>> for Rotor<T> { ... }
impl<T: Float> Mul<Motor<T>> for Motor<T> { ... }
```

Since versors are closed under multiplication (Rotor × Rotor = Rotor), the output type is well-defined and the Symbolica expansion will produce correct results.

### Phase 8: Products That Remain (Type-Safe)

| Product | Why Type-Safe |
|---------|---------------|
| `Wedge` | Output grade = gr(a) + gr(b) (single grade) |
| `Antiwedge` | Output antigrade = ag(a) + ag(b) (single grade) |
| `LeftContract` | Output grade = gr(b) - gr(a) (single grade) |
| `RightContract` | Output grade = gr(a) - gr(b) (single grade) |
| `BulkContract` | Reduces grade deterministically |
| `WeightContract` | Reduces grade deterministically |
| `BulkExpand` | Reduces antigrade deterministically |
| `WeightExpand` | Reduces antigrade deterministically |
| `Inner` | Output grade = \|gr(a) - gr(b)\| (single grade) |
| `Dot` | **NEW** - Returns scalar, same-grade only |
| `Antidot` | **NEW** - Returns scalar, same-antigrade only |
| `Sandwich` | Preserves operand type (versor action) |
| `Antisandwich` | Preserves operand type (versor action) |

### Phase 9: Update Documentation

Update:
- `CLAUDE.md` - Remove GeometricProduct, add Dot/Antidot
- `src/ops.rs` - Module docs
- PRD-23 - Mark GeometricProduct section as superseded

## Implementation Plan

### Step 1: Update `src/ops.rs`

1. Delete `GeometricProduct` trait definition (lines 39-57)
2. Delete `Antigeometric` trait definition (lines 315-333)
3. Add `Dot` trait definition (new)
4. Add `Antidot` trait definition (new)
5. Update module documentation

### Step 2: Update ProductTable in Codegen

In `crates/clifford-codegen/src/algebra/table.rs`:

1. Add `dot(a, b)` method - compute metric inner product (same-grade only)
2. Add `antidot(a, b)` method - compute antiproduct complement

### Step 3: Update Codegen Product Generation

In `crates/clifford-codegen/src/codegen/products.rs`:

1. Remove `generate_geometric_products()` and `generate_antigeometric_products()`
2. Add `generate_dot_products()` and `generate_antidot_products()`
3. Update `ProductKind` enum - remove Geometric/Antigeometric, add Dot/Antidot
4. Update sandwich product generation to use new decomposition

### Step 4: Update Codegen Trait Generation

In `crates/clifford-codegen/src/codegen/traits.rs`:

1. Remove `generate_geometric_product_trait()` function
2. Remove `generate_antigeometric_trait()` function
3. Add `generate_dot_trait()` function
4. Add `generate_antidot_trait()` function
5. Update `generate_mul_trait()` to use new product decomposition for versors

### Step 5: Update Discovery/Inference

In `crates/clifford-codegen/src/discovery/products.rs`:

1. Remove `ProductType::Geometric` and `ProductType::Antigeometric`
2. Add `ProductType::Dot` and `ProductType::Antidot`
3. Update product inference logic for dot/antidot

### Step 6: Update Symbolic Product Computation

In `crates/clifford-codegen/src/symbolic/product.rs`:

1. Remove `ProductKind::Geometric` handling
2. Add `ProductKind::Dot` and `ProductKind::Antidot` handling
3. Update sandwich computation to use decomposed products

### Step 7: Update TOML Specifications

In `crates/clifford-codegen/src/spec/`:

1. Remove `geometric` and `antigeometric` from `ProductsSpec`
2. Add `dot` and `antidot` to `ProductsSpec`
3. Update TOML parsing

### Step 8: Update Algebra TOML Files

In `algebras/*.toml`:

1. Remove `[products.geometric]` sections
2. Remove `[products.antigeometric]` sections
3. Add `[products.dot]` sections
4. Add `[products.antidot]` sections

### Step 9: Regenerate All Algebras

```bash
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done
```

### Step 10: Update Tests

1. Remove tests using `GeometricProduct` or `Antigeometric` traits
2. Remove tests using `geometric_*` or `antigeometric_*` free functions
3. Add tests for `Dot` and `Antidot` traits
4. Add tests for `dot_*` and `antidot_*` free functions
5. Verify sandwich products still work correctly

## What Stays

### Traits That Remain (with new additions)

```rust
// src/ops.rs - KEEP these traits
pub trait Wedge<Rhs = Self> { ... }
pub trait Antiwedge<Rhs = Self> { ... }
pub trait Inner<Rhs = Self> { ... }
pub trait LeftContract<Rhs = Self> { ... }
pub trait RightContract<Rhs = Self> { ... }
pub trait BulkContract<Rhs = Self> { ... }
pub trait WeightContract<Rhs = Self> { ... }
pub trait BulkExpand<Rhs = Self> { ... }
pub trait WeightExpand<Rhs = Self> { ... }
pub trait Sandwich<Operand> { ... }
pub trait Antisandwich<Operand> { ... }

// NEW traits
pub trait Dot<Rhs = Self> { ... }
pub trait Antidot<Rhs = Self> { ... }
```

### Operators That Remain

```rust
// std::ops impls - KEEP for versors (recomputed using new products)
impl<T: Float> Mul<Rotor<T>> for Rotor<T> { ... }
impl<T: Float> Mul<Motor<T>> for Motor<T> { ... }
impl<T: Float> Mul<Flector<T>> for Flector<T> { ... }
```

### Free Functions That Remain

```rust
// products.rs - KEEP grade-safe products
pub fn wedge_point_point(...) -> Line
pub fn antiwedge_plane_plane(...) -> Line
pub fn dot_vector_vector(...) -> T        // NEW
pub fn antidot_plane_plane(...) -> T      // NEW
pub fn sandwich_rotor_vector(...) -> Vector
// etc.
```

## Files Changed

| File | Action |
|------|--------|
| `src/ops.rs` | Remove GeometricProduct/Antigeometric, add Dot/Antidot |
| `crates/clifford-codegen/src/algebra/table.rs` | Add dot/antidot methods |
| `crates/clifford-codegen/src/codegen/products.rs` | Remove geometric, add dot products |
| `crates/clifford-codegen/src/codegen/traits.rs` | Remove geometric traits, add dot traits |
| `crates/clifford-codegen/src/discovery/products.rs` | Update ProductType enum |
| `crates/clifford-codegen/src/symbolic/product.rs` | Update ProductKind enum |
| `crates/clifford-codegen/src/spec/ir.rs` | Update ProductsSpec struct |
| `crates/clifford-codegen/src/spec/raw.rs` | Update TOML parsing |
| `algebras/euclidean2.toml` | Remove geometric, add dot sections |
| `algebras/euclidean3.toml` | Remove geometric, add dot sections |
| `algebras/projective2.toml` | Remove geometric, add dot sections |
| `algebras/projective3.toml` | Remove geometric, add dot sections |
| `src/specialized/*/generated/products.rs` | Regenerate (remove geometric_, add dot_) |
| `src/specialized/*/generated/traits.rs` | Regenerate (remove GeometricProduct, add Dot) |
| `docs/prd/prd-23-method-based-products.md` | Update to note traits removed |
| `CLAUDE.md` | Update product trait table |

## Testing Strategy

### Test New Dot Products

```rust
proptest! {
    #[test]
    fn dot_product_same_grade(v1 in arb_vector(), v2 in arb_vector()) {
        // Dot product of same-grade elements returns scalar
        let result: f64 = v1.dot(&v2);

        // Should be symmetric
        let reverse: f64 = v2.dot(&v1);
        prop_assert!(relative_eq!(result, reverse, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }

    #[test]
    fn dot_product_norm_squared(v in arb_vector()) {
        // v · v = |v|²
        let dot_self: f64 = v.dot(&v);
        let norm_sq = v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
        prop_assert!(relative_eq!(dot_self, norm_sq, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }
}
```

### Test Sandwich Products Still Work

```rust
proptest! {
    #[test]
    fn sandwich_preserves_grade(rotor in arb_unit_rotor(), v in arb_vector()) {
        // Sandwich must preserve grade - Vector in, Vector out
        let transformed: Vector<f64> = rotor.sandwich(&v);

        // Verify norm preservation (for unit rotor)
        let original_norm = v.dot(&v).sqrt();
        let transformed_norm = transformed.dot(&transformed).sqrt();
        prop_assert!(relative_eq!(original_norm, transformed_norm, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }
}
```

### Test Versor Composition

```rust
#[test]
fn rotor_mul_still_works() {
    let r1 = Rotor::from_angle_plane(0.5, Bivector::unit_xy());
    let r2 = Rotor::from_angle_plane(0.3, Bivector::unit_xy());
    let composed: Rotor<f64> = r1 * r2;  // Via Mul operator

    // Verify composition is correct
    let v = Vector::new(1.0, 0.0, 0.0);
    let via_composed = composed.sandwich(&v);
    let via_sequential = r1.sandwich(&r2.sandwich(&v));
    assert!(relative_eq!(via_composed, via_sequential, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
}
```

### Verify No Geometric Product Functions Exist

```rust
#[test]
fn geometric_product_removed() {
    // This should fail to compile if geometric products still exist
    // compile_fail test - geometric_rotor_rotor should not exist
}
```

## Verification

After implementation, run:
```bash
cargo fmt && cargo clippy && cargo doc --no-deps && cargo nextest run && cargo deny check
```

## Success Criteria

1. `GeometricProduct` and `Antigeometric` traits removed from `src/ops.rs`
2. `Dot` and `Antidot` traits added to `src/ops.rs`
3. No generated `impl GeometricProduct` or `impl Antigeometric` in any algebra
4. No generated `geometric_*` or `antigeometric_*` free functions
5. New `dot_*` and `antidot_*` free functions generated
6. All grade-preserving product traits (Wedge, Sandwich, etc.) still work
7. Versor `Mul` operators still work (recomputed via new products)
8. All tests pass
9. No compiler errors or warnings

## Migration Guide

### Before (Geometric Product Trait)
```rust
use clifford::ops::GeometricProduct;

let result = rotor1.geometric(&rotor2);
```

### After (Mul Operator)
```rust
// Use Mul operator for versor composition
let result = rotor1 * rotor2;  // Via std::ops::Mul
```

### Before (Geometric Free Function)
```rust
use clifford::specialized::euclidean::dim3::products::geometric_rotor_rotor;

let result = geometric_rotor_rotor(&rotor1, &rotor2);
```

### After (Mul Operator)
```rust
let result = rotor1 * rotor2;  // Via std::ops::Mul - same computation, cleaner API
```

### New Dot Product API
```rust
use clifford::ops::Dot;

// Compute angle between vectors
let cos_angle = v1.dot(&v2) / (v1.dot(&v1).sqrt() * v2.dot(&v2).sqrt());

// Compute magnitude squared
let magnitude_sq = bivector.dot(&bivector);
```

## References

- [RGA Dot Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products)
- [RGA Geometric Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_products)
- PRD-23: Method-Based Products (partially superseded)
