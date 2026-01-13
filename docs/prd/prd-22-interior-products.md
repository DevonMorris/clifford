# PRD-22: Interior Products and Product Naming Corrections

**Status**: Draft
**Goal**: Fix product naming to match RGA conventions and add the four true interior products

## Background

The current codegen conflates different products and uses non-standard naming. According to [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/wiki/), the terminology should be:

### RGA Symbol Conventions

| Symbol | Name | Description |
|--------|------|-------------|
| `∧` | **wedge** | Exterior product (grade-raising) |
| `∨` | **antiwedge** | Regressive product (antigrade-raising) |
| `★` | **dual** | Bulk dual (metric complement) |
| `☆` | **antidual** | Weight dual (antiproduct complement) |
| `ã` (tilde above) | **reverse** | Reverses order of basis vectors in each blade |
| `a̲` (tilde below) | **antireverse** | Antiproduct complement of reverse |
| `ā` (bar above) | **right complement** | Right complement (a ∧ ā = 𝟙) |
| `a̱` (bar below) | **left complement** | Left complement (a̱ ∧ a = 𝟙) |

### Exterior Products (what we currently have, need renaming)

| Current Name | RGA Name | Symbol | Definition |
|--------------|----------|--------|------------|
| `exterior` | **wedge** | `∧` | Grade-raising: gr(a ∧ b) = gr(a) + gr(b) |
| `regressive` | **antiwedge** | `∨` | Antigrade-raising: ag(a ∨ b) = ag(a) + ag(b) |

### Interior Products (what we're missing)

The current `interior` product in codegen is the **symmetric inner product** (Hestenes inner product), NOT what RGA calls "interior products". The RGA interior products are:

| Product | Symbol | Formula | Description |
|---------|--------|---------|-------------|
| **Bulk contraction** | `a ∨ b★` | antiwedge with bulk dual | Reduces grade by gr(a) - gr(b) |
| **Weight contraction** | `a ∨ b☆` | antiwedge with weight dual | Reduces grade by gr(a) - gr(b) |
| **Bulk expansion** | `a ∧ b★` | wedge with bulk dual | Reduces antigrade |
| **Weight expansion** | `a ∧ b☆` | wedge with weight dual | Reduces antigrade |

Where:
- `★` is the **bulk dual** (uses metric exomorphism)
- `☆` is the **weight dual** (uses metric antiexomorphism, also called "antidual")

### Duals Reference

From [RGA Duals](https://rigidgeometricalgebra.org/wiki/index.php?title=Duals):

- **Bulk dual** (`u★`): `u★ = ũ ⋙ 1` (geometric product with pseudoscalar on right)
- **Weight dual** (`u☆`): `u☆ = ũ ⋗ 1` (antigeometric product with antiscalar on right)

For non-degenerate metrics, the bulk and weight duals coincide. For degenerate metrics (PGA), they differ and both are needed.

## Problem Statement

### Issue 1: Non-standard naming

Current naming doesn't match RGA:
- `exterior_*` should be `wedge_*`
- `regressive_*` should be `antiwedge_*`

### Issue 2: Incorrect "interior" definition

The current `interior` product computes `|grade(a) - grade(b)|` (symmetric inner product), which is:
- Called the "Hestenes inner product" or "fat dot product" in GA literature
- NOT what RGA calls "interior products"

### Issue 3: Missing true interior products

The RGA interior products (contractions and expansions) are not implemented:
- `bulk_contraction`: `a ∨ b★`
- `weight_contraction`: `a ∨ b☆`
- `bulk_expansion`: `a ∧ b★`
- `weight_expansion`: `a ∧ b☆`

These are essential for projections, rejections, and geometric operations in PGA.

## Solution

### Phase 1: Rename Existing Products

#### 1.1 `exterior` → `wedge`

```rust
// Before
pub fn exterior_point_point<T: Float>(...) -> Line<T>

// After
pub fn wedge_point_point<T: Float>(...) -> Line<T>
```

TOML spec changes:
```toml
# Before
[products.exterior]
Point_Point = "Line"

# After
[products.wedge]
Point_Point = "Line"
```

#### 1.2 `regressive` → `antiwedge`

```rust
// Before
pub fn regressive_plane_plane<T: Float>(...) -> Line<T>

// After
pub fn antiwedge_plane_plane<T: Float>(...) -> Line<T>
```

TOML spec changes:
```toml
# Before
[products.regressive]
Plane_Plane = "Line"

# After
[products.antiwedge]
Plane_Plane = "Line"
```

#### 1.3 Decide what to do with `interior`

Options:
1. **Rename to `hestenes_inner`** - Technically accurate but verbose
2. **Keep as `inner`** - Common in GA libraries (ganja.js uses this)
3. **Remove** - Let users compose from left/right contractions

**Recommendation**: Rename to `inner` (shorter, widely understood in GA context). The left_contract and right_contract functions remain for asymmetric contractions.

### Phase 2: Add Dual Operations

Before implementing interior products, we need the dual operations in `ProductTable`:

#### 2.1 Bulk Dual (`★`)

```rust
impl ProductTable {
    /// Computes the bulk dual of a blade.
    ///
    /// The bulk dual is defined as: u★ = ũ × I⁻¹
    /// where I is the pseudoscalar (unit of highest grade).
    ///
    /// Returns (sign, result_blade).
    pub fn bulk_dual(&self, blade: usize) -> (i8, usize) {
        let pseudoscalar = (1 << self.dim()) - 1; // All bits set
        // ũ × I⁻¹ = reverse(blade) * inv(pseudoscalar)
        // For orthonormal bases, inv(I) = ±I
        let (rev_sign, _) = self.reverse_sign(blade);
        let (prod_sign, result) = self.geometric(blade, pseudoscalar);
        // Need to account for I⁻¹ = ±I based on signature
        let inv_sign = self.pseudoscalar_inverse_sign();
        (rev_sign * prod_sign * inv_sign, result)
    }
}
```

#### 2.2 Weight Dual (`☆`)

The weight dual already exists as `weight_dual` in `ProductTable` (used by projections).

### Phase 3: Add Interior Products

#### 3.1 Bulk Contraction

```rust
/// Bulk contraction: a ∨ b★
///
/// The bulk contraction reduces grade by grade(a) - grade(b).
/// Uses the bulk dual (★) which applies the metric complement.
impl ProductTable {
    pub fn bulk_contraction(&self, a: usize, b: usize) -> (i8, usize) {
        let (dual_sign, b_dual) = self.bulk_dual(b);
        if dual_sign == 0 {
            return (0, 0);
        }
        let (reg_sign, result) = self.regressive(a, b_dual);
        (dual_sign * reg_sign, result)
    }
}
```

#### 3.2 Weight Contraction

```rust
/// Weight contraction: a ∨ b☆
///
/// The weight contraction reduces grade by grade(a) - grade(b).
/// Uses the weight dual (☆) which applies the antiproduct complement.
impl ProductTable {
    pub fn weight_contraction(&self, a: usize, b: usize) -> (i8, usize) {
        let (dual_sign, b_dual) = self.weight_dual(b);
        if dual_sign == 0 {
            return (0, 0);
        }
        let (reg_sign, result) = self.regressive(a, b_dual);
        (dual_sign * reg_sign, result)
    }
}
```

#### 3.3 Bulk Expansion

```rust
/// Bulk expansion: a ∧ b★
///
/// The bulk expansion reduces antigrade.
/// Uses the bulk dual (★) with the exterior product.
impl ProductTable {
    pub fn bulk_expansion(&self, a: usize, b: usize) -> (i8, usize) {
        let (dual_sign, b_dual) = self.bulk_dual(b);
        if dual_sign == 0 {
            return (0, 0);
        }
        let (ext_sign, result) = self.exterior(a, b_dual);
        (dual_sign * ext_sign, result)
    }
}
```

#### 3.4 Weight Expansion

```rust
/// Weight expansion: a ∧ b☆
///
/// The weight expansion reduces antigrade.
/// Uses the weight dual (☆) with the exterior product.
impl ProductTable {
    pub fn weight_expansion(&self, a: usize, b: usize) -> (i8, usize) {
        let (dual_sign, b_dual) = self.weight_dual(b);
        if dual_sign == 0 {
            return (0, 0);
        }
        let (ext_sign, result) = self.exterior(a, b_dual);
        (dual_sign * ext_sign, result)
    }
}
```

### Phase 4: Expression Simplification with Symbolica

All generated product functions must use **Symbolica** for symbolic simplification before emitting Rust code. This applies to:
- All binary products (wedge, antiwedge, contractions, expansions)
- All sandwich products
- All unary operations (duals, reverse, etc.)

#### 4.1 Simplification Pipeline Using Symbolica

```rust
fn generate_product_expression(&self, ...) -> Vec<TokenStream> {
    // 1. Compute symbolic expression using Symbolica
    let symbolic_fields = symbolic_product.compute(...);

    // 2. Apply constraint substitutions (e.g., unit norm) via Symbolica
    let with_constraints = constraint_simplifier.apply(&field.expression);

    // 3. Simplify using Symbolica's expand() and collect()
    let simplified = field.expression.expand();  // Symbolica expansion
    let simplified = simplified.collect(...);     // Symbolica collection

    // 4. Convert to Rust code
    converter.convert(&simplified)
}
```

#### 4.2 Symbolica Simplification Operations

Use Symbolica's built-in simplification methods:

```rust
use symbolica::atom::Atom;

// Expand products: (a + b) * c → a*c + b*c
let expanded = expr.expand();

// Collect like terms around variables
let collected = expanded.collect::<Vec<_>>(&variables, None, None);

// Together for rational simplification
let simplified = expr.together();

// Cancel common factors
let cancelled = expr.cancel();
```

The key Symbolica operations:
1. **`expand()`**: Distributes products over sums
2. **`collect()`**: Groups terms by specified variables
3. **`together()`**: Combines fractions over common denominator
4. **`cancel()`**: Cancels common factors in numerator/denominator

#### 4.3 Symbolica for ALL Generated Products

**CRITICAL**: All product generation must use Symbolica for symbolic computation and simplification. This includes:

| Product Category | Examples | Simplification Benefit |
|-----------------|----------|----------------------|
| **Binary products** | `wedge_*`, `antiwedge_*`, `geometric_*` | Combine like terms, cancel zeros |
| **Interior products** | `bulk_contraction_*`, `weight_expansion_*` | Flatten composite operations |
| **Sandwich products** | `sandwich_motor_point` | Combine v*x*rev(v) terms |
| **Antisandwich products** | `antisandwich_motor_point` | Combine v⊛x⊛antirev(v) terms |
| **Scalar products** | `scalar_vector_vector` | Extract grade-0 component |
| **Unary operations** | `dual_*`, `reverse_*`, `complement_*` | Apply sign patterns |

**Interior products** especially benefit because they are composite:

```rust
// Bulk contraction: a ∨ b★
// Symbolica computes: antiwedge(a, bulk_dual(b))
// Then simplifies to direct formula:

pub fn bulk_contraction_point_plane<T: Float>(a: &Point<T>, b: &Plane<T>) -> Scalar<T> {
    // Simplified result after Symbolica computation of a ∨ b★
    Scalar::new(a.e1() * b.e023() + a.e2() * b.e031() + a.e3() * b.e012())
}
```

**Sandwich products** benefit from term combination:

```rust
// Sandwich: v * x * rev(v)
// Without Symbolica: many duplicate terms from v_i * x_j * v_k
// With Symbolica: terms are combined and simplified

pub fn sandwich_rotor_vector<T: Float>(v: &Rotor<T>, x: &Vector<T>) -> Vector<T> {
    // Symbolica combines terms like:
    // 2*v.s()*v.xy()*x.y() appears once, not as separate s*xy*y + xy*y*s
    Vector::new(
        v.s() * v.s() * x.x() + T::TWO * v.s() * v.xy() * x.y() - v.xy() * v.xy() * x.x(),
        // ... simplified expressions
    )
}
```

**No manual term computation**: The current `compute_sandwich_terms()` and `compute_terms()` methods that manually iterate over field combinations should be replaced with Symbolica-based symbolic computation for consistency and correctness.

### Phase 5: Update Code Generator

#### 5.1 New ProductKind Enum

```rust
pub enum ProductKind {
    // Fundamental products
    Geometric,
    Antigeometric,

    // Exterior products
    Wedge,           // renamed from Exterior
    Antiwedge,       // renamed from Regressive

    // Contractions (grade-lowering)
    LeftContraction,
    RightContraction,
    BulkContraction,    // NEW: a ∨ b★
    WeightContraction,  // NEW: a ∨ b☆

    // Expansions (antigrade-lowering)
    BulkExpansion,      // NEW: a ∧ b★
    WeightExpansion,    // NEW: a ∧ b☆

    // Inner products
    Inner,           // renamed from Interior (Hestenes inner)
    Scalar,          // grade-0 extraction
}
```

#### 5.2 Generated Function Names

| ProductKind | Function Name Pattern | Example |
|-------------|----------------------|---------|
| `Wedge` | `wedge_{a}_{b}` | `wedge_point_point` |
| `Antiwedge` | `antiwedge_{a}_{b}` | `antiwedge_plane_plane` |
| `BulkContraction` | `bulk_contraction_{a}_{b}` | `bulk_contraction_point_plane` |
| `WeightContraction` | `weight_contraction_{a}_{b}` | `weight_contraction_point_plane` |
| `BulkExpansion` | `bulk_expansion_{a}_{b}` | `bulk_expansion_line_point` |
| `WeightExpansion` | `weight_expansion_{a}_{b}` | `weight_expansion_line_point` |

### Phase 6: Update TOML Specification Format

```toml
[products]
# Fundamental
geometric = true
antigeometric = true

# Exterior (renamed)
wedge = true          # renamed from exterior
antiwedge = true      # renamed from regressive

# Contractions
left_contraction = true
right_contraction = true
bulk_contraction = true    # NEW
weight_contraction = true  # NEW

# Expansions
bulk_expansion = true      # NEW
weight_expansion = true    # NEW

# Inner products
inner = true          # renamed from interior
scalar = true
```

## Deliverables

### Core Changes

- [ ] Add `bulk_dual` to `ProductTable`
- [ ] Add `bulk_contraction` to `ProductTable`
- [ ] Add `weight_contraction` to `ProductTable`
- [ ] Add `bulk_expansion` to `ProductTable`
- [ ] Add `weight_expansion` to `ProductTable`

### Symbolica Integration

- [ ] Migrate all product generation to use Symbolica for symbolic computation
- [ ] Use Symbolica `expand()` and `collect()` for expression simplification
- [ ] Replace manual `compute_terms()` with Symbolica-based computation
- [ ] Replace manual `compute_sandwich_terms()` with Symbolica-based computation
- [ ] Ensure all generated code benefits from symbolic simplification

### Codegen Changes

- [ ] Rename `ProductKind::Exterior` → `ProductKind::Wedge`
- [ ] Rename `ProductKind::Regressive` → `ProductKind::Antiwedge`
- [ ] Rename `ProductKind::Interior` → `ProductKind::Inner`
- [ ] Add `ProductKind::BulkContraction`
- [ ] Add `ProductKind::WeightContraction`
- [ ] Add `ProductKind::BulkExpansion`
- [ ] Add `ProductKind::WeightExpansion`
- [ ] Update function name generation
- [ ] Update TOML parser for new product names

### TOML Updates

- [ ] Rename `[products.exterior]` → `[products.wedge]`
- [ ] Rename `[products.regressive]` → `[products.antiwedge]`
- [ ] Rename `[products.interior]` → `[products.inner]`
- [ ] Add `[products.bulk_contraction]`
- [ ] Add `[products.weight_contraction]`
- [ ] Add `[products.bulk_expansion]`
- [ ] Add `[products.weight_expansion]`

### Regeneration

- [ ] Update `algebras/euclidean2.toml`
- [ ] Update `algebras/euclidean3.toml`
- [ ] Update `algebras/projective3.toml`
- [ ] Regenerate all algebras

### Documentation

- [ ] Update generated products.rs header documentation
- [ ] Update CLAUDE.md with new product names and symbol conventions
- [ ] Add mathematical references to RGA wiki

### CLAUDE.md Updates

Add a new section on RGA product notation:

```markdown
### RGA Product Notation

This library follows [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/) conventions:

| Symbol | Name | Generated Function | Description |
|--------|------|-------------------|-------------|
| `∧` | **wedge** | `wedge_*` | Exterior product (grade-raising) |
| `∨` | **antiwedge** | `antiwedge_*` | Regressive product (antigrade-raising) |
| `★` | **dual** | `dual_*` | Bulk dual (metric complement) |
| `☆` | **antidual** | `antidual_*` | Weight dual (antiproduct complement) |
| `ã` (tilde above) | **reverse** | `reverse_*` | Reverses order of basis vectors |
| `a̲` (tilde below) | **antireverse** | `antireverse_*` | Antiproduct complement of reverse |
| `ā` (bar above) | **right complement** | `right_complement_*` | Right complement (a ∧ ā = 𝟙) |
| `a̱` (bar below) | **left complement** | `left_complement_*` | Left complement (a̱ ∧ a = 𝟙) |

**Interior Products** (contractions and expansions):

| Product | Formula | Generated Function |
|---------|---------|-------------------|
| Bulk contraction | `a ∨ b★` | `bulk_contraction_*` |
| Weight contraction | `a ∨ b☆` | `weight_contraction_*` |
| Bulk expansion | `a ∧ b★` | `bulk_expansion_*` |
| Weight expansion | `a ∧ b☆` | `weight_expansion_*` |

**Note**: The `inner_*` functions compute the Hestenes inner product (symmetric, grade `|ga - gb|`),
which is different from the RGA interior products defined above.

**Authoritative Reference**: The [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page)
is the authoritative reference for product definitions, formulas, and terminology used in this library.
When in doubt about GA operations, consult the RGA wiki first.
```

## Testing Strategy

### Unit Tests

```rust
#[test]
fn bulk_contraction_definition() {
    // Verify: bulk_contraction(a, b) = antiwedge(a, bulk_dual(b))
    let table = ProductTable::new(&Algebra::pga(3));

    for a in 0..16 {
        for b in 0..16 {
            let (bc_sign, bc_result) = table.bulk_contraction(a, b);

            let (dual_sign, b_dual) = table.bulk_dual(b);
            let (aw_sign, aw_result) = table.antiwedge(a, b_dual);

            assert_eq!(bc_sign, dual_sign * aw_sign);
            if bc_sign != 0 {
                assert_eq!(bc_result, aw_result);
            }
        }
    }
}
```

### Property Tests

```rust
proptest! {
    #[test]
    fn wedge_equals_exterior(a in any::<Point<f64>>(), b in any::<Point<f64>>()) {
        // New wedge should produce same result as old exterior
        let wedge_result = wedge_point_point(&a, &b);
        // Compare with Multivector wedge
        let mv_a = Multivector::from(a);
        let mv_b = Multivector::from(b);
        let mv_result = mv_a.wedge(&mv_b);

        prop_assert!(relative_eq!(
            Line::from(mv_result),
            wedge_result,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
```

### Verification Against RGA Examples

The RGA wiki provides Cayley tables for interior products in 4D. Verify our implementation matches.

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/algebra/table.rs` | Update - Add dual and interior methods |
| `crates/clifford-codegen/src/codegen/products.rs` | Update - Rename and add product kinds |
| `crates/clifford-codegen/src/spec/raw.rs` | Update - New TOML field names |
| `crates/clifford-codegen/src/spec/ir.rs` | Update - New product configuration |
| `algebras/euclidean2.toml` | Update - Rename product sections |
| `algebras/euclidean3.toml` | Update - Rename product sections |
| `algebras/projective3.toml` | Update - Rename product sections |
| `src/specialized/*/generated/products.rs` | Regenerate - New function names |
| `CLAUDE.md` | Update - Document new product names |

## References

**Authoritative Reference**:
- [Rigid Geometric Algebra Wiki (Main Page)](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page) - The primary reference for all GA operations in this library

**Specific Topics**:
- [RGA: Interior Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Interior_products)
- [RGA: Exterior Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Exterior_products)
- [RGA: Duals](https://rigidgeometricalgebra.org/wiki/index.php?title=Duals)
- [RGA: Cayley Tables for Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_products)

## Success Criteria

1. **Naming consistency**: Product names match RGA terminology (`wedge`, `antiwedge`, `bulk_contraction`, etc.)
2. **Mathematical correctness**: Interior products compute `a ∨ b★` and `a ∧ b★` correctly
3. **Symbolica simplification**: All generated products use Symbolica for symbolic computation and simplification
4. **Verification**: All products verified against Multivector implementation
5. **Documentation**: Clear mathematical references in generated code with RGA wiki links
6. **Performance**: Simplified expressions produce efficient generated code (no redundant operations)
