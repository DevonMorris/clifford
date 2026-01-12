# PRD-17: Codegen Product Completeness and Documentation

**Status**: Draft
**Goal**: Complete the codegen product suite with all standard GA products and fix documentation generation

## Sub-PRDs

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-17.1](prd-17.1-missing-products.md) | Draft | Missing Product Generation (interior, left/right contraction, scalar) |
| [PRD-17.2](prd-17.2-exterior-naming.md) | **Complete** | Exterior Product Naming (rename outer → exterior) |
| [PRD-17.3](prd-17.3-field-documentation.md) | **Complete** | Field Documentation Fix (use TOML field names) |
| [PRD-17.4](prd-17.4-guidance-updates.md) | **Complete** | Guidance Updates (CLAUDE.md, agents) |
| [PRD-17.5](prd-17.5-regenerate-algebras.md) | **Complete** | Regenerate All Algebras |
| [PRD-17.6](prd-17.6-antiproduct.md) | Draft | Geometric Antiproduct for PGA Transformations |
| [PRD-17.7](prd-17.7-public-products-api.md) | Draft | Expose Generated Products in Public API |

## Problem Statement

The current code generator (`clifford-codegen`) has several limitations that require manual workarounds in extension files:

### Issue 1: Missing Product Types

The generator currently produces:
- **Geometric product** (`geometric_*_*`)
- **Exterior product** (`outer_*_*`) - *incorrectly named "outer"*

But does NOT generate:
- **Interior product** (symmetric inner product)
- **Left contraction** (`A ⌋ B`)
- **Right contraction** (`A ⌊ B`)
- **Scalar product** (grade-0 extraction of geometric product)
- **Fat dot product** (Hestenes inner product)

This forces implementers to manually derive these products in extension files, which is error-prone and defeats the purpose of codegen.

### Issue 2: Inconsistent Naming ("outer" vs "exterior")

The standard mathematical terminology is **exterior product** (or **wedge product**), not "outer product". Current naming:
```rust
// Current (non-standard)
pub fn outer_point_point(...) -> Line

// Standard mathematical terminology
pub fn exterior_point_point(...) -> Line  // or wedge_point_point
```

### Issue 3: Incorrect Field Documentation

Generated type documentation uses canonical blade names from bitmask indices instead of the field names from the TOML spec:

```rust
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s`."]
    s: T,
    #[doc = "Coefficient of `e1e2`."]  // WRONG! Field is e23, should say "e2e3"
    e23: T,
    #[doc = "Coefficient of `e1e3`."]  // WRONG! Field is e31, should say "e3e1"
    e31: T,
    #[doc = "Coefficient of `e2e3`."]  // WRONG! Field is e12, should say "e1e2"
    e12: T,
    // ...
}
```

This happens because the codegen computes blade names from bitmask indices (`e1e2 = 0b0011`) instead of using the field names defined in the TOML specification.

### Issue 4: Extensions Using Manual Formulas

Due to missing products, extension files contain manually-derived formulas:

```rust
// From src/specialized/projective/dim3/extensions.rs
pub fn left_contract_plane(&self, plane: &Plane<T>) -> Point<T> {
    // Manual implementation because codegen doesn't generate left_contract_point_plane
    let e1 = -(plane.e023() * self.e03());
    // ...
}
```

These manual formulas are:
1. Error-prone (the Plane::meet bug was a parameter ordering issue)
2. Not verified against Multivector
3. Inconsistent with the "use codegen, not hand-rolled formulas" principle

## Solution

### Phase 1: Add Missing Products

#### 1.1 Product Type Definitions

Add to `ProductKind` enum:

```rust
pub enum ProductKind {
    Geometric,
    Exterior,        // Renamed from Outer
    Interior,        // New: symmetric inner product
    LeftContraction, // Exists but needs full generation
    RightContraction,// New
    Scalar,          // Exists but needs full generation
    FatDot,          // New: Hestenes inner product
}
```

#### 1.2 Grade Selection Rules

| Product | Output Grade | Condition |
|---------|--------------|-----------|
| Geometric | `|r - s|` to `r + s` | All grades |
| Exterior | `r + s` | Grade sum only |
| Interior | `|r - s|` | Grade diff (symmetric) |
| Left Contraction | `s - r` | `r <= s` only |
| Right Contraction | `r - s` | `s <= r` only |
| Scalar | `0` | Grade 0 only |
| Fat Dot | `|r - s|` | Like interior, but returns 0 if grades equal for vectors |

#### 1.3 Implementation

```rust
fn grade_filter(&self, grade_a: usize, grade_b: usize, result_grade: usize, kind: ProductKind) -> bool {
    match kind {
        ProductKind::Geometric => true,
        ProductKind::Exterior => result_grade == grade_a + grade_b,
        ProductKind::Interior => result_grade == grade_a.abs_diff(grade_b),
        ProductKind::LeftContraction => {
            grade_a <= grade_b && result_grade == grade_b - grade_a
        }
        ProductKind::RightContraction => {
            grade_b <= grade_a && result_grade == grade_a - grade_b
        }
        ProductKind::Scalar => result_grade == 0,
        ProductKind::FatDot => {
            // Hestenes inner: grade diff, but 0 for equal-grade vectors
            result_grade == grade_a.abs_diff(grade_b)
        }
    }
}
```

### Phase 2: Rename Outer to Exterior

#### 2.1 Code Changes

- Rename `ProductKind::Outer` to `ProductKind::Exterior`
- Rename generated functions: `outer_*_*` → `exterior_*_*`
- Update documentation to use "exterior product" terminology
- Update TOML product specifications

#### 2.2 Migration

Since this is a breaking change for anyone using generated code directly:

1. Keep `outer_*_*` as deprecated aliases initially
2. Document the change prominently
3. Remove aliases in next major version

```rust
#[deprecated(since = "0.2.0", note = "use exterior_* instead")]
#[inline]
pub fn outer_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    exterior_point_point(a, b)
}
```

### Phase 3: Fix Field Documentation

#### 3.1 Use Field Names, Not Blade Names

Change documentation generation to use the field name from the TOML spec:

```rust
fn generate_field_doc(field: &FieldSpec) -> String {
    // Use the field name directly, not the computed blade name
    format!("Coefficient of `{}`.", field.name)
}
```

Result:
```rust
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s`."]
    s: T,
    #[doc = "Coefficient of `e23`."]  // Correct!
    e23: T,
    // ...
}
```

#### 3.2 Optional: Add Blade Basis Expansion

For educational value, could add the basis blade expansion:

```rust
#[doc = "Coefficient of `e23` (basis: e2 ∧ e3)."]
e23: T,
```

### Phase 4: Update CLAUDE.md and Agents

#### 4.1 CLAUDE.md Updates

Add guidance on using generated products in extensions:

```markdown
### Using Generated Products in Extensions

When implementing domain-specific methods in extension files, **always use generated products** instead of manual formulas:

```rust
// CORRECT: Use generated product
use super::generated::products;

impl<T: Float> Point<T> {
    pub fn left_contract_plane(&self, plane: &Plane<T>) -> T {
        // Use generated left contraction
        products::left_contract_point_plane(self, plane)
    }
}

// WRONG: Manual formula (error-prone)
impl<T: Float> Point<T> {
    pub fn left_contract_plane(&self, plane: &Plane<T>) -> T {
        // Don't do this - prone to sign/ordering errors
        self.e1() * plane.e023() + self.e2() * plane.e031() + ...
    }
}
```

**Exceptions** (when manual formulas are acceptable):
1. The product combination isn't generated (file an issue!)
2. Performance-critical code where the generated form is suboptimal
3. Geometric shortcuts (e.g., `distance()` using specialized formula)

Even for exceptions, add a comment citing the source of the formula.
```

#### 4.2 Agent Updates

Update `implement.md`:

```markdown
### Using Generated Products

When implementing methods that involve algebraic products:

1. **Check if a generated product exists** in `generated/products.rs`
2. **Use the generated product** rather than deriving manually
3. **If no suitable product exists**, file an issue for codegen enhancement

```rust
// Check for existing products
use super::generated::products;

// Use generated left contraction
let result = products::left_contract_point_plane(&point, &plane);
```

**Red flags in code review**:
- Manual product formulas in extension files
- Sign corrections or coefficient adjustments
- Comments like "// manual implementation because codegen doesn't..."

These indicate missing codegen features that should be added.
```

Update `review.md`:

```markdown
### Generated Products Check

- [ ] **Uses generated products** - Extension methods use `products::*` functions
- [ ] **No manual formulas** - Algebraic products aren't hand-rolled
- [ ] **Missing products documented** - If manual formula needed, issue filed for codegen
```

## Deliverables

### Phase 1: Missing Products
- [ ] Add `Interior` product kind and generation
- [ ] Add `RightContraction` product kind and generation
- [ ] Complete `LeftContraction` generation for all type pairs
- [ ] Add `Scalar` product generation for all type pairs
- [ ] Add `FatDot` product kind (optional, based on need)
- [ ] Update TOML spec format to support new product types
- [ ] Generate verification tests for all new products

### Phase 2: Naming
- [ ] Rename `Outer` → `Exterior` in ProductKind
- [ ] Update function naming: `outer_*` → `exterior_*`
- [ ] Add deprecated aliases for backward compatibility
- [ ] Update all documentation to use "exterior" terminology

### Phase 3: Documentation
- [ ] Fix field documentation to use TOML field names
- [ ] Optionally add basis blade expansion in docs
- [ ] Regenerate all existing algebras

### Phase 4: Guidance
- [ ] Update CLAUDE.md with "Using Generated Products" section
- [ ] Update implement.md agent with generated product guidance
- [ ] Update review.md agent with generated product checklist

## TOML Specification Changes

### New Product Sections

```toml
[products]
# Existing
geometric = true
exterior = true    # renamed from outer

# New
interior = true
left_contraction = true
right_contraction = true
scalar = true

# Explicit product pairs (optional, inferred if not specified)
[products.left_contraction]
Point_Plane = "Scalar"
Point_Line = "Point"
Line_Plane = "Point"
```

## Testing Strategy

### Verification Against Multivector

All new products must have verification tests:

```rust
proptest! {
    #[test]
    fn left_contraction_matches_multivector(
        a in any::<Point<f64>>(),
        b in any::<Plane<f64>>()
    ) {
        let specialized = left_contract_point_plane(&a, &b);

        let mv_a = Multivector::from(a);
        let mv_b = Multivector::from(b);
        let mv_result = mv_a.left_contract(&mv_b);

        // Compare results
        prop_assert!(abs_diff_eq!(specialized, mv_result.scalar(), epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

### Algebraic Property Tests

```rust
proptest! {
    // Left contraction is not commutative
    #[test]
    fn left_contraction_not_commutative(a in any::<Line<f64>>(), b in any::<Plane<f64>>()) {
        let ab = left_contract_line_plane(&a, &b);
        // b ⌋ a has different type, so this tests grade behavior
    }

    // Exterior product is anticommutative for vectors
    #[test]
    fn exterior_anticommutative(a in any::<Point<f64>>(), b in any::<Point<f64>>()) {
        let ab = exterior_point_point(&a, &b);
        let ba = exterior_point_point(&b, &a);
        prop_assert!(abs_diff_eq!(ab, -ba, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Success Criteria

1. **Completeness**: All standard GA products are generated
2. **Correctness**: All products pass verification against Multivector
3. **Consistency**: Naming follows mathematical conventions (exterior, not outer)
4. **Documentation**: Field docs match TOML field names
5. **Guidance**: CLAUDE.md and agents direct users to use generated products

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/codegen/products.rs` | Update - Add new product kinds |
| `crates/clifford-codegen/src/codegen/types.rs` | Update - Fix field documentation |
| `crates/clifford-codegen/src/spec/raw.rs` | Update - New TOML fields |
| `crates/clifford-codegen/src/spec/ir.rs` | Update - New product configuration |
| `algebras/*.toml` | Update - Add new product sections |
| `src/generated/*/products.rs` | Regenerate - New products, renamed functions |
| `src/specialized/*/extensions.rs` | Update - Use generated products |
| `CLAUDE.md` | Update - Add codegen guidance |
| `.claude/agents/implement.md` | Update - Add product usage guidance |
| `.claude/agents/review.md` | Update - Add product review checklist |

## References

- [Geometric Algebra Primer](http://www.jaapsuter.com/geometric-algebra.pdf) - Product definitions
- [Clifford Algebra to Geometric Calculus](https://www.springer.com/gp/book/9789027725615) - Inner product variations
- [ganja.js Products](https://github.com/enkimute/ganja.js/blob/master/ganja.js) - Reference implementation
