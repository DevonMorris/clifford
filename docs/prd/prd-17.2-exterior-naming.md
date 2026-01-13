# PRD-17.2: Exterior Product Naming

**Status**: Draft
**Parent**: [PRD-17: Codegen Product Completeness and Documentation](prd-17-codegen-products.md)
**Goal**: Rename "outer product" to "exterior product" for mathematical consistency

## Problem Statement

The current codegen uses "outer" for the wedge product:

```rust
pub fn outer_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T>
```

However, the standard mathematical terminology is **exterior product** (or **wedge product**), not "outer product". The term "outer product" in mathematics typically refers to different operations:

1. **Tensor outer product**: `a ⊗ b` producing a matrix/tensor
2. **Cross product**: Sometimes called "outer" in older texts

Using "outer" creates confusion with these other meanings. The GA community and literature consistently use "exterior" or "wedge".

## Terminology

| Term | Symbol | GA Meaning | Other Meanings |
|------|--------|------------|----------------|
| **Exterior/Wedge** | `∧` | Grade-raising product | None (unambiguous) |
| Outer | various | Ambiguous | Tensor product, cross product |
| Inner | `·` | Grade-lowering product | Dot product (acceptable) |

**Standard GA texts use "exterior"**:
- Hestenes & Sobczyk, "Clifford Algebra to Geometric Calculus"
- Dorst, Fontijne & Mann, "Geometric Algebra for Computer Science"
- Doran & Lasenby, "Geometric Algebra for Physicists"

## Solution

### Rename ProductKind

```rust
// Before
pub enum ProductKind {
    Geometric,
    Outer,        // ← non-standard
    // ...
}

// After
pub enum ProductKind {
    Geometric,
    Exterior,     // ← standard mathematical term
    // ...
}
```

### Rename Generated Functions

```rust
// Before
pub fn outer_point_point<T: Float>(...) -> Line<T>
pub fn outer_line_point<T: Float>(...) -> Plane<T>

// After
pub fn exterior_point_point<T: Float>(...) -> Line<T>
pub fn exterior_line_point<T: Float>(...) -> Plane<T>
```

### Update Documentation

```rust
// Before
#[doc = "Outer product: Point ∧ Point -> Line"]

// After
#[doc = "Exterior product: Point ∧ Point -> Line"]
```

## Migration Strategy

Since this is a breaking change for users of generated code, provide a migration path:

### Phase 1: Add Deprecated Aliases

```rust
/// Exterior product: Point ∧ Point → Line
#[inline]
pub fn exterior_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    // ... implementation
}

/// Outer product (deprecated alias for exterior_point_point)
#[deprecated(since = "0.2.0", note = "use exterior_point_point instead")]
#[inline]
pub fn outer_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    exterior_point_point(a, b)
}
```

### Phase 2: Warn on Use (Next Minor Release)

The deprecation warning will alert users to migrate.

### Phase 3: Remove Aliases (Next Major Release)

Remove the `outer_*` functions entirely.

## Implementation

### ProductKind Rename

```rust
// In crates/clifford-codegen/src/codegen/products.rs

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProductKind {
    Geometric,
    Exterior,       // renamed from Outer
    Interior,
    LeftContraction,
    RightContraction,
    Scalar,
}

impl ProductKind {
    pub fn function_prefix(&self) -> &'static str {
        match self {
            ProductKind::Geometric => "geometric",
            ProductKind::Exterior => "exterior",  // changed from "outer"
            ProductKind::Interior => "interior",
            ProductKind::LeftContraction => "left_contract",
            ProductKind::RightContraction => "right_contract",
            ProductKind::Scalar => "scalar",
        }
    }

    pub fn doc_name(&self) -> &'static str {
        match self {
            ProductKind::Geometric => "Geometric product",
            ProductKind::Exterior => "Exterior product",  // changed from "Outer product"
            ProductKind::Interior => "Interior product",
            ProductKind::LeftContraction => "Left contraction",
            ProductKind::RightContraction => "Right contraction",
            ProductKind::Scalar => "Scalar product",
        }
    }

    pub fn symbol(&self) -> &'static str {
        match self {
            ProductKind::Geometric => "*",
            ProductKind::Exterior => "∧",
            ProductKind::Interior => "·",
            ProductKind::LeftContraction => "⌋",
            ProductKind::RightContraction => "⌊",
            ProductKind::Scalar => "*₀",
        }
    }
}
```

### TOML Spec Update

```toml
# Before
[products]
geometric = true
outer = true

# After
[products]
geometric = true
exterior = true

# Backward compatibility: accept both during transition
# Parser treats "outer" as alias for "exterior" with warning
```

### Parser Compatibility

```rust
// In crates/clifford-codegen/src/spec/parser.rs

fn parse_products(raw: &RawProducts) -> ProductsConfig {
    ProductsConfig {
        geometric: raw.geometric.unwrap_or(true),
        exterior: raw.exterior.or(raw.outer).unwrap_or(true),
        // If "outer" is used, emit deprecation warning
        // ...
    }
}
```

### Deprecated Alias Generation

```rust
fn generate_deprecated_aliases(&self) -> TokenStream {
    if !self.config.generate_deprecated_aliases {
        return quote! {};
    }

    let mut aliases = Vec::new();

    for product in &self.exterior_products {
        let exterior_name = &product.function_name;  // exterior_point_point
        let outer_name = format_ident!(
            "outer_{}",
            exterior_name.to_string().strip_prefix("exterior_").unwrap()
        );

        let a_type = &product.type_a;
        let b_type = &product.type_b;
        let output = &product.output_type;

        aliases.push(quote! {
            #[deprecated(since = "0.2.0", note = "use #exterior_name instead")]
            #[doc(hidden)]
            #[inline]
            pub fn #outer_name<T: Float>(a: &#a_type<T>, b: &#b_type<T>) -> #output<T> {
                #exterior_name(a, b)
            }
        });
    }

    quote! {
        // ============================================================
        // Deprecated Aliases (will be removed in next major version)
        // ============================================================
        #(#aliases)*
    }
}
```

## Extension File Updates

Update extension files to use the new naming:

```rust
// Before (in extensions.rs)
use super::generated::products::outer_point_point;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        outer_point_point(self, other)
    }
}

// After
use super::generated::products::exterior_point_point;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        exterior_point_point(self, other)
    }
}
```

## Documentation Updates

### Module Docs

```rust
//! # Product Operations
//!
//! This module provides the fundamental products of geometric algebra:
//!
//! | Product | Symbol | Function Prefix | Description |
//! |---------|--------|-----------------|-------------|
//! | Geometric | `*` | `geometric_` | Full product |
//! | Exterior | `∧` | `exterior_` | Grade-raising (wedge) |
//! | Interior | `·` | `interior_` | Grade-lowering (symmetric) |
//! | Left contraction | `⌋` | `left_contract_` | Asymmetric contraction |
//! | Right contraction | `⌊` | `right_contract_` | Asymmetric contraction |
//! | Scalar | `*₀` | `scalar_` | Grade-0 extraction |
```

### Function Docs

```rust
/// Exterior product: Point ∧ Point → Line
///
/// The exterior (wedge) product of two points yields a line through both points.
/// This is the fundamental operation for constructing higher-dimensional
/// geometric objects from lower-dimensional ones.
///
/// # Mathematical Definition
///
/// For blades A and B with grades r and s:
/// ```text
/// A ∧ B = ⟨AB⟩_{r+s}
/// ```
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim3::{Point, Line};
/// use clifford::specialized::projective::dim3::generated::products::exterior_point_point;
///
/// let p1 = Point::from_cartesian(0.0, 0.0, 0.0);
/// let p2 = Point::from_cartesian(1.0, 0.0, 0.0);
/// let line = exterior_point_point(&p1, &p2);
/// // line passes through origin in x-direction
/// ```
#[inline]
pub fn exterior_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    // ...
}
```

## Testing

### Verify Aliases Work

```rust
#[test]
#[allow(deprecated)]
fn deprecated_outer_alias_works() {
    let p1 = Point::from_cartesian(1.0, 0.0, 0.0);
    let p2 = Point::from_cartesian(0.0, 1.0, 0.0);

    let via_exterior = exterior_point_point(&p1, &p2);
    let via_outer = outer_point_point(&p1, &p2);  // deprecated

    assert_eq!(via_exterior, via_outer);
}
```

### Verify Deprecation Warning

```rust
#[test]
fn outer_produces_deprecation_warning() {
    // This test just ensures the code compiles with the warning
    // The actual warning is checked by CI/clippy
}
```

## Deliverables

- [ ] Rename `ProductKind::Outer` to `ProductKind::Exterior`
- [ ] Update function prefix from `outer_` to `exterior_`
- [ ] Update documentation to use "exterior" terminology
- [ ] Add deprecated `outer_*` aliases
- [ ] Update TOML spec to use `exterior` (with `outer` as deprecated alias)
- [ ] Update extension files to use `exterior_*` functions
- [ ] Update CLAUDE.md examples to use `exterior_*`
- [ ] Add deprecation tests
- [ ] Regenerate all algebras

## Success Criteria

1. **Consistent terminology**: All code uses "exterior", not "outer"
2. **Backward compatible**: Old `outer_*` functions still work (with warning)
3. **Documentation clear**: All docs use standard GA terminology

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/codegen/products.rs` | Update - Rename Outer → Exterior |
| `crates/clifford-codegen/src/spec/raw.rs` | Update - exterior field, outer deprecated |
| `algebras/*.toml` | Update - Use `exterior = true` |
| `src/generated/*/products.rs` | Regenerate - New function names |
| `src/specialized/*/extensions.rs` | Update - Use exterior_* imports |
| `CLAUDE.md` | Update - Use exterior in examples |
| `.claude/agents/*.md` | Update - Use exterior in examples |

## Timeline

1. **v0.2.0**: Add `exterior_*` functions, add deprecated `outer_*` aliases
2. **v0.3.0**: Emit compiler warnings for `outer_*` usage
3. **v1.0.0**: Remove `outer_*` aliases entirely
