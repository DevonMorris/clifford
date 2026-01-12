# PRD-17.7: Expose Generated Products in Public API

**Status**: Complete
**Parent**: [PRD-17](prd-17-codegen-products.md)
**Goal**: Remove dead code suppression by exposing generated products in the public API

## Problem Statement

The codegen generates many product functions (geometric, exterior, left_contraction, inner) for all type combinations. Currently, many of these are unused internally, leading to dead code warnings. We "fixed" this with:

```rust
#[allow(dead_code)] pub mod products;
```

This is a poor solution because:

1. **Hidden functionality**: Users can't discover these products
2. **Untested code**: Dead code may have bugs we never catch
3. **Wasted generation**: We generate code nobody can use
4. **Poor API design**: Products are implementation details, not public features

## Solution

Expose all generated products as part of the public API, making them discoverable and usable.

### Current Structure

```
src/specialized/euclidean/dim2/
├── mod.rs              # Re-exports types only
├── generated/
│   ├── mod.rs
│   ├── types.rs        # Scalar, Vector, Bivector, Rotor
│   ├── products.rs     # geometric_*, exterior_*, left_contract_*, inner_*
│   ├── traits.rs       # Operator impls, verification tests
│   └── conversions.rs  # From/Into Multivector
└── extensions.rs       # Domain-specific methods
```

### Proposed Structure

```
src/specialized/euclidean/dim2/
├── mod.rs              # Re-exports types AND products module
├── generated/
│   ├── mod.rs          # Remove #[allow(dead_code)]
│   ├── types.rs
│   ├── products.rs     # All products are public API
│   ├── traits.rs
│   └── conversions.rs
└── extensions.rs
```

### API Design

#### Option A: Re-export Products Module

```rust
// src/specialized/euclidean/dim2/mod.rs
pub use generated::types::{Bivector, Rotor, Scalar, Vector};
pub use generated::products;  // NEW: expose products module

// Usage:
use clifford::specialized::euclidean::dim2::{Vector, products};
let result = products::geometric_vector_vector(&v1, &v2);
```

#### Option B: Re-export Individual Functions (Too Verbose)

```rust
// Not recommended - too many functions
pub use generated::products::{
    geometric_vector_vector,
    geometric_vector_bivector,
    // ... 50+ functions
};
```

#### Option C: Prelude with Common Products

```rust
// src/specialized/euclidean/dim2/mod.rs
pub mod products {
    pub use super::generated::products::*;
}

// Also add a prelude for the most common ones
pub mod prelude {
    pub use super::{Vector, Bivector, Rotor, Scalar};
    pub use super::products::{
        geometric_vector_vector,
        exterior_vector_vector,
        // ... most commonly used
    };
}
```

**Recommendation**: Option A is simplest and most discoverable.

## Implementation

### Phase 1: Remove Dead Code Suppression

Update codegen to not emit `#[allow(dead_code)]`:

```rust
// crates/clifford-codegen/src/main.rs
fn generate_mod_file(...) -> TokenStream {
    let mut mods = vec![
        quote! { pub mod types; },
        quote! { pub mod products; },  // Remove #[allow(dead_code)]
        quote! { pub mod traits; },
        quote! { pub mod conversions; },
    ];
    // ...
}
```

### Phase 2: Update Module Re-exports

Update each specialized module to re-export products:

```rust
// src/specialized/euclidean/dim2/mod.rs

// Generated code (do not edit manually)
mod generated;

// Domain-specific extensions
mod extensions;

// Re-export generated types
pub use generated::types::{Bivector, Rotor, Scalar, Vector};

// Re-export products module for direct access
pub use generated::products;

// Re-export Even alias from extensions
pub use extensions::Even;
```

### Phase 3: Documentation

Add module-level documentation for products:

```rust
// In generated products.rs
//! Product functions for euclidean2 algebra.
//!
//! This module provides all algebraic products between types in the algebra.
//! Each function is named `{product}_{lhs}_{rhs}` where:
//! - `product` is one of: `geometric`, `exterior`, `left_contract`, `inner`
//! - `lhs` and `rhs` are the input type names in lowercase
//!
//! # Example
//!
//! ```
//! use clifford::specialized::euclidean::dim2::{Vector, products};
//!
//! let v1 = Vector::new(1.0, 0.0);
//! let v2 = Vector::new(0.0, 1.0);
//!
//! // Geometric product: Vector × Vector → Rotor
//! let rotor = products::geometric_vector_vector(&v1, &v2);
//!
//! // Exterior product: Vector ∧ Vector → Bivector
//! let bivector = products::exterior_vector_vector(&v1, &v2);
//!
//! // Inner product: Vector · Vector → Scalar
//! let scalar = products::inner_vector_vector(&v1, &v2);
//! ```
//!
//! # Available Products
//!
//! | Product | Symbol | Description |
//! |---------|--------|-------------|
//! | `geometric_*` | `×` | Full geometric product |
//! | `exterior_*` | `∧` | Wedge/exterior product |
//! | `left_contract_*` | `⌋` | Left contraction |
//! | `inner_*` | `·` | Symmetric inner product |
```

### Phase 4: Regenerate All Algebras

```bash
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" -o "appropriate/path" --force
done
```

## Deliverables

- [ ] Remove `#[allow(dead_code)]` from codegen mod.rs generation
- [ ] Update euclidean2 mod.rs to re-export products
- [ ] Update euclidean3 mod.rs to re-export products
- [ ] Update projective3 mod.rs to re-export products
- [ ] Add products module documentation to codegen
- [ ] Regenerate all algebras
- [ ] Verify no dead code warnings
- [ ] Verify clippy passes
- [ ] Add usage examples to module docs

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/main.rs` | Remove `#[allow(dead_code)]` |
| `crates/clifford-codegen/src/codegen/products.rs` | Add module documentation |
| `src/specialized/euclidean/dim2/mod.rs` | Add `pub use generated::products` |
| `src/specialized/euclidean/dim3/mod.rs` | Add `pub use generated::products` |
| `src/specialized/projective/dim3/mod.rs` | Add `pub use generated::products` |
| `src/specialized/*/generated/*` | Regenerate |

## Testing Strategy

1. **Compilation**: Verify no dead code warnings after changes
2. **Documentation**: Run `cargo doc` and verify products are documented
3. **Usage**: Add doc tests showing product usage
4. **Clippy**: Verify `cargo clippy` passes

## Success Criteria

1. No `#[allow(dead_code)]` attributes in generated code
2. All products accessible via `module::products::*`
3. Products module has comprehensive documentation
4. Doc tests demonstrate usage patterns
5. Clippy and fmt pass

## API Examples After Implementation

```rust
use clifford::specialized::euclidean::dim3::{Vector, Bivector, products};

// Geometric product
let v1 = Vector::new(1.0, 2.0, 3.0);
let v2 = Vector::new(4.0, 5.0, 6.0);
let rotor = products::geometric_vector_vector(&v1, &v2);

// Exterior product (wedge)
let bivector = products::exterior_vector_vector(&v1, &v2);

// Left contraction
let scalar = products::left_contract_vector_vector(&v1, &v2);

// Chaining products
let b = Bivector::new(1.0, 0.0, 0.0);
let result = products::geometric_bivector_vector(&b, &v1);
```

```rust
use clifford::specialized::projective::dim3::{Point, Line, Plane, products};

// Join points to make a line (exterior product)
let p1 = Point::from_cartesian(0.0, 0.0, 0.0);
let p2 = Point::from_cartesian(1.0, 0.0, 0.0);
let line = products::exterior_point_point(&p1, &p2);

// Join line and point to make a plane
let p3 = Point::from_cartesian(0.0, 1.0, 0.0);
let plane = products::exterior_line_point(&line, &p3);
```
