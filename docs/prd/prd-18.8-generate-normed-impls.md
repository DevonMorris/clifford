# PRD-18.8: Generate Normed Trait Implementations

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Update codegen to generate `impl Normed for Type<T>` blocks for all generated types

## Overview

The `Normed` trait hierarchy is defined in `src/norm.rs` but the generated types don't implement these traits. Each generated type has standalone methods (`norm_squared()`, `norm()`, `try_normalize()`, `normalize()`, `scale()`) but doesn't implement the `Normed` trait.

This PRD adds codegen support to generate proper trait implementations, enabling use of the wrapper types (`Unit<T>`, `Bulk<T>`, etc.).

## Problem

Currently, generated types have:
```rust
impl<T: Float> Vector<T> {
    pub fn norm_squared(&self) -> T { ... }
    pub fn norm(&self) -> T { ... }
    pub fn try_normalize(&self) -> Option<Self> { ... }
    pub fn normalize(&self) -> Self { ... }
    pub fn scale(&self, s: T) -> Self { ... }
}
```

But they don't implement `Normed`:
```rust
// This doesn't exist yet!
impl<T: Float> Normed for Vector<T> {
    type Scalar = T;
    fn norm_squared(&self) -> T { ... }
    fn try_normalize(&self) -> Option<Self> { ... }
    fn scale(&self, factor: T) -> Self { ... }
}
```

This prevents using wrapper types:
```rust
// ERROR: Unit requires T: Normed
let unit: Unit<Vector<f64>> = Unit::new_normalize(v);
```

## Deliverables

### 1. Update Codegen to Generate `Normed` Implementation

**File:** `crates/clifford-codegen/src/codegen/traits.rs`

Add a new function to generate `Normed` trait implementations:

```rust
/// Generates `impl Normed for Type<T>` for all types.
pub fn generate_normed_impl(type_spec: &TypeSpec) -> TokenStream {
    let type_name = format_ident!("{}", type_spec.name);
    let fields = &type_spec.fields;

    // Generate norm_squared: sum of squares of all field components
    let norm_squared_terms: Vec<_> = fields.iter().map(|f| {
        let field_name = format_ident!("{}", f.name);
        quote! { self.#field_name * self.#field_name }
    }).collect();

    // Generate scale: multiply each field by factor
    let scale_fields: Vec<_> = fields.iter().map(|f| {
        let field_name = format_ident!("{}", f.name);
        quote! { self.#field_name * factor }
    }).collect();

    let field_names: Vec<_> = fields.iter().map(|f| format_ident!("{}", f.name)).collect();

    quote! {
        impl<T: Float> crate::norm::Normed for #type_name<T> {
            type Scalar = T;

            #[inline]
            fn norm_squared(&self) -> T {
                #(#norm_squared_terms)+*
            }

            fn try_normalize(&self) -> Option<Self> {
                let n = self.norm();
                if n < T::epsilon() {
                    None
                } else {
                    Some(self.scale(T::one() / n))
                }
            }

            #[inline]
            fn scale(&self, factor: T) -> Self {
                Self::new(#(#scale_fields),*)
            }
        }
    }
}
```

### 2. Generate `DegenerateNormed` for PGA Types

For PGA types, also generate `DegenerateNormed`:

```rust
/// Generates `impl DegenerateNormed for Type<T>` for PGA types.
pub fn generate_degenerate_normed_impl(
    type_spec: &TypeSpec,
    algebra: &AlgebraSpec,
) -> Option<TokenStream> {
    // Only generate for algebras with degenerate basis (zero signature)
    if algebra.signature.zero.is_empty() {
        return None;
    }

    let type_name = format_ident!("{}", type_spec.name);

    // Bulk norm: components NOT involving e0 (degenerate basis)
    // Weight norm: components involving e0 via antidot product

    // This requires knowledge of which fields involve e0
    let bulk_fields = get_bulk_fields(type_spec, algebra);
    let weight_norm_expr = compute_weight_norm_expr(type_spec, algebra);

    let bulk_norm_terms: Vec<_> = bulk_fields.iter().map(|f| {
        let field_name = format_ident!("{}", f);
        quote! { self.#field_name * self.#field_name }
    }).collect();

    Some(quote! {
        impl<T: Float> crate::norm::DegenerateNormed for #type_name<T> {
            fn bulk_norm_squared(&self) -> T {
                #(#bulk_norm_terms)+*
            }

            fn weight_norm_squared(&self) -> T {
                #weight_norm_expr
            }

            fn try_unitize(&self) -> Option<Self> {
                let w = self.weight_norm();
                if w < T::epsilon() {
                    None
                } else {
                    Some(self.scale(T::one() / w))
                }
            }
        }
    })
}

/// Identifies fields that don't involve the degenerate basis (bulk part).
fn get_bulk_fields(type_spec: &TypeSpec, algebra: &AlgebraSpec) -> Vec<String> {
    let degenerate_basis = &algebra.signature.zero;
    type_spec.fields.iter()
        .filter(|f| {
            // Field is bulk if its blade doesn't contain degenerate basis
            !f.blade_indices.iter().any(|idx| degenerate_basis.contains(idx))
        })
        .map(|f| f.name.clone())
        .collect()
}
```

### 3. Update `generate_traits` to Include Normed

**File:** `crates/clifford-codegen/src/codegen/traits.rs`

```rust
pub fn generate_traits(spec: &AlgebraSpec) -> TokenStream {
    let mut output = TokenStream::new();

    for type_spec in &spec.types {
        // Existing trait generation...
        output.extend(generate_approx_traits(type_spec));
        output.extend(generate_arbitrary_impl(type_spec));

        // NEW: Generate Normed implementations
        output.extend(generate_normed_impl(type_spec));

        // NEW: Generate DegenerateNormed for PGA types
        if let Some(degenerate_impl) = generate_degenerate_normed_impl(type_spec, spec) {
            output.extend(degenerate_impl);
        }
    }

    output
}
```

### 4. Remove Standalone Methods from Types

Once `Normed` is implemented, the standalone methods can delegate to the trait:

```rust
impl<T: Float> Vector<T> {
    /// Returns the squared Euclidean norm.
    #[inline]
    pub fn norm_squared(&self) -> T {
        <Self as crate::norm::Normed>::norm_squared(self)
    }

    // ... other methods delegate similarly
}
```

Or keep both (standalone for direct use, trait for generic code).

### 5. Add Required Imports to Generated Code

**File:** `crates/clifford-codegen/src/codegen/mod.rs`

Ensure the generated file includes necessary imports:

```rust
fn generate_traits_file(spec: &AlgebraSpec) -> String {
    let traits = generate_traits(spec);

    quote! {
        //! Trait implementations for the algebra types.

        use crate::scalar::Float;
        use crate::norm::{Normed, DegenerateNormed};  // NEW

        #traits
    }.to_string()
}
```

## Testing

### Unit Tests

```rust
#[test]
fn vector_implements_normed() {
    use crate::norm::Normed;
    use crate::specialized::euclidean::dim3::Vector;

    let v = Vector::new(3.0, 4.0, 0.0);

    // Via trait
    assert!(relative_eq!(Normed::norm(&v), 5.0, epsilon = 1e-10, max_relative = 1e-10));

    // Via standalone method
    assert!(relative_eq!(v.norm(), 5.0, epsilon = 1e-10, max_relative = 1e-10));
}

#[test]
fn motor_implements_degenerate_normed() {
    use crate::norm::{Normed, DegenerateNormed};
    use crate::specialized::projective::dim3::Motor;

    let m = Motor::from_rotation_z(0.5);

    // Bulk norm (rotor part)
    assert!(relative_eq!(m.bulk_norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));

    // Pure rotation has zero weight
    assert!(relative_eq!(m.weight_norm(), 0.0, epsilon = 1e-10, max_relative = 1e-10));
}

#[test]
fn unit_wrapper_works_with_normed_types() {
    use crate::wrappers::Unit;
    use crate::specialized::euclidean::dim3::Vector;

    let v = Vector::new(3.0, 4.0, 0.0);
    let unit: Unit<Vector<f64>> = Unit::new_normalize(v);

    assert!(relative_eq!(unit.norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
}
```

### Property-Based Tests

```rust
proptest! {
    #[test]
    fn normed_trait_consistent_with_methods(v in any::<Vector<f64>>()) {
        use crate::norm::Normed;

        // Trait and method should return same result
        prop_assert!(relative_eq!(
            Normed::norm(&v),
            v.norm(),
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn unit_wrapper_has_unit_norm(v in any::<NonZeroVector<f64>>()) {
        use crate::wrappers::Unit;

        let unit = Unit::new_normalize(v.0);
        prop_assert!(relative_eq!(
            unit.norm(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
```

## Implementation Notes

### Bulk vs Weight Norm for PGA

For PGA types, the bulk and weight norms depend on the grade structure:

| Type | Bulk Fields | Weight Fields |
|------|-------------|---------------|
| Point | e1, e2, e3 | e0 |
| Line | e23, e31, e12 (direction) | e01, e02, e03 (moment) |
| Plane | e123 | e023, e031, e012 |
| Motor | s, e23, e31, e12 (rotor) | e01, e02, e03, e0123 (translator) |
| Flector | e1, e2, e3, e123 | e0, e023, e031, e012 |

The weight norm specifically uses the antidot product formula from Rigid GA.

### Algebra Type Detection

To determine which norm traits to implement:

```rust
fn get_norm_traits_for_algebra(algebra: &AlgebraSpec) -> Vec<&'static str> {
    let mut traits = vec!["Normed"];

    if !algebra.signature.zero.is_empty() {
        // Has degenerate basis -> PGA
        traits.push("DegenerateNormed");
    }

    if !algebra.signature.negative.is_empty() && !algebra.signature.positive.is_empty() {
        // Has mixed signature -> Minkowski
        traits.push("IndefiniteNormed");
    }

    // CGA detection would need additional markers

    traits
}
```

## Success Criteria

1. All generated types implement `Normed` trait
2. PGA types additionally implement `DegenerateNormed`
3. `Unit<T>` wrapper works with all Euclidean types
4. `Bulk<T>` wrapper works with all PGA versor types
5. `Ideal<T>` wrapper works with all PGA homogeneous types
6. Standalone methods remain available (backward compatibility)
7. No clippy warnings in generated code
8. All tests pass

## Dependencies

- PRD-18.1 (Normed trait hierarchy) - Must be implemented first
- PRD-18.2 (Wrapper types) - Must be implemented first

## References

- [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
- `src/norm.rs` - Trait definitions
- `src/wrappers.rs` - Wrapper type definitions
