# PRD-18.9: Generate Wrapper Type Aliases

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Update codegen to generate type aliases for wrapper types (e.g., `type UnitVector<T> = Unit<Vector<T>>`)

## Overview

The wrapper types (`Unit<T>`, `Bulk<T>`, `Ideal<T>`, `Proper<T>`) are defined in `src/wrappers.rs`, but users must manually construct the wrapped types:

```rust
// Current: verbose
let unit: Unit<Vector<f64>> = Unit::new_normalize(v);
let rigid: Bulk<Motor<f64>> = Bulk::new_normalize(m);
```

This PRD adds codegen support to generate convenient type aliases:

```rust
// With aliases: cleaner
let unit: UnitVector<f64> = Unit::new_normalize(v);
let rigid: BulkMotor<f64> = Bulk::new_normalize(m);
```

## Deliverables

### 1. Define Wrapper Alias Configuration

**File:** `crates/clifford-codegen/src/spec/ir.rs`

Add configuration for which wrappers apply to each type:

```rust
/// Specifies which wrapper type aliases to generate for a type.
#[derive(Debug, Clone)]
pub struct WrapperConfig {
    /// Generate `Unit<Type>` alias (for Euclidean types)
    pub unit: bool,
    /// Generate `Bulk<Type>` alias (for PGA versors)
    pub bulk: bool,
    /// Generate `Ideal<Type>` alias (for PGA homogeneous types)
    pub ideal: bool,
    /// Generate `Proper<Type>` alias (for Minkowski timelike)
    pub proper: bool,
}

impl WrapperConfig {
    /// Infer wrapper config based on algebra and type characteristics.
    pub fn infer(type_spec: &TypeSpec, algebra: &AlgebraSpec) -> Self {
        let is_pga = !algebra.signature.zero.is_empty();
        let is_minkowski = !algebra.signature.negative.is_empty()
            && !algebra.signature.positive.is_empty();
        let is_versor = type_spec.versor;
        let is_homogeneous = type_spec.grades.iter().any(|&g| g == 1 || g == 3);

        Self {
            // Euclidean types get Unit<T>
            unit: !is_pga && !is_minkowski,
            // PGA versors get Bulk<T>
            bulk: is_pga && is_versor,
            // PGA homogeneous types (points, planes) get Ideal<T>
            ideal: is_pga && is_homogeneous && !is_versor,
            // Minkowski vectors get Proper<T>
            proper: is_minkowski && type_spec.grades == vec![1],
        }
    }
}
```

### 2. Generate Type Aliases

**File:** `crates/clifford-codegen/src/codegen/types.rs`

Add function to generate wrapper type aliases:

```rust
/// Generates wrapper type aliases for a type.
pub fn generate_wrapper_aliases(type_spec: &TypeSpec, algebra: &AlgebraSpec) -> TokenStream {
    let config = WrapperConfig::infer(type_spec, algebra);
    let type_name = format_ident!("{}", type_spec.name);
    let mut aliases = TokenStream::new();

    if config.unit {
        let alias_name = format_ident!("Unit{}", type_spec.name);
        aliases.extend(quote! {
            #[doc = concat!("A unit ", stringify!(#type_name), " (Euclidean norm = 1).")]
            #[doc = ""]
            #[doc = concat!("This is an alias for `Unit<", stringify!(#type_name), "<T>>`.")]
            #[doc = ""]
            #[doc = "# Example"]
            #[doc = ""]
            #[doc = "```"]
            #[doc = concat!("use clifford::specialized::euclidean::dim3::{", stringify!(#type_name), ", ", stringify!(#alias_name), "};")]
            #[doc = "use clifford::wrappers::Unit;"]
            #[doc = ""]
            #[doc = concat!("let v = ", stringify!(#type_name), "::new(3.0, 4.0, 0.0);")]
            #[doc = concat!("let unit: ", stringify!(#alias_name), "<f64> = Unit::new_normalize(v);")]
            #[doc = "assert!((unit.norm() - 1.0).abs() < 1e-10);"]
            #[doc = "```"]
            pub type #alias_name<T> = crate::wrappers::Unit<#type_name<T>>;
        });
    }

    if config.bulk {
        let alias_name = format_ident!("Bulk{}", type_spec.name);
        aliases.extend(quote! {
            #[doc = concat!("A bulk-normalized ", stringify!(#type_name), " (bulk norm = 1).")]
            #[doc = ""]
            #[doc = concat!("This is an alias for `Bulk<", stringify!(#type_name), "<T>>`.")]
            #[doc = ""]
            #[doc = "For PGA versors, bulk normalization ensures rigid body transformations"]
            #[doc = "(rotation + translation, no scaling)."]
            pub type #alias_name<T> = crate::wrappers::Bulk<#type_name<T>>;
        });
    }

    if config.ideal {
        let alias_name = format_ident!("Ideal{}", type_spec.name);
        aliases.extend(quote! {
            #[doc = concat!("An ", stringify!(#type_name), " in standard homogeneous form (weight norm = 1).")]
            #[doc = ""]
            #[doc = concat!("This is an alias for `Ideal<", stringify!(#type_name), "<T>>`.")]
            #[doc = ""]
            #[doc = "For PGA homogeneous coordinates, weight normalization gives the standard"]
            #[doc = "form where the projective weight is 1."]
            pub type #alias_name<T> = crate::wrappers::Ideal<#type_name<T>>;
        });
    }

    if config.proper {
        let alias_name = format_ident!("Proper{}", type_spec.name);
        aliases.extend(quote! {
            #[doc = concat!("A proper (timelike) ", stringify!(#type_name), ".")]
            #[doc = ""]
            #[doc = concat!("This is an alias for `Proper<", stringify!(#type_name), "<T>>`.")]
            #[doc = ""]
            #[doc = "In Minkowski spacetime, timelike vectors have positive squared norm"]
            #[doc = "and represent valid 4-velocities."]
            pub type #alias_name<T> = crate::wrappers::Proper<#type_name<T>>;
        });
    }

    aliases
}
```

### 3. Include Aliases in Generated Types File

**File:** `crates/clifford-codegen/src/codegen/types.rs`

Update the main generation function:

```rust
pub fn generate_types_file(spec: &AlgebraSpec) -> TokenStream {
    let mut output = TokenStream::new();

    // File header
    output.extend(generate_types_header(spec));

    // Type definitions
    for type_spec in &spec.types {
        output.extend(generate_type_definition(type_spec));
        output.extend(generate_type_impl(type_spec));
    }

    // NEW: Wrapper type aliases at the end of the file
    output.extend(quote! {
        // ============================================================================
        // Wrapper Type Aliases
        // ============================================================================
    });

    for type_spec in &spec.types {
        output.extend(generate_wrapper_aliases(type_spec, spec));
    }

    output
}
```

### 4. Expected Generated Output

#### Euclidean 3D (`src/specialized/euclidean/dim3/generated/types.rs`)

```rust
// ... type definitions ...

// ============================================================================
// Wrapper Type Aliases
// ============================================================================

/// A unit Vector (Euclidean norm = 1).
///
/// This is an alias for `Unit<Vector<T>>`.
///
/// # Example
///
/// ```
/// use clifford::specialized::euclidean::dim3::{Vector, UnitVector};
/// use clifford::wrappers::Unit;
///
/// let v = Vector::new(3.0, 4.0, 0.0);
/// let unit: UnitVector<f64> = Unit::new_normalize(v);
/// assert!((unit.norm() - 1.0).abs() < 1e-10);
/// ```
pub type UnitVector<T> = crate::wrappers::Unit<Vector<T>>;

/// A unit Bivector (Euclidean norm = 1).
pub type UnitBivector<T> = crate::wrappers::Unit<Bivector<T>>;

/// A unit Trivector (Euclidean norm = 1).
pub type UnitTrivector<T> = crate::wrappers::Unit<Trivector<T>>;

/// A unit Rotor (Euclidean norm = 1).
pub type UnitRotor<T> = crate::wrappers::Unit<Rotor<T>>;
```

#### Projective 3D (`src/specialized/projective/dim3/generated/types.rs`)

```rust
// ... type definitions ...

// ============================================================================
// Wrapper Type Aliases
// ============================================================================

/// A Point in standard homogeneous form (weight norm = 1).
///
/// This is an alias for `Ideal<Point<T>>`.
pub type IdealPoint<T> = crate::wrappers::Ideal<Point<T>>;

/// A Line in standard form.
///
/// Note: Lines can use either Bulk (direction normalized) or Ideal (moment normalized).
pub type BulkLine<T> = crate::wrappers::Bulk<Line<T>>;

/// A Plane in standard homogeneous form (weight norm = 1).
pub type IdealPlane<T> = crate::wrappers::Ideal<Plane<T>>;

/// A bulk-normalized Motor (bulk norm = 1).
///
/// For PGA motors, bulk normalization ensures rigid body transformations
/// (rotation + translation, no scaling).
pub type BulkMotor<T> = crate::wrappers::Bulk<Motor<T>>;

/// A bulk-normalized Flector (bulk norm = 1).
pub type BulkFlector<T> = crate::wrappers::Bulk<Flector<T>>;
```

### 5. Re-export Aliases from Module

**File:** Generated `mod.rs` or ensure aliases are accessible

The aliases should be re-exported from the main module:

```rust
// In src/specialized/euclidean/dim3/mod.rs
pub use generated::types::{
    Vector, Bivector, Trivector, Rotor, Scalar,
    // Wrapper aliases
    UnitVector, UnitBivector, UnitTrivector, UnitRotor,
};
```

## Testing

### Compile Tests

```rust
#[test]
fn type_aliases_compile() {
    use crate::specialized::euclidean::dim3::{Vector, UnitVector};
    use crate::wrappers::Unit;

    // Both should work
    let _unit1: UnitVector<f64> = Unit::new_normalize(Vector::new(1.0, 0.0, 0.0));
    let _unit2: Unit<Vector<f64>> = Unit::new_normalize(Vector::new(1.0, 0.0, 0.0));

    // Should be the same type
    fn takes_unit_vector(_: UnitVector<f64>) {}
    takes_unit_vector(_unit1);
}

#[test]
fn pga_aliases_compile() {
    use crate::specialized::projective::dim3::{Motor, BulkMotor, Point, IdealPoint};
    use crate::wrappers::{Bulk, Ideal};

    let motor = Motor::identity();
    let _bulk: BulkMotor<f64> = Bulk::new_normalize(motor);

    let point = Point::from_cartesian(1.0, 2.0, 3.0);
    let _ideal: IdealPoint<f64> = Ideal::new_normalize(point);
}
```

### Documentation Tests

The generated doc examples should compile and run:

```rust
/// ```
/// use clifford::specialized::euclidean::dim3::{Vector, UnitVector};
/// use clifford::wrappers::Unit;
///
/// let v = Vector::new(3.0, 4.0, 0.0);
/// let unit: UnitVector<f64> = Unit::new_normalize(v);
/// assert!((unit.norm() - 1.0).abs() < 1e-10);
/// ```
```

## Alias Naming Convention

| Type | Euclidean Alias | PGA Bulk Alias | PGA Ideal Alias |
|------|-----------------|----------------|-----------------|
| Vector | UnitVector | - | - |
| Bivector | UnitBivector | - | - |
| Trivector | UnitTrivector | - | - |
| Rotor | UnitRotor | - | - |
| Point | - | - | IdealPoint |
| Line | - | BulkLine | - |
| Plane | - | - | IdealPlane |
| Motor | - | BulkMotor | - |
| Flector | - | BulkFlector | - |

## Implementation Notes

### Lines Have Two Normalizations

Lines in PGA can be normalized two ways:
- **Bulk normalization**: Direction has unit length (useful for ray tracing)
- **Ideal/Weight normalization**: Standard Plücker coordinates

We generate `BulkLine` by default since direction normalization is more common. Users can manually use `Ideal<Line<T>>` if needed.

### Scalars Don't Get Aliases

Scalar types don't need wrapper aliases since:
- Unit scalar is just 1.0 or -1.0
- Normalization is trivial

### Future: Minkowski Types

When Minkowski algebras are added, generate:
```rust
pub type ProperFourVector<T> = crate::wrappers::Proper<FourVector<T>>;
```

## Success Criteria

1. Type aliases generated for all applicable types
2. Euclidean types get `Unit*` aliases
3. PGA versors get `Bulk*` aliases
4. PGA homogeneous types get `Ideal*` aliases
5. Aliases exported from module's public API
6. Documentation includes usage examples
7. Doc tests compile and pass
8. No naming conflicts with existing types

## Dependencies

- PRD-18.1 (Normed trait hierarchy)
- PRD-18.2 (Wrapper types)
- PRD-18.8 (Generate Normed impls) - Wrappers require types to implement `Normed`

## References

- `src/wrappers.rs` - Wrapper type definitions
- nalgebra's `Unit<T>` pattern for inspiration
