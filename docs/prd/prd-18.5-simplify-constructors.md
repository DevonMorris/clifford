# PRD-18.5: Simplify Constructors

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Remove `new()` constructor, keep only `new_unchecked()` and `new_checked()`

## Overview

Current generated constructors are confusing:
- `new(params)` - Solves constraint, returns `Option<Self>` for quadratic
- `new_checked(all_fields, tolerance)` - Validates constraint
- `new_unchecked(all_fields)` - Raw construction

The `new()` constructor is problematic because:
1. It solves constraints, which is counterintuitive
2. Most users want raw construction
3. The name suggests it's the "normal" constructor

This PRD simplifies to just two constructors:
- `new_unchecked()` - Raw construction, no validation
- `new_checked()` - Validates inferred constraints

## Deliverables

### Modified Files (clifford-codegen)

#### `src/codegen/types.rs`

Remove `new()` generation, update remaining constructors:

```rust
fn generate_constructors(type_spec: &TypeSpec) -> TokenStream {
    let fields = &type_spec.fields;
    let field_names: Vec<_> = fields.iter().map(|f| &f.name).collect();
    let field_types: Vec<_> = fields.iter().map(|_| quote!(T)).collect();

    let new_unchecked = quote! {
        /// Creates a new instance without any validation.
        ///
        /// Use this for:
        /// - Product outputs (algebraically guaranteed valid)
        /// - Factory methods that construct valid instances
        /// - Performance-critical code where you've pre-validated
        #[inline]
        pub fn new_unchecked(#(#field_names: #field_types),*) -> Self {
            Self { #(#field_names),* }
        }
    };

    let new_checked = generate_new_checked(type_spec);

    quote! {
        #new_unchecked
        #new_checked
    }
}

fn generate_new_checked(type_spec: &TypeSpec) -> TokenStream {
    let fields = &type_spec.fields;
    let field_names: Vec<_> = fields.iter().map(|f| &f.name).collect();
    let field_types: Vec<_> = fields.iter().map(|_| quote!(T)).collect();
    let type_name = &type_spec.name;

    if let Some(constraint) = &type_spec.inferred_constraint {
        let residual_expr = parse_constraint_to_residual(&constraint.expression);
        let constraint_name = constraint.name;

        quote! {
            /// Creates a new instance with constraint validation.
            ///
            /// Validates that the inferred constraint is satisfied within
            /// the given tolerance. Returns an error if violated.
            ///
            /// # Constraint
            #[doc = #constraint.expression]
            pub fn new_checked(
                #(#field_names: #field_types,)*
                tolerance: T,
            ) -> Result<Self, crate::ConstraintError> {
                let residual = #residual_expr;
                if residual.abs() > tolerance {
                    return Err(crate::ConstraintError::new(
                        #type_name,
                        #constraint_name,
                        residual.to_f64().unwrap_or(0.0),
                    ));
                }
                Ok(Self { #(#field_names),* })
            }
        }
    } else {
        // No constraint - new_checked is trivially Ok
        quote! {
            /// Creates a new instance with constraint validation.
            ///
            /// This type has no algebraic constraint, so this always succeeds.
            #[inline]
            pub fn new_checked(
                #(#field_names: #field_types,)*
                _tolerance: T,
            ) -> Result<Self, crate::ConstraintError> {
                Ok(Self { #(#field_names),* })
            }
        }
    }
}
```

#### Remove from `src/codegen/types.rs`

```rust
// REMOVE: generate_new() function
// REMOVE: generate_constraint_solving_new() function
// REMOVE: SolutionType handling for new()
```

### Generated Code Example

**Before:**
```rust
impl<T: Float> Motor<T> {
    // Confusing: solves for e0123
    pub fn new(s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T) -> Self {
        let e0123 = -(e23 * e01 + e31 * e02 + e12 * e03) / s;
        Self { s, e23, e31, e12, e01, e02, e03, e0123 }
    }

    pub fn new_checked(..., tolerance: T) -> Result<Self, ConstraintError> { ... }
    pub fn new_unchecked(...) -> Self { ... }
}
```

**After:**
```rust
impl<T: Float> Motor<T> {
    /// Creates a Motor without any validation.
    #[inline]
    pub fn new_unchecked(
        s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T, e0123: T
    ) -> Self {
        Self { s, e23, e31, e12, e01, e02, e03, e0123 }
    }

    /// Creates a Motor with Study condition validation.
    ///
    /// # Constraint
    /// `s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0`
    pub fn new_checked(
        s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T, e0123: T,
        tolerance: T,
    ) -> Result<Self, ConstraintError> {
        let residual = s * e0123 + e23 * e01 + e31 * e02 + e12 * e03;
        if residual.abs() > tolerance {
            return Err(ConstraintError::new("Motor", "study", residual.to_f64().unwrap_or(0.0)));
        }
        Ok(Self { s, e23, e31, e12, e01, e02, e03, e0123 })
    }
}
```

## Migration

Users who relied on `new()` solving constraints should:
1. Use factory methods that guarantee validity (e.g., `Motor::from_rotation_z()`)
2. Construct with `new_unchecked()` and normalize with wrapper types
3. Validate external data with `new_checked()`

## Testing

```rust
#[test]
fn new_unchecked_accepts_any_values() {
    // Even invalid values are accepted
    let m = Motor::new_unchecked(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 999.0);
    assert_eq!(m.e0123(), 999.0);
}

#[test]
fn new_checked_validates_constraint() {
    // Valid motor (Study condition satisfied)
    let valid = Motor::new_checked(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-10);
    assert!(valid.is_ok());

    // Invalid motor (Study condition violated)
    let invalid = Motor::new_checked(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1e-10);
    assert!(invalid.is_err());
}
```

## Success Criteria

1. `new()` constructor removed from generation
2. `new_unchecked()` takes all fields, no validation
3. `new_checked()` validates inferred constraint
4. Generated docs explain each constructor's purpose
5. All existing tests updated or passing
