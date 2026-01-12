# PRD-17.3: Field Documentation Fix

**Status**: Draft
**Parent**: [PRD-17: Codegen Product Completeness and Documentation](prd-17-codegen-products.md)
**Goal**: Fix generated field documentation to use TOML field names instead of canonical blade names

## Problem Statement

Generated type documentation uses canonical blade names from bitmask indices instead of the field names defined in the TOML specification:

```rust
// Current (WRONG)
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s`."]
    s: T,
    #[doc = "Coefficient of `e1e2`."]  // WRONG! Field is e23
    e23: T,
    #[doc = "Coefficient of `e1e3`."]  // WRONG! Field is e31
    e31: T,
    #[doc = "Coefficient of `e2e3`."]  // WRONG! Field is e12
    e12: T,
    #[doc = "Coefficient of `e1e4`."]  // WRONG! Should be e0e1 or just e01
    e01: T,
    // ...
}
```

This happens because the codegen computes blade names from bitmask indices:
- `e12` has bitmask `0b0011` → interpreted as `e1 ∧ e2` → doc says `e1e2`
- But the TOML defines field name as `e23` (meaning `e2 ∧ e3`)

The mismatch makes code review and debugging extremely difficult.

## Root Cause

In `types.rs` generation, the blade name is computed from the bitmask:

```rust
// Current (problematic) code
fn generate_field_doc(field: &FieldSpec, algebra: &Algebra) -> String {
    let blade = Blade::from_index(field.blade_index);
    let blade_name = blade.to_string(algebra);  // Computes from bitmask
    format!("Coefficient of `{}`.", blade_name)
}
```

The `blade.to_string()` computes the canonical name from indices, ignoring the user-specified field name.

## Solution

### Use Field Name Directly

Change documentation generation to use the field name from the TOML spec:

```rust
// Fixed code
fn generate_field_doc(field: &FieldSpec) -> String {
    format!("Coefficient of `{}`.", field.name)
}
```

Result:
```rust
// Correct
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s`."]
    s: T,
    #[doc = "Coefficient of `e23`."]  // Correct!
    e23: T,
    #[doc = "Coefficient of `e31`."]  // Correct!
    e31: T,
    #[doc = "Coefficient of `e12`."]  // Correct!
    e12: T,
    #[doc = "Coefficient of `e01`."]  // Correct!
    e01: T,
    // ...
}
```

### Optional: Add Basis Expansion

For educational value, optionally include the basis blade expansion:

```rust
fn generate_field_doc(field: &FieldSpec, include_basis: bool) -> String {
    if include_basis && !field.is_scalar() {
        let basis = compute_basis_expansion(&field.name);
        format!("Coefficient of `{}` (basis: {}).", field.name, basis)
    } else {
        format!("Coefficient of `{}`.", field.name)
    }
}

fn compute_basis_expansion(name: &str) -> String {
    // e23 → "e₂ ∧ e₃"
    // e012 → "e₀ ∧ e₁ ∧ e₂"
    // ...
}
```

Result:
```rust
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s`."]
    s: T,
    #[doc = "Coefficient of `e23` (basis: e₂ ∧ e₃)."]
    e23: T,
    #[doc = "Coefficient of `e31` (basis: e₃ ∧ e₁)."]
    e31: T,
    // ...
}
```

## Implementation

### Update TypeGenerator

```rust
// In crates/clifford-codegen/src/codegen/types.rs

impl<'a> TypeGenerator<'a> {
    fn generate_field(&self, field: &FieldSpec) -> TokenStream {
        let field_name = format_ident!("{}", field.name);
        let field_type = quote! { T };

        // Use field name directly, not computed blade name
        let doc = format!("Coefficient of `{}`.", field.name);

        quote! {
            #[doc = #doc]
            #field_name: #field_type,
        }
    }
}
```

### Update Accessor Documentation

Also fix accessor method documentation:

```rust
// Current (may have same issue)
impl<T: Float> Motor<T> {
    #[doc = "Returns the `e1e2` coefficient."]  // WRONG
    pub fn e23(&self) -> T { self.e23 }
}

// Fixed
impl<T: Float> Motor<T> {
    #[doc = "Returns the `e23` coefficient."]  // Correct
    pub fn e23(&self) -> T { self.e23 }
}
```

```rust
fn generate_accessor(&self, field: &FieldSpec) -> TokenStream {
    let field_name = format_ident!("{}", field.name);
    let doc = format!("Returns the `{}` coefficient.", field.name);

    quote! {
        #[doc = #doc]
        #[inline]
        pub fn #field_name(&self) -> T {
            self.#field_name
        }
    }
}
```

### Configuration Option

Add TOML option for verbose documentation:

```toml
[codegen]
verbose_field_docs = true  # Include basis expansion in docs
```

```rust
// In parser
struct CodegenConfig {
    verbose_field_docs: bool,
}

// In generator
fn generate_field_doc(&self, field: &FieldSpec) -> String {
    if self.config.verbose_field_docs {
        let basis = self.compute_basis_expansion(field);
        format!("Coefficient of `{}` (basis: {}).", field.name, basis)
    } else {
        format!("Coefficient of `{}`.", field.name)
    }
}
```

## Basis Expansion Algorithm

For the optional verbose mode:

```rust
fn compute_basis_expansion(field_name: &str) -> String {
    // Handle scalar
    if field_name == "s" || field_name == "scalar" {
        return "1".to_string();
    }

    // Parse field name like "e23", "e012", "e0123"
    let indices: Vec<char> = field_name
        .strip_prefix('e')
        .unwrap_or(field_name)
        .chars()
        .collect();

    // Build basis string with subscripts and wedge symbols
    indices.iter()
        .map(|&c| format!("e{}", subscript(c)))
        .collect::<Vec<_>>()
        .join(" ∧ ")
}

fn subscript(c: char) -> &'static str {
    match c {
        '0' => "₀",
        '1' => "₁",
        '2' => "₂",
        '3' => "₃",
        '4' => "₄",
        _ => "?",
    }
}
```

Examples:
- `e23` → `"e₂ ∧ e₃"`
- `e012` → `"e₀ ∧ e₁ ∧ e₂"`
- `e0123` → `"e₀ ∧ e₁ ∧ e₂ ∧ e₃"`
- `s` → `"1"`

## Testing

### Verify Documentation Correctness

```rust
#[test]
fn field_docs_match_field_names() {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();
    let generator = TypeGenerator::new(&spec);

    for ty in &spec.types {
        let generated = generator.generate_type(ty);
        let code = generated.to_string();

        for field in &ty.fields {
            // Check that doc contains field name, not computed blade name
            let expected_doc = format!("Coefficient of `{}`", field.name);
            assert!(
                code.contains(&expected_doc),
                "Type {} field {} has incorrect doc",
                ty.name, field.name
            );
        }
    }
}
```

### Verify Accessor Documentation

```rust
#[test]
fn accessor_docs_match_field_names() {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();
    let generator = TypeGenerator::new(&spec);

    for ty in &spec.types {
        let generated = generator.generate_type(ty);
        let code = generated.to_string();

        for field in &ty.fields {
            let expected_doc = format!("Returns the `{}` coefficient", field.name);
            assert!(
                code.contains(&expected_doc),
                "Type {} accessor {} has incorrect doc",
                ty.name, field.name
            );
        }
    }
}
```

## Generated Output Example

### Before (Current)

```rust
/// Motor represents a rigid transformation in 3D PGA.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s`."]
    s: T,
    #[doc = "Coefficient of `e1e2`."]
    e23: T,
    #[doc = "Coefficient of `e1e3`."]
    e31: T,
    #[doc = "Coefficient of `e2e3`."]
    e12: T,
    #[doc = "Coefficient of `e1e4`."]
    e01: T,
    #[doc = "Coefficient of `e2e4`."]
    e02: T,
    #[doc = "Coefficient of `e3e4`."]
    e03: T,
    #[doc = "Coefficient of `e1e2e3e4`."]
    e0123: T,
}

impl<T: Float> Motor<T> {
    #[doc = "Returns the `e1e2` coefficient."]
    pub fn e23(&self) -> T { self.e23 }
    // ...
}
```

### After (Fixed)

```rust
/// Motor represents a rigid transformation in 3D PGA.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s`."]
    s: T,
    #[doc = "Coefficient of `e23`."]
    e23: T,
    #[doc = "Coefficient of `e31`."]
    e31: T,
    #[doc = "Coefficient of `e12`."]
    e12: T,
    #[doc = "Coefficient of `e01`."]
    e01: T,
    #[doc = "Coefficient of `e02`."]
    e02: T,
    #[doc = "Coefficient of `e03`."]
    e03: T,
    #[doc = "Coefficient of `e0123`."]
    e0123: T,
}

impl<T: Float> Motor<T> {
    #[doc = "Returns the `e23` coefficient."]
    pub fn e23(&self) -> T { self.e23 }
    // ...
}
```

### After (Verbose Mode)

```rust
pub struct Motor<T: Float> {
    #[doc = "Coefficient of `s` (basis: 1)."]
    s: T,
    #[doc = "Coefficient of `e23` (basis: e₂ ∧ e₃)."]
    e23: T,
    #[doc = "Coefficient of `e31` (basis: e₃ ∧ e₁)."]
    e31: T,
    #[doc = "Coefficient of `e12` (basis: e₁ ∧ e₂)."]
    e12: T,
    #[doc = "Coefficient of `e01` (basis: e₀ ∧ e₁)."]
    e01: T,
    // ...
}
```

## Deliverables

- [ ] Update `generate_field()` to use field name from spec
- [ ] Update `generate_accessor()` to use field name from spec
- [ ] Add optional `verbose_field_docs` config
- [ ] Implement `compute_basis_expansion()` for verbose mode
- [ ] Add tests verifying doc correctness
- [ ] Regenerate all existing algebras
- [ ] Verify docs render correctly with `cargo doc`

## Success Criteria

1. **Field docs correct**: All field docs show TOML field name
2. **Accessor docs correct**: All accessor docs show field name
3. **Cargo doc builds**: Documentation renders correctly
4. **No regressions**: Existing tests still pass

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/codegen/types.rs` | Update - Use field.name in docs |
| `crates/clifford-codegen/src/spec/raw.rs` | Update - Add verbose_field_docs option |
| `src/generated/*/types.rs` | Regenerate - Corrected documentation |
| `src/specialized/*/generated/types.rs` | Regenerate - Corrected documentation |

## Estimated Scope

- Field doc fix: ~10 lines
- Accessor doc fix: ~10 lines
- Verbose mode: ~50 lines
- Tests: ~50 lines
- Total: ~120 lines of code

This is a small but high-impact fix that significantly improves code readability and reviewability.
