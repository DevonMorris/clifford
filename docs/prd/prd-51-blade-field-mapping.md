# PRD-51: Per-Type Blade-to-Field Mapping

**Status**: Draft
**Goal**: Enable explicit, user-defined blade-to-field mappings in algebra TOML files, supporting alternative blade conventions (e.g., RGA's e20 vs lexicographic e02)

## Background

Different geometric algebra references use different blade ordering conventions:

1. **Lexicographic (ascending)**: Always write indices in ascending order
   - e01, e02, e12, e13, e23, etc.

2. **RGA/Cyclic**: Use conventions that make geometric formulas cleaner
   - e41, e42, e43 (projective basis first)
   - e23, e31, e12 (cyclic, matching right-hand rule)

These conventions differ by **sign**: `e20 = -e02`, `e31 = -e13`, etc.

### Current Behavior

The codegen assumes fields are listed in lexicographic blade order:

```toml
[types.Line]
grades = [2]
fields = ["dist", "normal_x", "normal_y"]
# Implicitly maps: dist -> e12, normal_x -> e01, normal_y -> e02
```

This implicit mapping causes problems when the geometric interpretation doesn't match lexicographic order. For example, in 2D PGA:

- The wedge of two points produces bivector coefficients (e12, e01, e02)
- The line incidence formula is: `e01*y - e02*x + e12 = 0`
- If we want `nx*x + ny*y + d = 0`, we need: `nx = -e02`, `ny = e01`, `d = e12`

The current implicit mapping gives wrong results because it assumes `normal_x = e01` when geometrically `normal_x = -e02`.

## Problem Statement

1. **Sign errors in geometric formulas**: Field values don't match expected geometric semantics
2. **No way to use alternative conventions**: Can't follow RGA conventions even when they make code cleaner
3. **Silent failures**: Wrong results with no indication of the underlying convention mismatch

## Solution

Add support for **explicit blade-to-field mappings** in the TOML specification.

### New TOML Syntax

**All types must use explicit `field_map`** - no more implicit lexicographic ordering:

```toml
[types.Line]
grades = [2]
field_map = [
  { name = "nx", blade = "e20" },   # e20 = -e02, so this stores -e02
  { name = "ny", blade = "e01" },   # e01 directly
  { name = "d", blade = "e12" }     # e12 directly
]

[types.Point]
grades = [1]
field_map = [
  { name = "x", blade = "e1" },
  { name = "y", blade = "e2" },
  { name = "w", blade = "e0" }
]
```

The old `fields = [...]` syntax is **removed**. This eliminates ambiguity and makes every blade mapping explicit.

### Blade Name Parsing

The codegen will parse blade names and determine:
1. **Which basis vectors** are involved (e.g., "e20" -> {e2, e0})
2. **The sign relative to canonical form** (e.g., e20 = -e02, so sign = -1)

```rust
struct BladeSpec {
    name: String,           // User-provided name, e.g., "e20"
    canonical_index: usize, // Bitmask index in canonical order
    sign: i8,               // +1 or -1 relative to canonical
}

fn parse_blade_name(name: &str, basis_names: &[String]) -> Result<BladeSpec, Error> {
    // Parse "e20" -> indices [2, 0]
    // Compute canonical order [0, 2] with sign from permutation parity
    // Return canonical bitmask index with sign correction
}
```

### Sign Correction in Products

When generating products, apply sign corrections:

```rust
// In product generation
let canonical_coeff = compute_canonical_product(lhs_blade, rhs_blade);
let output_coeff = canonical_coeff * output_field.sign;
```

## Detailed Design

### TOML Schema Changes

```toml
# Option 1: Simple fields (lexicographic order, positive sign)
fields = ["x", "y", "z"]

# Option 2: Explicit mapping with blade names
field_map = [
  { name = "field_name", blade = "e123" },
  ...
]

# Option 3: Explicit mapping with sign override (if needed)
field_map = [
  { name = "field_name", blade = "e12", sign = -1 },
  ...
]
```

### Parser Changes

In `spec/parser.rs`:

```rust
#[derive(Debug, Clone, Deserialize)]
#[serde(untagged)]
enum FieldSpec {
    Simple(Vec<String>),
    Mapped(Vec<FieldMapping>),
}

#[derive(Debug, Clone, Deserialize)]
struct FieldMapping {
    name: String,
    blade: String,
    #[serde(default = "default_sign")]
    sign: i8,
}

fn default_sign() -> i8 { 1 }
```

### Blade Name Parser

```rust
/// Parses a blade name like "e20", "e123", "e31" into a canonical form.
fn parse_blade_name(
    blade_str: &str,
    basis_order: &[String],
) -> Result<(usize, i8), ParseError> {
    // 1. Extract indices from blade name (e.g., "e20" -> [2, 0])
    // 2. Map to canonical indices based on basis_order
    // 3. Sort indices to canonical order, counting swaps for sign
    // 4. Compute bitmask from canonical indices
    // 5. Return (bitmask, sign)
}

fn permutation_sign(indices: &[usize]) -> i8 {
    // Count inversions (pairs where i < j but indices[i] > indices[j])
    let mut inversions = 0;
    for i in 0..indices.len() {
        for j in (i + 1)..indices.len() {
            if indices[i] > indices[j] {
                inversions += 1;
            }
        }
    }
    if inversions % 2 == 0 { 1 } else { -1 }
}
```

### Updated Field Structure

```rust
#[derive(Debug, Clone)]
pub struct Field {
    pub name: String,
    pub blade_index: usize,  // Canonical bitmask
    pub grade: usize,
    pub sign: i8,            // NEW: +1 or -1 relative to canonical
}
```

### Product Generation Changes

In `generation/products.rs`, when accumulating terms:

```rust
fn generate_product_term(
    lhs_field: &Field,
    rhs_field: &Field,
    output_field: &Field,
    canonical_sign: i8,
) -> Term {
    // Total sign = canonical_sign * lhs.sign * rhs.sign * output.sign
    // (output.sign because we're storing in a field with that convention)
    let total_sign = canonical_sign * lhs_field.sign * rhs_field.sign * output_field.sign;

    Term {
        coefficient: total_sign,
        lhs: lhs_field.name.clone(),
        rhs: rhs_field.name.clone(),
    }
}
```

## Example: Fixing 2D PGA Line

### Before (broken)

```toml
[types.Line]
grades = [2]
fields = ["dist", "normal_x", "normal_y"]
# Implicit: dist -> e12, normal_x -> e01, normal_y -> e02
# Line equation interpreted as: normal_x*x + normal_y*y + dist = 0
# But PGA gives: e01*y - e02*x + e12 = 0
# MISMATCH!
```

### After (fixed)

```toml
[types.Line]
grades = [2]
description = "2D line: nx*x + ny*y + d = 0"
field_map = [
  { name = "nx", blade = "e20" },   # e20 = -e02, stores coefficient of x
  { name = "ny", blade = "e01" },   # e01 directly, stores coefficient of y
  { name = "d", blade = "e12" }     # e12 directly, stores constant term
]
```

Now when Point ∧ Point produces a Line:
- The e01 coefficient goes to `ny` (correct: coeff of y)
- The e02 coefficient is negated and goes to `nx` (correct: coeff of x)
- The e12 coefficient goes to `d` (correct: constant term)

## Implementation Plan

### Phase 1: Core Parser Changes
- [ ] Add `FieldMapping` struct and `FieldSpec` enum
- [ ] Implement `parse_blade_name()` function
- [ ] Add `sign` field to `Field` struct
- [ ] Update `build_type()` to handle both field syntaxes

### Phase 2: Product Generation Changes
- [ ] Update `generate_product_term()` to apply sign corrections
- [ ] Update all product generators (wedge, geometric, etc.)
- [ ] Add tests for sign-corrected products

### Phase 3: Type Generation Changes
- [ ] Update struct generation to use field signs correctly
- [ ] Update accessor method generation
- [ ] Update `From` conversions for Multivector

### Phase 4: Migration
- [ ] Migrate ALL algebra TOML files to use `field_map`
- [ ] Update `projective2.toml` Line type to use correct blade conventions
- [ ] Regenerate all algebras

### Phase 5: Testing
- [ ] Add unit tests for blade name parsing
- [ ] Add integration tests for sign-corrected products
- [ ] Add property-based tests comparing to Symbolica

## Backwards Compatibility

**This is a breaking change.** The `fields = [...]` syntax is removed. All existing algebra TOML files must be migrated to use `field_map`.

Migration is straightforward: convert each field to an explicit blade mapping with the canonical blade name (ascending indices, sign = +1).

## Testing Strategy

### Unit Tests

```rust
#[test]
fn parse_canonical_blade() {
    let (idx, sign) = parse_blade_name("e12", &["e0", "e1", "e2"]).unwrap();
    assert_eq!(idx, 0b110);  // bits 1 and 2
    assert_eq!(sign, 1);     // already canonical
}

#[test]
fn parse_swapped_blade() {
    let (idx, sign) = parse_blade_name("e20", &["e0", "e1", "e2"]).unwrap();
    assert_eq!(idx, 0b101);  // bits 0 and 2 (e0 and e2)
    assert_eq!(sign, -1);    // one swap needed: e20 = -e02
}

#[test]
fn parse_triple_blade() {
    let (idx, sign) = parse_blade_name("e312", &["e1", "e2", "e3"]).unwrap();
    assert_eq!(idx, 0b111);  // all three bits
    // e312 -> e123 requires: 312 -> 132 -> 123 = 2 swaps = even = +1
    assert_eq!(sign, 1);
}
```

### Integration Tests

```rust
#[test]
fn line_from_point_wedge_has_correct_semantics() {
    let p1 = Point::from_cartesian(0.0, 0.0);
    let p2 = Point::from_cartesian(1.0, 0.0);

    let line = p1.wedge(&p2);  // Should be y = 0 (x-axis)

    // Line equation: nx*x + ny*y + d = 0
    // For y = 0: nx = 0, ny = 1, d = 0
    assert!(line.nx().abs() < 1e-10);
    assert!((line.ny() - 1.0).abs() < 1e-10 || (line.ny() + 1.0).abs() < 1e-10);
    assert!(line.d().abs() < 1e-10);
}
```

## Success Criteria

1. [ ] `field_map` syntax parses correctly in TOML
2. [ ] Blade names with non-canonical order produce correct sign corrections
3. [ ] Products generate correct signs when output type uses `field_map`
4. [ ] 2D PGA Line wedge product gives geometrically correct results
5. [ ] All existing tests pass (backwards compatibility)
6. [ ] Demo visualization shows correct line-through-points

## References

- [RGA Wiki - Line](https://rigidgeometricalgebra.org/wiki/index.php?title=Line) (uses e41, e42, e43, e23, e31, e12)
- [CGA Wiki - Line](https://conformalgeometricalgebra.org/wiki/index.php?title=Line)
- PRD-16.1: Blade Ordering Audit (related but different focus)
- PRD-29: Semantic Field Names (complements this PRD)
