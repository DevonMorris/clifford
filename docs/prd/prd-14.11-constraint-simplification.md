# PRD-14.11: Constraint Simplification

## Status: Draft

## Problem Statement

The current constraint system has unnecessary complexity:

1. **Bespoke constraint kinds**: `Unit`, `NonZero`, `Blade` require special-case code generation
2. **Wrapper types**: `UnitRotor`, `NonZeroVector` duplicate the base type with validation
3. **Mixed approaches**: Some constraints use `norm = "euclidean"`, others use `condition = "expr"`

This complexity is unnecessary because:
- All types already satisfy geometric constraints (PRD-14.8 ensures `u * ũ = scalar`)
- Field-based expressions can represent any constraint uniformly
- Wrapper types add API surface without adding value

## Solution

### 1. Remove Bespoke Constraint Kinds

Replace special constraint kinds with field expressions:

**Before:**
```toml
[types.Rotor.constraints.unit]
norm = "euclidean"
```

**After:**
```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy"]
constraint = "s * s + xy * xy = 1"  # Direct field expression
```

### 2. Remove Wrapper Types

Instead of generating `UnitRotor` wrapper around `Rotor`, constraints become:
- Documentation on the type
- Optional runtime validation methods on the base type
- Type aliases if needed for API clarity

**Before:**
```rust
pub struct Rotor<T> { s: T, xy: T }
pub struct UnitRotor<T>(Rotor<T>);  // wrapper

impl<T: Float> UnitRotor<T> {
    pub fn new_checked(inner: Rotor<T>) -> Self { ... }
}
```

**After:**
```rust
/// Rotor (even subalgebra element).
///
/// # Constraint
///
/// Unit rotors satisfy: `s² + xy² = 1`
pub struct Rotor<T> { s: T, xy: T }

impl<T: Float> Rotor<T> {
    /// Returns true if this rotor satisfies the unit constraint.
    pub fn is_unit(&self, tolerance: T) -> bool {
        (self.s * self.s + self.xy * self.xy - T::one()).abs() < tolerance
    }

    /// Normalizes to unit norm.
    pub fn normalize(&self) -> Option<Self> { ... }
}
```

### 3. Constraint Specification Format

Single field expression in TOML:

```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy"]
constraint = "s * s + xy * xy = 1"  # equality constraint
description = "2D rotor (even subalgebra element)"

[types.Vector]
grades = [1]
fields = ["x", "y"]
# No constraint - general vectors don't have one
```

Constraint expressions support:
- Field references: `s`, `xy`, `x`, `y`, etc.
- Arithmetic: `+`, `-`, `*`, `/`
- Equality: `=` (with tolerance in generated code)
- Comparisons: `>`, `<` for non-zero constraints

### 4. Generated Code

For types with constraints:

```rust
impl<T: Float> Rotor<T> {
    /// Checks if this element satisfies its constraint.
    ///
    /// Constraint: `s * s + xy * xy = 1`
    #[inline]
    pub fn satisfies_constraint(&self, tolerance: T) -> bool {
        (self.s * self.s + self.xy * self.xy - T::one()).abs() < tolerance
    }

    /// Asserts the constraint holds, panicking if not.
    #[inline]
    pub fn assert_constraint(&self, tolerance: T) {
        assert!(
            self.satisfies_constraint(tolerance),
            "Constraint violated: s * s + xy * xy = 1"
        );
    }
}
```

For types without constraints, no constraint methods are generated.

## Changes Required

### Files to Modify

1. **`src/spec/raw.rs`**
   - Remove `RawConstraint` struct
   - Add `constraint: Option<String>` to `RawTypeSpec`

2. **`src/spec/ir.rs`**
   - Remove `ConstraintSpec`, `ConstraintKind`
   - Add `constraint: Option<String>` to `TypeSpec`

3. **`src/spec/parser.rs`**
   - Remove constraint parsing logic
   - Parse simple constraint string

4. **`src/codegen/constraints.rs`**
   - Remove entirely (wrapper type generation)

5. **`src/codegen/types.rs`**
   - Generate `satisfies_constraint()` method when constraint present
   - Add constraint to type documentation

6. **`src/codegen/mod.rs`**
   - Remove `ConstraintGenerator`
   - Remove `constrained.rs` output

### Files to Remove

- `src/codegen/constraints.rs` - No longer needed

## Example: Before and After

### Before (euclidean2.toml)
```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy"]

[types.Rotor.constraints.unit]
norm = "euclidean"
```

### After (euclidean2.toml)
```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy"]
constraint = "s * s + xy * xy = 1"
```

## Migration

Existing TOML files with `[types.X.constraints.Y]` sections will need updating:
- `unit` with `norm = "euclidean"` → `constraint = "<field>² + ... = 1"`
- `nonzero` → `constraint = "<field>² + ... > 0"`
- Custom `condition = "expr"` → `constraint = "expr"`

## Success Criteria

1. All constraint logic removed from codegen
2. No wrapper types generated
3. Simple `constraint = "expr"` in TOML
4. `satisfies_constraint()` method on types with constraints
5. All existing tests pass
6. Generated code compiles and works correctly

## Non-Goals

- Runtime constraint enforcement (users call `satisfies_constraint()` as needed)
- Compile-time constraint checking (beyond documentation)
- Automatic constraint inference from algebra properties
