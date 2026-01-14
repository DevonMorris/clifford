# PRD-27: Auto-Identify Constraints and Generate Checked Constructors

**Status**: Complete
**Goal**: Remove vestigial constraint infrastructure and generate constraint-checking constructors

## Background

The algebra TOML specification files historically required manual constraint configuration:

```toml
[types.Motor]
grades = [0, 2, 4]
constraint = "s*s + ... = 1"  # Was required manually
solve_for = "e0123"           # Field to compute
```

This has been obviated by:
1. **Automatic constraint derivation** via `ConstraintDeriver::derive_geometric_constraint()`
2. **Wrapper types** (`Unit<T>`, `Unitized<T>`) that enforce normalization
3. **Constraint solver** that computes dependent fields from independent ones

## Problem Statement

### Issue 1: Dead Code
The `UserConstraint` struct and related code in `ir.rs` is never populated:
- Parser sets `constraints = Vec::new()` (line 283)
- `SignConvention` enum is unused
- `TypeSpec` constraint methods are dead code

### Issue 2: Missing Validation Constructors
Types with geometric constraints (Motor, Line, Flector) have:
- `new()` - accepts any values, no validation
- `new_unchecked()` - just an alias for `new()`
- **No `new_checked()`** - documented but doesn't exist

### Issue 3: No Ergonomic Reduced-DOF Constructors
Motor has 8 fields but only 7 degrees of freedom (Study condition).
Users must compute `e0123` manually instead of having a constructor that does it.

## Solution

### Part 1: Remove Dead Code
Remove from `ir.rs`:
- `UserConstraint` struct
- `SignConvention` enum
- `TypeSpec::constraints` field
- Related methods (`has_constraint`, `solve_for_fields`, etc.)

### Part 2: Generate new_checked() Constructors
For types with derived constraints, generate:
```rust
pub fn new_checked(..., tolerance: T) -> Result<Self, ConstraintError> {
    let expected = /* computed from constraint */;
    if (actual - expected).abs() > tolerance {
        return Err(ConstraintError::new("Motor", "Study condition"));
    }
    Ok(Self::new_unchecked(...))
}
```

### Part 3: Generate from_components() Constructors
For ergonomics, generate reduced-DOF constructors:
```rust
pub fn from_components(/* N-1 free params */) -> Option<Self> {
    // Compute constrained field from free params
    // Return None if computation would be unstable
}
```

## Implementation Plan

### Step 1: Remove dead code from IR
- Remove `UserConstraint`, `SignConvention` from `ir.rs`
- Remove `constraints` field from `TypeSpec`
- Update `mod.rs` exports

### Step 2: Add constraint-checking constructor generation
- In `codegen/types.rs`, detect types with derived constraints
- Generate `new_checked()` that validates the constraint
- Use `ConstraintDeriver` to get the constraint expression

### Step 3: Add reduced-DOF constructor generation
- Generate `from_components()` for constrained types
- Use `ConstraintSolver` to compute dependent field

## Affected Types

| Type | Algebra | Constraint | Free DOF |
|------|---------|------------|----------|
| Motor | PGA 3D | Study condition on e0123 | 7 |
| Line | PGA 3D | Plücker condition on e12 | 5 |
| Flector | PGA 3D | Constraint on e123 | 7 |

## Testing

1. Build codegen: `cargo build -p clifford-codegen`
2. Run tests: `cargo nextest run`
3. Verify `Motor::new_checked()` rejects invalid values
4. Verify `Motor::from_components()` computes correct e0123

## Success Criteria

1. ✅ No `UserConstraint` or `SignConvention` in codebase
2. ✅ `Motor::new_checked()` exists and validates Study condition
3. ✅ `Motor::from_components()` exists and computes e0123
4. ✅ All existing tests pass

## Completed Work

### Part 1: Dead Code Removal

**`ir.rs`**:
- Removed `UserConstraint` struct
- Removed `SignConvention` enum
- Removed `constraints: Vec<UserConstraint>` field from `TypeSpec`
- Removed methods: `has_constraint()`, `solve_for_fields()`, `has_domain_restrictions()`

**`mod.rs`**:
- Removed `SignConvention` and `UserConstraint` from exports

**`parser.rs`**:
- Removed `constraints` field from `TypeSpec` construction

**`constraint_simplify.rs`**:
- Updated to return empty simplifier (constraints are auto-derived)
- Removed dead helper functions: `prefix_expression`, `replace_identifier`, `parse_atom`

**Codegen files** (`types.rs`, `traits.rs`, `projections.rs`, `unary.rs`, `conversions.rs`):
- Simplified constructor selection (removed `solve_for_fields()` checks)
- All files now use `new()` directly since `new` and `new_unchecked` are identical

### Part 2: new_checked() Generation

Added `generate_constraint_constructors()` method to `TypeGenerator` in `types.rs`:
- Uses `ConstraintDeriver` to detect types with geometric constraints
- Uses `ConstraintSolver` to solve for the constrained field
- Generates `new_checked(fields..., tolerance) -> Result<Self, &str>`
- Validates that the constrained field matches the expected value within tolerance

Example generated code for Motor:
```rust
pub fn new_checked(
    s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T, e0123: T,
    tolerance: T,
) -> Result<Self, &'static str> {
    let expected = (e01 * e23 - e02 * e31 + e03 * e12) / (s);
    let actual = e0123;
    if (actual - expected).abs() > tolerance {
        return Err("Motor constraint");
    }
    Ok(Self::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123))
}
```

### Part 3: from_components() Generation

Added `generate_from_components()` method to `TypeGenerator`:
- Takes N-1 free parameters (all except the constrained field)
- Computes the constrained field from the solution
- Returns `Option<Self>` (None if divisor is zero)

Example generated code for Motor:
```rust
pub fn from_components(
    s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T,
) -> Option<Self> {
    if (s).abs() < T::epsilon() {
        return None;
    }
    Some(Self::new_unchecked(
        s, e23, e31, e12, e01, e02, e03,
        (e01 * e23 - e02 * e31 + e03 * e12) / (s),
    ))
}
```

### Part 4: Remove `new()` for Constrained Types

Constrained types should only have `new_unchecked()`, `new_checked()`, and `from_components()`:

**types.rs**:
- Modified `generate_constructor()` to check if type has constraint
- Constrained types: only generate `new_unchecked()`
- Unconstrained types: generate both `new()` and `new_unchecked()` (alias)

**All codegen files** (`types.rs`, `traits.rs`, `projections.rs`, `unary.rs`, `conversions.rs`):
- Updated all generated code to use `new_unchecked()` instead of `new()`
- This ensures generated operations work for both constrained and unconstrained types

### Bug Fixes During Implementation

**Parser auto-versor identification** (`parser.rs`):
- Fixed: Single-grade types (Point, Line) were incorrectly marked as versors
- Solution: Only mark types with multiple grades as versors

**Sandwich target inference** (`traits.rs`):
- Fixed: `infer_sandwich_targets()` was doing worst-case grade analysis
- Solution: Versors preserve grades, so all types are valid targets
