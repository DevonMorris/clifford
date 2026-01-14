# PRD-27: Auto-Identify Constraints and Generate Checked Constructors

**Status**: Partial (Part 1 complete, Parts 2-3 deferred)
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
2. ⏳ `Motor::new_checked()` exists and validates Study condition (deferred)
3. ⏳ `Motor::from_components()` exists and computes e0123 (deferred)
4. ✅ All existing tests pass

## Completed Work

### Part 1: Dead Code Removal (Complete)

The following dead code was removed:

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

## Deferred Work

### Part 2: new_checked() Generation (Future)

Generating `new_checked()` constructors requires:
1. Converting symbolic constraint expressions to Rust code
2. Implementing expression-to-code conversion in codegen
3. Handling different constraint forms (Study, Plücker, etc.)

This is tracked for future implementation.

### Part 3: from_components() Generation (Future)

Generating `from_components()` constructors requires:
1. Solving constraints symbolically for the dependent field
2. Generating code to compute the dependent field
3. Handling edge cases (division by zero, etc.)

This is tracked for future implementation.
