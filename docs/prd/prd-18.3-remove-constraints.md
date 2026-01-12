# PRD-18.3: Remove User-Configurable Constraints

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Remove user-configurable constraints from TOML files and codegen

## Overview

Currently, users must specify constraints in TOML files:

```toml
[[types.Motor.constraints]]
name = "study"
expression = "s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0"
solve_for = "e0123"
```

This is redundant because constraints are fully determined by the algebra's signature and grade structure. This PRD removes user-configurable constraints.

## Deliverables

### Modified Files (clifford-codegen)

#### `src/spec/raw.rs`

Remove constraint-related fields from TOML parsing:

```rust
// REMOVE these fields from RawTypeSpec:
// - constraints: Vec<RawUserConstraint>
// - geometric_constraint: Option<String>
// - geometric_solve_for: Option<String>
// - antiproduct_constraint: Option<String>
// - antiproduct_solve_for: Option<String>

// REMOVE this struct entirely:
// pub struct RawUserConstraint { ... }
```

#### `src/spec/ir.rs`

Remove constraint-related types from IR:

```rust
// REMOVE:
// pub struct UserConstraint { ... }

// REMOVE from TypeSpec:
// pub constraints: Vec<UserConstraint>
```

#### `src/spec/parser.rs`

Remove constraint parsing logic:

```rust
// REMOVE:
// - parse_constraints() function
// - constraint validation logic
// - constraint dependency detection
```

### Modified Files (algebras)

#### `algebras/euclidean2.toml`

**Before:**
```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy"]
constraints = [
    { name = "unit", expression = "s*s + xy*xy = 1", solve_for = "s" }
]
```

**After:**
```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy"]
```

#### `algebras/euclidean3.toml`

**Before:**
```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy", "xz", "yz"]
constraints = [
    { name = "unit", expression = "s*s + xy*xy + xz*xz + yz*yz = 1", solve_for = "s" }
]
```

**After:**
```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy", "xz", "yz"]
```

#### `algebras/projective3.toml`

**Before:**
```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03", "e0123"]
versor = true

[[types.Motor.constraints]]
name = "study"
expression = "s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0"
solve_for = "e0123"

[types.Line]
grades = [2]
fields = ["e01", "e02", "e03", "e23", "e31", "e12"]

[[types.Line.constraints]]
name = "plucker"
expression = "e01*e23 + e02*e31 + e03*e12 = 0"
solve_for = "e03"
```

**After:**
```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03", "e0123"]
versor = true

[types.Line]
grades = [2]
fields = ["e01", "e02", "e03", "e23", "e31", "e12"]
```

## Testing

After removal, verify:

```bash
# Should still parse without errors
cargo run --package clifford-codegen -- generate algebras/euclidean3.toml --force

# Existing tests should pass (constraint inference handles validation)
cargo test
```

## Success Criteria

1. No `constraints` field in TOML parsing
2. No `UserConstraint` struct in IR
3. All algebra TOML files updated
4. Codegen still works (constraints will be inferred in PRD-18.4)
5. No breaking changes to existing tests
