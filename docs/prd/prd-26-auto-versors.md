# PRD-26: Auto-Identify Versors and Simplify TOML

**Status**: Draft
**Goal**: Simplify algebra TOML specs by auto-identifying versors and their sandwich targets

## Background

The algebra TOML specification files require manual configuration for versors:

```toml
[types.Rotor]
grades = [0, 2]
versor = true  # Must manually mark as versor

[types.Motor]
grades = [0, 2, 4]
versor = true
sandwich = { targets = ["Point", "Line", "Plane", "Motor", "Flector"] }
```

This is unnecessary because:
1. **Versors can be identified automatically** - they have uniform grade parity (all even or all odd)
2. **Sandwich targets can be inferred** - already implemented in codegen
3. **Products are already auto-inferred** - `[products]` section is ignored

## Problem Statement

### Issue 1: Manual Versor Marking
Users must remember to add `versor = true` for rotor/motor/flector types. This is error-prone and adds complexity to the TOML format.

### Issue 2: Manual Sandwich Targets
The `sandwich = { targets = [...] }` configuration is redundant because the codegen already has `infer_sandwich_targets()` that computes valid targets automatically.

### Issue 3: Unnecessary Complexity
New users must learn about versor configuration when defining algebras, even though the information can be derived from the grade structure.

## Solution

### Auto-Identify Versors

A type is automatically identified as a versor if all its grades have the same parity:
- **Even versors**: grades `[0, 2]`, `[0, 2, 4]` (rotors, motors)
- **Odd versors**: grades `[1, 3]`, `[1, 3, 5]` (flectors, reflectors)

Use the existing `versor_parity()` function from `algebra/versor.rs`.

### Auto-Infer Sandwich Targets

Already implemented - `infer_sandwich_targets()` computes valid targets where `V * X * rev(V)` preserves grades. Simply use empty targets to trigger auto-inference.

## Implementation Plan

### Step 1: Remove versor fields from TOML parsing
- Remove `versor: bool` and `sandwich: Option<RawSandwichConfig>` from `RawTypeSpec`
- Remove `RawSandwichConfig` struct entirely

### Step 2: Add auto-versor identification
- Add `identify_versors()` function in parser
- Use `versor_parity()` to check grade parity
- Set `versor: Some(VersorSpec { is_unit: false, sandwich_targets: vec![] })`

### Step 3: Update TOML files
Remove from all algebras/*.toml:
- `versor = true`
- `sandwich = { targets = [...] }`

### Step 4: Update template generation
Remove versor-related output from discovery templates.

## Testing

1. Build codegen: `cargo build -p clifford-codegen`
2. Run all tests: `cargo nextest run`
3. Verify versors are correctly identified (Rotor, Motor, Flector)
4. Verify sandwich/antisandwich traits are generated

## Success Criteria

1. No `versor` or `sandwich` fields in TOML files
2. All existing versors (Rotor, Motor, Flector) still identified
3. All tests pass
4. Generated code unchanged
