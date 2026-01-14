# PRD-25: Clean Up TOML Codegen Options

**Status**: Done
**Goal**: Remove unused options from algebra TOML specs and simplify defaults

## Background

The algebra TOML specification files contain several options that control code generation:

```toml
[options]
generate_serde = true
generate_arbitrary = true
generate_nalgebra = true
```

However, analysis of the codegen reveals:
1. `generate_nalgebra` is never parsed - it only exists in template generation but not in the actual parser
2. `generate_tests` is parsed but never used anywhere in code generation (dead code)
3. Default values don't match intended behavior: `generate_serde` defaults to `false` but should default to `true`

## Problem Statement

### Issue 1: Dead Code - `generate_tests`

The `generate_tests` option is:
- Defined in `RawOptions` struct
- Parsed from TOML
- Converted to `GenerationOptions`
- **Never actually used** in any code generation logic

### Issue 2: Orphaned Option - `generate_nalgebra`

The `generate_nalgebra` option:
- Appears in TOML files
- Appears in template generation (`template.rs`)
- **Is not defined in** `RawOptions` or `GenerationOptions`
- Cannot actually control anything since it's not parsed

### Issue 3: Incorrect Defaults

Current defaults:
- `generate_serde`: `false` (should be `true` - we always want serde)
- `generate_arbitrary`: `true` (correct)

## Solution

### Phase 1: Remove Dead Code

Remove `generate_tests` from:
- `RawOptions` struct in `raw.rs`
- `GenerationOptions` struct in `ir.rs`
- Parser mapping in `parser.rs`

### Phase 2: Remove Orphaned Option

Remove `generate_nalgebra` from:
- Template generation in `template.rs`
- All TOML files

### Phase 3: Fix Defaults

Update `raw.rs` to default `generate_serde` to `true`:
```rust
#[serde(default = "default_true")]
pub generate_serde: bool,
```

### Phase 4: Clean Up TOML Files

Since both `generate_serde` and `generate_arbitrary` will default to `true`, and `generate_nalgebra` is removed, the `[options]` section becomes unnecessary and can be removed from all TOML files.

## Implementation Plan

### Step 1: Update `crates/clifford-codegen/src/spec/raw.rs`
- Remove `generate_tests` field
- Add `#[serde(default = "default_true")]` to `generate_serde`

### Step 2: Update `crates/clifford-codegen/src/spec/ir.rs`
- Remove `generate_tests` field from `GenerationOptions`

### Step 3: Update `crates/clifford-codegen/src/spec/parser.rs`
- Remove `generate_tests` mapping

### Step 4: Update `crates/clifford-codegen/src/discovery/template.rs`
- Remove `generate_nalgebra` from template generation

### Step 5: Update TOML Files
Remove `[options]` section from:
- `algebras/euclidean2.toml`
- `algebras/euclidean3.toml`
- `algebras/projective2.toml`
- `algebras/projective3.toml`

## Testing

1. Build codegen: `cargo build -p clifford-codegen`
2. Run all tests: `cargo nextest run`
3. Verify TOML parsing works without explicit options

## Success Criteria

1. No `generate_tests` references remain in codebase
2. No `generate_nalgebra` references remain in codebase
3. TOML files have no `[options]` section
4. All tests pass
5. Generated code unchanged (since options were already `true`)
