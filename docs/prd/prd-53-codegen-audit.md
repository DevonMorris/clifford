# PRD-53: Code Generation Pipeline Audit and Cleanup

**Status**: Planning
**Goal**: Comprehensive audit and cleanup of the clifford-codegen pipeline to eliminate technical debt, standardize patterns, and create a formal TOML specification schema

## Executive Summary

This PRD documents a thorough audit of the `clifford-codegen` crate, identifying deprecated features, inconsistencies, technical debt, and areas requiring cleanup. It also proposes creating a formal TOML specification schema to replace the current example-based documentation.

---

## Table of Contents

1. [Audit Findings Summary](#audit-findings-summary)
2. [Critical Issues](#critical-issues)
3. [Technical Debt Inventory](#technical-debt-inventory)
4. [TOML Specification Schema](#toml-specification-schema)
5. [Implementation Plan](#implementation-plan)
6. [Effort Estimates](#effort-estimates)

---

## Audit Findings Summary

### Overview Statistics

| Category | Count | Severity Distribution |
|----------|-------|----------------------|
| TODO comments in code | 10 | Low |
| `#[allow(...)]` suppressions | 31 | Medium (violates CLAUDE.md) |
| Unused/dead code | 6 functions | Medium |
| Code duplication | 12+ locations | High |
| TOML spec inconsistencies | 15 files | High |
| Missing required sections | 4 algebras | High |
| Deprecated/legacy patterns | 3 modules | Medium |

---

## Critical Issues

### Issue 1: Missing `[norm]` Section in Core Algebras

**Severity**: HIGH
**Files Affected**:
- `algebras/euclidean2.toml`
- `algebras/euclidean3.toml`
- `algebras/projective2.toml`
- `algebras/projective3.toml`

**Problem**: These four core algebras lack the `[norm]` section entirely. The parser defaults to `"reverse"` which happens to be correct for these algebras, but the missing section:
1. Makes the intent unclear
2. Creates inconsistency with other algebras
3. Could cause silent bugs if defaults change

**Fix**: Add explicit `[norm]` sections to all four files.

### Issue 2: Extensive `#[allow(...)]` Suppressions

**Severity**: HIGH (violates CLAUDE.md principle: "No warning suppression - fix root causes")
**Location**: `src/codegen/traits.rs` (25 instances)

**Problematic Patterns**:
```rust
#[allow(unused_variables)]  // 18 instances - for wrapper constraints
#[allow(unused_imports)]    // 4 instances - for generated code
#[allow(clippy::vec_init_then_push)]  // 1 instance
#[allow(clippy::type_complexity)]     // 1 instance
#[allow(clippy::missing_docs_in_private_items)]  // 2 instances
```

**Root Cause**: Generated code sometimes produces unused variables when constraints simplify products to constants. Instead of suppressing, the generator should:
1. Conditionally generate variable bindings only when used
2. Use `_` prefix for intentionally unused bindings
3. Generate cleaner code paths for trivial products

### Issue 3: Identical Conditional Branches (Dead Code Pattern)

**Severity**: HIGH
**Location**: `src/codegen/traits.rs` (12+ locations)

**Example** (lines 2665-2669):
```rust
let constructor_call = if operand.versor.is_some() {
    quote! { #operand_name::new_unchecked(#(#field_exprs),*) }
} else {
    quote! { #operand_name::new_unchecked(#(#field_exprs),*) }  // IDENTICAL!
};
```

**Affected Functions**:
- `generate_sandwich_trait_from_versor()` (line 2665)
- `generate_antisandwich_trait_from_versor()` (line 2699)
- `generate_project_trait()` (line 2299)
- `generate_antiproject_trait()` (line 2334)
- `generate_wrapper_product_trait()` (line 2382)
- `generate_wrapper_sandwich_trait()` (line 2863)
- `generate_wrapper_antisandwich_trait()` (line 2943)
- And 5+ more locations

**Fix**: Remove dead conditional branches or implement the intended differentiation.

### Issue 4: Code Duplication in Traits Generator

**Severity**: HIGH
**Location**: `src/codegen/traits.rs`

**Duplicated Patterns**:

| Pattern | Locations | Lines Each |
|---------|-----------|------------|
| `compute_sandwich_expressions` / `compute_antisandwich_expressions` | Lines 715-756 | ~40 |
| `compute_project_expressions` / `compute_antiproject_expressions` | Lines 893-921 | ~30 |
| `generate_sandwich_trait` / `generate_antisandwich_trait` | Lines 2653, 2687 | ~50 |
| Expression string building | 4 locations in constraints.rs | ~15 |
| Symbol creation | constraint_derive.rs:111, product.rs:83 | ~10 |
| Degenerate index checking | constraint_derive.rs: 391, 428, 464 | ~20 |

**Fix**: Extract parametric helper functions that take product kind as argument.

### Issue 5: Legacy/Stub Module

**Severity**: MEDIUM
**Location**: `src/symbolic/constraint_simplify.rs`

**Problem**: This module is effectively dead code:
- Constructor returns `Vec::new()` (line 45)
- `apply()` method returns input unchanged when empty (line 54)
- Comment states: "can be reimplemented using ConstraintDeriver if needed"
- Despite being instantiated in `traits.rs`, it does nothing

**Decision Required**: Remove entirely or implement properly.

---

## Technical Debt Inventory

### 3.1 TODO Comments

| Location | Line | Content | Priority |
|----------|------|---------|----------|
| `src/spec/parser.rs` | 883 | PRD-45 blade-level inference for sparse types | Medium |
| `src/codegen/types.rs` | 672 | Handle non-Euclidean metrics for norm computation | High |
| `src/discovery/template.rs` | 54, 55, 178, 184, 190, 193 | Template placeholders (acceptable) | Low |
| `README.md` | 47, 48 | Template placeholders (acceptable) | Low |

### 3.2 Unused Functions

| Function | Location | Status |
|----------|----------|--------|
| `suggest_required_constraints()` | `discovery/mod.rs:215-246` | Public but not exported, never called |
| `normalize_fields()` | `symbolic/parser.rs:111` | Parameter `_field_names` unused |
| `ExpressionSimplifier::simplify()` | `symbolic/simplify.rs:36` | Just calls `expand()` - adds no value |
| Error variants in `error.rs` | Lines 70-78, 150-173, 130-137 | Defined but never triggered |

### 3.3 Inconsistent Patterns

| Pattern | Issue | Files |
|---------|-------|-------|
| HashMap import | `std::collections::HashMap::new()` vs imported | traits.rs vs projections.rs |
| Error handling | Some use `Result`, some `unwrap()`, some `unwrap_or_else()` | Multiple |
| Clone usage | Unnecessary `.clone()` in 15+ locations | traits.rs, spec/parser.rs |
| Function params | 4-5 parameter functions that should use builder pattern | 8+ functions in traits.rs |

### 3.4 File Size Imbalance

| File | Lines | Functions | Recommendation |
|------|-------|-----------|----------------|
| `traits.rs` | 5,558 | 201 | Split into: operators.rs, products.rs, wrappers.rs, normed.rs |
| `projections.rs` | 514 | ~20 | Fine |
| `conversions.rs` | ~200 | 13 | Fine |
| `types.rs` | ~800 | 55 | Fine |

### 3.5 Hardcoded Values

| Value | Location | Recommendation |
|-------|----------|----------------|
| `"--edition=2024"` | format.rs:38 | Extract to constant |
| `REL_EPSILON: f64 = 1e-10` | traits.rs:4641 | Make configurable or export |
| Signature-to-name mapping | conversions.rs:125-145 | Use registry pattern |
| Grade names 1-6 | template.rs:131-140 | Use formula for higher grades |

---

## TOML Specification Schema

### Current State

The TOML specification format is **defined only through examples**. The implicit schema is embedded across:
- `src/spec/raw.rs` - Serde deserialization structs
- `src/spec/ir.rs` - Intermediate representation
- `src/spec/parser.rs` - Validation logic
- 15 example algebra files in `algebras/`

### Proposed Formal Schema

Create a formal specification document at `docs/toml-schema.md`:

```toml
# CLIFFORD CODEGEN TOML SPECIFICATION v1.0

#=============================================================================
# [algebra] - REQUIRED
# Core metadata about the algebra
#=============================================================================
[algebra]
name = "string"           # REQUIRED: Unique identifier (e.g., "euclidean3")
module_path = "string"    # OPTIONAL: Rust module path (e.g., "euclidean::dim3")
description = "string"    # OPTIONAL: Documentation for the algebra
complete = bool           # OPTIONAL: Default true. If true, validates all products have matching output types

#=============================================================================
# [signature] - REQUIRED (at least one of positive/negative/zero must be non-empty)
# Defines the metric signature of the algebra
#=============================================================================
[signature]
positive = ["e1", "e2", ...]  # CONDITIONAL: Basis vectors squaring to +1
negative = ["e3", ...]        # OPTIONAL: Basis vectors squaring to -1
zero = ["e0", ...]            # OPTIONAL: Degenerate basis vectors squaring to 0

# Indexing Rules:
# - Numeric names (e1, e2, e3): Index = number - 1 (e1 -> 0, e2 -> 1)
# - Special case: e0 -> dim - 1 (PGA convention for degenerate basis)
# - Non-numeric names (ep, em): Positional indexing (positive first, then negative, then zero)
# - Total dimension (p + q + r) must be 1-6

#=============================================================================
# [norm] - OPTIONAL
# Controls which involution is used for norm computation
#=============================================================================
[norm]
primary_involution = "reverse"  # OPTIONAL: Default "reverse"
                                # Options: "reverse", "grade_involution", "clifford_conjugate"

#=============================================================================
# [blades] - OPTIONAL
# Custom display names for blades
#=============================================================================
[blades]
e12 = "xy"    # Maps internal blade name to display name
e123 = "I"    # Pseudoscalar display name

#=============================================================================
# [types.<TypeName>] - REQUIRED (at least one type)
# Defines algebraic types (elements of the algebra)
#=============================================================================
[types.Vector]
grades = [1]              # REQUIRED: Array of grades this type contains
description = "string"    # OPTIONAL: Documentation
field_map = [             # CONDITIONAL: Required unless alias_of is set
  { name = "x", blade = "e1" },
  { name = "y", blade = "e2" },
  { name = "z", blade = "e3" }
]
alias_of = "OtherType"           # OPTIONAL: Makes this type an alias
inverse_sandwich_targets = ["T"] # OPTIONAL: For non-versor types that support inverse sandwich

# Field Map Rules:
# - name: Unique field name within type, should be semantic (what it DOES, not blade index)
# - blade: Valid blade reference
#   - "s" for scalar (grade 0)
#   - "e1", "e2", etc. for basis vectors
#   - "e12", "e123", etc. for higher-grade blades
#   - Non-canonical orderings supported: "e21" = -e12 (sign auto-computed)

# Versor Detection (automatic):
# - Type is a versor if all grades have same parity (all even OR all odd)
# - Versors can transform other elements via sandwich product

# Sparse Type Detection (automatic):
# - Type is sparse if field count < expected blade count for its grades
# - Example: CGA Line has grade [3] but only 6 of 10 grade-3 blades
```

### Schema Implementation Approach

> **Note**: TOML has no native schema language. While JSON Schema can be used with tools like Taplo,
> we prefer documentation + enhanced CLI validation to avoid maintaining a separate schema format.

**Phase 1: Document Existing Behavior** (Effort: 2-3 days)
1. Create `docs/toml-schema.md` with complete field reference
2. Document all validation rules from `parser.rs`
3. Document all defaults from `raw.rs`
4. Add examples for each field

**Phase 2: Enhanced `validate` CLI Command** (Effort: 3-4 days)
1. Expand existing `validate` subcommand with comprehensive checks
2. Validate all fields against documented schema rules
3. Report all issues with helpful, actionable messages
4. Suggest fixes for common problems (e.g., "Missing [norm] section, add: primary_involution = \"reverse\"")
5. Add `--strict` mode that fails on warnings (missing optional sections, etc.)

**Phase 3: Migration Tool** (Effort: 1-2 days)
1. Add `migrate` subcommand
2. Updates old specs to current schema
3. Adds missing sections with correct defaults

---

## Implementation Plan

### Phase 1: Critical Fixes (Week 1)

| Task | Effort | Files |
|------|--------|-------|
| Add missing `[norm]` sections | 1h | 4 TOML files |
| Remove dead conditional branches | 4h | traits.rs |
| Remove `constraint_simplify.rs` stub | 1h | symbolic/ |
| Fix TODO in types.rs:672 | 2h | types.rs |

### Phase 2: Code Deduplication (Week 2)

| Task | Effort | Files |
|------|--------|-------|
| Extract sandwich/antisandwich helper | 4h | traits.rs |
| Extract projection/antiprojection helper | 3h | traits.rs, projections.rs |
| Extract symbol creation utility | 2h | symbolic/*.rs |
| Extract expression string builder | 1h | discovery/constraints.rs |

### Phase 3: Warning Suppression Removal (Week 2-3)

| Task | Effort | Files |
|------|--------|-------|
| Conditional variable generation | 8h | traits.rs |
| Clean up unused imports in generated code | 2h | traits.rs |
| Fix clippy warnings at source | 2h | table.rs, traits.rs |

### Phase 4: TOML Schema & Validation (Week 3-4)

| Task | Effort | Files |
|------|--------|-------|
| Write schema documentation | 16h | docs/toml-schema.md |
| Enhance `validate` CLI command | 12h | main.rs, spec/validate.rs |
| Add `migrate` CLI command | 6h | main.rs, spec/migrate.rs |

### Phase 5: File Reorganization (Week 4-5)

| Task | Effort | Files |
|------|--------|-------|
| Split traits.rs into modules | 8h | codegen/*.rs |
| Standardize error handling | 4h | Multiple |
| Remove unnecessary clones | 2h | traits.rs, parser.rs |
| Extract constants | 1h | format.rs, conversions.rs |

### Phase 6: TOML Standardization (Week 5-6)

| Task | Effort | Files |
|------|--------|-------|
| Standardize field naming | 4h | 15 TOML files |
| Add missing `[blades]` to conformal3 | 1h | conformal3.toml |
| Standardize descriptions | 4h | 15 TOML files |
| Document `complete` flag usage | 1h | All incomplete algebras |

---

## Effort Estimates

### Summary by Category

| Category | Estimated Effort | Priority |
|----------|-----------------|----------|
| Critical fixes | 8 hours | P0 |
| Code deduplication | 10 hours | P1 |
| Warning suppression removal | 12 hours | P1 |
| TOML schema documentation | 16 hours | P1 |
| Enhanced validate CLI | 12 hours | P1 |
| Migrate CLI command | 6 hours | P2 |
| File reorganization | 15 hours | P2 |
| TOML standardization | 10 hours | P2 |
| **Total** | **89 hours** | - |

### Recommended Phasing

**Sprint 1 (P0 + P1 Critical)**: 46 hours
- Critical fixes
- Code deduplication
- TOML schema documentation
- Enhanced validate CLI

**Sprint 2 (P1 Remaining)**: 12 hours
- Warning suppression removal

**Sprint 3 (P2)**: 31 hours
- Migrate CLI command
- File reorganization
- TOML standardization

---

## TOML Specification Inconsistencies (Detailed)

### Field Naming Patterns

| Algebra | Pattern | Example |
|---------|---------|---------|
| euclidean3 | Semantic | `rz, ry, rx` (rotation axes) |
| projective3 | Semantic | `tz, ty, tx` (translation), `rx, ry, rz` (rotation) |
| conformal2 | Mixed | `x, y, ep, em` (coords + basis names) |
| complex | Mathematical | `real, imag` |
| quaternion | Standard | `w, i, j, k` |
| euclidean2 | Generic | `b` (too vague) |

**Recommendation**: Document when to use each pattern in schema.

### Missing Features by Algebra

| Feature | Present In | Missing From |
|---------|-----------|--------------|
| `[norm]` section | 11 algebras | euclidean2, euclidean3, projective2, projective3 |
| `[blades]` section | conformal2, projective2, projective3 | conformal3 (inconsistent) |
| `inverse_sandwich_targets` | conformal2 | conformal3 (should it have it?) |
| Extensive documentation | conformal2, projective2/3 | Most others |

### Type Naming Inconsistencies

| Concept | Names Used | Files |
|---------|-----------|-------|
| Pure imaginary | "Imaginary" vs "ImagUnit" | quaternion vs complex |
| Full algebra element | "Spacetime" vs "Multivector" | minkowski2 vs - |
| Even subalgebra | "Rotor" vs "Motor" vs "Eventor" | euclidean vs projective vs minkowski |

---

## Success Criteria

1. **Zero `#[allow(...)]` suppressions** in codegen crate
2. **No dead conditional branches** (identical if/else)
3. **All algebras have explicit `[norm]` section**
4. **TOML schema documented** in `docs/toml-schema.md`
5. **Enhanced `validate` CLI** catches all schema violations with helpful messages
6. **No code duplication** > 20 lines
7. **traits.rs split** into logical modules (< 1500 lines each)
8. **All TODO comments** resolved or tracked in PRDs
9. **Consistent TOML field naming** documented and enforced

---

## Appendix: Complete Findings by Module

### A. algebra/ Module
- `table.rs:214`: `#[allow(clippy::type_complexity)]` for return type - consider type alias

### B. spec/ Module
- `parser.rs:415`: Legacy compatibility comment for blade name parsing
- `parser.rs:883`: TODO for PRD-45 blade-level inference
- Unused error variants in `error.rs`

### C. codegen/ Module
- `traits.rs`: 5,558 lines - needs splitting
- `traits.rs`: 25 `#[allow(...)]` suppressions
- `traits.rs`: 12+ identical conditional branches
- `traits.rs`: 4+ pairs of nearly-identical functions
- `types.rs:672`: TODO for non-Euclidean metrics
- `types.rs:142, 447, 483`: Unwraps without error context

### D. symbolic/ Module
- `constraint_simplify.rs`: Entire module is a stub
- `parser.rs:111`: Unused parameter
- `simplify.rs:36`: Trivial wrapper function
- `verify.rs:33-66`: Incomplete verifier stub
- `constraint_derive.rs`: Creates new simplifier on every equivalence check (inefficient)

### E. discovery/ Module
- `mod.rs:215-246`: Unused public function
- `naming.rs:26`: Unused `_dim` parameter
- `products.rs:1146-1148`: Disabled test validation
- `constraints.rs`: Expression string building duplicated 4x

### F. CLI (main.rs)
- Lines 175-177: CLI args accepted but not implemented (`serde`, `no_tests`, `no_arbitrary`)

---

## Related PRDs

- PRD-25: TOML Cleanup (completed - removed dead options)
- PRD-45: Blade-Level Inference (referenced TODO)
- PRD-46: Groebner Constraints
- PRD-47: Wrapper Constraint Simplification
