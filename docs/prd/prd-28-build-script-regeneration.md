# PRD-28: Build Script Auto-Regeneration

**Status**: Complete
**Goal**: Automatically regenerate all algebra code on every build to prevent manual edits to generated files

## Background

The `clifford-codegen` tool generates Rust code for specialized geometric algebra types from TOML specifications. Previously, regeneration required manual CLI invocation:

```bash
cargo run -p clifford-codegen -- generate algebras/euclidean2.toml --output src/specialized/euclidean/dim2/generated --force
```

This created a risk: LLMs or developers could manually edit generated files in `src/specialized/*/generated/` instead of fixing the codegen tool when issues arise.

## Problem Statement

### Issue 1: Manual Edits to Generated Code
Without enforcement, generated code can be manually modified, leading to:
- Divergence between TOML specs and generated output
- Bugs that disappear when code is regenerated
- Violation of the principle: "fix the codegen tool, not the generated code"

### Issue 2: Stale Generated Code
Generated files could become out of sync with:
- Changes to the codegen tool
- Updates to TOML specifications
- New features added to code generation

### Issue 3: LLM Trust Boundary
LLMs assisting with development might take shortcuts by patching generated code directly rather than fixing root causes in the codegen system.

## Solution

Add a `build.rs` script that regenerates all algebras on every build:

```rust
use clifford_codegen::algebra::{Algebra, ProductTable};
use clifford_codegen::codegen::{ConversionsGenerator, TraitsGenerator, TypeGenerator, format_tokens};
use clifford_codegen::spec::parse_spec;

const ALGEBRAS: &[(&str, &str, &str)] = &[
    ("euclidean2", "algebras/euclidean2.toml", "src/specialized/euclidean/dim2/generated"),
    ("euclidean3", "algebras/euclidean3.toml", "src/specialized/euclidean/dim3/generated"),
    ("projective2", "algebras/projective2.toml", "src/specialized/projective/dim2/generated"),
    ("projective3", "algebras/projective3.toml", "src/specialized/projective/dim3/generated"),
];

fn main() {
    for (name, toml_path, output_dir) in ALGEBRAS {
        generate_algebra(name, toml_path, output_dir);
    }
}
```

### Key Design Decisions

1. **Regenerate on every build** - No caching or change detection; always overwrites
2. **Write to `src/`** - Files are committed to git and visible for review
3. **Use codegen as build-dependency** - `clifford-codegen` compiles first, then generates code

### Trade-offs

| Benefit | Cost |
|---------|------|
| Guarantees generated code matches specs | Slower incremental builds |
| Prevents manual edits from persisting | File lock contention with rust-analyzer |
| Enforces codegen-first development | Must wait for codegen on every build |

## Implementation

### Cargo.toml Changes

```toml
[build-dependencies]
clifford-codegen = { path = "crates/clifford-codegen" }
proc-macro2 = "1"
quote = "1"
```

### build.rs

The build script:
1. Iterates over all algebra configurations
2. Parses each TOML specification
3. Builds the algebra and product table
4. Generates types, traits, and conversions
5. Writes formatted output to the generated directories

### Generated Files

For each algebra, generates:
- `mod.rs` - Module declarations
- `types.rs` - Type definitions and constructors
- `traits.rs` - Product implementations and tests
- `conversions.rs` - Type conversions

## Testing

1. `cargo clean && cargo build` - Verify clean build succeeds
2. Manually edit a generated file, then `cargo build` - Verify edit is overwritten
3. `cargo nextest run` - Verify all tests pass with regenerated code

## Success Criteria

1. **Build script runs on every build** - Verified by build output
2. **Manual edits are overwritten** - Any changes to `src/specialized/*/generated/` are lost on next build
3. **All existing tests pass** - No regressions from regeneration
4. **Codegen changes propagate automatically** - Updating codegen tool affects all algebras

## Completed Work

### build.rs Implementation
- Created `build.rs` at project root
- Added `clifford-codegen`, `proc-macro2`, `quote` as build-dependencies
- Implemented `generate_algebra()` function using codegen library directly
- Added `generate_mod_file()` helper for module declarations
- Configured all four algebras: euclidean2, euclidean3, projective2, projective3

### Code Formatting
- Added `rustfmt` call after writing each generated file
- Uses `--edition 2024` flag for consistent formatting
- Gracefully handles rustfmt failures (logs warning, doesn't fail build)

### License Exceptions (deny.toml)
- Added `LicenseRef-Symbolica` to allowed licenses with clarification
- Added LGPL-3.0 exceptions for `gmp-mpfr-sys` and `rug`
- These are build-time only dependencies; generated code is not a derivative work

### Git Configuration
- Generated directories added to `.gitignore`
- Generated files removed from git tracking
- Users must build to generate code (enforces regeneration)

### Verification
- Clean build succeeds (~32s including codegen compilation)
- Incremental builds succeed (~2s)
- All tests pass with regenerated code
- `cargo deny check` passes with license exceptions
