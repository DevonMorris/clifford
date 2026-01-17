# PRD-49: Workspace Restructure

**Status**: Draft
**Goal**: Move all crates into `crates/` directory for consistent workspace structure

## Motivation

Currently the workspace has an inconsistent structure:
```
clifford/              # Main crate at root
  src/
  Cargo.toml
  build.rs
  crates/
    clifford-codegen/  # Sub-crate in crates/
```

This is awkward when adding more crates. A cleaner structure places all crates in `crates/`:
```
clifford/              # Workspace root only
  Cargo.toml           # Workspace manifest
  crates/
    clifford/          # Main library
    clifford-codegen/  # Code generation
    clifford-viz/      # Visualization (PRD-48)
```

This pattern is used by bevy, rust-analyzer, tokio, and other major Rust projects.

## Migration Steps

### 1. Create New Directory Structure

```bash
mkdir -p crates/clifford
```

### 2. Move Main Crate Files

```bash
git mv src crates/clifford/
git mv build.rs crates/clifford/
git mv Cargo.toml crates/clifford/  # Will recreate workspace manifest at root
```

### 3. Create Workspace Manifest

Root `Cargo.toml` becomes workspace-only:

```toml
[workspace]
resolver = "2"
members = [
    "crates/clifford",
    "crates/clifford-codegen",
]

[workspace.package]
version = "0.1.0"
edition = "2024"
license = "MIT OR Apache-2.0"
repository = "https://github.com/DevonMorris/clifford"

[workspace.dependencies]
# Shared dependencies across crates
proptest = "1.6"
approx = "0.5"
num-traits = "0.2"
```

### 4. Update Main Crate Cargo.toml

`crates/clifford/Cargo.toml`:

```toml
[package]
name = "clifford"
version.workspace = true
edition.workspace = true
license.workspace = true
repository.workspace = true
description = "Geometric Algebra (Clifford Algebra) library for Rust"

[dependencies]
num-traits.workspace = true
# ... rest of dependencies

[build-dependencies]
clifford-codegen = { path = "../clifford-codegen" }

[dev-dependencies]
proptest.workspace = true
approx.workspace = true
```

### 5. Update clifford-codegen Path

`crates/clifford-codegen/Cargo.toml` - update any paths if needed.

### 6. Update build.rs Paths

`crates/clifford/build.rs` - update TOML spec paths:

```rust
// Before
const ALGEBRAS: ... = &[
    AlgebraConfig::new(
        "euclidean2",
        "algebras/euclidean2.toml",  // Relative to workspace root
        "src/specialized/euclidean/dim2/generated",
    ),
    // ...
];

// After
const ALGEBRAS: ... = &[
    AlgebraConfig::new(
        "euclidean2",
        "../../algebras/euclidean2.toml",  // Relative to crates/clifford/
        "src/specialized/euclidean/dim2/generated",
    ),
    // ...
];
```

**Alternative**: Use `CARGO_MANIFEST_DIR` for robustness:

```rust
fn get_workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()  // crates/
        .unwrap()
        .parent()  // workspace root
        .unwrap()
        .to_path_buf()
}

fn main() {
    let root = get_workspace_root();
    let toml_path = root.join("algebras/euclidean2.toml");
    // ...
}
```

### 7. Update CI Workflows

`.github/workflows/*.yml` - update any hardcoded paths:

```yaml
# Before
- run: cargo test

# After (should work the same, but verify)
- run: cargo test -p clifford
```

### 8. Update Documentation

Files to update:
- `CLAUDE.md` - module structure paths
- `README.md` - any path references
- `docs/prd/*.md` - path references in PRDs

### 9. Update .gitignore

Ensure patterns still work:
```gitignore
/target/           # Workspace target (unchanged)
```

### 10. Verify rerun-if-changed Paths

In `build.rs`, the `cargo::rerun-if-changed` paths need updating:

```rust
// Before
println!("cargo::rerun-if-changed=algebras/euclidean2.toml");
println!("cargo::rerun-if-changed=crates/clifford-codegen/src");

// After
println!("cargo::rerun-if-changed=../../algebras/euclidean2.toml");
println!("cargo::rerun-if-changed=../clifford-codegen/src");
```

## Files to Move

| From | To |
|------|-----|
| `src/` | `crates/clifford/src/` |
| `build.rs` | `crates/clifford/build.rs` |
| `Cargo.toml` | `crates/clifford/Cargo.toml` |

## Files to Create

| File | Purpose |
|------|---------|
| `Cargo.toml` (root) | Workspace manifest |

## Files to Update

| File | Changes |
|------|---------|
| `crates/clifford/build.rs` | Update paths to algebras/ and codegen |
| `crates/clifford/Cargo.toml` | Update path to clifford-codegen |
| `crates/clifford-codegen/Cargo.toml` | Use workspace dependencies |
| `CLAUDE.md` | Update module structure docs |
| `.github/workflows/*.yml` | Verify/update paths |
| `.gitignore` | Verify patterns still match (e.g., `/target/`) |

## Verification Checklist

```bash
# From workspace root
cargo build                    # Builds all crates
cargo build -p clifford        # Builds main crate
cargo build -p clifford-codegen # Builds codegen
cargo test -p clifford         # Tests main crate
cargo doc --workspace          # Docs for all crates
cargo clippy --workspace       # Lint all crates
cargo fmt --all                # Format all crates
```

- [ ] `cargo build` succeeds
- [ ] `cargo test` passes
- [ ] `cargo doc` generates docs
- [ ] Generated code paths are correct
- [ ] CI workflows pass
- [ ] `algebras/*.toml` files are found by build.rs

## Rollback Plan

If issues arise:
```bash
git mv crates/clifford/src ./
git mv crates/clifford/build.rs ./
git mv crates/clifford/Cargo.toml ./
rmdir crates/clifford
# Restore original root Cargo.toml from git
```

## Future Crates

After this restructure, adding new crates is trivial:
```bash
cargo new crates/clifford-viz --lib
# Add to workspace members
```

## Dependencies

None - this is foundational.

## Dependent PRDs

- PRD-48: Visualization Showcase (will use `crates/clifford-viz/`)
