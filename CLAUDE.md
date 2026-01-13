# Clifford - Geometric Algebra Library

A Rust library for Geometric Algebra (Clifford Algebra).

## Project Principles

### 1. Educational Focus
- Code should be readable and well-documented to facilitate learning GA
- Include mathematical explanations in doc comments
- Link to resources and papers where appropriate

### 2. Documentation First
- Every public API must have comprehensive rustdoc
- Include examples in doc comments
- Explain the geometric intuition, not just the code

### 3. Clean Git History
- **Never push directly to `main`** - all changes go through pull requests
- **Always branch from latest `origin/main`**:
  ```bash
  git fetch origin main
  git checkout -b feat/your-feature origin/main
  ```
- Work in feature branches: `feat/multivector`, `fix/product-sign`
- **Small, logical commits** - each commit should be a single logical change
- Use conventional commit format: `type(scope): description`
- **Confirm before creating PRs** - always ask for user confirmation before running `gh pr create`
- **Review before merging** - wait for CI and Greptile, address feedback, then merge with `--squash --delete-branch`

### 4. Performance & Benchmarking
- Use SIMD instructions where beneficial (via `std::arch` or `portable_simd`)
- Benchmark critical paths with Criterion
- Profile before optimizing
- See the **optimize agent** for detailed benchmarking workflows

### 5. Minimal Dependencies
- Prefer std library where possible
- Only add dependencies when they provide significant value
- Dependencies must use permissive licenses (run `cargo deny check`)

### 6. Rust API Guidelines
- Follow the official [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- Implement standard traits (Debug, Clone, PartialEq, etc.)
- **No warning suppression** - never use `#[allow(...)]`, fix root causes instead
- **Private fields with public accessors** for specialized types
- **Use Clifford types in APIs**, not tuples

### 7. Testing
- **Property-based testing is mandatory** with `proptest`
- Use `prop_assert!` inside `proptest!` blocks
- Use `relative_eq!` with both `epsilon` AND `max_relative` parameters
- Use `RELATIVE_EQ_EPS` constant from `crate::test_utils`
- **Symbolica tests** in `clifford-codegen` must be prefixed with `symbolica_`
- See the **test agent** for detailed testing patterns

### 8. Code Review
- PRs are reviewed by Greptile (AI-powered review)
- Address feedback before merging; comment `@greptileai review` after fixes
- **No `todo!()` macros** - code must be complete

### 9. Claude Code Agents

Specialized agents in `.claude/agents/` handle different tasks:

| Agent | Purpose |
|-------|---------|
| **implement** | Implementing features, code changes, codegen |
| **test** | Writing property-based tests |
| **document** | Writing rustdoc and documentation |
| **review** | Code review checklist |
| **explain** | Teaching GA concepts |
| **optimize** | Performance optimization, benchmarking |
| **devops** | CI/CD and infrastructure |
| **release** | Version bumps, publishing |

### 10. Code Generation

**CRITICAL: Do NOT manually derive algebraic formulas.**

Use the `clifford-codegen` tool for all algebraic operations:
- Never manually derive geometric products, sandwich products, or transformations
- If codegen produces wrong results, fix the codegen tool, not the generated code
- Regenerate all algebras after codegen changes

See the **implement agent** for detailed codegen usage.

## Development Commands

```bash
cargo build           # Build the library
cargo nextest run     # Run tests (recommended - handles Symbolica correctly)
cargo test            # Run tests (fallback)
cargo doc --open      # Generate and view documentation
cargo bench           # Run benchmarks
cargo clippy          # Run linter
cargo fmt             # Format code
cargo deny check      # Check licenses and advisories
```

### Verification Workflow

**Before every commit**:
```bash
cargo fmt && cargo clippy && cargo doc --no-deps && cargo nextest run && cargo deny check
```

## Module Structure

```
specialized/
  euclidean/
    dim2/             # 2D: Vector, Bivector, Rotor
    dim3/             # 3D: Vector, Bivector, Trivector, Rotor
  projective/
    dim3/             # 3D PGA: Point, Line, Plane, Motor, Flector
```

Naming: Use full words (`Vector` not `Vec`), don't include dimension in type names.

## Resources

### Authoritative Reference

**The [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page) is the authoritative reference for all GA operations.**

### RGA Product Notation

| Symbol | Name | Function Pattern |
|--------|------|------------------|
| `∧` | wedge | `wedge_*` |
| `∨` | antiwedge | `antiwedge_*` |
| `★` | dual | `dual_*` |
| `☆` | antidual | `antidual_*` |

**Interior Products**: `bulk_contraction_*`, `weight_contraction_*`, `bulk_expansion_*`, `weight_expansion_*`

### Learning Resources
- [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page)
- [Look, Ma, No Matrices!](https://enkimute.github.io/LookMaNoMatrices/)
