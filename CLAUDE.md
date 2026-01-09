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
- Work in feature branches, e.g. `feat/multivector`, `fix/product-sign`
- Atomic commits with clear messages
- Each commit should be buildable and pass tests
- Use conventional commit format: `type(scope): description`
- **Review before merging** - after creating a PR:
  ```bash
  gh pr create --title "..." --body "..."
  gh pr checks <PR_NUMBER> --watch  # Wait for CI and Greptile
  # Review and address Greptile comments
  gh pr merge <PR_NUMBER> --squash --delete-branch
  ```
  Do not auto-merge. Always review Greptile feedback before merging.
- **Wait for PRs to merge before starting new work** - don't begin the next task until the current PR has merged. This prevents cascading issues if a PR fails.

### 4. Performance
- Use SIMD instructions where beneficial (via `std::arch` or `portable_simd`)
- Benchmark critical paths
- Profile before optimizing

### 5. Minimal Dependencies
- Prefer std library where possible
- Only add dependencies when they provide significant value
- Document why each dependency is needed

### 6. Rust API Guidelines
- Follow the official Rust API Guidelines (https://rust-lang.github.io/api-guidelines/)
- Ergonomic builder patterns where appropriate
- Implement standard traits (Debug, Clone, PartialEq, etc.)
- **Prefer structs with associated methods over free functions**
  - Avoid primitive obsession (e.g., using `usize` instead of a proper `Blade` type)
  - Encapsulate internal details; don't leak implementation (e.g., bitmask indices)
  - Methods on types are more discoverable and provide better IDE support
  - Example: `blade.grade()` is better than `grade_of_blade(index)`

### 7. Testing
- **Property-based testing is mandatory**: Use `proptest` for all tests where possible. Tests that only pass for hardcoded inputs are insufficient—correctness must hold across the full input domain.
- Unit tests with specific examples are acceptable only when property-based testing is not feasible
- Doc tests for examples
- All tests run automatically via GitHub Actions CI on every push and PR

### 8. Code Review
- PRs are reviewed by Greptile (AI-powered review)
- Address or acknowledge Greptile feedback before merging
- **Retrigger Greptile after pushing new commits**: Comment `@greptileai` on the PR to request a new review

## Development Commands

```bash
cargo build           # Build the library
cargo test            # Run all tests
cargo doc --open      # Generate and view documentation
cargo bench           # Run benchmarks
cargo clippy          # Run linter
cargo fmt             # Format code
```

## Verification Workflow

**Before every commit**, run these commands to ensure CI will pass:

```bash
cargo fmt             # Format code (CI checks this!)
cargo clippy          # Lint check
cargo test            # Run all tests including doc tests
```

CI will reject PRs that fail any of these checks. Always run `cargo fmt` before committing.

## Architecture Notes

(To be expanded as the library develops)

## Status

### Completed
- [x] Initial project setup with Cargo.toml (edition 2024)
- [x] GitHub Actions CI workflow (check, test, fmt, clippy, audit)
- [x] Strict lint configuration: deny all warnings, require docs on all items
- [x] GitHub repository created: https://github.com/DevonMorris/clifford
- [x] Branch protection: PRs required, CI must pass, no direct pushes to main
- [x] Minimal README
- [x] Claude Code agents for specialized workflows (implement, test, document, review, explain, devops)
- [x] Greptile AI code review integration
- [x] proptest dependency for property-based testing
- [x] crates.io publish workflow (triggered by version tags)
- [x] Implementation PRDs (see `docs/prd/`)

### Next Steps
- [x] **PRD-1: Foundation** - Float trait, Signature trait, Blade type
- [ ] **PRD-2: Core Multivector** - Multivector type, geometric product
- [ ] **PRD-3: Products** - inner, outer, regressive products, grade operations
- [ ] **PRD-4: Specialized** - optimized 2D/3D types
- [ ] **PRD-5: PGA** - Projective GA, motors
- [ ] **PRD-6: CGA** - Conformal GA, polish
