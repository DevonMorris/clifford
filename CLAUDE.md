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
  This prevents merge conflicts and rebase issues.
- Work in feature branches, e.g. `feat/multivector`, `fix/product-sign`
- **Small, logical commits** - each commit should be a single logical change:
  - Separate documentation updates from implementation
  - Separate different modules (e.g., ga2d and ga3d in different commits)
  - Separate refactoring from new features
  - Each commit should be independently reviewable
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

### 4. Performance & Benchmarking

**NOTE:** When updating benchmarking rules in this section, also update the corresponding agents in `.claude/agents/` (especially `implement.md` and `review.md`).

- Use SIMD instructions where beneficial (via `std::arch` or `portable_simd`)
- Benchmark critical paths
- Profile before optimizing

#### Benchmark Workflow

1. **Run benchmarks regularly** - Run `cargo bench` to verify performance hasn't regressed
2. **Update benchmarks when changing code** - If you modify an operation that's benchmarked, run benchmarks before and after to check for regressions
3. **Add new features to benchmarks** - When adding new operations:
   - Add generic operations to `benches/generic.rs`
   - Add specialized 2D/3D operations to `benches/specialized.rs`

#### Capturing Benchmark Reports

After running `cargo bench`, capture and commit SVG plots:

```bash
# Run all benchmarks
cargo bench

# Copy all SVG reports to benches/reports/
for svg in target/criterion/*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  cp "$svg" "benches/reports/${bench}_pdf.svg"
done

# Update benches/README.md with new timing data if needed
# Commit the updated SVGs and README
```

#### Benchmark File Structure

- `benches/generic.rs` - Benchmarks for generic `Multivector` operations
- `benches/specialized.rs` - Benchmarks for specialized 2D/3D types (`Vec2`, `Vec3`, `Rotor2`, `Rotor3`, etc.)
- `benches/reports/` - SVG timing distribution plots (committed to repo)
- `benches/README.md` - Performance summary with embedded SVG plots

### 5. Minimal Dependencies
- Prefer std library where possible
- Only add dependencies when they provide significant value
- Document why each dependency is needed

### 6. Rust API Guidelines

**NOTE:** When updating style rules in this section, also update the corresponding agents in `.claude/agents/` (especially `implement.md` and `review.md`).

- Follow the official Rust API Guidelines (https://rust-lang.github.io/api-guidelines/)
- Ergonomic builder patterns where appropriate
- Implement standard traits (Debug, Clone, PartialEq, etc.)
- **Prefer structs with associated methods over free functions**
  - Avoid primitive obsession (e.g., using `usize` instead of a proper `Blade` type)
  - Encapsulate internal details; don't leak implementation (e.g., bitmask indices)
  - Methods on types are more discoverable and provide better IDE support
  - Example: `blade.grade()` is better than `grade_of_blade(index)`
- **Avoid fully-qualified syntax** - Prefer `Type::method()` over `<Type as Trait>::method()` when possible. Add helper methods or type aliases to make simpler syntax work.
- **Don't expose foreign traits in public API** - When our public API depends on a foreign trait (e.g., `typenum::Unsigned`), either:
  - Re-export the trait in our prelude so users don't need to import the dependency directly, or
  - Add helper methods that encapsulate the foreign trait usage (preferred when feasible)

### 7. Testing
- **Property-based testing is mandatory**: Use `proptest` for all tests where possible. Tests that only pass for hardcoded inputs are insufficient—correctness must hold across the full input domain.
- **Implement `Arbitrary` trait** for types instead of writing free functions:
  ```rust
  // Good: implement Arbitrary trait, use any::<Type>()
  impl Arbitrary for Vec3<f64> {
      type Parameters = ();
      type Strategy = BoxedStrategy<Self>;
      fn arbitrary_with(_: Self::Parameters) -> Self::Strategy { ... }
  }
  // Then use: any::<Vec3<f64>>()

  // Avoid: free functions
  fn arb_vec3() -> impl Strategy<Value = Vec3<f64>> { ... }
  ```
- Unit tests with specific examples are acceptable only when property-based testing is not feasible
- Doc tests for examples
- All tests run automatically via GitHub Actions CI on every push and PR

### 8. Code Review
- PRs are reviewed by Greptile (AI-powered review)
- Address or acknowledge Greptile feedback before merging
- **IMPORTANT: After pushing commits that address Greptile feedback**, always comment `@greptileai review` on the PR to request a re-review. This ensures the fixes are verified before merging.

### 9. Claude Code Agents

This project has specialized agents in `.claude/agents/` for different tasks. **Use these agents proactively** when your work matches their purpose:

| Agent | When to Use |
|-------|-------------|
| **implement** | Implementing new features or making code changes. Use when adding functionality, writing new modules, or modifying existing code. |
| **test** | Writing or improving tests. Use when adding property-based tests, creating test strategies, or verifying algebraic properties. |
| **document** | Writing or improving documentation. Use when adding rustdoc, writing module docs, or improving mathematical explanations. |
| **review** | Reviewing code before merging. Use to check code against the project's quality checklist (docs, correctness, tests, style, performance). |
| **explain** | Teaching GA concepts. Use when explaining geometric algebra ideas, bridging from linear algebra, or helping users understand the math. |
| **optimize** | Performance optimization. Use when profiling, benchmarking, adding SIMD, or improving algorithmic efficiency. |
| **devops** | CI/CD and infrastructure tasks. Use for GitHub Actions, publishing, branch protection, or build/release automation. |

#### How to Invoke Agents

Run agents using the Task tool with the appropriate agent file:

```
Use the implement agent to add the new rotor type
Use the test agent to write property-based tests for the geometric product
Use the review agent to check this PR before merging
```

#### Agent Selection Guide

- **Writing code?** → `implement`
- **Writing tests?** → `test`
- **Writing docs?** → `document`
- **Checking quality?** → `review`
- **Teaching/explaining?** → `explain`
- **Improving speed?** → `optimize`
- **CI/CD/infra?** → `devops`

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
- [x] **PRD-2: Core Multivector** - Multivector type, geometric product
- [x] **PRD-3: Products** - inner, outer, regressive products, grade operations
- [x] **PRD-4: Specialized** - optimized 2D/3D types (Vec2, Vec3, Rotor2, Rotor3, etc.)
- [ ] **PRD-5: PGA** - Projective GA, motors
- [ ] **PRD-6: CGA** - Conformal GA, polish
