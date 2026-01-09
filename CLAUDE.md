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
  - Separate different modules (e.g., euclidean::dim2 and euclidean::dim3 in different commits)
  - Separate refactoring from new features
  - Each commit should be independently reviewable
- Each commit should be buildable and pass tests
- Use conventional commit format: `type(scope): description`
- **Confirm before creating PRs** - always ask for user confirmation before running `gh pr create`. This minimizes churn on GitHub and ensures the user has reviewed the changes.
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

# Copy generic benchmarks
for svg in target/criterion/generic_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#generic_}
  cp "$svg" "benches/reports/generic/${name}.svg"
done

# Copy euclidean/dim2 benchmarks
for svg in target/criterion/euclidean_dim2_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#euclidean_dim2_}
  cp "$svg" "benches/reports/euclidean/dim2/${name}.svg"
done

# Copy euclidean/dim3 benchmarks
for svg in target/criterion/euclidean_dim3_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#euclidean_dim3_}
  cp "$svg" "benches/reports/euclidean/dim3/${name}.svg"
done

# Update benches/README.md with new timing data if needed
# Commit the updated SVGs and README
```

#### Benchmark File Structure

- `benches/generic.rs` - Benchmarks for generic `Multivector` operations
- `benches/specialized.rs` - Benchmarks for specialized 2D/3D euclidean types
- `benches/reports/` - SVG timing distribution plots (committed to repo)
  - `generic/` - Generic multivector benchmarks
  - `euclidean/dim2/` - Specialized 2D benchmarks
  - `euclidean/dim3/` - Specialized 3D benchmarks
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

#### Module Structure and Naming

The `specialized` module is organized by algebra type, then dimension:

```
specialized/
  euclidean/          # Euclidean GA (standard geometry)
    dim2/             # 2D: Vector, Bivector, Rotor
    dim3/             # 3D: Vector, Bivector, Trivector, Rotor, Even
  projective/         # (future) Projective GA
  conformal/          # (future) Conformal GA
```

**Naming conventions for specialized types:**
- **Don't include dimension in type names** - The module path provides context:
  - Good: `euclidean::dim3::Vector` (clear from path)
  - Avoid: `euclidean::dim3::Vec3` (redundant "3")
- **Use full words for clarity**:
  - `Vector` not `Vec` (avoids confusion with `std::vec::Vec`)
  - `Bivector` not `Bivec`
  - `Trivector` not `Trivec`
- **Same names across dimensions** - `dim2::Vector` and `dim3::Vector` are both just `Vector`
- **Use module prefixes when importing from multiple dimensions**:
  ```rust
  use clifford::specialized::euclidean::{dim2, dim3};

  let v2 = dim2::Vector::new(1.0, 2.0);
  let v3 = dim3::Vector::new(1.0, 2.0, 3.0);
  ```

### 7. Testing
- **Property-based testing is mandatory**: Use `proptest` for all tests where possible. Tests that only pass for hardcoded inputs are insufficient—correctness must hold across the full input domain.
- **Implement `Arbitrary` trait** for types instead of writing free functions:
  ```rust
  // Good: implement Arbitrary trait, use any::<Type>()
  impl Arbitrary for Vector<f64> {
      type Parameters = ();
      type Strategy = BoxedStrategy<Self>;
      fn arbitrary_with(_: Self::Parameters) -> Self::Strategy { ... }
  }
  // Then use: any::<Vector<f64>>()

  // Avoid: free functions
  fn arb_vec3() -> impl Strategy<Value = Vector<f64>> { ... }
  ```
- **Arbitrary modules**: Each module with types has an `arbitrary` submodule containing:
  - Generic `Arbitrary` implementations for all types using `Float::from_f64()` for conversion
  - Generic wrapper types for constrained values (e.g., `NonZeroVector<T>`, `UnitVector<T>`, `UnitRotor<T>`)
  - Compile with `#[cfg(any(test, feature = "proptest-support"))]`
  ```rust
  // Import wrapper types from the arbitrary module
  use crate::specialized::euclidean::dim3::arbitrary::{NonZeroVector, UnitVector, UnitRotor};

  // Use any::<Type<f64>>() - always specify the float type explicitly
  // Internal tests should use f64 for consistency
  proptest! {
      #[test]
      fn rotor_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
          let rotated = r.rotate(v);  // Deref allows direct method access
          prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));
      }
  }
  ```
- **Generic Arbitrary pattern**: All types use generic `Arbitrary` impls with `Float::from_f64()`:
  ```rust
  // Base types generate f64 values and convert
  impl<T: Float + Debug> Arbitrary for Vector<T> {
      fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
          (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
              .prop_map(|(x, y, z)| Vector::new(T::from_f64(x), T::from_f64(y), T::from_f64(z)))
              .boxed()
      }
  }

  // Wrapper types use where clauses requiring the inner type to be Arbitrary
  impl<T> Arbitrary for NonZeroVector<T>
  where
      T: Float + Debug,
      Vector<T>: Arbitrary + Debug,
      <Vector<T> as Arbitrary>::Strategy: 'static,
  { ... }
  ```
- **proptest-support feature**: External consumers enable `proptest-support` feature to access arbitrary modules
- **Arbitrary wrapper ergonomics**: All wrapper types implement `Deref`, `AsRef`, `From`, and `into_inner()` for easy access to the inner value
- **Use `approx` crate for comparisons**: Never hand-roll floating-point comparisons. Use `abs_diff_eq!`, `relative_eq!`, or `ulps_eq!` macros from the `approx` crate.
- **Use `ABS_DIFF_EQ_EPS` constant**: Don't use magic numbers like `1e-10` for epsilon values. Use the standard constant `ABS_DIFF_EQ_EPS` defined in `src/lib.rs::test_utils`:
  ```rust
  use crate::test_utils::ABS_DIFF_EQ_EPS;
  use approx::abs_diff_eq;

  // Good: use the standard constant
  assert!(abs_diff_eq!(a.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));

  // Avoid: magic numbers
  assert!(abs_diff_eq!(a.norm(), 1.0, epsilon = 1e-10));
  ```
  For integration tests (`tests/` directory), define the constant locally since `test_utils` is `pub(crate)`.
- **Use `prop_assert!` in proptest blocks**: Inside `proptest!` blocks, always use `prop_assert!` instead of `assert!`. This provides better error reporting with counterexamples:
  ```rust
  proptest! {
      #[test]
      fn rotor_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
          let rotated = r.rotate(v);
          // Good: prop_assert! with standard epsilon constant
          prop_assert!(abs_diff_eq!(v.norm(), rotated.norm(), epsilon = ABS_DIFF_EQ_EPS));

          // Avoid: assert! loses proptest's shrinking and reporting benefits
          // Avoid: magic numbers like 1e-9
      }
  }
  ```
- All types implement `AbsDiffEq`, `RelativeEq`, and `UlpsEq` traits for f32 and f64 variants
- The `approx` traits are re-exported from `clifford::prelude`
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
cargo build --all-features    # Build the library with all features
cargo test --all-features     # Run all tests
cargo doc --open              # Generate and view documentation
cargo bench                   # Run benchmarks
cargo clippy --all-features   # Run linter
cargo fmt                     # Format code
```

## Verification Workflow

**Before every commit**, run these commands to ensure CI will pass:

```bash
cargo fmt                         # Format code (CI checks this!)
cargo clippy --all-features       # Lint check
cargo doc --all-features --no-deps # Documentation build (CI checks this!)
cargo test --all-features         # Run all tests including doc tests
```

CI will reject PRs that fail any of these checks. Always run `cargo fmt` before committing.

**Important**: Always use `--all-features` to ensure all code paths are tested, including feature-gated modules like `proptest-support`.

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
- [x] **PRD-4: Specialized** - optimized 2D/3D types (Vector, Vector, Rotor, Rotor, etc.) with conversions to generic Multivector
- [ ] **PRD-5: PGA** - Projective GA, motors
- [ ] **PRD-6: CGA** - Conformal GA, polish
