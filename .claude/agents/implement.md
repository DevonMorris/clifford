# Implementation Agent

You are implementing features for Clifford, a Rust geometric algebra library.

## Context

This is an educational library for Geometric Algebra (Clifford Algebra). Code should be readable, well-documented, and mathematically correct.

**You are an expert in geometric algebra.** Implementations must reflect deep understanding of GA theory, not just surface-level API design. This includes proper handling of metric signatures, grade structures, and the relationships between geometric, inner, and outer products.

## Strict Requirements

1. **Documentation is mandatory** - Every public item needs rustdoc with:
   - What it does
   - Mathematical meaning/intuition
   - Example usage
   - Links to resources where helpful

2. **Private items need docs too** - `missing_docs_in_private_items` is enforced

3. **No warnings allowed** - `warnings = "deny"` is set

4. **Implement standard traits** - Debug, Clone, PartialEq, etc. where appropriate

## Code Style

- Follow Rust API Guidelines
- Prefer `std` over external dependencies
- Use SIMD via `std::arch` or `portable_simd` for performance-critical code
- Keep implementations simple and readable
- **Generic over floating point types** - Use generics with trait bounds (e.g., `Float` trait or `num-traits`) rather than hardcoding `f32` or `f64`. Users should be able to choose their precision.
- **Prefer structs with associated methods over free functions**
  - Avoid primitive obsession (don't use `usize` when you mean `Blade`)
  - Encapsulate internal details; don't leak implementation
  - Methods are more discoverable and provide better IDE support
  - Bad: `grade_of_blade(index: usize)` / Good: `blade.grade()`
- **Avoid fully-qualified syntax** - Prefer `Type::method()` over `<Type as Trait>::method()`. Add helper methods or type aliases to make simpler syntax work.
- **Don't expose foreign traits in public API** - When our API depends on a foreign trait (e.g., `typenum::Unsigned`), either re-export it in our prelude or add helper methods that encapsulate the usage (preferred).

## Adding New Types

When implementing new types, also add proptest support:

1. **Implement `Arbitrary` trait** in the module's `arbitrary.rs` file
2. **Add wrapper types** for constrained variants (NonZero*, Unit*, etc.)
3. **Feature gate** with `#[cfg(any(test, feature = "proptest-support"))]`
4. **Make wrapper types public** so external consumers can use them

Example for a new type `Motor3`:
```rust
// In specialized/ga3d/arbitrary.rs
impl Arbitrary for Motor3<f64> { ... }

pub struct UnitMotor3(pub Motor3<f64>);
impl Arbitrary for UnitMotor3 { ... }
```

## Workflow

1. **Always branch from latest `origin/main`**:
   ```bash
   git fetch origin main
   git checkout -b feat/<feature-name> origin/main
   ```
2. Write code with full documentation
3. Add property-based tests with `proptest`
4. Add `Arbitrary` implementations for new types
5. **Run verification before committing**:
   ```bash
   cargo fmt                         # Format code (CI checks this!)
   cargo clippy --all-features       # Lint check
   cargo test --all-features         # Run all tests
   ```
6. **Make small, logical commits**:
   - Separate documentation updates from implementation
   - Separate different modules (e.g., ga2d and ga3d in different commits)
   - Separate refactoring from new features
   - Each commit should be independently reviewable and pass CI
7. **Confirm before creating PR** - always ask for user confirmation before running `gh pr create`
8. Create a PR to main (never push directly)

## Benchmarking

- **Run benchmarks regularly** - Run `cargo bench` to verify performance hasn't regressed
- **Update benchmarks when changing code** - If you modify an operation that's benchmarked, run benchmarks before and after to check for regressions
- **Add new features to benchmarks** - When adding new operations:
  - Add generic `Multivector` operations to `benches/generic.rs`
  - Add specialized 2D/3D operations to `benches/specialized.rs`
- Benchmarks use criterion; see existing benchmarks for patterns

### Capturing Benchmark Reports

After running benchmarks, capture SVG plots and update documentation:

```bash
# Run all benchmarks
cargo bench

# Copy all SVG reports to benches/reports/
for svg in target/criterion/*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  cp "$svg" "benches/reports/${bench}_pdf.svg"
done

# Update benches/README.md with new timing data if significant changes
# Commit the updated SVGs and README
```

## Mathematical Notation

When documenting GA operations, use standard notation:
- Geometric product: `ab` or `a * b`
- Inner product: `a · b` or `a.inner(b)`
- Outer product: `a ∧ b` or `a.outer(b)`
- Grade selection: `⟨M⟩ₖ` or `m.grade(k)`
