# PRD-19: Symbolica Testing Strategy

## Overview

This PRD defines a testing strategy for the clifford-codegen crate that accommodates Symbolica's global state requirements while maintaining fast parallel execution for the majority of tests.

## Problem Statement

Symbolica uses global state internally that causes conflicts when multiple tests run in parallel. Currently, we have:

1. **26+ Symbolica tests** in `crates/clifford-codegen/src/symbolic/` marked with `#[ignore = "Symbolica global state conflicts"]`
2. **400+ tests** in the main library that can run in parallel
3. **Property-based tests** (proptest) that are computationally expensive and benefit significantly from parallel execution

The current workaround (`#[ignore]`) means Symbolica tests are never run in CI, leaving a significant portion of the codegen logic unverified. We need a strategy that:

- Runs Symbolica tests (serially) to ensure correctness
- Runs other tests in parallel for speed
- Works consistently in both local development and CI
- Is simple to understand and maintain

## Goals

1. **Run all tests in CI** - No tests should be permanently ignored
2. **Maintain test speed** - Proptests and non-Symbolica tests should run in parallel
3. **Simple developer experience** - `cargo test` should "just work" locally
4. **Deterministic results** - Tests should pass consistently, not flaky due to race conditions

## Non-Goals

- Fixing Symbolica's global state (upstream issue, out of our control)
- Rewriting codegen to avoid Symbolica (too much effort, Symbolica is valuable)
- Supporting test parallelism within Symbolica tests (not possible due to global state)

## Analysis

### Current Test Distribution

| Location | Test Count | Parallelizable | Notes |
|----------|------------|----------------|-------|
| `clifford` (main lib) | ~400 | Yes | Proptests, unit tests |
| `clifford-codegen` (non-symbolic) | ~50 | Yes | Parser, algebra, codegen |
| `clifford-codegen` (symbolic) | ~26 | **No** | Symbolica global state |

### Approaches Considered

#### 1. `serial_test` Crate
Use the `serial_test` crate to mark Symbolica tests with `#[serial]`.

**Pros:**
- Simple attribute-based approach
- Works with standard `cargo test`

**Cons:**
- Requires adding dependency
- All `#[serial]` tests across all crates serialize together, potentially slowing unrelated tests
- Doesn't integrate well with `cargo nextest`

#### 2. `cargo nextest` with Test Groups
Use nextest's test groups feature to run Symbolica tests in a single-threaded group.

**Pros:**
- No code changes needed (configuration only)
- Fine-grained control over which tests serialize
- Fast parallel execution for everything else
- Good CI integration

**Cons:**
- Requires nextest (additional tool)
- Configuration file to maintain
- Developers need to know to use `cargo nextest run`

#### 3. Separate Integration Test Binary
Move Symbolica tests to a separate integration test file that runs serially.

**Pros:**
- Clear separation of concerns
- Works with standard `cargo test`
- Each test binary can have different parallelism settings

**Cons:**
- Requires restructuring test code
- May need to expose internal APIs for testing
- More complex project structure

#### 4. Process-Level Isolation
Run Symbolica tests in separate processes using `#[test]` with environment checks.

**Pros:**
- Complete isolation
- Works with any test runner

**Cons:**
- Complex implementation
- Slower due to process overhead
- Not a standard pattern

#### 5. Feature Flag for Symbolica Tests
Use a feature flag to conditionally compile Symbolica tests.

**Pros:**
- Clear separation
- Can run different test sets easily

**Cons:**
- Tests still need to be serialized when enabled
- Doesn't solve the core parallelism problem

## Recommended Approach

**Combination of approaches 2 and 3**: Use `cargo nextest` with test groups as the primary strategy, with clear naming conventions for Symbolica tests.

### Implementation Plan

#### Phase 1: Test Organization

1. **Rename Symbolica tests with prefix** - All Symbolica tests get a `symbolica_` prefix:
   ```rust
   #[test]
   fn symbolica_parser_handles_constraint() { ... }
   ```

2. **Remove `#[ignore]` attributes** - Tests should run, not be skipped:
   ```rust
   // Before
   #[ignore = "Symbolica global state conflicts"]
   #[test]
   fn verify_motor_constraint() { ... }

   // After
   #[test]
   fn symbolica_verify_motor_constraint() { ... }
   ```

3. **Add module-level documentation** explaining the naming convention:
   ```rust
   //! # Test Naming Convention
   //!
   //! Tests that use Symbolica's global state are prefixed with `symbolica_`.
   //! These tests must run serially. Use `cargo nextest run` which automatically
   //! handles this via `.config/nextest.toml`.
   ```

#### Phase 2: Nextest Configuration

Create `.config/nextest.toml`:

```toml
[profile.default]
# Default test threads (use all available)
test-threads = "num-cpus"

# Symbolica tests must run serially due to global state
[[profile.default.overrides]]
filter = "test(symbolica_)"
test-threads = 1
```

This configuration:
- Runs tests matching `symbolica_*` with a single thread
- Runs all other tests with full parallelism
- Applies automatically without developer intervention

#### Phase 3: CI Integration

Update `.github/workflows/ci.yml`:

```yaml
test:
  name: Test
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@stable
    - uses: Swatinem/rust-cache@v2
    - uses: taiki-e/install-action@nextest
    - run: cargo nextest run --features proptest-support,serde

test-codegen:
  name: Test Codegen
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@stable
    - uses: Swatinem/rust-cache@v2
    - uses: taiki-e/install-action@nextest
    - run: cargo nextest run --package clifford-codegen
```

#### Phase 4: Local Development Support

Add to `CLAUDE.md` and project documentation:

```markdown
## Running Tests

We use `cargo nextest` for test execution to properly handle Symbolica's global state:

```bash
# Run all tests (recommended)
cargo nextest run

# Run only clifford-codegen tests
cargo nextest run --package clifford-codegen

# Run with standard cargo test (works but Symbolica tests may conflict)
cargo test
```

### Why nextest?

The codegen crate uses Symbolica for symbolic algebra, which has global state
that conflicts when tests run in parallel. Nextest's test groups feature allows
us to run Symbolica tests serially while keeping everything else parallel.
```

### Fallback for `cargo test`

For developers who prefer `cargo test`, we can add a compile-time check:

```rust
#[cfg(test)]
mod tests {
    use std::sync::Mutex;

    /// Global lock for Symbolica tests to prevent parallel execution.
    /// This is a fallback for when tests are run with `cargo test` instead of nextest.
    static SYMBOLICA_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn symbolica_example_test() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        // ... test code using Symbolica ...
    }
}
```

This provides a safety net but with performance cost (mutex overhead and potential deadlocks if tests panic).

## Migration Checklist

- [ ] Create `.config/nextest.toml` with Symbolica test group
- [ ] Rename all Symbolica tests with `symbolica_` prefix
- [ ] Remove all `#[ignore = "Symbolica global state conflicts"]` attributes
- [ ] Add `SYMBOLICA_LOCK` mutex as fallback for `cargo test`
- [ ] Update CI to use `cargo nextest run`
- [ ] Update `CLAUDE.md` with testing instructions
- [ ] Update project README with testing instructions
- [ ] Verify all tests pass locally with `cargo nextest run`
- [ ] Verify all tests pass in CI

## Files to Create/Modify

### New Files
- `.config/nextest.toml` - Nextest configuration

### Modified Files
- `.github/workflows/ci.yml` - Use nextest
- `CLAUDE.md` - Add testing instructions
- `crates/clifford-codegen/src/symbolic/parser.rs` - Rename tests, remove ignore
- `crates/clifford-codegen/src/symbolic/product.rs` - Rename tests, remove ignore
- `crates/clifford-codegen/src/symbolic/verify.rs` - Rename tests, remove ignore
- `crates/clifford-codegen/src/symbolic/simplify.rs` - Rename tests, remove ignore
- `crates/clifford-codegen/src/symbolic/to_rust.rs` - Rename tests, remove ignore
- `crates/clifford-codegen/src/symbolic/constraint_simplify.rs` - Rename tests, remove ignore
- `crates/clifford-codegen/src/symbolic/constraint_derive.rs` - Add prefix if needed

## Success Criteria

1. `cargo nextest run` passes all tests (including Symbolica tests)
2. CI runs all tests and passes
3. Test execution time remains reasonable (< 2 minutes for full suite)
4. `cargo test` works (with mutex fallback) for developers who prefer it
5. No tests are marked `#[ignore]` without a good reason

## Testing the Testing Strategy

Before full rollout:

1. Create a test branch with the nextest configuration
2. Rename one Symbolica test as a pilot
3. Run `cargo nextest run` locally multiple times to verify no flakiness
4. Push to CI and verify it passes
5. Proceed with remaining tests

## References

- [cargo-nextest documentation](https://nexte.st/)
- [nextest test groups](https://nexte.st/book/test-groups.html)
- [Symbolica issue on global state](https://github.com/benruijl/symbolica/issues) (if applicable)
- [serial_test crate](https://crates.io/crates/serial_test) (alternative approach)
