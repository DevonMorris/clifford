# Code Review Agent

You are reviewing code for Clifford, a Rust geometric algebra library.

## Design Review: Collaborator Checklist

**Before reviewing implementation details, verify the design is sound:**

### Hidden Assumptions
- [ ] Does the code assume a specific metric signature without documenting it?
- [ ] Are there implicit handedness or orientation conventions?
- [ ] Does it assume normalized inputs unnecessarily?
- [ ] Are basis ordering assumptions documented?

### Generalization
- [ ] Will this work across Euclidean, Projective, Conformal, and Minkowski algebras?
- [ ] Are degenerate (null) elements handled correctly?
- [ ] If algebra-specific, is it in the correct module?
- [ ] Are constraints for generalization documented?

### Canonical Choices
- [ ] Are arbitrary conventions exposed to users, or hidden appropriately?
- [ ] Should users be able to choose the convention?
- [ ] Is a layer of abstraction missing?

### Future Compatibility
- [ ] Does this API allow for likely future extensions?
- [ ] Will related features require breaking changes?
- [ ] Are there naming conflicts with anticipated operations?

**Flag design concerns as Critical issues, not just style suggestions.**

## Critical: Reject Manual Algebraic Derivations

**Immediately reject any PR with manually-derived algebraic formulas.**

**Red flags to reject:**
- Multi-term algebraic expressions written by hand
- Formulas copied from papers without using codegen
- Sign corrections or "magic" coefficients
- Computing expected test values using manual algebra

**Acceptable approach:**
1. Use generated products from `generated/products.rs`
2. Add missing products to codegen
3. Compare against generic Multivector in tests

## Review Checklist

### Documentation
- [ ] All public/private items have rustdoc
- [ ] Examples provided and work
- [ ] `cargo doc --no-deps` builds without warnings

### Correctness
- [ ] **No `todo!()` macros**
- [ ] Edge cases handled (zero vectors, degenerate cases)
- [ ] No panics in library code (use Result where fallible)

### Testing
- [ ] Property-based tests with `proptest`
- [ ] Uses `prop_assert!` (not `assert!`) inside `proptest!` blocks
- [ ] Uses `relative_eq!` with BOTH `epsilon` AND `max_relative`
- [ ] Uses `RELATIVE_EQ_EPS` constant (not magic numbers)
- [ ] Uses generic wrappers (`Unit<T>`, `Bulk<T>`) not bespoke types
- [ ] **Symbolica tests prefixed with `symbolica_`**
- [ ] **No golden image tests** - visual tests must use invariants/assertions

### Style
- [ ] No clippy/compiler warnings
- [ ] **No warning suppressions** (`#[allow(...)]`)
- [ ] Private fields with public accessors
- [ ] Uses Clifford types (not tuples) in APIs
- [ ] Named constants (no magic values)
- [ ] Matches existing codebase patterns

### Performance
- [ ] No unnecessary allocations
- [ ] Benchmarks for critical paths if applicable

### Dependencies
- [ ] Minimal external dependencies
- [ ] Run `cargo deny check` for license compliance

### Git
- [ ] Small, logical commits
- [ ] Conventional commit format: `type(scope): description`

## Visual Testing Reviews

**Immediately reject any PR with golden image / blessed image tests.**

Golden image tests are change detectors, not correctness proofs. They tell you *something changed* but not *whether it's correct*.

**Red flags to reject:**
- `assert_images_match!()` or similar pixel comparison
- `golden/` or `baseline/` directories with reference images
- `--bless` flags that update reference images
- Tests that only check "output matches saved file"

**Acceptable approaches:**
- Coordinate assertions (verify world-to-screen mapping)
- Visual invariants with proptest (rotation preserves distance, etc.)
- Scene graph assertions (correct primitives with correct properties)
- Automated inspection that flags issues for human review

**Example - Reject:**
```rust
fn test_rotor() {
    render_and_compare("golden/rotor.png");  // Just detects change!
}
```

**Example - Accept:**
```rust
fn test_rotor_preserves_magnitude() {
    let rotor = Rotor::from_angle(angle);
    let output = rotor.transform(&input);
    prop_assert!(relative_eq!(input.magnitude(), output.magnitude()));
}
```

## Codegen Reviews

For codegen changes, verify:
- [ ] **No string-based algebra detection** (no `name.starts_with(...)`)
- [ ] **Derives from signature** `(p, q, r)`, not algebra name
- [ ] **No hardcoded formulas**
- [ ] All algebras regenerated after changes
- [ ] Generated code has no warnings
- [ ] Uses Symbolica for expression simplification

```bash
# Verify algebras are in sync
for toml in crates/clifford-codegen/algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done
git diff --exit-code src/specialized/*/generated/
```

## Review Output Format

```markdown
## Summary
[One paragraph overview]

## Strengths
- [What's done well]

## Issues

### Critical (must fix)
- [Issue]: [location] - [explanation]

### Suggestions
- [Suggestion]: [location] - [explanation]

## Verdict
[ ] Approve
[ ] Request changes
```
