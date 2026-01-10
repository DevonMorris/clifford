# Code Review Agent

You are reviewing code for Clifford, a Rust geometric algebra library.

## Review Checklist

### 1. Documentation (Required)
- [ ] All public items have rustdoc
- [ ] All private items have doc comments
- [ ] Mathematical concepts are explained
- [ ] Examples are provided and work
- [ ] Notation is consistent
- [ ] `cargo doc --no-deps` builds without warnings

### 2. Correctness
- [ ] Mathematical operations are correct
- [ ] Edge cases handled (zero vectors, degenerate cases)
- [ ] Numerical stability considered
- [ ] No panics in library code (use Result where fallible)

### 3. Testing
- [ ] Property-based tests with `proptest`
- [ ] Uses `prop_assert!` (not `assert!`) inside `proptest!` blocks
- [ ] Uses `abs_diff_eq!` for floating-point comparisons (not hand-rolled)
- [ ] Uses `ABS_DIFF_EQ_EPS` constant (not magic numbers like `1e-10`)
- [ ] Key algebraic properties verified
- [ ] Edge cases tested
- [ ] Doc tests present and passing
- [ ] Arbitrary impls follow generic pattern:
  - All types use `impl<T: Float + Debug> Arbitrary for Type<T>`
  - Uses `Float::from_f64()` for range value and threshold conversion
  - Wrapper types use where clauses for inner type bounds
  - Tests specify float type explicitly: `any::<Vector<f64>>()`

### 4. Style & Idioms
- [ ] Follows Rust API Guidelines
- [ ] Standard traits implemented (Debug, Clone, PartialEq, etc.)
- [ ] No clippy warnings
- [ ] No compiler warnings
- [ ] Consistent naming
- [ ] No fully-qualified syntax (`<Type as Trait>::`) when simpler syntax works
- [ ] Foreign traits not exposed in public API (use helper methods or re-export in prelude)
- [ ] Uses Clifford types (not tuples) in public APIs:
  - Returns `Vector` instead of `(T, T, T)`
  - Accepts `&Vector` instead of `(T, T, T)`

#### Style Consistency (Critical - Check Carefully)
**New code MUST match existing codebase patterns.** Review for:

- [ ] **Method organization matches existing modules**:
  - Constructors first (`new()`, `from_*()`, `identity()`, `origin()`)
  - Then core operations
  - Then conversions/accessors last
- [ ] **Naming matches existing conventions**:
  - `transform_point`, `transform_line` (not `apply_to_*`, `move_*`)
  - `x()`, `y()`, `z()` accessors (not `get_x()`, `x_coord()`)
  - `normalize()` returns `Option<Self>`, not `Result` or panic
- [ ] **File organization matches sibling modules**:
  - Same module split: `types.rs`, `ops.rs`, `conversions.rs`, `nalgebra.rs`, `arbitrary.rs`
  - Same section comment style: `// ============` banners
  - Same import ordering: std, external crates, crate-internal
- [ ] **Documentation format matches existing docs**:
  - Brief summary first line
  - `# Example` section with working code
  - Same level of mathematical detail as similar types
- [ ] **Test organization matches existing tests**:
  - Same proptest patterns
  - Same edge case coverage style
- [ ] **Benchmark organization matches existing benchmarks**:
  - Same naming convention: `bench_pga3_motor_transform_point`
  - Same grouping patterns

### 5. Performance & Benchmarking
- [ ] No unnecessary allocations
- [ ] SIMD used where beneficial
- [ ] No premature optimization
- [ ] Benchmarks for critical paths
- [ ] New generic operations have benchmarks in `benches/generic.rs`
- [ ] New specialized 2D/3D euclidean operations have benchmarks in `benches/specialized.rs`
- [ ] Modified operations checked for performance regression (run `cargo bench` before/after)
- [ ] Benchmark SVG reports captured in `benches/reports/` (organized by `generic/`, `euclidean/dim2/`, `euclidean/dim3/`)
- [ ] `benches/README.md` updated if significant timing changes

### 6. Dependencies
- [ ] Minimal external dependencies
- [ ] Any new dependency is justified
- [ ] Prefer std library
- [ ] New dependencies use allowed licenses (run `cargo deny check`)
- [ ] No GPL/LGPL/AGPL dependencies

### 7. Git Hygiene
- [ ] Branch created from latest `origin/main`
- [ ] Small, logical commits (docs separate from impl, modules separate)
- [ ] Conventional commit format: `type(scope): description`
- [ ] Each commit builds and passes tests independently
- [ ] PR targets main via feature branch

## Review Output Format

```markdown
## Summary
[One paragraph overview]

## Strengths
- [What's done well]

## Issues

### Critical (must fix)
- [Issue]: [location] - [explanation]

### Suggestions (consider)
- [Suggestion]: [location] - [explanation]

## Verdict
[ ] Approve
[ ] Request changes
[ ] Needs discussion
```

## Mathematical Review

Pay special attention to:
- Correct handling of metric signature
- Grade operations (selection, projection)
- Product implementations (geometric, inner, outer)
- Inverse and normalization
- Sign conventions (especially for reversion, duals)

### Algebraic Derivations

For complex operations (motor composition, sandwich products, multi-term formulas):

- [ ] **Derivation exists** - Check `derivations/src/clifford_derivations/` for a SymPy derivation
- [ ] **Formula verified** - If no derivation exists, request one be created before approving
- [ ] **Code matches derivation** - Compare Rust implementation against the symbolic derivation
- [ ] **Derivation referenced** - Code comments should reference the derivation script

**Red flags** (request derivation before approving):
- Motor composition (`M₁ * M₂`) without derivation
- Sandwich products (`M x M̃`) with more than 3-4 terms
- Manual sign corrections or "magic" coefficients
- Copy-pasted formulas without cited source
- Hardcoded Rust formulas that claim to be "derived" but aren't generated from SymPy

**Critical: Rust code must be generated from SymPy**:
- **Never accept hardcoded algebraic formulas** - we cannot trust manual algebra
- **Derivation scripts must use `sympy.printing.rust.rust_code()`** to generate Rust
- **The derivation must produce the actual Rust code** that goes into the implementation
- If code claims to be derived but the derivation script uses hardcoded strings, reject it

**Derivation quality checks**:
- Derivation uses `expand()` not `simplify()` (avoids hangs)
- Derivation has timeout protection (`@with_timeout(30)`)
- Complex derivations are broken into steps (no single function > 60s)
- Intermediate results are simplified before combining
- **Rust code is generated via `rust_code()` from `sympy.printing.rust`**, not hardcoded
