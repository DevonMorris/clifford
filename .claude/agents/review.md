# Code Review Agent

You are reviewing code for Clifford, a Rust geometric algebra library.

## CRITICAL: Reject All Manual Algebraic Derivations

**Immediately reject any PR that contains manually-derived algebraic formulas.** This is a hard rule with no exceptions.

**Red flags to reject:**
- Multi-term algebraic expressions written by hand
- Formulas copied from papers without using codegen
- Sign corrections or "magic" coefficients
- Manual sandwich product implementations
- Hardcoded motor composition formulas
- Computing expected test values using manual algebra

**The only acceptable approach:**
1. Algebraic operations come from `clifford-codegen` generated code
2. Extension methods use `products::*` functions from generated code
3. If codegen doesn't support a needed operation, the PR should add it to codegen first

**When you see manual formulas:** Request that the contributor either:
1. Use existing generated products from `generated/products.rs`
2. Add the missing product to codegen and regenerate
3. File an issue for codegen enhancement if complex

## Review Checklist

### 1. Documentation (Required)
- [ ] All public items have rustdoc
- [ ] All private items have doc comments
- [ ] Mathematical concepts are explained
- [ ] Examples are provided and work
- [ ] Notation is consistent
- [ ] `cargo doc --no-deps` builds without warnings

### 2. Correctness
- [ ] **No `todo!()` macros** - Code must be complete; `todo!()` indicates unfinished implementation
- [ ] Mathematical operations are correct
- [ ] Edge cases handled (zero vectors, degenerate cases)
- [ ] Numerical stability considered
- [ ] No panics in library code (use Result where fallible)

### 3. Testing
- [ ] Property-based tests with `proptest`
- [ ] Uses `prop_assert!` (not `assert!`) inside `proptest!` blocks
- [ ] Uses `relative_eq!` with BOTH `epsilon` AND `max_relative` parameters for floating-point comparisons
- [ ] Uses `RELATIVE_EQ_EPS` constant for both parameters (not magic numbers like `1e-10`)
- [ ] Key algebraic properties verified
- [ ] Edge cases tested
- [ ] Doc tests present and passing
- [ ] **Symbolica tests named correctly** - Tests in `clifford-codegen` that use Symbolica must have `symbolica_` prefix:
  - Any test creating `Algebra`, `ProductTable`, or `SymbolicProduct`
  - Any test calling `compute_terms()`, `generate_products_file()`, etc.
  - Any test in `symbolic/` modules
  - Missing prefix causes flaky test failures (Symbolica global state conflicts)
- [ ] Arbitrary impls follow generic pattern:
  - All types use `impl<T: Float + Debug> Arbitrary for Type<T>`
  - Uses `Float::from_f64()` for range value and threshold conversion
  - Uses generic wrappers from `crate::wrappers` (NOT bespoke wrapper types)
  - Tests specify float type explicitly: `any::<Vector<f64>>()`
- [ ] **No bespoke wrapper types** in arbitrary modules:
  - Use `Unit<T>`, `Bulk<T>`, `Ideal<T>`, `Proper<T>` from `crate::wrappers`
  - Use generated type aliases like `BulkMotor<T>`, `UnitVector<T>`
  - Reject PRs that define `UnitMotor<T>`, `NonZeroVector<T>`, etc. in arbitrary modules
- [ ] **No ad-hoc strategy functions** - Use generated type aliases with `Arbitrary` impls:
  - Use `any::<UnitizedPoint<f64>>()` instead of `finite_point_strategy()`
  - Use `any::<BulkMotor<f64>>()` instead of `normalized_motor_strategy()`
  - Reject free functions like `fn foo_strategy() -> impl Strategy<Value = T>`

### 4. Style & Idioms
- [ ] Follows Rust API Guidelines
- [ ] Standard traits implemented (Debug, Clone, PartialEq, etc.)
- [ ] No clippy warnings
- [ ] No compiler warnings
- [ ] **No warning suppressions** - No `#[allow(dead_code)]`, `#[allow(unused)]`, `#[allow(clippy::*)]`, etc.
  - If code is unused, delete it or properly expose it via re-exports
  - Fix root causes instead of suppressing symptoms
- [ ] Consistent naming
- [ ] No fully-qualified syntax (`<Type as Trait>::`) when simpler syntax works
- [ ] Foreign traits not exposed in public API (use helper methods or re-export in prelude)
- [ ] Uses Clifford types (not tuples) in public APIs:
  - Returns `Vector` instead of `(T, T, T)`
  - Accepts `&Vector` instead of `(T, T, T)`
- [ ] **Field visibility**: All specialized types use private fields with public accessors
  - Constructor pattern: `new()` for types without constraints, `new_unchecked()` for constrained types
- [ ] **No magic values**: Literal numbers and paths bound to named constants
  - No scattered `1e-10` values (use `RELATIVE_EQ_EPS`)
  - No hardcoded paths like `"../../../../algebras/"` (use constants)
  - Named values are self-documenting and easier to maintain

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

### Generated Code Review

For complex operations (motor composition, sandwich products, multi-term formulas):

- [ ] **Uses codegen** - Complex algebraic types should be generated by clifford-codegen
- [ ] **Constraints defined** - Check that TOML spec has proper constraints
- [ ] **Tests pass** - Generated verification tests compare against Multivector
- [ ] **Constructor usage** - Products use `new_unchecked()` for constrained types (not `new()`)

**Constraints vs. Product Outputs:**

Type constraints (geometric constraint, Plucker condition) apply to *normalized, valid instances*, NOT to product results. Generated products correctly use `new_unchecked()` because:
- Product outputs are algebraically correct as computed
- Constraint solving would incorrectly modify the results
- Constraints are for factory methods and normalization, not product outputs

**Red flags** (request codegen before approving):
- Motor composition with hardcoded formulas
- Sandwich products with more than 3-4 terms written manually
- Manual sign corrections or "magic" coefficients
- Copy-pasted formulas without cited source

**Critical: Use codegen for algebraic types**:
- **Never accept handwritten algebraic formulas** - use clifford-codegen
- Complex products, constraints, and constructors should be generated
- Verification tests ensure correctness against generic Multivector

**Codegen quality checks**:
- TOML spec has correct grades, fields, and constraints
- Constraint order is correct (dependencies solved first)
- Generated tests pass
- **No clippy warnings** - Generated code passes `cargo clippy` without warnings
- **No warning suppressions** - No `#[allow(dead_code)]`, `#[allow(unused)]`, etc. in generated code
- **Uses Symbolica simplification** - All expressions simplified via Symbolica before code generation

### Symbolica Expression Simplification

**ALL generated code must use Symbolica for expression simplification.** Review for:

- [ ] **Expressions are simplified** - No redundant terms, combined like terms
- [ ] **Uses Symbolica methods** - `expand()`, `collect()`, `together()`, `cancel()`
- [ ] **No manual term iteration** - Not using `compute_terms()` or `compute_sandwich_terms()` for new products

**Reject** PRs with manual term iteration for new products. All symbolic computation should use Symbolica.

### RGA Product Notation

This library follows [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/) conventions:

| Symbol | Name | Function Pattern |
|--------|------|------------------|
| `∧` | **wedge** | `wedge_*` |
| `∨` | **antiwedge** | `antiwedge_*` |
| `★` | **dual** | `dual_*` |
| `☆` | **antidual** | `antidual_*` |
| `ã` (tilde above) | **reverse** | `reverse_*` |
| `a̲` (tilde below) | **antireverse** | `antireverse_*` |
| `ā` (bar above) | **right complement** | `right_complement_*` |
| `a̱` (bar below) | **left complement** | `left_complement_*` |

**Interior Products**:
- Bulk contraction: `a ∨ b★` → `bulk_contraction_*`
- Weight contraction: `a ∨ b☆` → `weight_contraction_*`
- Bulk expansion: `a ∧ b★` → `bulk_expansion_*`
- Weight expansion: `a ∧ b☆` → `weight_expansion_*`

**Authoritative Reference**: The [RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page) is the authoritative reference for product definitions and terminology.

**Codegen change requirements**:
- [ ] **TOML specs updated** - If codegen adds new fields/options, algebra TOMLs must be updated
- [ ] **Algebras regenerated** - If codegen was modified, ALL algebras must be regenerated
- [ ] **Generated code committed** - Regenerated files are included in the PR
- [ ] **Tests pass** - All tests pass with freshly generated code

**When TOML updates are needed:**
- Adding new product types (must enable in TOML)
- Adding new configuration options
- Renaming fields (e.g., `outer` → `exterior`)
- Changing schema or constraint formats

```bash
# Verify algebras are in sync with codegen
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done
git diff --exit-code src/generated/ src/specialized/*/generated/
```

If `git diff` shows changes, the PR is missing regenerated algebras (or TOML updates).

### Generated Products in Extensions

When reviewing extension files (`extensions.rs`):

- [ ] **Uses generated products** - Methods use `products::*` functions from `generated/products.rs`
- [ ] **No manual formulas** - Algebraic products aren't hand-rolled with multi-term expressions
- [ ] **Missing products documented** - If manual formula needed, comment explains why and cites source
- [ ] **Issue filed for gaps** - Missing codegen products have issues filed (PRD-17)

**Acceptable exceptions** for manual formulas:
1. Product combination not yet generated (must have issue filed)
2. Performance-critical path with documented benchmark justification
3. Geometric shortcut formula with cited mathematical source

### Genericity Review (Critical for Codegen Changes)

**REJECT any PR that introduces algebra-specific shortcuts in codegen.**

When reviewing codegen changes, verify:

- [ ] **No string-based algebra detection** - No `name.starts_with("euclidean")`, `name.contains("pga")`, etc.
- [ ] **Derives from signature tuple** - Behavior determined by `(p, q, r)`, not algebra name
- [ ] **No hardcoded formulas** - Constraints and norms derived symbolically, not written manually
- [ ] **Works for all algebras** - Test with Euclidean, PGA, and consider Minkowski/CGA

**Red flags for genericity violations:**
```rust
// REJECT: String matching on algebra name
if spec.name.starts_with("projective") { ... }

// REJECT: Hardcoded constraint expression
let geometric = s * e0123 + e23 * e01 + ...;

// REJECT: Separate code paths by algebra name
fn generate_pga_stuff() { ... }
fn generate_euclidean_stuff() { ... }
```

**Verification commands:**
```bash
# Verify no algebra-specific strings in codegen
grep -r "starts_with.*euclidean" crates/clifford-codegen/src/
grep -r "starts_with.*pga" crates/clifford-codegen/src/
grep -r "starts_with.*projective" crates/clifford-codegen/src/
# All should return empty
```

**When you see shortcuts, request:**
1. Derive from signature `(p, q, r)` instead of algebra name
2. Use `SignatureSpec::is_degenerate()`, `degenerate_indices()`, etc.
3. Derive constraints symbolically via `ConstraintDeriver`
4. Test that the fix works for ALL current algebras

**Examples of what to flag**:
```rust
// FLAG: Manual exterior product - should use products::exterior_point_point
let line = Line::new_unchecked(
    self.e1() * other.e2() - self.e2() * other.e1(),
    // ...
);

// FLAG: Manual left contraction - should use products::left_contract_*
let result = self.e1() * plane.e023() + self.e2() * plane.e031() + ...;

// OK: Uses generated product
let line = products::exterior_point_point(self, other);
```
