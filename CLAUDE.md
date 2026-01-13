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

#### Style Consistency (Critical)

**New code must match the style of existing code in the repository.** Before writing new code:

1. **Study existing patterns** - Look at similar modules/files to understand:
   - How types are structured (field ordering, visibility)
   - How methods are organized (constructors first, then operations, then conversions)
   - How documentation is formatted (what sections, what level of detail)
   - How tests are organized (property-based tests, edge cases)
   - How benchmarks are structured

2. **Match naming conventions exactly**:
   - Method names: `transform_point`, `transform_line` (not `apply_to_point`, `move_point`)
   - Constructor patterns: `new()`, `from_*()`, `identity()`, `origin()`
   - Accessor patterns: `x()`, `y()`, `z()` (not `get_x()`, `x_coord()`)
   - Normalization: `normalize()` returns `Option<Self>`, `normalized()` panics or returns copy

3. **Match file organization**:
   - Module structure: `types.rs`, `ops.rs`, `conversions.rs`, `nalgebra.rs`, `arbitrary.rs`
   - Import ordering: std first, then external crates, then crate-internal
   - Section comments: `// ============` banners for major sections

4. **Match documentation style**:
   - First line is a brief summary
   - `# Example` section with working code
   - `# Panics` / `# Errors` sections where applicable
   - Mathematical notation in backticks: `` `e₁ ∧ e₂` ``

5. **When in doubt, grep the codebase** - Search for similar functionality and follow that pattern exactly
- **Prefer structs with associated methods over free functions**
  - Avoid primitive obsession (e.g., using `usize` instead of a proper `Blade` type)
  - Encapsulate internal details; don't leak implementation (e.g., bitmask indices)
  - Methods on types are more discoverable and provide better IDE support
  - Example: `blade.grade()` is better than `grade_of_blade(index)`
- **Use Clifford types in APIs, not tuples**
  - Return Clifford types (e.g., `euclidean::dim3::Vector`) instead of tuples `(T, T, T)`
  - Accept Clifford types as parameters instead of tuples
  - Examples:
    - `line.direction() -> Vector` (not `-> (T, T, T)`)
    - `line.normal() -> Vector` (not `-> (T, T)`)
    - `Motor::from_axis_angle(&Vector, T)` (not `(T, T, T), T`)
  - This provides better IDE support, type safety, and method chaining
- **Avoid fully-qualified syntax** - Prefer `Type::method()` over `<Type as Trait>::method()` when possible. Add helper methods or type aliases to make simpler syntax work.
- **Don't expose foreign traits in public API** - When our public API depends on a foreign trait (e.g., `typenum::Unsigned`), either:
  - Re-export the trait in our prelude so users don't need to import the dependency directly, or
  - Add helper methods that encapsulate the foreign trait usage (preferred when feasible)

#### No Warning Suppression

**Never use `#[allow(...)]` attributes to suppress warnings.** Always fix the root cause:

- **`#[allow(dead_code)]`** - If code is unused, either delete it or properly expose it through the public API (e.g., add re-exports in `mod.rs`)
- **`#[allow(unused_imports)]`** - Remove the unused import
- **`#[allow(unused_variables)]`** - Use `_` prefix for intentionally unused variables, or remove them
- **`#[allow(clippy::*)]`** - Fix the code to satisfy clippy, or if truly a false positive, document why with a comment

**Underscore prefix on definitions is NOT a workaround:**
- **Don't use `_field_name` on struct fields** to suppress dead_code warnings - this is equivalent to `#[allow(dead_code)]`. Either use the field or delete it.
- **Don't use `_VariantName` on enum variants** - either construct the variant somewhere or remove it.
- **Exception: Function parameters** - `_param` on unused function parameters is acceptable when the signature is fixed (trait implementations, API compatibility).
- **Exception: PhantomData** - `_marker: PhantomData<T>` is legitimate when you need the type parameter for type system reasons but don't store data of that type.

**Why this matters:**
- Warnings exist to catch real problems - suppressing them hides bugs
- Unused code becomes maintenance burden and confuses readers
- CI enforces `warnings = "deny"` - suppressions work around safety checks

**The right approach:**
```rust
// BAD: Suppressing dead code warning
#[allow(dead_code)]
pub type UnitVector<T> = Unit<Vector<T>>;

// ALSO BAD: Using underscore prefix to suppress warning
struct Options {
    _unused_field: bool,  // Don't do this - delete it instead
}

// GOOD: Properly expose through module's public API
// In mod.rs:
pub use generated::types::UnitVector;

// GOOD: PhantomData when you need the type parameter
struct Container<T> {
    data: Vec<u8>,
    _marker: PhantomData<T>,  // OK - needed for type variance
}

// GOOD: Unused parameter in trait implementation
impl Visitor for MyVisitor {
    fn visit(&self, _node: &Node) {  // OK - signature is fixed by trait
        // ...
    }
}
```

#### Field Visibility

**All specialized types use private fields with public accessors** for consistency:

```rust
pub struct Vector<T: Float> {
    x: T,  // private
    y: T,  // private
    z: T,  // private
}

impl<T: Float> Vector<T> {
    // Accessors
    pub fn x(&self) -> T { self.x }
    pub fn y(&self) -> T { self.y }
    pub fn z(&self) -> T { self.z }

    // Constructor
    pub fn new(x: T, y: T, z: T) -> Self { ... }
}
```

**For types with constraints** (Motor, Line), add validated constructors:
```rust
impl<T: Float> Motor<T> {
    pub fn new_unchecked(...) -> Self { ... }               // raw construction
    pub fn try_from_components(...) -> Result<Self, E> { }  // validated
    pub fn from_rotation_z(angle: T) -> Self { ... }        // factory (safe)
}
```

**Constructor naming conventions:**
- `new()` - Standard constructor for types without constraints.
- `new_unchecked()` - Raw construction for constrained types. For internal operations and AD.
- `try_from_components()` - Validates constraint, returns `Result`.
- `from_*()` - Factory methods that guarantee validity by construction.

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

#### Scalar and Float Types

Our `Float` trait extends `num_traits::Float`, which provides standard floating-point operations.

**Accessing constants via methods (not associated constants):**
```rust
// CORRECT: Use methods from num_traits::Float
let zero: T = T::zero();
let one: T = T::one();
let eps: T = T::epsilon();

// WRONG: Don't use associated constants (these don't exist)
// let zero: T = T::ZERO;  // Compile error!
```

**Using approx crate for comparisons:**
```rust
use approx::relative_eq;

// CORRECT: Use relative_eq! with both epsilon and max_relative
// - epsilon: for near-zero comparisons (absolute tolerance)
// - max_relative: for larger values (relative tolerance)
assert!(relative_eq!(a, b, epsilon = 1e-10, max_relative = 1e-10));
prop_assert!(relative_eq!(result, expected, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));

// WRONG: Using only max_relative fails for values near zero
// assert!(relative_eq!(a, 0.0, max_relative = 1e-10));  // Will fail unless a == 0

// WRONG: Don't use custom comparison methods
// assert!(a.approx_eq(&b, 1e-10));  // This method doesn't exist
```

**Custom constants on our Float trait:**
```rust
// Additional constants provided by clifford::scalar::Float
let two: T = T::TWO;     // 2.0
let pi: T = T::PI;       // π

// Conversions
let from_int: T = T::from_i8(1);
let from_idx: T = T::from_usize(42);
let from_f64: T = T::from_f64(3.14);
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
- **Use generic wrappers for normalized types**: The `crate::wrappers` module provides generic wrappers with `Arbitrary` implementations:
  - `Unit<T>` - Euclidean norm = 1 (requires `Normed` trait)
  - `Bulk<T>` - PGA bulk norm = 1 (requires `DegenerateNormed` trait)
  - `Ideal<T>` - PGA weight norm = 1 (requires `DegenerateNormed` trait)

  **Do NOT create bespoke wrapper types** like `UnitMotor<T>` in arbitrary modules. Use the generic wrappers instead:
  ```rust
  // GOOD: Use generic wrappers from crate::wrappers
  use clifford::specialized::projective::dim3::{BulkMotor, Motor, Point};

  proptest! {
      #[test]
      fn motor_preserves_distance(motor in any::<BulkMotor<f64>>(), p in any::<Point<f64>>()) {
          // BulkMotor<T> = Bulk<Motor<T>> with Arbitrary impl
      }
  }

  // BAD: Bespoke wrapper in arbitrary module
  pub struct UnitMotor<T>(pub Motor<T>);  // Don't do this!
  impl<T> Arbitrary for UnitMotor<T> { ... }  // Don't do this!
  ```
- **proptest-support feature**: External consumers enable `proptest-support` feature to access arbitrary modules
- **Arbitrary wrapper ergonomics**: All wrapper types implement `Deref`, `AsRef`, `From`, and `into_inner()` for easy access to the inner value
- **Use generated type aliases, not ad-hoc strategy functions**: The codegen generates type aliases like `UnitizedPoint<T>`, `BulkMotor<T>`, etc. that have `Arbitrary` impls. Use these instead of writing ad-hoc strategy functions:
  ```rust
  // GOOD: Use generated type aliases with Arbitrary impls
  use clifford::specialized::projective::dim3::UnitizedPoint;

  proptest! {
      #[test]
      fn test_finite_point(p in any::<UnitizedPoint<f64>>()) {
          // UnitizedPoint guarantees weight_norm = 1 (finite point)
          // Access inner Point via Deref: &*p or p.as_inner()
          let na_p: na::Point3<f64> = (*p).try_into().unwrap();
      }
  }

  // BAD: Ad-hoc strategy function
  fn finite_point_strategy() -> impl Strategy<Value = Point<f64>> {
      (-100.0..100.0, -100.0..100.0, -100.0..100.0)
          .prop_map(|(x, y, z)| Point::from_cartesian(x, y, z))
  }
  ```
  **Available type aliases** (generated by codegen):
  - `UnitizedPoint<T>` - Finite point (weight = 1)
  - `UnitizedPlane<T>` - Finite plane
  - `UnitizedLine<T>` - Finite line
  - `BulkMotor<T>` - Normalized motor (bulk_norm = 1)
  - `BulkFlector<T>` - Normalized flector
  - `IdealPoint<T>` - Point at infinity (weight ≈ 0)
  - `IdealLine<T>` - Line at infinity
- **Use `approx` crate with `relative_eq!`**: Never hand-roll floating-point comparisons. Use `relative_eq!` with both `epsilon` and `max_relative` parameters:
  ```rust
  use approx::relative_eq;

  // Use BOTH epsilon (for near-zero) and max_relative (for larger values)
  assert!(relative_eq!(a, b, epsilon = 1e-10, max_relative = 1e-10));
  ```
- **Avoid magic values (numbers and paths)**: Don't use literal numbers or paths scattered throughout code. Bind them to named constants or variables:
  ```rust
  // GOOD: Named constants for numbers
  const RELATIVE_EQ_EPS: f64 = 1e-10;
  const MAX_ITERATIONS: usize = 100;
  const DEFAULT_TOLERANCE: f64 = 1e-6;

  assert!(relative_eq!(a, b, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));

  // BAD: Magic numbers
  assert!(relative_eq!(a, b, epsilon = 1e-10, max_relative = 1e-10));  // What is 1e-10?
  for _ in 0..100 { ... }  // Why 100?
  ```
  ```rust
  // GOOD: Named constants for paths (especially in tests)
  const ALGEBRAS_DIR: &str = "algebras";
  const EUCLIDEAN3_TOML: &str = "algebras/euclidean3.toml";

  let spec = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/", ALGEBRAS_DIR, "/euclidean3.toml"));

  // BAD: Magic paths scattered in code
  let spec = include_str!("../../../../algebras/euclidean3.toml");  // Fragile, hard to update
  ```
  **Why this matters:**
  - Named values are self-documenting
  - Changing a value requires updating only one place
  - Reduces bugs from inconsistent values across files
  - Makes code review easier (reviewers can verify the constant definition once)

- **Use `RELATIVE_EQ_EPS` constant**: Use the standard constant `RELATIVE_EQ_EPS` defined in `src/lib.rs::test_utils`:
  ```rust
  use crate::test_utils::RELATIVE_EQ_EPS;
  use approx::relative_eq;

  // Good: use the standard constant with both parameters
  assert!(relative_eq!(a.norm(), 1.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));

  // Avoid: magic numbers
  assert!(relative_eq!(a.norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));

  // WRONG: using only max_relative fails for values near zero
  // assert!(relative_eq!(a, 0.0, max_relative = RELATIVE_EQ_EPS));
  ```
  For integration tests (`tests/` directory), define the constant locally since `test_utils` is `pub(crate)`.
- **Always use both `epsilon` AND `max_relative`**: The `epsilon` parameter provides absolute tolerance for values near zero, while `max_relative` provides relative tolerance for larger values. Using both ensures robust comparisons across all magnitudes.
- **Use `prop_assert!` in proptest blocks**: Inside `proptest!` blocks, always use `prop_assert!` instead of `assert!`. This provides better error reporting with counterexamples:
  ```rust
  proptest! {
      #[test]
      fn rotor_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
          let rotated = r.rotate(v);
          // Good: prop_assert! with relative_eq! using both epsilon and max_relative
          prop_assert!(relative_eq!(v.norm(), rotated.norm(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));

          // Avoid: assert! loses proptest's shrinking and reporting benefits
          // Avoid: magic numbers like 1e-9
          // Avoid: using only max_relative (fails for values near zero)
      }
  }
  ```
- All types implement `AbsDiffEq`, `RelativeEq`, and `UlpsEq` traits for f32 and f64 variants
- The `approx` traits are re-exported from `clifford::prelude`
- Unit tests with specific examples are acceptable only when property-based testing is not feasible
- Doc tests for examples
- All tests run automatically via GitHub Actions CI on every push and PR
- **Symbolica test naming convention**: Tests in `clifford-codegen` that use Symbolica for symbolic computation **must** be prefixed with `symbolica_`. This enables nextest to run them serially (avoiding global state conflicts).
  ```rust
  // GOOD: Test uses Symbolica, has symbolica_ prefix
  #[test]
  fn symbolica_generates_geometric_product() {
      let algebra = Algebra::euclidean(3);
      // ... uses Symbolica for symbolic computation
  }

  // BAD: Test uses Symbolica but lacks prefix (will cause flaky failures)
  #[test]
  fn generates_geometric_product() {  // Missing symbolica_ prefix!
      let algebra = Algebra::euclidean(3);
      // ...
  }

  // GOOD: Test doesn't use Symbolica, no prefix needed
  #[test]
  fn parse_spec_handles_empty_types() {
      // Pure parsing test, no Symbolica
  }
  ```
  **When to add the prefix:**
  - Any test that creates `Algebra`, `ProductTable`, or `SymbolicProduct`
  - Any test that calls `compute_terms()`, `generate_products_file()`, or similar
  - Any test in `symbolic/` modules

  See `.config/nextest.toml` for the test group configuration.

### 8. Code Review
- PRs are reviewed by Greptile (AI-powered review)
- Address or acknowledge Greptile feedback before merging
- **IMPORTANT: After pushing commits that address Greptile feedback**, always comment `@greptileai review` on the PR to request a re-review. This ensures the fixes are verified before merging.
- **No `todo!()` macros** - Code containing `todo!()` must not pass review. The macro indicates incomplete implementation and will panic at runtime. Either implement the functionality or remove the code path entirely.

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
| **devops** | CI/CD and infrastructure tasks. Use for GitHub Actions, branch protection, or build automation. |
| **release** | Release and publishing tasks. Use for version bumps, CHANGELOG updates, crates.io publishing, and creating GitHub releases. |

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
- **Publishing a release?** → `release`

### 10. Code Generation for Algebraic Types

**CRITICAL: Do NOT manually derive or reason about algebraic formulas.**

Geometric algebra formulas are complex and error-prone. Signs, orderings, metric contractions, and grade projections are easy to get wrong. This project uses automated tools to ensure correctness:

1. **clifford-codegen** - Generates Rust code for products, transformations, and constraints
2. **Symbolica** - The codegen tool uses Symbolica for symbolic computation

**What this means for you:**
- **Never manually derive** geometric product formulas, sandwich products, or transformations
- **Never manually compute** signs, metric contractions, or blade orderings
- **Never try to "fix"** algebraic formulas by adjusting signs or coefficients
- **Always use** the codegen tool to generate correct implementations
- **If codegen produces wrong results**, fix the codegen tool, not the generated code

**When you encounter algebraic issues:**
1. Check if the product/operation is already generated - use the generated version
2. If not generated, add it to the codegen tool (see `crates/clifford-codegen/`)
3. If codegen seems wrong, debug and fix the codegen tool itself
4. Regenerate all algebras after fixing codegen

**The codegen tool uses Symbolica** for symbolic algebra, which:
- Computes geometric products symbolically with correct signs
- Handles metric signatures (positive, negative, zero basis vectors)
- Solves constraints (geometric constraint, Plücker condition, etc.)
- Generates optimized Rust code with verified formulas

#### Genericity Principle: No Shortcuts

**CRITICAL: The code generator must be completely generic. No algebra-specific shortcuts.**

When implementing codegen features, always ask: **"How does this generalize to arbitrary signatures?"**

**What this means:**
- **Never check algebra names** - No `if name.starts_with("euclidean")` or similar string matching
- **Derive everything from signature** - The tuple `(p, q, r)` (positive, negative, degenerate basis count) determines all behavior
- **No hardcoded formulas for specific algebras** - Constraint expressions, norm formulas, and products must be derived symbolically
- **Test with multiple signatures** - If it works for Euclidean but not PGA, fix the root cause

**Examples of violations to avoid:**
```rust
// BAD: String-based algebra detection
if name.starts_with("pga") { use_pga_formula() }

// BAD: Hardcoded constraint for one algebra
let geometric_residual = s * e0123 + e23 * e01 + ...;  // Motor-specific

// BAD: Separate code paths for "PGA" vs "Euclidean"
if is_pga { generate_pga_wrappers() } else { generate_euclidean_wrappers() }
```

**Correct approach:**
```rust
// GOOD: Derive from signature tuple
match (sig.p, sig.q, sig.r) {
    (n, 0, 0) => /* Euclidean behavior */,
    (n, 0, 1) => /* Degenerate metric behavior */,
    _ => /* Generic fallback */,
}

// GOOD: Derive constraints symbolically from u * rev(u)
let constraint = deriver.derive_geometric_constraint(ty);

// GOOD: Partition fields by degenerate basis involvement
let bulk_fields = fields.filter(|f| !contains_degenerate_index(f));
```

**Root cause fixes over workarounds:**
- If a formula doesn't work for one algebra, the fix should make it work for ALL algebras
- If you're tempted to add a special case, step back and find the general solution
- Workarounds accumulate tech debt; generic solutions scale to future algebras (CGA, Minkowski, etc.)

Complex algebraic operations in geometric algebra (products, motor compositions, transformations) are handled by the **clifford-codegen** tool. This ensures correctness through:

- **Automatic constraint solving**: Constraints like geometric constraints and Plücker conditions are solved symbolically
- **Type-safe code generation**: Products, conversions, and traits are generated with correct formulas
- **Property-based testing**: Generated verification tests compare specialized types against generic Multivector

#### Algebra TOML Files

**All algebra specifications live in the top-level `algebras/` directory.** This is the single source of truth for algebra definitions.

```
algebras/
  euclidean2.toml    # 2D Euclidean GA
  euclidean3.toml    # 3D Euclidean GA
  projective3.toml   # 3D Projective GA (PGA)
```

**Important guidelines:**
- **Never duplicate** algebra TOMLs elsewhere in the repo
- **Always regenerate** to the correct specialized module path:
  - `euclidean::dim2` → `src/specialized/euclidean/dim2/generated/`
  - `euclidean::dim3` → `src/specialized/euclidean/dim3/generated/`
  - `projective::dim3` → `src/specialized/projective/dim3/generated/`
- The `module_path` field in the TOML must match the target location

#### Usage

```bash
# Generate code for an algebra (specify output path)
cargo run --package clifford-codegen -- generate algebras/euclidean3.toml -o src/specialized/euclidean/dim3/generated --force

# Discover valid entities for a signature
cargo run --package clifford-codegen -- discover 3 0 1

# List blades for an algebra
cargo run --package clifford-codegen -- blades algebras/projective3.toml

# Regenerate all algebras after codegen changes
cargo run --package clifford-codegen -- generate algebras/euclidean2.toml -o src/specialized/euclidean/dim2/generated --force
cargo run --package clifford-codegen -- generate algebras/euclidean3.toml -o src/specialized/euclidean/dim3/generated --force
cargo run --package clifford-codegen -- generate algebras/projective3.toml -o src/specialized/projective/dim3/generated --force
```

#### Generated Code Quality

**The code generator must produce clean, warning-free Rust code:**

- **No clippy warnings** - Generated code must pass `cargo clippy` without warnings
- **No warning suppressions** - Never add `#[allow(dead_code)]`, `#[allow(unused)]`, or similar attributes to hide warnings
- **If generated code produces warnings**, fix the generator to produce correct code instead

**When you see warnings in generated code:**
1. Identify the root cause in the codegen tool
2. Fix the generator to avoid producing that pattern
3. Regenerate all algebras
4. Verify no warnings remain

#### When to Use Codegen

- **Adding new algebras**: Define types and constraints in TOML, generate Rust code
- **Complex products**: Motor composition, sandwich products, etc. are generated correctly
- **Constraint handling**: Geometric constraints, Plücker conditions, unit norm - all solved automatically

**Never hardcode algebraic formulas.** Use the codegen tool to generate correct implementations.

#### Constraints vs. Product Outputs

**Important distinction:** Type constraints (geometric constraint, Plucker condition, etc.) apply to *normalized, valid instances* of a type, NOT to algebraic product results.

Generated products use `new_unchecked()` for constrained types because:
1. **Product outputs are algebraically correct** - The formulas are derived from the multiplication table
2. **Constraint solving would modify results** - e.g., computing `e0123` via geometric constraint when the product naturally produces `e0123 = 0`
3. **Constraints apply to normalized instances** - A `Motor` from a product is "motor-shaped", not necessarily constraint-satisfying

**When constraints matter:**
- `Motor::normalize()` - explicit re-normalization when drift accumulates
- Factory methods like `Motor::from_translation()` - guarantee validity by construction
- `Motor::try_from_components()` - validate user inputs

**When constraints don't apply:**
- Product outputs (`geometric_motor_motor`, etc.)
- Intermediate computations
- Any operation where the algebraic result is what's needed

See [PRD-17.8](docs/prd/prd-17.8-product-normalization.md) for detailed investigation.

#### Using Generated Products in Extensions

When implementing domain-specific methods in extension files (`extensions.rs`), **always use generated products** instead of manual formulas:

```rust
// CORRECT: Use generated product from generated/products.rs
use super::generated::products;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Use generated exterior product
        products::exterior_point_point(self, other)
    }
}

// WRONG: Manual formula (error-prone, not verified)
impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Don't do this - prone to sign/ordering errors
        Line::new_unchecked(
            self.e1() * other.e2() - self.e2() * other.e1(),  // might be wrong!
            // ...
        )
    }
}
```

**Acceptable exceptions** (when manual formulas are okay):
1. The product combination isn't generated yet (file an issue to add it!)
2. Performance-critical code where the generated form is suboptimal
3. Geometric shortcuts (e.g., `distance()` using a specialized formula)

Even for exceptions, add a comment citing the mathematical source of the formula.

**Red flags in code review**:
- Manual product formulas in extension files
- Sign corrections or coefficient adjustments
- Comments like "// manual implementation because codegen doesn't..."

These indicate missing codegen features that should be added via PRD-17.

#### Regenerating Algebras After Codegen Changes

**Whenever you modify the code generator (`clifford-codegen`), you MUST regenerate all algebras.** This ensures the generated code stays in sync with codegen changes.

**Step 1: Update TOML specs if needed**

Some codegen changes require updating the algebra TOML files:
- New configuration fields (e.g., adding `interior = true` to products)
- Renamed fields (e.g., `outer` → `exterior`)
- New constraint formats
- Schema changes

```bash
# Example: enabling a new product type in all algebras
for toml in algebras/*.toml; do
    # Add new configuration (manually or via script)
done
```

**Step 2: Regenerate all algebras**

```bash
# After any codegen modification, regenerate all algebras
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done

# Then verify everything still works
cargo fmt
cargo clippy
cargo test
cargo doc --no-deps
```

**When TOML updates are needed:**
- Adding new product types (must enable in TOML)
- Adding new type-level options
- Changing field naming conventions
- Adding new constraint syntax

**When only regeneration is needed:**
- Fixing bugs in existing generation
- Improving generated code quality
- Fixing documentation generation
- Optimizing generated expressions

**Commit pattern for codegen changes:**
1. First commit: The codegen fix/feature
2. Second commit: TOML updates (if needed) + regenerated algebras
3. Or: Single commit with all if they're tightly coupled

This prevents situations where codegen is fixed but generated code is stale.

## Development Commands

```bash
cargo build           # Build the library (default features include serde, proptest-support, nalgebra-0_33)
cargo nextest run     # Run all tests (recommended - handles Symbolica correctly)
cargo test            # Run all tests (fallback - uses mutex for Symbolica tests)
cargo doc --open      # Generate and view documentation
cargo bench           # Run benchmarks
cargo clippy          # Run linter
cargo fmt             # Format code
cargo deny check      # Check licenses and advisories
```

### Why `cargo nextest`?

The `clifford-codegen` crate uses [Symbolica](https://symbolica.io) for symbolic algebra, which has global state that conflicts when tests run in parallel. We use `cargo nextest` with test groups to:

- Run Symbolica tests (prefixed with `symbolica_`) serially
- Run all other tests in parallel for speed

The configuration is in `.config/nextest.toml`. Tests also have a mutex fallback for `cargo test` users, but `nextest` is preferred for reliability.

Install nextest with:
```bash
cargo install cargo-nextest
```

## Verification Workflow

**Before every commit**, run these commands to ensure CI will pass:

```bash
cargo fmt             # Format code (CI checks this!)
cargo clippy          # Lint check
cargo doc --no-deps   # Documentation build (CI checks this!)
cargo nextest run     # Run all tests (recommended)
cargo deny check      # License and security audit (CI checks this!)
```

CI will reject PRs that fail any of these checks. Always run `cargo fmt` before committing.

**Note on features**: Use default features for local development. CI handles testing the full feature matrix:
- Default features: `serde`, `proptest-support`, `nalgebra-0_33`
- CI tests each nalgebra version separately (`nalgebra-0_32`, `nalgebra-0_33`, `nalgebra-0_34`)
- The nalgebra features are **mutually exclusive** - do not use `--all-features`

### License Policy

Dependencies must use permissive licenses. The `deny.toml` configuration:
- **Allowed**: MIT, Apache-2.0, BSD-2-Clause, BSD-3-Clause, ISC, Zlib, CC0-1.0, Unlicense, MPL-2.0, Unicode-3.0
- **Denied**: GPL, LGPL, AGPL (copyleft licenses)

Run `cargo deny check` to verify license compliance before adding new dependencies.

## Architecture Notes

(To be expanded as the library develops)

## Status

### Completed
- [x] Initial project setup with Cargo.toml (edition 2024)
- [x] GitHub Actions CI workflow (check, test, fmt, clippy, cargo-deny)
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

## Resources

### Projective Geometric Algebra (PGA)
- [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page) - Comprehensive reference for 3D PGA formulas, motor transformations, and the geometric antiproduct
- [Look, Ma, No Matrices!](https://enkimute.github.io/LookMaNoMatrices/) - Optimized PGA formulas with shader implementations
