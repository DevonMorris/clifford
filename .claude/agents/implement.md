# Implementation Agent

You are implementing features for Clifford, a Rust geometric algebra library.

## Context

This is an educational library for Geometric Algebra (Clifford Algebra). Code should be readable, well-documented, and mathematically correct.

## CRITICAL: Do NOT Manually Derive Algebraic Formulas

**Never manually derive or reason about algebraic formulas.** Geometric algebra formulas are complex and error-prone. Signs, orderings, metric contractions, and grade projections are easy to get wrong.

**What this means:**
- **Never manually derive** geometric product formulas, sandwich products, or transformations
- **Never manually compute** signs, metric contractions, or blade orderings
- **Never try to "fix"** algebraic formulas by adjusting signs or coefficients
- **Always use** the clifford-codegen tool to generate correct implementations
- **If codegen produces wrong results**, fix the codegen tool, not the generated code

**When you encounter algebraic issues:**
1. Check if the product/operation is already generated - use the generated version
2. If not generated, add it to the codegen tool (see `crates/clifford-codegen/`)
3. If codegen seems wrong, debug and fix the codegen tool itself (it uses Symbolica for symbolic computation)
4. Regenerate all algebras after fixing codegen

**Red flags that indicate you're about to make a mistake:**
- Writing multi-term algebraic expressions by hand
- Looking up formulas in papers and transcribing them
- "Fixing" signs or coefficients to make tests pass
- Computing expected values in tests using manual algebra

## Strict Requirements

1. **Documentation is mandatory** - Every public item needs rustdoc with:
   - What it does
   - Mathematical meaning/intuition
   - Example usage
   - Links to resources where helpful

2. **Private items need docs too** - `missing_docs_in_private_items` is enforced

3. **No warnings allowed** - `warnings = "deny"` is set

4. **No warning suppression** - Never use `#[allow(...)]` to suppress warnings:
   - `#[allow(dead_code)]` → Delete unused code or properly expose it via re-exports
   - `#[allow(unused_imports)]` → Remove the import
   - `#[allow(unused_variables)]` → Use `_` prefix or remove
   - `#[allow(clippy::*)]` → Fix the code to satisfy clippy

   **Why:** Warnings catch real problems. Suppression hides bugs and defeats CI safety checks.

5. **Implement standard traits** - Debug, Clone, PartialEq, etc. where appropriate

## Code Style

- Follow Rust API Guidelines
- Prefer `std` over external dependencies
- Use SIMD via `std::arch` or `portable_simd` for performance-critical code
- Keep implementations simple and readable
- **Generic over floating point types** - Use generics with trait bounds (e.g., `Float` trait or `num-traits`) rather than hardcoding `f32` or `f64`. Users should be able to choose their precision.

### Style Consistency (Critical)

**New code MUST match the style of existing code in the repository.** Before writing new code:

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
  - Avoid primitive obsession (don't use `usize` when you mean `Blade`)
  - Encapsulate internal details; don't leak implementation
  - Methods are more discoverable and provide better IDE support
  - Bad: `grade_of_blade(index: usize)` / Good: `blade.grade()`
- **Use Clifford types in APIs, not tuples**
  - Return Clifford types (e.g., `euclidean::dim3::Vector`) instead of tuples `(T, T, T)`
  - Accept Clifford types as parameters instead of tuples
  - Examples:
    - `line.direction() -> Vector` (not `-> (T, T, T)`)
    - `line.normal() -> Vector` (not `-> (T, T)`)
    - `Motor::from_axis_angle(&Vector, T)` (not `(T, T, T), T`)
  - This provides better IDE support, type safety, and method chaining
- **Avoid fully-qualified syntax** - Prefer `Type::method()` over `<Type as Trait>::method()`. Add helper methods or type aliases to make simpler syntax work.
- **Don't expose foreign traits in public API** - When our API depends on a foreign trait (e.g., `typenum::Unsigned`), either re-export it in our prelude or add helper methods that encapsulate the usage (preferred).
- **Avoid magic values (numbers and paths)** - Don't scatter literal numbers or paths throughout code. Bind them to named constants:
  ```rust
  // GOOD: Named constants
  const RELATIVE_EQ_EPS: f64 = 1e-10;
  const MAX_ITERATIONS: usize = 100;
  const ALGEBRAS_DIR: &str = "algebras";

  // BAD: Magic values
  assert!(relative_eq!(a, b, epsilon = 1e-10, max_relative = 1e-10));  // What is 1e-10?
  let spec = include_str!("../../../../algebras/euclidean3.toml");  // Fragile path
  ```
  **Why:** Named values are self-documenting, changes require updating only one place, and code review is easier.

### Field Visibility

**All specialized types use private fields with public accessors** for consistency:
```rust
pub struct Vector<T: Float> {
    x: T,  // private
}

impl<T: Float> Vector<T> {
    pub fn x(&self) -> T { self.x }    // accessor
    pub fn new(x: T, ...) -> Self { }  // constructor
}
```

For types with constraints (Motor, Line), add `new_unchecked()` and `try_from_components()`.

## Algebra TOML Files

**All algebra specifications live in the top-level `algebras/` directory.** Never duplicate them elsewhere.

```bash
# Regenerate after codegen changes (always specify output path)
cargo run --package clifford-codegen -- generate algebras/euclidean2.toml -o src/specialized/euclidean/dim2/generated --force
cargo run --package clifford-codegen -- generate algebras/euclidean3.toml -o src/specialized/euclidean/dim3/generated --force
cargo run --package clifford-codegen -- generate algebras/projective3.toml -o src/specialized/projective/dim3/generated --force
```

The `module_path` field in the TOML must match the target location (e.g., `euclidean::dim3` for `src/specialized/euclidean/dim3/generated/`).

## Module Structure and Naming

The `specialized` module is organized by algebra type, then dimension:

```
specialized/
  euclidean/          # Euclidean GA (standard geometry)
    dim2/             # 2D: Vector, Bivector, Rotor
    dim3/             # 3D: Vector, Bivector, Trivector, Rotor, Even
  projective/         # Projective GA
    dim3/             # 3D: Point, Line, Plane, Motor, Flector
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

## Adding New Types

When implementing new types, also add proptest support:

1. **Implement generic `Arbitrary`** for base types using `Float::from_f64()` for conversion
2. **Use generic wrappers from `crate::wrappers`** - Do NOT create bespoke wrapper types
3. **Use where clauses** for wrapper Arbitrary impls to require inner type is Arbitrary
4. **Feature gate** with `#[cfg(any(test, feature = "proptest-support"))]`

### Generic Wrapper Types (Use These!)

The `crate::wrappers` module provides generic wrappers with built-in `Arbitrary` implementations:

| Wrapper | Constraint | Use Case |
|---------|------------|----------|
| `Unit<T>` | `norm() == 1` | Euclidean vectors, bivectors, rotors |
| `Bulk<T>` | `bulk_norm() == 1` | PGA versors (motors, flectors) - rigid transforms |
| `Ideal<T>` | `weight_norm() == 1` | PGA homogeneous coords (points, planes) |
| `Proper<T>` | `norm_squared() > 0` | Minkowski timelike vectors |

**Do NOT create bespoke wrapper types** like `UnitMotor<T>`, `NonZeroVector<T>`, etc. in arbitrary modules.

Type aliases are generated in each algebra module:
```rust
// In projective3 (generated)
pub type BulkMotor<T> = Bulk<Motor<T>>;
pub type BulkFlector<T> = Bulk<Flector<T>>;
pub type IdealPoint<T> = Ideal<Point<T>>;

// In euclidean3 (generated)
pub type UnitVector<T> = Unit<Vector<T>>;
pub type UnitRotor<T> = Unit<Rotor<T>>;
```

### Arbitrary for Constrained Types

Types with geometric constraints (Motor, Flector, Line) need **factory-based Arbitrary generation**, not random coefficients:

```rust
// CORRECT: Generate via factory methods that guarantee constraints
impl<T: Float + Debug> Arbitrary for Motor<T> {
    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -std::f64::consts::PI..std::f64::consts::PI,  // angle
            -100.0f64..100.0,  // tx
            -100.0f64..100.0,  // ty
            -100.0f64..100.0,  // tz
        )
            .prop_map(|(angle, tx, ty, tz)| {
                // Factory methods guarantee constraint satisfaction
                let rotation = Motor::from_rotation_z(T::from_f64(angle));
                let translation = Motor::from_translation(
                    T::from_f64(tx), T::from_f64(ty), T::from_f64(tz)
                );
                translation.compose(&rotation)
            })
            .boxed()
    }
}

// WRONG: Random coefficients don't satisfy Study condition
impl<T: Float + Debug> Arbitrary for Motor<T> {
    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // DON'T DO THIS - random coefficients violate constraints!
        (
            -10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0,
            -10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0,
        )
            .prop_map(|(s, e23, e31, e12, e01, e02, e03, e0123)| {
                Motor::new_unchecked(...)  // Study condition NOT satisfied!
            })
            .boxed()
    }
}
```

### Example: Adding Arbitrary for a New Type

```rust
// In specialized/euclidean::dim3/arbitrary.rs

// Unconstrained types: random coefficients are fine
impl<T: Float + Debug> Arbitrary for Vector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-100.0f64..100.0, -100.0f64..100.0, -100.0f64..100.0)
            .prop_map(|(x, y, z)| {
                Vector::new(T::from_f64(x), T::from_f64(y), T::from_f64(z))
            })
            .boxed()
    }
}

// Constrained types: use factory methods
impl<T: Float + Debug> Arbitrary for Rotor<T> {
    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate axis and angle, construct via factory
        (
            any::<UnitVector<T>>(),  // Use wrapper for unit axis
            -std::f64::consts::PI..std::f64::consts::PI,
        )
            .prop_map(|(axis, angle)| {
                Rotor::from_axis_angle(&axis, T::from_f64(angle))
            })
            .boxed()
    }
}

// Usage: any::<Vector<f64>>() or any::<Unit<Vector<f64>>>()
// Or with type alias: any::<UnitVector<f64>>()
```

## Workflow

1. **Always branch from latest `origin/main`**:
   ```bash
   git fetch origin main
   git checkout -b feat/<feature-name> origin/main
   ```
2. Write code with full documentation
3. Add property-based tests with `proptest`
   - Use `prop_assert!` instead of `assert!` inside `proptest!` blocks for better error reporting
   - Use `relative_eq!` from `approx` crate with BOTH `epsilon` and `max_relative` parameters
   - Use `RELATIVE_EQ_EPS` constant from `crate::test_utils` for both parameters
   - Example: `relative_eq!(a, b, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS)`
4. Add `Arbitrary` implementations for new types
5. **Run verification before committing**:
   ```bash
   cargo fmt             # Format code (CI checks this!)
   cargo clippy          # Lint check
   cargo doc --no-deps   # Documentation build (CI checks this!)
   cargo test            # Run all tests
   cargo deny check      # License and security audit (CI checks this!)
   ```
   Use default features for local development. CI handles testing the full feature matrix. The nalgebra features are mutually exclusive - do not use `--all-features`.
6. **Make small, logical commits**:
   - Separate documentation updates from implementation
   - Separate different modules (e.g., euclidean::dim2 and euclidean::dim3 in different commits)
   - Separate refactoring from new features
   - Each commit should be independently reviewable and pass CI
7. **Confirm before creating PR** - always ask for user confirmation before running `gh pr create`
8. Create a PR to main (never push directly)

## Benchmarking

- **Run benchmarks regularly** - Run `cargo bench` to verify performance hasn't regressed
- **Update benchmarks when changing code** - If you modify an operation that's benchmarked, run benchmarks before and after to check for regressions
- **Add new features to benchmarks** - When adding new operations:
  - Add generic `Multivector` operations to `benches/generic.rs`
  - Add specialized 2D/3D euclidean operations to `benches/specialized.rs`
- Benchmarks use criterion; see existing benchmarks for patterns

### Capturing Benchmark Reports

After running benchmarks, capture SVG plots and update documentation:

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

# Update benches/README.md with new timing data if significant changes
# Commit the updated SVGs and README
```

## Mathematical Notation

When documenting GA operations, use standard notation:
- Geometric product: `ab` or `a * b`
- Inner product: `a · b` or `a.inner(b)`
- Outer product: `a ∧ b` or `a.outer(b)`
- Grade selection: `⟨M⟩ₖ` or `m.grade(k)`

## Code Generation for Algebraic Types

**Use the clifford-codegen tool to generate correct implementations.** Manual algebraic derivations are error-prone.

### Using the Codegen Tool

```bash
# Generate code for an algebra
cargo run --package clifford-codegen -- generate algebras/projective3.toml --force

# Discover valid entities for a signature
cargo run --package clifford-codegen -- discover 3 0 1

# List blades for an algebra
cargo run --package clifford-codegen -- blades algebras/projective3.toml
```

### What Codegen Handles

- **Products**: Geometric, exterior, interior, left/right contraction - all generated correctly
- **Constraints**: Study conditions, Plucker conditions, unit norm - solved automatically
- **Constructors**: `new()`, `new_checked()`, `new_unchecked()` with proper constraint solving
- **Verification tests**: Property-based tests comparing against generic Multivector

### Constraints vs. Product Outputs

**Important:** Type constraints apply to *normalized, valid instances*, NOT to algebraic product results.

Generated products use `new_unchecked()` for constrained types because:
1. Product outputs are algebraically correct as computed
2. Constraint solving would incorrectly modify the results
3. A `Motor` from a product is "motor-shaped", not necessarily constraint-satisfying

**When constraints matter:**
- `Motor::normalize()` - explicit re-normalization
- Factory methods like `Motor::from_translation()` - valid by construction
- `Motor::try_from_components()` - validate user inputs

**When constraints don't apply:**
- Product outputs (they are algebraically correct)
- Intermediate computations

### When Adding New Algebras

1. Create a TOML spec in `algebras/` (see `projective3.toml` for example)
2. Define types with grades, fields, and constraints
3. Run `cargo run --package clifford-codegen -- generate algebras/your_algebra.toml`
4. Review generated code in `src/generated/`

**Never hardcode algebraic formulas.** Use the codegen tool to generate correct implementations.

### Generated Code Quality

**The code generator must produce clean, warning-free Rust code:**

- **No clippy warnings** - Generated code must pass `cargo clippy` without warnings
- **No warning suppressions** - Never add `#[allow(dead_code)]`, `#[allow(unused)]`, or similar attributes to hide warnings
- **If generated code produces warnings**, fix the generator to produce correct code instead

**When you see warnings in generated code:**
1. Identify the root cause in the codegen tool
2. Fix the generator to avoid producing that pattern
3. Regenerate all algebras
4. Verify no warnings remain

### Using Generated Products in Extensions

When implementing domain-specific methods in extension files:

1. **Check if a generated product exists** in `generated/products.rs`
2. **Use the generated product** rather than deriving formulas manually
3. **If no suitable product exists**, file an issue for codegen enhancement (PRD-17)

```rust
// CORRECT: Use generated products
use super::generated::products;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::exterior_point_point(self, other)
    }

    pub fn left_contract_plane(&self, plane: &Plane<T>) -> T {
        products::left_contract_point_plane(self, plane)
    }
}

// WRONG: Manual formulas
impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Don't hand-roll - use products::exterior_point_point
        Line::new_unchecked(
            self.e1() * other.e2() - self.e2() * other.e1(),
            // ...
        )
    }
}
```

**Red flags during implementation**:
- Writing multi-term algebraic expressions manually
- Sign corrections or "magic" coefficients
- Copy-pasting formulas from papers without codegen

These indicate missing codegen features - file an issue rather than working around.

### Regenerating Algebras After Codegen Changes

**Whenever you modify the code generator, you MUST regenerate all algebras.**

**Step 1: Update TOML specs if needed**

Some codegen changes require updating the algebra TOML files first:
- New configuration fields (e.g., adding `interior = true`)
- Renamed fields (e.g., `outer` → `exterior`)
- New constraint formats or schema changes

**Step 2: Regenerate all algebras**

```bash
# After any codegen modification
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done

# Verify everything works
cargo fmt
cargo clippy
cargo test
```

**When TOML updates are needed:**
- Adding new product types (must enable in TOML)
- Adding new type-level options
- Changing naming conventions

**When only regeneration is needed:**
- Fixing bugs in existing generation
- Improving generated code quality
- Fixing documentation generation

**Commit pattern:**
1. First commit: The codegen fix/feature
2. Second commit: TOML updates (if needed) + regenerated algebras

Never leave generated code out of sync with codegen.
