# Clifford - Geometric Algebra Library

A Rust library for Geometric Algebra (Clifford Algebra).

## Claude as Collaborator

**Claude is a collaborator on this project, not a task runner.**

When implementing new features, think deeply about the design. Don't just implement what's asked—question assumptions and consider the broader implications:

### A. Challenge Hidden Assumptions

Does the feature operate under false or unstated assumptions? Common traps:
- Assuming a particular metric signature (Euclidean vs Minkowski vs degenerate)
- Assuming handedness or orientation conventions
- Assuming specific basis ordering (e01 vs e10)
- Assuming normalized inputs when the math works for any magnitude

### B. Verify Generalization

Does the feature generalize across all algebras and implementations?
- Will this work for Euclidean, Projective, Conformal, and Minkowski algebras?
- Does it handle degenerate (null) elements correctly?
- What constraints must hold for it to generalize? Document them.
- If it's algebra-specific, is it in the right module?

### C. Question Canonical Choices

Does the feature assume something canonical about ordering or invariants?
- Is there a "natural" ordering, or is it an arbitrary convention?
- Should the user control the convention, or is one objectively correct?
- Is a layer of abstraction missing that would let users choose?
- Example: `to_array()` implies an ordering—should it be `to_array_grade_order()` vs `to_array_blade_order()`?

### D. Anticipate Future Features

Think 3-4 steps ahead:
- What features might build on this one?
- Does this design paint us into a corner?
- Will future features require breaking changes to this API?
- Example: A `normalize()` method today might conflict with `unitize()` vs `bulk_normalize()` distinctions needed for PGA later.

**When in doubt, raise concerns before implementing.** A conversation about design is more valuable than code that needs to be rewritten.

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
- Work in feature branches: `feat/multivector`, `fix/product-sign`
- **Small, logical commits** - each commit should be a single logical change
- Use conventional commit format: `type(scope): description`
- **Confirm before creating PRs** - always ask for user confirmation before running `gh pr create`
- **Review before merging** - wait for CI and Greptile, address feedback, then merge with `--squash --delete-branch`

### 4. Performance & Benchmarking
- Use SIMD instructions where beneficial (via `std::arch` or `portable_simd`)
- Benchmark critical paths with Criterion
- Profile before optimizing
- See the **optimize agent** for detailed benchmarking workflows

### 5. Minimal Dependencies
- Prefer std library where possible
- Only add dependencies when they provide significant value
- Dependencies must use permissive licenses (run `cargo deny check`)

### 6. Rust API Guidelines
- Follow the official [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- Implement standard traits (Debug, Clone, PartialEq, etc.)
- **No warning suppression** - never use `#[allow(...)]`, fix root causes instead
- **Private fields with public accessors** for specialized types
- **Use Clifford types in APIs**, not tuples

### 7. Testing
- **Property-based testing is mandatory** with `proptest`
- Use `prop_assert!` inside `proptest!` blocks
- Use `relative_eq!` with both `epsilon` AND `max_relative` parameters
- Use `RELATIVE_EQ_EPS` constant from `crate::test_utils`
- See the **test agent** for detailed testing patterns

#### Visual Testing Philosophy

**CRITICAL: Never use golden image ("blessed") tests for visual code.**

Golden image tests are change detectors, not correctness proofs. They tell you *something changed* but not *whether it's correct*. For visualization code, use:

1. **Coordinate assertions** - Verify world-to-screen coordinate mappings
2. **Visual invariants** - Property-based tests for geometric properties (rotation preserves distance, joined points are incident with line, etc.)
3. **Scene graph assertions** - Test that primitives exist with correct properties
4. **Automated inspection** - Flag blank frames, NaN artifacts, clipping issues for human review

**Wrong:**
```rust
// Just detects change, doesn't prove correctness
fn test_rotor_visual() {
    render_frame();
    assert_images_match!("golden/rotor_45deg.png");
}
```

**Right:**
```rust
// Proves the rotation is actually correct
fn test_rotor_rotates_correctly() {
    let rotor = Rotor::from_angle(FRAC_PI_4);
    let input = Vector::new(1.0, 0.0);
    let output = rotor.transform(&input);

    // Verify the actual geometric property
    relative_eq!(output.x(), FRAC_1_SQRT_2, epsilon = 1e-6);
    relative_eq!(output.y(), FRAC_1_SQRT_2, epsilon = 1e-6);
}
```

See PRD-48.10 for detailed visual testing patterns.

### 8. Code Review
- PRs are reviewed by Greptile (AI-powered review)
- Address feedback before merging; comment `@greptileai review` after fixes
- **No `todo!()` macros** - code must be complete

### 9. Claude Code Agents

Specialized agents in `.claude/agents/` handle different tasks:

| Agent | Purpose |
|-------|---------|
| **implement** | Implementing features, code changes, codegen |
| **test** | Writing property-based tests |
| **document** | Writing rustdoc and documentation |
| **review** | Code review checklist |
| **explain** | Teaching GA concepts |
| **optimize** | Performance optimization, benchmarking |
| **devops** | CI/CD and infrastructure |
| **release** | Version bumps, publishing |

### 10. Code Navigation with ctags

Use the ctags index (`.claude/tags`) for efficient code navigation:
- When looking up definitions (classes, functions, traits), use `/ctags-usage` skill
- Prefer ctags lookups over broad grep searches for specific symbol definitions
- Run `/refresh` to regenerate the index after major code changes

### 11. Code Generation

**CRITICAL: Do NOT manually derive algebraic formulas.**

Use the `clifford-codegen` tool for all algebraic operations:
- Never manually derive geometric products, sandwich products, or transformations
- If codegen produces wrong results, fix the codegen tool, not the generated code
- Regenerate all algebras after codegen changes

See the **implement agent** for detailed codegen usage.

### 12. Semantic Field Naming

**CRITICAL: Field names must reflect what the field DOES, not which blade it corresponds to.**

When naming fields in algebra TOML files:
- Name fields based on their **geometric/physical meaning**, not their blade indices
- Test by reading the code: `from_translation(dx, dy, dz)` should use `tx`, `ty`, `tz` fields
- If a field controls x-rotation, name it `rx`, not based on which e_ij blade it is

**Example - 3D PGA Motor:**
The Motor has 8 components. In point-based PGA with antisandwich:
- Positions 1,2,3 (e12, e13, e14) control **translation** → name as `tz`, `ty`, `tx`
- Positions 4,5,6 (e23, e24, e34) control **rotation** → name as `rx`, `ry`, `rz`
- Position 7 (e0123) is the pseudoscalar → name as `ps`

**Wrong:** Naming e12 as `bz` because it's the xy-plane bivector
**Right:** Naming e12 as `tz` because it controls z-translation in this formulation

Always verify by checking how factory methods use the fields.

### 13. Extension Methods

**CRITICAL: Prefer generated traits over extension methods.**

The `extensions.rs` files in each specialized module should be minimal:

1. **NEVER shadow a trait with an extension method** - If a trait like `Wedge`, `Antiwedge`, `Transform`, or `Sandwich` already provides the operation, do NOT add a method like `meet()` or `join()` that just calls the trait.

2. **Extensions are primarily for constructors** - Factory methods like `from_cartesian()`, `from_angle_plane()`, `from_translation()`, `origin()`, `identity()` belong in extensions.

3. **Do NOT add trivial wrapper methods** - Methods like `meet_plane(&self, plane)` that just call `self.antiwedge(plane)` add no value. Users should use the trait directly.

**What belongs in extensions.rs:**
- Constructors and factory methods
- Coordinate extraction (e.g., `cartesian_x()`, `to_cartesian()`)
- Geometric queries that require multiple operations (e.g., `is_parallel()`, `distance()`)
- Methods with domain-specific validation or normalization

**What does NOT belong in extensions.rs:**
- `join()` → use `Wedge::wedge()`
- `meet()` → use `Antiwedge::antiwedge()`
- `transform()` → use `Transform::transform()`
- `meet_plane()`, `meet_line()` → use `Antiwedge::antiwedge()`
- Any method that is just `self.some_trait_method(other)`

**Example - Wrong:**
```rust
impl Line {
    pub fn meet(&self, other: &Line) -> Point {
        self.antiwedge(other)  // Just shadows the trait!
    }
}
```

**Example - Right:**
```rust
// Let users call the trait directly:
use clifford::ops::Antiwedge;
let intersection = line1.antiwedge(&line2);
```

## Development Commands

```bash
cargo build           # Build the library
cargo nextest run     # Run tests (recommended - handles Symbolica correctly)
cargo test            # Run tests (fallback)
cargo doc --open      # Generate and view documentation
cargo bench           # Run benchmarks
cargo clippy          # Run linter
cargo fmt             # Format code
cargo deny check      # Check licenses and advisories
```

### Verification Workflow

**Before every commit**:
```bash
cargo fmt && cargo clippy && cargo doc --no-deps && cargo nextest run && cargo deny check
```

## Module Structure

```
crates/
  clifford/           # Main library crate
    src/
      specialized/
        euclidean/
          dim2/       # 2D: Vector, Bivector, Rotor
          dim3/       # 3D: Vector, Bivector, Trivector, Rotor
        projective/
          dim3/       # 3D PGA: Point, Line, Plane, Motor, Flector
  clifford-codegen/   # Code generation tool
```

Naming: Use full words (`Vector` not `Vec`), don't include dimension in type names.

## Resources

### Authoritative Reference

**The [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page) is the authoritative reference for all GA operations.**

### RGA Product Notation

| Symbol | Name | Function Pattern |
|--------|------|------------------|
| `∧` | wedge | `wedge_*` |
| `∨` | antiwedge | `antiwedge_*` |
| `★` | dual | `dual_*` |
| `☆` | antidual | `antidual_*` |

**Interior Products**: `bulk_contraction_*`, `weight_contraction_*`, `bulk_expansion_*`, `weight_expansion_*`

### Learning Resources
- [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page)
- [Look, Ma, No Matrices!](https://enkimute.github.io/LookMaNoMatrices/)
