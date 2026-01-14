# Implementation Agent

You are implementing features for Clifford, a Rust geometric algebra library.

## Critical Rules

### Do NOT Manually Derive Algebraic Formulas

**Never manually derive or reason about algebraic formulas.** Use `clifford-codegen` for all products, transformations, and constraints.

**Red flags:**
- Writing multi-term algebraic expressions by hand
- Looking up formulas in papers and transcribing them
- "Fixing" signs or coefficients to make tests pass

**When you encounter algebraic issues:**
1. Check if the product is already generated - use `generated/products.rs`
2. If not generated, add it to codegen (`crates/clifford-codegen/`)
3. If codegen seems wrong, fix the codegen tool itself
4. Regenerate all algebras after fixing

### Codegen Genericity

**No algebra-specific shortcuts in codegen.** Derive everything from signature `(p, q, r)`:
- Never check algebra names with string matching
- Never hardcode formulas for specific algebras
- Test with multiple signatures (Euclidean, PGA)

## Requirements

1. **Documentation** - Every public and private item needs rustdoc
2. **No warnings** - `warnings = "deny"` is enforced
3. **No `#[allow(...)]`** - Fix root causes, don't suppress warnings
4. **Standard traits** - Implement Debug, Clone, PartialEq where appropriate
5. **Generic over Float** - Use `T: Float` bounds, not hardcoded f32/f64

## Code Style

### Match Existing Patterns

Before writing new code, study similar modules for:
- Type structure (field ordering, visibility)
- Method organization (constructors → operations → conversions)
- Documentation format
- Test organization

### Naming Conventions

- Methods: `transform_point`, `normalize`, not `apply_to_point`, `get_normalized`
- Constructors: `new()`, `from_*()`, `identity()`, `origin()`
- Accessors: `x()`, `y()`, `z()`, not `get_x()`
- Normalization: `normalize()` returns `Option<Self>`

### Field Visibility

Private fields with public accessors:
```rust
pub struct Vector<T: Float> {
    x: T,  // private
}

impl<T: Float> Vector<T> {
    pub fn x(&self) -> T { self.x }
    pub fn new(x: T, ...) -> Self { ... }
}
```

For constrained types, add `new_unchecked()` and `try_from_components()`.

### Constants Over Magic Values

```rust
// Good
const RELATIVE_EQ_EPS: f64 = 1e-10;
assert!(relative_eq!(a, b, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));

// Bad
assert!(relative_eq!(a, b, epsilon = 1e-10, max_relative = 1e-10));
```

## Code Generation

### Usage

```bash
# Generate for one algebra
cargo run --package clifford-codegen -- generate algebras/projective3.toml --force

# Regenerate all algebras
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done

# Discover entities for a signature
cargo run --package clifford-codegen -- discover 3 0 1
```

### Using Generated Products

Always use generated products in extension files:
```rust
use super::generated::products;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::exterior_point_point(self, other)  // Correct
    }
}
```

Never write manual formulas like:
```rust
Line::new_unchecked(
    self.e1() * other.e2() - self.e2() * other.e1(),  // Wrong - use codegen
    ...
)
```

### Constraints vs Product Outputs

Generated products use `new_unchecked()` for constrained types because:
- Product outputs are algebraically correct as computed
- Constraints apply to normalized instances, not product results

Constraints matter for: `normalize()`, factory methods, `try_from_components()`

### Versor Parity Rules and the Antiproduct

**Critical for PGA algebras:** Versors must be defined with correct grade parity for the antiproduct to work.

The **antiproduct** is defined as: `a ⊛ b = ∁(∁a × ∁b)` where `∁` is the complement.

The complement maps grade `k` → grade `(n-k)`. This has different effects based on dimension:

| Dimension | Complement Effect | Antiproduct Parity |
|-----------|-------------------|-------------------|
| Even (n=4, 3D PGA) | Preserves parity | Same as geometric product |
| Odd (n=3, 2D PGA) | Flips parity | Opposite of geometric product |

**Geometric product parity rules:**
- Even × Even = Even
- Odd × Odd = Even
- Even × Odd = Odd

**Antiproduct parity rules (odd dimension):**
- Even ⊛ Even = **Odd** (flipped!)
- Odd ⊛ Odd = **Odd** (closed!)
- Even ⊛ Odd = **Even** (flipped!)

**Consequence for versor definitions:**

For versors to be closed under antiproduct composition:

| Algebra | Motor grades | Flector grades |
|---------|--------------|----------------|
| 3D PGA (n=4, even) | [0, 2, 4] (even) | [1, 3] (odd) |
| 2D PGA (n=3, odd) | [1, 3] (odd) | [0, 2] (even) |

This ensures `Motor ⊛ Motor = Motor` in all PGA algebras.

**Why antisandwich for transformations:**

In PGA, the regular sandwich `V × X × rev(V)` fails for translations because `e0² = 0` causes terms to vanish. The antisandwich `V ⊛ X ⊛ antirev(V)` goes through complements, avoiding the degenerate metric issue. The `Transform` trait uses antisandwich for this reason.

### Adding New Operations to Codegen

When adding a new product or unary operation to codegen, you must:

1. **Add the trait to `src/ops.rs`** (clifford main crate)
   ```rust
   pub trait MyNewProduct<Rhs = Self> {
       type Output;
       fn my_new_product(&self, rhs: &Rhs) -> Self::Output;
   }
   ```

2. **Export from `src/prelude.rs`**
   ```rust
   pub use crate::ops::{..., MyNewProduct};
   ```

3. **Generate free functions** in `crates/clifford-codegen/src/codegen/products.rs`

4. **Generate trait impls** in `crates/clifford-codegen/src/codegen/traits.rs`:
   - Add import for the new trait in `generate_imports()`
   - Add generation loop in `generate_all_product_traits()`
   - Add trait impl generator method (e.g., `generate_my_new_product_trait()`)

5. **Regenerate all algebras**
   ```bash
   for toml in algebras/*.toml; do
       cargo run --package clifford-codegen -- generate "$toml" --force
   done
   ```

**Existing traits in `clifford::ops`:**
- Binary: `GeometricProduct`, `Wedge`, `Antiwedge`, `Inner`, `LeftContract`, `RightContract`, `Sandwich`, `Antisandwich`, `ScalarProduct`, `BulkContract`, `WeightContract`, `BulkExpand`, `WeightExpand`, `Antigeometric`
- Unary: `Reverse`, `Antireverse`, `LeftComplement`, `RightComplement`, `BulkDual`, `WeightDual`

## Testing

### Symbolica Test Naming

Tests in `clifford-codegen` that use Symbolica **must** be prefixed with `symbolica_`:
```rust
#[test]
fn symbolica_generates_geometric_product() { ... }  // Uses Symbolica - needs prefix

#[test]
fn parse_spec_handles_empty_types() { ... }  // No Symbolica - no prefix needed
```

**When to add the prefix:**
- Any test creating `Algebra`, `ProductTable`, `SymbolicProduct`
- Any test calling `compute_terms()`, `generate_products_file()`
- Any test in `symbolic/` modules

### Property-Based Testing

```rust
proptest! {
    #[test]
    fn rotor_preserves_norm(r in any::<UnitRotor<f64>>(), v in any::<Vector<f64>>()) {
        let rotated = r.rotate(v);
        prop_assert!(relative_eq!(v.norm(), rotated.norm(),
            epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }
}
```

### Wrapper Types

Use generic wrappers from `crate::wrappers`:
- `Unit<T>` - Euclidean norm = 1
- `Bulk<T>` - PGA bulk norm = 1
- `Ideal<T>` - PGA weight norm = 1

Type aliases are generated: `BulkMotor<T>`, `UnitVector<T>`, etc.

## Workflow

1. Branch from latest `origin/main`
2. Write code with full documentation
3. Add property-based tests
4. **Verify before commit**: `cargo fmt && cargo clippy && cargo doc --no-deps && cargo test`
5. Make small, logical commits
6. Ask for confirmation before creating PR

## Benchmarking

When adding or modifying benchmarked operations:
```bash
cargo bench -- --save-baseline before
# make changes
cargo bench -- --baseline before
```

Add new benchmarks to:
- `benches/generic.rs` for generic Multivector operations
- `benches/specialized.rs` for specialized 2D/3D operations
