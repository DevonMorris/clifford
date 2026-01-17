# PRD-47: Wrapper Type Constraint Simplification

**Status**: In Progress
**Depends on**: [PRD-46 Groebner Constraints](prd-46-groebner-constraints.md)
**Goal**: Leverage wrapper type constraints (norm=1, weight=0, etc.) to simplify generated code via Groebner basis reduction

## Implementation Progress

### Completed
- [x] Remove `Deref` from all 7 wrapper types
- [x] Add `WrapperKind` enum to codegen
- [x] Add `WrapperPosition` enum for binary product positions (Lhs/Rhs/Both)
- [x] Extend `AtomToRust` with `.as_inner()` support for wrapper field access
- [x] Add `create_groebner_simplifier_with_wrappers()` for constraint-based reduction
- [x] Generate `Wedge` trait impls for `Unit<T>` / `Unitized<T>`
- [x] Generate `Antiwedge` trait impls for `Unit<T>` / `Unitized<T>`
- [x] Add wrapper equivalence tests for positive-definite norm algebras
- [x] Fix test generation to exclude indefinite-norm algebras (hyperbolic, conformal, etc.)

### Not Yet Implemented
- [ ] **Sandwich/Antisandwich products** ← *This is where real speedup comes from*
- [ ] **Geometric products** (for versors like Rotor, Motor)
- [ ] **Interior products** (contractions, expansions)
- [ ] **Transform trait** for wrapper versors
- [ ] Bulk<Motor> optimized sandwich for point transformation
- [ ] Accessor delegation methods on wrappers
- [ ] Benchmarks comparing wrapper vs bare type performance

## Problem Statement

### Current State

The codebase has 7 wrapper types that enforce constraints at the type level:

| Wrapper | Constraint | Polynomial Form |
|---------|------------|-----------------|
| `Unit<T>` | `norm() == 1` | `x² + y² + z² - 1 = 0` |
| `Bulk<T>` | `bulk_norm() == 1` | `s² + r² + ... - 1 = 0` |
| `Unitized<T>` | `weight_norm() == 1` | `w² - 1 = 0` (simplified) |
| `Ideal<T>` | `weight_norm() == 0` | `w = 0` |
| `Proper<T>` | timelike, `|norm²| == 1` | `t² - x² - y² - z² - 1 = 0` |
| `Spacelike<T>` | spacelike, `|norm²| == 1` | `t² - x² - y² - z² + 1 = 0` |
| `Null<T>` | `norm² == 0` | `t² - x² - y² - z² = 0` |

**Problem 1: Deref Leaks Constraints**

All wrapper types implement `Deref<Target = T>`, which allows implicit access to the inner value. This creates issues:

```rust
let unit_v: Unit<Vector> = Unit::new_normalize(v);

// Deref allows this - but the result is NOT Unit<Vector>!
let scaled = unit_v.scale(2.0);  // Returns Vector, not Unit<Vector>

// User might expect this to still be normalized, but it isn't
```

More critically, `Deref` makes it impossible to implement specialized operations on wrappers that leverage their constraints, because method resolution always finds the inner type's methods first.

**Problem 2: Wrapper Constraints Unused in Codegen**

PRD-46 introduced Groebner basis constraint simplification for algebraic constraints (e.g., Plücker condition for Lines). However, wrapper type constraints are NOT used during code generation:

```rust
// When computing Unit<Rotor> * Unit<Rotor>, we could use:
//   r₁.s² + r₁.xy² + r₁.yz² + r₁.xz² = 1
//   r₂.s² + r₂.xy² + r₂.yz² + r₂.xz² = 1
//
// But currently, this constraint information is lost.
```

**Problem 3: No Wrapper-Specific Trait Implementations**

Without removing `Deref`, we cannot implement traits like `GeometricProduct` directly on wrapper types with optimized code paths:

```rust
// This doesn't work well with Deref because Rotor's impl takes precedence
impl GeometricProduct<Unit<Rotor>> for Unit<Rotor> {
    type Output = Unit<Rotor>;  // Result is guaranteed normalized!

    fn geometric_product(&self, rhs: &Unit<Rotor>) -> Self::Output {
        // Could use constraint-simplified code here
    }
}
```

### Example: Unit Rotor Composition

**Current generated code for `Rotor * Rotor`:**
```rust
// All 16 term products computed naively
Rotor::new(
    self.s() * rhs.s() - self.xy() * rhs.xy() - self.yz() * rhs.yz() - self.xz() * rhs.xz(),
    self.s() * rhs.xy() + self.xy() * rhs.s() + self.yz() * rhs.xz() - self.xz() * rhs.yz(),
    // ... more terms
)
```

**With Unit constraint (norm=1) applied:**
```rust
// Groebner reduction can eliminate redundant terms using:
//   s² + xy² + yz² + xz² = 1
// This substitution can simplify cross-terms and eliminate
// expressions that are algebraically equivalent under the constraint.
```

## Solution

### Phase 1: Remove Deref, Add Explicit Access

Remove `Deref` from all wrapper types. Users must use explicit methods:

```rust
// Before (with Deref)
let unit_v: Unit<Vector> = Unit::new_normalize(v);
let x = unit_v.x();  // Implicit deref to Vector

// After (without Deref)
let unit_v: Unit<Vector> = Unit::new_normalize(v);
let x = unit_v.as_inner().x();  // Explicit access
// OR via trait if implemented on Unit<Vector>
let x = unit_v.x();  // If Unit<Vector> has x() method
```

**Keep `AsRef<T>` for borrowing scenarios** - this is explicit and appropriate.

### Phase 2: Define Wrapper Constraint Polynomials

Add a new trait for wrapper types that provides their constraint polynomials:

```rust
/// Trait for wrapper types that impose algebraic constraints.
///
/// The constraints are expressed as polynomial equations that
/// equal zero for valid wrapped values.
pub trait WrapperConstraint {
    /// Returns constraint polynomials as Symbolica Atoms.
    ///
    /// Each returned polynomial `p` represents `p = 0`.
    /// These constraints can be used by the Groebner basis
    /// simplification system during code generation.
    fn constraint_polynomials(field_prefix: &str) -> Vec<Atom>;

    /// Human-readable name for the constraint.
    fn constraint_name() -> &'static str;
}
```

### Phase 3: Implement Constraints for Each Wrapper

```rust
impl<T: Normed> WrapperConstraint for Unit<T> {
    fn constraint_polynomials(prefix: &str) -> Vec<Atom> {
        // norm² = 1  →  norm² - 1 = 0
        // For a 3D vector: x² + y² + z² - 1 = 0
        T::norm_squared_polynomial(prefix) - Atom::num(1)
    }

    fn constraint_name() -> &'static str {
        "unit_norm"
    }
}

impl<T: DegenerateNormed> WrapperConstraint for Unitized<T> {
    fn constraint_polynomials(prefix: &str) -> Vec<Atom> {
        // weight_norm² = 1  →  weight_norm² - 1 = 0
        T::weight_norm_squared_polynomial(prefix) - Atom::num(1)
    }

    fn constraint_name() -> &'static str {
        "unitized"
    }
}

impl<T: DegenerateNormed> WrapperConstraint for Ideal<T> {
    fn constraint_polynomials(prefix: &str) -> Vec<Atom> {
        // weight_norm = 0  →  all weight components = 0
        T::weight_field_polynomials(prefix)  // Returns each weight field = 0
    }

    fn constraint_name() -> &'static str {
        "ideal"
    }
}

impl<T: IndefiniteNormed> WrapperConstraint for Null<T> {
    fn constraint_polynomials(prefix: &str) -> Vec<Atom> {
        // norm² = 0 (Minkowski)
        T::norm_squared_polynomial(prefix)  // The polynomial itself = 0
    }

    fn constraint_name() -> &'static str {
        "null"
    }
}
```

### Phase 4: Extend Codegen to Use Wrapper Constraints

#### 4.1 Add Wrapper Type Specification to Codegen

Extend the type specification system to understand wrapped types:

```rust
/// Represents a possibly-wrapped type for code generation.
#[derive(Debug, Clone)]
pub enum WrappedTypeSpec {
    /// Bare type with no wrapper constraints
    Bare(TypeSpec),
    /// Type wrapped with constraint
    Wrapped {
        wrapper: WrapperKind,
        inner: TypeSpec,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WrapperKind {
    Unit,      // norm = 1
    Bulk,      // bulk_norm = 1
    Unitized,  // weight_norm = 1
    Ideal,     // weight_norm = 0
    Proper,    // timelike, |norm²| = 1
    Spacelike, // spacelike, |norm²| = 1
    Null,      // norm² = 0
}
```

#### 4.2 Extend ProductConstraintCollector

```rust
impl<'a> ProductConstraintCollector<'a> {
    /// Collects constraints for a wrapped type.
    pub fn collect_wrapped_constraints(
        &self,
        wrapped: &WrappedTypeSpec,
        prefix: &str,
    ) -> Vec<Atom> {
        let mut constraints = Vec::new();

        match wrapped {
            WrappedTypeSpec::Bare(ty) => {
                // Existing constraint collection (Plücker, Study, etc.)
                constraints.extend(self.collect_constraints(ty, prefix));
            }
            WrappedTypeSpec::Wrapped { wrapper, inner } => {
                // 1. Inner type constraints
                constraints.extend(self.collect_constraints(inner, prefix));

                // 2. Wrapper constraints
                constraints.extend(self.wrapper_constraints(*wrapper, inner, prefix));
            }
        }

        constraints
    }

    fn wrapper_constraints(
        &self,
        wrapper: WrapperKind,
        inner: &TypeSpec,
        prefix: &str,
    ) -> Vec<Atom> {
        match wrapper {
            WrapperKind::Unit => {
                // norm² - 1 = 0
                let norm_sq = self.deriver.derive_norm_squared(inner, prefix);
                vec![norm_sq - Atom::num(1)]
            }
            WrapperKind::Unitized => {
                // weight_norm² - 1 = 0
                let weight_sq = self.deriver.derive_weight_norm_squared(inner, prefix);
                vec![weight_sq - Atom::num(1)]
            }
            WrapperKind::Ideal => {
                // Each weight component = 0
                self.deriver.derive_weight_components(inner, prefix)
            }
            WrapperKind::Null => {
                // norm² = 0 (the polynomial itself)
                let norm_sq = self.deriver.derive_norm_squared(inner, prefix);
                vec![norm_sq]
            }
            // ... other wrappers
        }
    }
}
```

### Phase 5: Generate Wrapper Trait Implementations

For each wrapper type combination, generate optimized trait implementations:

```rust
// Generated code for Unit<Rotor<f64>>
impl GeometricProduct<Unit<Rotor<f64>>> for Unit<Rotor<f64>> {
    type Output = Unit<Rotor<f64>>;

    #[inline]
    fn geometric_product(&self, rhs: &Unit<Rotor<f64>>) -> Self::Output {
        // Groebner-simplified expression using:
        //   self.s² + self.xy² + self.yz² + self.xz² = 1
        //   rhs.s² + rhs.xy² + rhs.yz² + rhs.xz² = 1
        Unit::new_unchecked(Rotor::new(
            // Simplified terms here
        ))
    }
}
```

### Phase 6: Accessor Delegation (Optional)

To maintain ergonomics after removing `Deref`, generate accessor methods on wrappers:

```rust
// Generated for Unit<Vector<S>>
impl<S: Scalar> Unit<Vector<S>> {
    #[inline]
    pub fn x(&self) -> S { self.as_inner().x() }
    #[inline]
    pub fn y(&self) -> S { self.as_inner().y() }
    #[inline]
    pub fn z(&self) -> S { self.as_inner().z() }
}
```

This is optional - users can always use `as_inner()` explicitly.

## Configuration

### TOML Specification (Optional)

For types that should have wrapper-specific codegen:

```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy", "yz", "xz"]

# Optional: enable wrapper-specific trait generation
[types.Rotor.wrappers]
unit = true  # Generate Unit<Rotor> trait impls
```

If not specified, wrapper codegen is disabled for that type (fallback to generic wrapper behavior).

### CLI Flags

```bash
# Generate wrapper trait implementations
clifford-codegen --wrapper-traits algebras/euclidean3.toml

# Disable wrapper codegen (use generic wrappers only)
clifford-codegen --no-wrapper-traits algebras/euclidean3.toml
```

## Migration Guide

### Breaking Change: Deref Removal

**Before:**
```rust
let unit: Unit<Vector<f64>> = Unit::new_normalize(v);
let x = unit.x();           // Via Deref
let norm = unit.norm();     // Via Deref
let scaled = unit.scale(2.0);  // Via Deref - returns Vector!
```

**After:**
```rust
let unit: Unit<Vector<f64>> = Unit::new_normalize(v);

// Option 1: Explicit access
let x = unit.as_inner().x();
let norm = unit.as_inner().norm();

// Option 2: Generated accessors (if enabled)
let x = unit.x();           // Delegated to inner

// Option 3: Use wrapper trait impls
let norm = unit.norm();     // Implemented on Unit<T> directly

// Scaling now requires explicit unwrap
let scaled = unit.into_inner().scale(2.0);
```

### Semver Implications

This is a **breaking change** requiring a major version bump:
- `Deref` removal breaks method resolution
- Some trait bounds may need updating
- Pattern matching on wrapper types unchanged

## Implementation Phases

### Phase 1: Deref Removal ✅ COMPLETE
- [x] Remove `impl Deref for Unit<T>`
- [x] Remove `impl Deref for Bulk<T>`
- [x] Remove `impl Deref for Unitized<T>`
- [x] Remove `impl Deref for Ideal<T>`
- [x] Remove `impl Deref for Proper<T>`
- [x] Remove `impl Deref for Spacelike<T>`
- [x] Remove `impl Deref for Null<T>`
- [ ] Add explicit accessor delegation methods
- [x] Update all internal usages

### Phase 2: Constraint Infrastructure ✅ COMPLETE
- [x] Add `WrapperKind` enum to `spec/ir.rs`
- [x] Add `WrapperPosition` enum to `codegen/traits.rs`
- [x] Add norm polynomial derivation to `ConstraintDeriver`
- [x] Add `collect_wrapper_constraints()` to `ProductConstraintCollector`
- [x] Add wrapper-aware `AtomToRust` converter with `.as_inner()` support

### Phase 3: Wrapper Product Trait Generation ⚠️ PARTIAL

**Approach**: For every `impl Trait<B> for A` we generate, also generate:
- `impl Trait<B> for Unit<A>` (LHS wrapped)
- `impl Trait<Unit<B>> for A` (RHS wrapped)
- `impl Trait<Unit<B>> for Unit<A>` (Both wrapped)

Same pattern for `Unitized<T>` in PGA algebras.

| Product | Status | Constraint Simplification? |
|---------|--------|---------------------------|
| `Wedge` | ✅ Done | No (exterior, no metric) |
| `Antiwedge` | ✅ Done | No (exterior, no metric) |
| `LeftContract` | ❌ TODO | Yes (uses metric) |
| `RightContract` | ❌ TODO | Yes (uses metric) |
| `ScalarProduct` | ❌ TODO | Yes (uses metric) |
| `Dot` | ❌ TODO | Yes (uses metric) |
| `Antidot` | ❌ TODO | Yes (uses metric) |
| `BulkContract` | ❌ TODO | Yes (PGA metric) |
| `WeightContract` | ❌ TODO | Yes (PGA metric) |
| `BulkExpand` | ❌ TODO | Yes (PGA metric) |
| `WeightExpand` | ❌ TODO | Yes (PGA metric) |
| `Project` | ❌ TODO | Yes (uses metric) |
| `Antiproject` | ❌ TODO | Yes (uses metric) |
| `Sandwich` | ❌ TODO | Yes (versor, big win) |
| `Antisandwich` | ❌ TODO | Yes (versor, big win) |
| `Transform` | ❌ TODO | Yes (delegates to above) |
| `Versor` | ❌ TODO | Yes (composition) |

Implementation tasks:
- [x] Add `generate_wrapper_product_trait()` helper
- [x] Generate wrapper variants for `Wedge`
- [x] Generate wrapper variants for `Antiwedge`
- [ ] Generate wrapper variants for `LeftContract` / `RightContract`
- [ ] Generate wrapper variants for `ScalarProduct`
- [ ] Generate wrapper variants for `Dot` / `Antidot`
- [ ] Generate wrapper variants for `BulkContract` / `WeightContract`
- [ ] Generate wrapper variants for `BulkExpand` / `WeightExpand`
- [ ] Generate wrapper variants for `Project` / `Antiproject`
- [ ] Generate wrapper variants for `Sandwich` / `Antisandwich`
- [ ] Generate wrapper variants for `Transform`
- [ ] Generate wrapper variants for `Versor`

**Note**: For products that use the metric (all except Wedge/Antiwedge), the Groebner basis with wrapper constraints (e.g., `norm²=1`) can simplify the generated expressions. For exterior products, wrapper impls just delegate with no simplification benefit.

### Phase 6: Testing ⚠️ PARTIAL
- [x] Property tests: Unit/Bulk wrapper norm equivalence
- [ ] Property tests: wrapped products equal unwrapped
- [ ] Verify constraint polynomials reduce correctly
- [ ] Benchmark: wrapped vs unwrapped performance
- [x] Fix test generation for indefinite-norm algebras

### Phase 7: Documentation
- [ ] Update wrapper rustdoc for new API
- [ ] Add migration guide
- [ ] Document constraint simplification benefits

## Expected Impact

### Term Reduction Examples

| Product | Wrapper Constraint | Before | After | Reduction |
|---------|-------------------|--------|-------|-----------|
| `Bulk<Motor>.antisandwich(Point)` | bulk_norm²=1 | 48+ | TBD | TBD |
| `Unit<Rotor>.sandwich(Vector)` | norm²=1 | 27 | TBD | TBD |
| `Unitized<Point>.antiwedge(Plane)` | weight=1 | 12 | ~8 | 33% |
| `Unit<Vector>.dot(Vector)` | norm²=1 | 3 | TBD | TBD |

**Note**: Actual term reduction needs to be measured after implementation.
Exterior products (Wedge/Antiwedge) don't use the metric, so constraint simplification provides no benefit there.
The real benefit comes from interior products and sandwich products which involve the metric.

### Type Safety Benefits

1. **Constraint preservation**: Versor composition preserves normalization
2. **Compile-time guarantees**: No accidental denormalization
3. **API clarity**: Explicit access via `as_inner()` prevents confusion

### Performance Benefits

1. **Fewer operations**: Groebner reduction eliminates redundant terms
2. **No runtime checks**: Constraints verified at type level, not runtime
3. **Better inlining**: Smaller functions inline more aggressively

## Risks and Mitigations

### Risk: Breaking Change Disruption

**Mitigation:**
- Provide detailed migration guide
- Consider deprecation period with both APIs
- Generate accessor delegation methods for common patterns

### Risk: Codegen Complexity

**Mitigation:**
- Start with `Unit<T>` only, expand incrementally
- Use feature flag for wrapper codegen
- Thorough testing before enabling by default

### Risk: Combinatorial Explosion

For N types and 7 wrappers, naive approach generates 7N² trait implementations.

**Mitigation:**
- Only generate for types marked in TOML
- Focus on high-value combinations (rotors, motors)
- Use macro generation for accessor delegation

### Risk: Constraint Polynomial Derivation Complexity

Different algebras have different norm definitions.

**Mitigation:**
- Extend existing `ConstraintDeriver` infrastructure
- Norm computation already exists in algebra spec
- Leverage `derive_norm_squared()` pattern

## Files Changed

| File | Action | Description |
|------|--------|-------------|
| `src/wrappers.rs` | Modified | Remove Deref impls, add accessor methods |
| `crates/clifford-codegen/src/symbolic/constraint_derive.rs` | Modified | Add norm polynomial derivation |
| `crates/clifford-codegen/src/symbolic/groebner.rs` | Modified | Add wrapper constraint collection |
| `crates/clifford-codegen/src/codegen/traits.rs` | Modified | Generate wrapper trait implementations |
| `crates/clifford-codegen/src/codegen/wrappers.rs` | Created | Wrapper-specific codegen logic |
| `tests/wrapper_constraints.rs` | Created | Property tests for wrapper constraints |

## Success Criteria

1. **Deref removed**: All 7 wrapper types no longer implement `Deref`
2. **AsRef preserved**: `AsRef<T>` still available for explicit borrowing
3. **Constraints derived**: `ConstraintDeriver` can produce norm polynomials
4. **Groebner integration**: Wrapper constraints feed into simplification
5. **Trait generation**: At least `Unit<Rotor>` has optimized geometric product
6. **Tests pass**: All existing tests updated and passing
7. **Performance**: Measurable term reduction for constrained products

## References

- [PRD-46: Groebner Basis Constraint Simplification](prd-46-groebner-constraints.md)
- [PRD-18.2: Geometry-Specific Wrapper Types](prd-18.2-wrappers.md)
- [Rust API Guidelines: Smart pointers](https://rust-lang.github.io/api-guidelines/predictability.html#c-smart-ptr)
- [Deref polymorphism antipattern](https://rust-lang.github.io/api-guidelines/predictability.html#c-deref)
