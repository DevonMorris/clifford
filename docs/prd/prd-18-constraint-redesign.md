# PRD-18: Constraint System Redesign with Normed Traits

**Status**: Complete
**Goal**: Simplify the constraint system by removing user configuration, inferring constraints from algebra structure, and introducing geometry-specific wrapper types for normalized entities.

## Overview

The current constraint system exposes unnecessary complexity to users via configurable constraints in TOML files. The signature, metric, and antimetric of an algebra fully determine valid constraints - users shouldn't need to specify them manually.

Additionally, some entities are useful both normalized and unnormalized:
- `Motor` can represent scaling transformations
- `Bulk<Motor>` should represent _only_ rigid body transformations

This PRD redesigns the constraint system to:
1. Remove user-configurable constraints from TOML
2. Infer constraints automatically during code generation
3. Introduce a `Normed` trait hierarchy for different algebra types
4. Introduce geometry-specific wrappers (`Unit`, `Bulk`, `Ideal`, `Proper`)
5. Simplify constructors to just `new_unchecked()` and `new_checked()`

## Problem Statement

### Current Complexity

Users must specify constraints in TOML like this:

```toml
[[types.Motor.constraints]]
name = "study"
description = "Study condition for valid motors"
expression = "s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0"
solve_for = "e0123"
```

This has several issues:
1. **Redundant**: Constraints are derivable from the algebra's structure
2. **Error-prone**: Users might specify incorrect expressions
3. **Confusing**: `new()` vs `new_unchecked()` vs `new_checked()` is unclear
4. **Incomplete**: No way to represent "normalized" vs "unnormalized" variants

### Missing Unit Type Pattern

nalgebra has `Unit<T>` which wraps any type and guarantees unit norm. Clifford has no equivalent - there's no way to distinguish a `Motor` that might have any norm from a unit `Motor` representing only rigid transformations.

### Constructor Confusion

Current generated constructors:
- `new(params)` - Solves constraint, returns `Option<Self>` for quadratic
- `new_checked(all_fields, tolerance)` - Validates constraint
- `new_unchecked(all_fields)` - Raw construction

The distinction between `new` and `new_unchecked` is confusing - most users want raw construction, not constraint solving.

## Proposed Solution

### 1. Remove User-Configurable Constraints from TOML

**Before:**
```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03", "e0123"]
versor = true

[[types.Motor.constraints]]
name = "study"
expression = "s*e0123 + e23*e01 + e31*e02 + e12*e03 = 0"
solve_for = "e0123"
```

**After:**
```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03", "e0123"]
versor = true
# Constraints inferred automatically from algebra structure
```

### 2. Compute Constraints During Generation

Extend existing `satisfies_geometric_constraint()` and `satisfies_antiproduct_constraint()` to:

1. **Detect constraint type** based on grade structure:
   - Grade [0, 2] in 3D Euclidean → Unit norm constraint
   - Grade [0, 2, 4] in 4D PGA → Study condition
   - Grade [2] in 4D PGA → Plücker condition
   - Grade [1, 3] in 4D PGA → Flector geometric constraint

2. **Derive constraint expressions** symbolically using Symbolica

3. **Generate appropriate validation code** in `new_checked()`

### 3. Simplify Constructors

Remove `new()` in favor of just two constructors:

```rust
impl<T: Float> Motor<T> {
    /// Creates a Motor without any validation.
    pub fn new_unchecked(s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T, e0123: T) -> Self {
        Self { s, e23, e31, e12, e01, e02, e03, e0123 }
    }

    /// Creates a Motor with constraint validation.
    pub fn new_checked(
        s: T, e23: T, e31: T, e12: T, e01: T, e02: T, e03: T, e0123: T,
        tolerance: T,
    ) -> Result<Self, ConstraintError> {
        let residual = s * e0123 + e23 * e01 + e31 * e02 + e12 * e03;
        if residual.abs() > tolerance {
            return Err(ConstraintError::new("Motor", "study", residual.to_f64().unwrap_or(0.0)));
        }
        Ok(Self { s, e23, e31, e12, e01, e02, e03, e0123 })
    }
}
```

### 4. Introduce `Normed` Trait Hierarchy

Different geometric algebras have different norm semantics:

| Algebra | Traits | Notes |
|---------|--------|-------|
| Euclidean (Cl(n,0,0)) | `Normed` | Standard positive-definite norm |
| PGA (Cl(n,0,1)) | `Normed`, `DegenerateNormed` | Bulk/weight norms for degenerate basis |
| CGA (Cl(n+1,1,0)) | `Normed`, `ConformalNormed` | Null vector handling |
| Minkowski (Cl(1,3,0)) | `Normed`, `IndefiniteNormed` | Timelike/spacelike classification |

#### Base `Normed` Trait

```rust
pub trait Normed {
    type Scalar: Float;

    fn norm_squared(&self) -> Self::Scalar;
    fn norm(&self) -> Self::Scalar { self.norm_squared().abs().sqrt() }
    fn try_normalize(&self) -> Option<Self> where Self: Sized;
    fn normalize(&self) -> Self where Self: Sized;
    fn scale(&self, factor: Self::Scalar) -> Self where Self: Sized;
}
```

#### `DegenerateNormed` Trait (PGA)

```rust
pub trait DegenerateNormed: Normed {
    fn bulk_norm_squared(&self) -> Self::Scalar;
    fn bulk_norm(&self) -> Self::Scalar;
    fn weight_norm_squared(&self) -> Self::Scalar;
    fn weight_norm(&self) -> Self::Scalar;
    fn try_unitize(&self) -> Option<Self> where Self: Sized;
    fn unitize(&self) -> Self where Self: Sized;
}
```

#### `IndefiniteNormed` Trait (Minkowski)

```rust
pub trait IndefiniteNormed: Normed {
    fn is_timelike(&self) -> bool;
    fn is_spacelike(&self) -> bool;
    fn is_lightlike(&self) -> bool;
    fn causal_character(&self) -> CausalCharacter;
}
```

#### `ConformalNormed` Trait (CGA)

```rust
pub trait ConformalNormed: Normed {
    fn einf_coefficient(&self) -> Self::Scalar;
    fn is_null(&self) -> bool;
    fn try_normalize_cga(&self) -> Option<Self> where Self: Sized;
}
```

### 5. Introduce Geometry-Specific Wrapper Types

Each geometry type has its own unambiguous wrapper:

| Wrapper | Constraint | Use Case |
|---------|------------|----------|
| `Unit<T>` | `norm() == 1` | Euclidean vectors, bivectors, rotors |
| `Bulk<T>` | `bulk_norm() == 1` | PGA versors (motors, flectors) - rigid transforms |
| `Ideal<T>` | `weight_norm() == 1` | PGA homogeneous coords (points, planes) |
| `Proper<T>` | `is_timelike()` | Minkowski timelike vectors |

**Rationale for names:**
- `Unit<T>` - standard mathematical term for norm = 1
- `Bulk<T>` - PGA terminology for the non-degenerate part
- `Ideal<T>` - PGA terminology for the degenerate/projective part
- `Proper<T>` - physics terminology for proper time/proper length

### 6. Generate Type Aliases

```rust
// Euclidean
pub type UnitVector<T> = Unit<Vector<T>>;
pub type UnitRotor<T> = Unit<Rotor<T>>;

// PGA
pub type BulkMotor<T> = Bulk<Motor<T>>;
pub type BulkFlector<T> = Bulk<Flector<T>>;
pub type IdealPoint<T> = Ideal<Point<T>>;
pub type IdealPlane<T> = Ideal<Plane<T>>;

// Minkowski (future)
pub type ProperFourVelocity<T> = Proper<FourVector<T>>;
```

## Implementation Plan

### Phase 1: Add `Normed` Trait Hierarchy

**Files:** `src/norm.rs` (new), `src/lib.rs`

1. Define base `Normed` trait
2. Define `DegenerateNormed` for PGA
3. Define `IndefiniteNormed` for Minkowski
4. Define `ConformalNormed` for CGA
5. Define `CausalCharacter` enum

### Phase 2: Add Geometry-Specific Wrappers

**Files:** `src/wrappers.rs` (new), `src/lib.rs`

1. Implement `Unit<T>` for Euclidean types
2. Implement `Bulk<T>` for PGA versors
3. Implement `Ideal<T>` for PGA homogeneous types
4. Implement `Proper<T>` for Minkowski timelike types

### Phase 3: Update Codegen - Remove User Constraints

**Files:** `crates/clifford-codegen/src/spec/*.rs`, `algebras/*.toml`

1. Remove `constraints` field from TOML parsing
2. Remove `UserConstraint` from IR
3. Update all TOML files

### Phase 4: Update Codegen - Infer Constraints

**Files:** `crates/clifford-codegen/src/algebra/constraints.rs`, `src/codegen/types.rs`

1. Add `infer_constraint()` function
2. Use Symbolica for symbolic computation
3. Generate `new_checked()` with inferred constraints

### Phase 5: Update Codegen - Simplify Constructors

**Files:** `crates/clifford-codegen/src/codegen/types.rs`, `src/codegen/traits.rs`

1. Remove `new()` constructor generation
2. Generate `Normed` trait implementations
3. Generate type aliases

### Phase 6: Update Codegen - Bespoke Arbitrary

**Files:** `crates/clifford-codegen/src/codegen/traits.rs`

1. Generate `Arbitrary` that solves inferred constraints
2. Generate `Arbitrary` for wrapper types

### Phase 7: Regenerate All Algebras

1. Update TOML files (remove constraints)
2. Regenerate euclidean2, euclidean3, projective3
3. Verify with `cargo fmt && cargo clippy && cargo test && cargo doc --no-deps`

## Success Criteria

1. **No user-configurable constraints** - TOML files contain only type definitions
2. **Constraints inferred** - Codegen derives Study, Plücker, etc. from algebra structure
3. **Simple constructors** - Only `new_unchecked()` and `new_checked()`
4. **`Normed` trait hierarchy** - Base trait plus algebra-specific sub-traits
5. **Geometry-specific wrappers** - `Unit`, `Bulk`, `Ideal`, `Proper`
6. **Type aliases** - `UnitVector`, `BulkMotor`, `IdealPoint`, etc.
7. **Property tests pass** - All constraints verified via proptest
8. **Generalizes across algebras** - Works for Euclidean, PGA, CGA, and Minkowski

## References

- [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
- [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
- [CGA Documentation](https://clifford.readthedocs.io/en/latest/tutorials/cga/index.html)
- [PRD-17.8](prd-17.8-product-normalization.md) - Why products don't normalize

## Sub-PRDs

This PRD is broken into the following sub-PRDs:

| Sub-PRD | Title | Status | Description |
|---------|-------|--------|-------------|
| [PRD-18.1](prd-18.1-normed-traits.md) | Normed Trait Hierarchy | Complete | Define `Normed`, `DegenerateNormed`, `IndefiniteNormed`, `ConformalNormed` traits |
| [PRD-18.2](prd-18.2-wrappers.md) | Wrapper Types | Complete | Define `Unit<T>`, `Bulk<T>`, `Ideal<T>`, `Proper<T>` wrappers |
| [PRD-18.3](prd-18.3-remove-constraints.md) | Remove User Constraints | Complete | Remove constraint sections from TOML parsing |
| [PRD-18.4](prd-18.4-infer-constraints.md) | Infer Constraints | Complete | Add constraint inference from algebra structure |
| [PRD-18.5](prd-18.5-simplify-constructors.md) | Simplify Constructors | Complete | Remove `new()`, keep only `new_unchecked()` and `new_checked()` |
| [PRD-18.6](prd-18.6-arbitrary.md) | Bespoke Arbitrary | Complete | Generate `Arbitrary` implementations that satisfy constraints |
| [PRD-18.7](prd-18.7-regenerate.md) | Regenerate Algebras | Complete | Regenerate all algebras with new constraint system |
| [PRD-18.8](prd-18.8-generate-normed-impls.md) | Generate Normed Impls | Superseded | Already implemented in codegen |
| [PRD-18.9](prd-18.9-generate-wrapper-aliases.md) | Generate Wrapper Aliases | Superseded | Already implemented in codegen |

## Dependencies

- PRD-14 (Code Generator) - Base codegen infrastructure
