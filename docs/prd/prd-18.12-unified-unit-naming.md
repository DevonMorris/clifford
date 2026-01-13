# PRD-18.12: Type-Safe Normalization and Constraint Wrappers

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Provide correct wrapper types for PGA normalization and ideal/finite distinction

## Problem Statement

In degenerate algebras like PGA, elements have **two independent norms**:
- **Bulk norm**: Norm of the non-degenerate (Euclidean) part
- **Weight norm**: Norm of the degenerate (projective) part

These lead to different wrapper types with distinct semantics:

| Wrapper | Constraint | Semantics |
|---------|------------|-----------|
| `Bulk<T>` | bulk_norm = 1 | **Normalization** - scales to unit bulk |
| `Unitized<T>` | weight_norm = 1 | **Normalization** - scales to unit weight (finite, standard form) |
| `Ideal<T>` | weight_norm = 0 | **Constraint** - verifies element is ideal (at infinity) |

**Key insight**: `Ideal<T>` is NOT a normalization wrapper. It's a constraint wrapper that guarantees the element has zero weight (is at infinity). You cannot "normalize" an ideal element by weight because weight = 0.

### Current Issues

1. Our current `Ideal<T>` wrapper divides by weight_norm, which is wrong - that should be called `Unitized<T>`
2. We need a separate `Ideal<T>` wrapper that *constrains* weight = 0 (not normalizes)
3. Terminology confusion: "ideal" in PGA means weight = 0, not weight = 1

## Proposed Solution

### 1. Fix Wrapper Semantics

Rename and fix the wrapper types:

```rust
// Normalization wrappers (scale to unit norm)
pub struct Unit<T> { ... }      // Euclidean norm = 1
pub struct Bulk<T> { ... }      // Bulk norm = 1 (PGA normalized)
pub struct Unitized<T> { ... }  // Weight norm = 1 (PGA unitized, finite)

// Constraint wrappers (verify property, don't scale)
pub struct Ideal<T> { ... }     // Weight norm = 0 (at infinity)
pub struct Proper<T> { ... }    // Timelike (Minkowski)
```

### 2. Wrapper Implementations

```rust
impl<T: DegenerateNormed> Unitized<T> {
    /// Unitizes by dividing by weight norm. Returns None for ideal elements.
    pub fn try_new(inner: T) -> Option<Self> {
        let w = inner.weight_norm();
        if w < T::Scalar::epsilon() {
            None  // Can't unitize ideal elements (weight = 0)
        } else {
            Some(Self { inner: inner.scale(T::Scalar::one() / w) })
        }
    }
}

impl<T: DegenerateNormed> Ideal<T> {
    /// Verifies element is ideal (weight ≈ 0). Does NOT normalize.
    pub fn try_new(inner: T) -> Option<Self> {
        if inner.weight_norm() < T::Scalar::epsilon() {
            Some(Self { inner })
        } else {
            None  // Not ideal - has finite weight
        }
    }
}
```

### 3. Type Aliases by Algebra

#### Euclidean 3D

```rust
pub type UnitVector<T> = crate::wrappers::Unit<Vector<T>>;
pub type UnitBivector<T> = crate::wrappers::Unit<Bivector<T>>;
pub type UnitRotor<T> = crate::wrappers::Unit<Rotor<T>>;
```

#### Projective 3D (PGA)

```rust
// Versors - bulk normalization for rigid transforms
pub type BulkMotor<T> = crate::wrappers::Bulk<Motor<T>>;      // Normalized motor
pub type UnitizedMotor<T> = crate::wrappers::Unitized<Motor<T>>; // Unitized motor

pub type BulkFlector<T> = crate::wrappers::Bulk<Flector<T>>;

// Points - finite vs ideal distinction
pub type UnitizedPoint<T> = crate::wrappers::Unitized<Point<T>>; // Finite point (w=1)
pub type IdealPoint<T> = crate::wrappers::Ideal<Point<T>>;       // Point at infinity (w=0)

// Lines
pub type BulkLine<T> = crate::wrappers::Bulk<Line<T>>;        // Direction normalized
pub type UnitizedLine<T> = crate::wrappers::Unitized<Line<T>>; // Moment normalized

// Planes
pub type UnitizedPlane<T> = crate::wrappers::Unitized<Plane<T>>; // Standard form
pub type IdealPlane<T> = crate::wrappers::Ideal<Plane<T>>;       // Plane at infinity
```

## Deliverables

### 1. Update Wrapper Inference Logic

**File:** `crates/clifford-codegen/src/codegen/types.rs`

```rust
/// Determines which wrappers are applicable for a type based on signature.
fn applicable_wrappers(sig: &Signature, type_spec: &TypeSpec) -> Vec<WrapperKind> {
    let (p, q, r) = (sig.p, sig.q, sig.r);
    let mut wrappers = Vec::new();

    match (p, q, r) {
        // Euclidean: Unit<T> only
        (_, 0, 0) => {
            wrappers.push(WrapperKind::Unit);
        }

        // Projective (PGA): Both Bulk<T> and Ideal<T>
        (_, 0, r) if r > 0 => {
            wrappers.push(WrapperKind::Bulk);
            wrappers.push(WrapperKind::Ideal);
        }

        // Minkowski: Unit<T> and Proper<T>
        (p, q, 0) if p > 0 && q > 0 => {
            wrappers.push(WrapperKind::Unit);
            if type_spec.grades == vec![1] {
                wrappers.push(WrapperKind::Proper);
            }
        }

        // Mixed degenerate + indefinite
        _ => {
            wrappers.push(WrapperKind::Bulk);
            wrappers.push(WrapperKind::Ideal);
        }
    }

    wrappers
}

#[derive(Debug, Clone, Copy)]
enum WrapperKind {
    Unit,   // Euclidean norm = 1
    Bulk,   // Bulk norm = 1 (non-degenerate part)
    Ideal,  // Weight norm = 1 (degenerate part)
    Proper, // Timelike, normalized
}
```

### 2. Generate All Applicable Aliases

```rust
fn generate_wrapper_aliases(type_spec: &TypeSpec, sig: &Signature) -> TokenStream {
    let wrappers = applicable_wrappers(sig, type_spec);
    let type_name = format_ident!("{}", type_spec.name);
    let mut aliases = TokenStream::new();

    for wrapper in wrappers {
        let (prefix, wrapper_path, doc) = match wrapper {
            WrapperKind::Unit => (
                "Unit",
                quote!(crate::wrappers::Unit),
                "Euclidean norm = 1",
            ),
            WrapperKind::Bulk => (
                "Bulk",
                quote!(crate::wrappers::Bulk),
                "bulk norm = 1 (non-degenerate part normalized)",
            ),
            WrapperKind::Ideal => (
                "Ideal",
                quote!(crate::wrappers::Ideal),
                "weight norm = 1 (degenerate part normalized)",
            ),
            WrapperKind::Proper => (
                "Proper",
                quote!(crate::wrappers::Proper),
                "timelike and normalized",
            ),
        };

        let alias_name = format_ident!("{}{}", prefix, type_spec.name);
        aliases.extend(quote! {
            #[doc = concat!("A ", stringify!(#type_name), " with ", #doc, ".")]
            pub type #alias_name<T> = #wrapper_path<#type_name<T>>;
        });
    }

    aliases
}
```

### 3. Expected Output for PGA Motor

```rust
/// A Motor with bulk norm = 1 (non-degenerate part normalized).
///
/// Use this for rigid body transformations where the rotation
/// component should have unit magnitude.
pub type BulkMotor<T> = crate::wrappers::Bulk<Motor<T>>;

/// A Motor with weight norm = 1 (degenerate part normalized).
///
/// Use this when the translation component should have unit magnitude.
pub type IdealMotor<T> = crate::wrappers::Ideal<Motor<T>>;
```

### 4. "Unit" as Canonical Alias (Optional)

For user convenience, we could also generate a `Unit*` alias that maps to the "most common" normalization:

```rust
// In PGA, UnitMotor = BulkMotor (rotation-normalized is more common)
pub type UnitMotor<T> = BulkMotor<T>;

// In Euclidean, UnitRotor = Unit<Rotor> (only one option)
pub type UnitRotor<T> = crate::wrappers::Unit<Rotor<T>>;
```

This gives users a consistent `Unit*` naming while still exposing the explicit `Bulk*`/`Ideal*` options.

## PGA Terminology Correspondence

The Rigid GA literature uses specific terms:

| Our Wrapper | PGA Term | Constraint | Reference |
|-------------|----------|------------|-----------|
| `Bulk<T>` | **Normalized** | bulk_norm = 1 | [Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm) |
| `Unitized<T>` | **Unitized** | weight_norm = 1 | [Unitization](https://rigidgeometricalgebra.org/wiki/index.php?title=Unitization) |
| `Ideal<T>` | **Ideal** | weight_norm = 0 | Points/planes at infinity |

**Key distinction**:
- `Unitized<T>` is a **normalization** (scales weight to 1)
- `Ideal<T>` is a **constraint** (verifies weight is 0, no scaling)

## Semantic Meaning by Type

### Motors (Versors)

| Alias | PGA Term | Constraint | Geometric Meaning |
|-------|----------|------------|-------------------|
| `BulkMotor` | Normalized | bulk = 1 | Standard rigid transformation |
| `UnitizedMotor` | Unitized | weight = 1 | Translation-unitized motor |

### Points

| Alias | Constraint | Geometric Meaning |
|-------|------------|-------------------|
| `UnitizedPoint` | weight = 1 | Finite point in standard form (w=1) |
| `IdealPoint` | weight = 0 | Point at infinity (direction only) |

Note: `BulkPoint` doesn't make sense - points in PGA have only weight components (grade 1 involves e₀).

### Lines

| Alias | Constraint | Geometric Meaning |
|-------|------------|-------------------|
| `BulkLine` | bulk = 1 | Unit direction vector |
| `UnitizedLine` | weight = 1 | Standard Plücker coordinates |

### Planes

| Alias | Constraint | Geometric Meaning |
|-------|------------|-------------------|
| `UnitizedPlane` | weight = 1 | Standard form (unit normal) |
| `IdealPlane` | weight = 0 | Plane at infinity |

## Implementation Notes

### Not All Combinations Make Sense

Some wrapper combinations are mathematically meaningless:
- `BulkPoint` - Points have no bulk components in PGA (all components involve e₀)
- `IdealMotor` - Motors can't be "ideal" (they always have bulk from the scalar/bivector parts)
- `UnitizedRotor` - Rotors in Euclidean have no weight components

The codegen should skip generating aliases for types where the wrapper doesn't apply:

```rust
fn wrapper_applies(wrapper: WrapperKind, type_spec: &TypeSpec, sig: &Signature) -> bool {
    match wrapper {
        WrapperKind::Unit => sig.is_euclidean(),
        WrapperKind::Bulk => type_spec.has_bulk_components(sig),
        WrapperKind::Unitized => type_spec.has_weight_components(sig),
        WrapperKind::Ideal => type_spec.can_be_ideal(sig),  // Has weight AND can be zero
        WrapperKind::Proper => sig.is_minkowski() && type_spec.is_vector(),
    }
}
```

### Trait Requirements

The wrapper types require specific traits:
- `Unit<T>` requires `T: Normed`
- `Bulk<T>` requires `T: DegenerateNormed`
- `Unitized<T>` requires `T: DegenerateNormed`
- `Ideal<T>` requires `T: DegenerateNormed`
- `Proper<T>` requires `T: IndefiniteNormed`

PRD-18.8 handles generating these trait implementations.

### Trait Hierarchy for All Signature Types

The trait hierarchy must handle all combinations of signature components:

| Signature | Bulk Metric | Weight Metric | Example |
|-----------|-------------|---------------|---------|
| Cl(n,0,0) | Positive-definite | N/A | Euclidean |
| Cl(p,q,0) | Indefinite | N/A | Minkowski |
| Cl(n,0,r) | Positive-definite | Positive-definite* | PGA |
| Cl(p,q,r) | Indefinite | Possibly indefinite | Conformal Minkowski |

*The weight metric (anti-inner product) can also be indefinite in some algebras.

**Key insight:** Use `_squared()` methods as the fundamental operations, since they work for all metric types:

```rust
/// Base trait for types with a degenerate metric component.
/// Works for ALL degenerate algebras, including those with indefinite bulk.
pub trait DegenerateNormed {
    type Scalar: Float;

    /// Bulk norm squared: (self • self) where • is the inner product.
    /// May be positive, negative, or zero depending on the bulk metric.
    fn bulk_norm_squared(&self) -> Self::Scalar;

    /// Weight norm squared: (self ⊛ self) where ⊛ is the anti-inner product.
    /// May be negative if the anti-metric is indefinite.
    fn weight_norm_squared(&self) -> Self::Scalar;

    /// Bulk magnitude: √|bulk_norm_squared|, always non-negative.
    fn bulk_magnitude(&self) -> Self::Scalar {
        self.bulk_norm_squared().abs().sqrt()
    }

    /// Weight magnitude: √|weight_norm_squared|, always non-negative.
    fn weight_magnitude(&self) -> Self::Scalar {
        self.weight_norm_squared().abs().sqrt()
    }

    /// Full geometric norm as magnitudes: (bulk_magnitude, weight_magnitude).
    fn geometric_norm(&self) -> (Self::Scalar, Self::Scalar) {
        (self.bulk_magnitude(), self.weight_magnitude())
    }

    /// Unitized geometric norm: bulk_magnitude / weight_magnitude.
    /// Returns None for ideal elements (weight ≈ 0).
    fn unitized_geometric_norm(&self) -> Option<Self::Scalar> {
        let w = self.weight_magnitude();
        if w < Self::Scalar::epsilon() {
            None
        } else {
            Some(self.bulk_magnitude() / w)
        }
    }
}

/// Extension trait for algebras with positive-definite bulk (e.g., PGA).
/// Provides `bulk_norm()` which is only valid when bulk·bulk ≥ 0.
pub trait PositiveDefiniteBulk: DegenerateNormed {
    /// Bulk norm: √(bulk·bulk). Only valid for positive-definite bulk metrics.
    /// Panics or returns NaN if bulk_norm_squared < 0.
    fn bulk_norm(&self) -> Self::Scalar {
        let sq = self.bulk_norm_squared();
        debug_assert!(sq >= Self::Scalar::zero(), "bulk_norm requires positive-definite bulk");
        sq.sqrt()
    }
}

/// Extension trait for algebras with positive-definite weight (most PGA variants).
pub trait PositiveDefiniteWeight: DegenerateNormed {
    /// Weight norm: √(weight·weight). Only valid for positive-definite weight.
    fn weight_norm(&self) -> Self::Scalar {
        let sq = self.weight_norm_squared();
        debug_assert!(sq >= Self::Scalar::zero(), "weight_norm requires positive-definite weight");
        sq.sqrt()
    }
}

/// Causal character of the bulk part (for indefinite bulk metrics).
pub trait BulkCausalCharacter: DegenerateNormed {
    fn is_bulk_timelike(&self) -> bool {
        self.bulk_norm_squared() > Self::Scalar::epsilon()
    }

    fn is_bulk_spacelike(&self) -> bool {
        self.bulk_norm_squared() < -Self::Scalar::epsilon()
    }

    fn is_bulk_null(&self) -> bool {
        self.bulk_norm_squared().abs() < Self::Scalar::epsilon()
    }
}
```

### Wrapper Trait Requirements (Revised)

```rust
// Euclidean (positive-definite, no degenerate)
impl<T: Normed> Unit<T> { ... }

// PGA with positive-definite bulk (e.g., Cl(3,0,1))
impl<T: PositiveDefiniteBulk> Bulk<T> { ... }

// Any degenerate algebra
impl<T: DegenerateNormed> Unitized<T> { ... }  // Uses weight_magnitude()
impl<T: DegenerateNormed> Ideal<T> { ... }     // Checks weight_magnitude() ≈ 0

// Indefinite algebras (no degenerate)
impl<T: IndefiniteNormed> Proper<T> { ... }
impl<T: IndefiniteNormed> Spacelike<T> { ... }
impl<T: IndefiniteNormed> Null<T> { ... }

// Mixed indefinite + degenerate (e.g., Cl(3,1,1))
impl<T: BulkCausalCharacter> ProperBulk<T> { ... }     // Timelike bulk, normalized
impl<T: BulkCausalCharacter> SpacelikeBulk<T> { ... }  // Spacelike bulk, normalized
impl<T: BulkCausalCharacter> NullBulk<T> { ... }       // Null bulk (constraint)
```

### Type Safety: Mutually Exclusive Trait Families

**CRITICAL:** Types must implement exactly ONE trait family based on their algebra's signature. This prevents misuse of inappropriate wrappers.

```rust
/// The norm trait families are MUTUALLY EXCLUSIVE.
/// Codegen selects exactly one based on signature (p, q, r).
enum NormTraitFamily {
    /// Cl(n,0,0): Euclidean - implements `Normed` only
    Euclidean,

    /// Cl(n,0,r): Degenerate with positive-definite bulk (PGA)
    /// Implements: `DegenerateNormed`, `PositiveDefiniteBulk`, `PositiveDefiniteWeight`
    DegeneratePositive,

    /// Cl(p,q,0): Indefinite, no degenerate (Minkowski)
    /// Implements: `IndefiniteNormed`, `CausalCharacter`
    Indefinite,

    /// Cl(p,q,r): Degenerate with indefinite bulk (Conformal Minkowski)
    /// Implements: `DegenerateNormed`, `BulkCausalCharacter`
    DegenerateIndefinite,
}

fn select_norm_traits(sig: &Signature) -> NormTraitFamily {
    match (sig.p, sig.q, sig.r) {
        (_, 0, 0) => NormTraitFamily::Euclidean,
        (_, 0, r) if r > 0 => NormTraitFamily::DegeneratePositive,
        (p, q, 0) if p > 0 && q > 0 => NormTraitFamily::Indefinite,
        (p, q, r) if q > 0 && r > 0 => NormTraitFamily::DegenerateIndefinite,
        _ => panic!("unsupported signature combination"),
    }
}
```

**What this prevents:**

| Algebra | Forbidden Wrappers | Why |
|---------|-------------------|-----|
| Cl(3,0,0) Euclidean | `Bulk<T>`, `Proper<T>` | No degenerate, no indefinite |
| Cl(3,0,1) PGA | `Unit<T>`, `Proper<T>` | Has degenerate (use Bulk), no indefinite |
| Cl(3,1,0) Minkowski | `Unit<T>`, `Bulk<T>` | Indefinite (use Proper), no degenerate |
| Cl(3,1,1) Mixed | `Unit<T>`, `Bulk<T>`, `Proper<T>` | Use `ProperBulk<T>` instead |

**Compile-time enforcement:**

```rust
// In Cl(3,0,1) - PGA
// Codegen generates:
impl<T: Float> DegenerateNormed for Motor<T> { ... }
impl<T: Float> PositiveDefiniteBulk for Motor<T> { ... }
// Does NOT generate: impl Normed, impl IndefiniteNormed, impl BulkCausalCharacter

// This means:
let m: Motor<f64> = ...;
let _ = Bulk::new(m);        // ✓ Compiles - PositiveDefiniteBulk implemented
let _ = Unit::new(m);        // ✗ Error - Normed not implemented
let _ = Proper::new(m);      // ✗ Error - IndefiniteNormed not implemented
let _ = ProperBulk::new(m);  // ✗ Error - BulkCausalCharacter not implemented
```

```rust
// In Cl(3,1,1) - Conformal Minkowski
// Codegen generates:
impl<T: Float> DegenerateNormed for Motor<T> { ... }
impl<T: Float> BulkCausalCharacter for Motor<T> { ... }
// Does NOT generate: impl Normed, impl PositiveDefiniteBulk, impl IndefiniteNormed

// This means:
let m: Motor<f64> = ...;
let _ = ProperBulk::try_new(m);  // ✓ Compiles - BulkCausalCharacter implemented
let _ = Bulk::new(m);            // ✗ Error - PositiveDefiniteBulk not implemented
let _ = Unit::new(m);            // ✗ Error - Normed not implemented
let _ = Proper::new(m);          // ✗ Error - IndefiniteNormed not implemented
```

This ensures **compile-time rejection** of semantically incorrect normalizations.

### Example: Cl(3,1,1) - Conformal Minkowski

```rust
// In Cl(3,1,1):
// - e₁² = e₂² = e₃² = +1 (spatial)
// - e₄² = -1 (temporal)
// - e₅² = 0 (conformal/degenerate)

// A vector could be:
let v = Vector::new(1.0, 0.0, 0.0, 2.0, 0.5);  // (x, y, z, t, w)

// Bulk part: (1, 0, 0, 2) in Cl(3,1,0)
// bulk·bulk = 1 + 0 + 0 - 4 = -3 (SPACELIKE!)

// Weight part: involves e₅
// weight·weight uses anti-metric

// With new traits:
assert!(v.is_bulk_spacelike());  // bulk_norm_squared < 0
assert!(v.bulk_magnitude() > 0.0);  // √|-3| = √3

// Can normalize by bulk magnitude (not bulk norm!)
let normalized = SpacelikeBulk::try_new(v)?;
```

This design ensures Cl(3,1,1) and similar algebras work correctly.

### Ideal vs Unitized - The Key Difference

```rust
// Unitized: NORMALIZES by dividing by weight
let finite_point = Unitized::try_new(point)?;  // Scales so w = 1
assert!((finite_point.weight_norm() - 1.0).abs() < 1e-10);

// Ideal: CONSTRAINS to zero weight (no scaling)
let infinite_point = Ideal::try_new(direction)?;  // Verifies w ≈ 0
assert!(infinite_point.weight_norm() < 1e-10);
```

## Testing

```rust
#[test]
fn pga_motor_normalizations() {
    use crate::specialized::projective::dim3::{Motor, BulkMotor, UnitizedMotor};
    use crate::wrappers::{Bulk, Unitized};

    let m = Motor::from_translation(1.0, 2.0, 3.0);

    // Both normalizations work
    let bulk: BulkMotor<f64> = Bulk::new_normalize(m);
    let unitized: UnitizedMotor<f64> = Unitized::new_normalize(m);

    // They produce different results
    assert!((bulk.bulk_norm() - 1.0).abs() < 1e-10);
    assert!((unitized.weight_norm() - 1.0).abs() < 1e-10);
}

#[test]
fn pga_point_finite_vs_ideal() {
    use crate::specialized::projective::dim3::{Point, UnitizedPoint, IdealPoint};
    use crate::wrappers::{Unitized, Ideal};

    // Finite point can be unitized
    let finite = Point::from_cartesian(1.0, 2.0, 3.0);  // w = 1
    let unitized: UnitizedPoint<f64> = Unitized::try_new(finite).unwrap();
    assert!((unitized.weight_norm() - 1.0).abs() < 1e-10);

    // Ideal point (direction) cannot be unitized
    let direction = Point::new(1.0, 0.0, 0.0, 0.0);  // w = 0
    assert!(Unitized::<Point<f64>>::try_new(direction).is_none());

    // But it CAN be wrapped as Ideal
    let ideal: IdealPoint<f64> = Ideal::try_new(direction).unwrap();
    assert!(ideal.weight_norm() < 1e-10);
}

#[test]
fn euclidean_rotor_has_unit_only() {
    use crate::specialized::euclidean::dim3::{Rotor, UnitRotor, Vector};
    use crate::wrappers::Unit;

    let r = Rotor::from_axis_angle(&Vector::unit_z(), 0.5);
    let unit: UnitRotor<f64> = Unit::new_normalize(r);

    assert!((unit.norm() - 1.0).abs() < 1e-10);
}
```

## Migration

### Breaking Changes

**Rename `Ideal<T>` to `Unitized<T>`** for the weight-normalization wrapper.

This is a breaking change if users are using `Ideal<T>` directly. The type aliases (e.g., `IdealPoint`) will need to be updated.

### Migration Path

```rust
// Before (incorrect semantics)
let p: Ideal<Point<f64>> = Ideal::new_normalize(point);

// After (correct semantics)
let p: Unitized<Point<f64>> = Unitized::new_normalize(point);  // Finite point
let p: Ideal<Point<f64>> = Ideal::try_new(direction).unwrap(); // Point at infinity
```

### New Types

- `Unitized<T>` - new wrapper for weight_norm = 1
- `UnitizedPoint<T>`, `UnitizedMotor<T>`, etc. - new aliases

### Renamed Aliases

- `IdealPoint<T>` now means weight = 0 (point at infinity), not weight = 1

## Success Criteria

### Degenerate Algebras Cl(n,0,r) - PGA

1. `Unitized<T>` normalizes to weight_magnitude = 1 (finite elements)
2. `Ideal<T>` constrains weight_magnitude = 0 (elements at infinity)
3. PGA versors get `Bulk*` and `Unitized*` aliases
4. PGA geometric entities get `Unitized*` (finite) and `Ideal*` (infinite) aliases
5. `PositiveDefiniteBulk` trait provides `bulk_norm()` for PGA types

### Indefinite Algebras Cl(p,q,0) - Minkowski

6. `Proper<T>` normalizes timelike elements to |v·v| = 1
7. `Spacelike<T>` normalizes spacelike elements to |v·v| = 1
8. `Null<T>` constrains elements to v·v ≈ 0 (no normalization)
9. `IndefiniteNormed` trait provides `bulk_norm_squared()` for signed norms
10. `CausalCharacter` trait provides `is_timelike()`, `is_spacelike()`, `is_lightlike()`

### Mixed Algebras Cl(p,q,r) - Conformal Minkowski

11. `DegenerateNormed` uses `bulk_norm_squared()` as primitive (handles indefinite bulk)
12. `BulkCausalCharacter` trait provides `is_bulk_timelike()`, `is_bulk_spacelike()`, `is_bulk_null()`
13. `ProperBulk<T>`, `SpacelikeBulk<T>`, `NullBulk<T>` wrappers for mixed algebras
14. Both bulk and weight can have indefinite metrics (use `_squared()` primitives)

### All Algebras

15. Euclidean types get `Unit*` aliases
16. Invalid combinations are not generated (e.g., no `BulkPoint` in PGA)
17. Documentation clearly explains the distinction between normalization and constraint wrappers
18. All behavior derived from signature (p, q, r), not hardcoded algebra names

### Type Safety (Critical)

19. **Mutually exclusive trait families** - each type implements exactly ONE norm trait family
20. `Unit<T>` only compiles for Euclidean types (Cl(n,0,0))
21. `Bulk<T>` only compiles for positive-definite degenerate types (Cl(n,0,r))
22. `Proper<T>` only compiles for indefinite non-degenerate types (Cl(p,q,0))
23. `ProperBulk<T>` only compiles for mixed indefinite+degenerate types (Cl(p,q,r))
24. Attempting wrong wrapper produces compile-time error, not runtime panic

## Dependencies

- PRD-18.8 (Generate Normed trait impls)
- PRD-18.2 (Wrapper types)
- New: `IndefiniteNormed` trait definition
- New: `CausalCharacter` trait definition
- New: `Proper<T>`, `Spacelike<T>`, `Null<T>` wrapper implementations

## Indefinite Signatures (Minkowski, STA)

For algebras with negative signature components (q > 0 in Cl(p,q,r)), the norm squared can be negative. This requires additional wrappers and terminology consistent with physics literature.

### The Problem with Indefinite Metrics

In Minkowski space Cl(3,1,0) or Cl(1,3,0), the "norm" of a vector depends on its causal character:

```rust
// In Cl(3,1,0) with metric (+,+,+,-):
let v = Vector::new(1.0, 0.0, 0.0, 2.0);  // (x, y, z, t)

// v · v = x² + y² + z² - t² = 1 + 0 + 0 - 4 = -3
// The norm squared is NEGATIVE, so √(v·v) would be imaginary!
```

Physics literature distinguishes three cases:

| Causal Character | Condition | Physics Term | Example |
|-----------------|-----------|--------------|---------|
| **Timelike** | v·v > 0 | "Proper" | 4-velocity, worldlines |
| **Spacelike** | v·v < 0 | (no standard term) | Spatial separations |
| **Lightlike** | v·v = 0 | "Null" | Photon worldlines |

Note: The sign convention varies by author. Some use (+,-,-,-) where timelike is *negative*. We'll support both via the signature tuple.

### IndefiniteNormed Trait

For indefinite metrics, we provide `bulk_norm_squared()` instead of assuming a real square root:

```rust
pub trait IndefiniteNormed {
    type Scalar: Float;

    /// Returns the norm squared, which may be positive, negative, or zero.
    /// - Positive: timelike (or spacelike, depending on convention)
    /// - Negative: spacelike (or timelike, depending on convention)
    /// - Zero: lightlike/null
    fn bulk_norm_squared(&self) -> Self::Scalar;

    /// Returns the magnitude, always non-negative.
    /// This is |√|v·v|| = √|v·v|
    fn bulk_magnitude(&self) -> Self::Scalar {
        self.bulk_norm_squared().abs().sqrt()
    }

    /// True if v·v > ε (positive-definite direction)
    fn is_positive_definite(&self) -> bool {
        self.bulk_norm_squared() > Self::Scalar::epsilon()
    }

    /// True if v·v < -ε (negative-definite direction)
    fn is_negative_definite(&self) -> bool {
        self.bulk_norm_squared() < -Self::Scalar::epsilon()
    }

    /// True if |v·v| < ε (null/lightlike)
    fn is_null(&self) -> bool {
        self.bulk_norm_squared().abs() < Self::Scalar::epsilon()
    }
}
```

### Causal Character Helper Trait

For Minkowski-specific applications, we provide physics-friendly terminology:

```rust
/// Trait for types with Minkowski-style causal structure.
/// Assumes the conventional physics metric where timelike is v·v > 0.
pub trait CausalCharacter: IndefiniteNormed {
    /// True if this element is timelike (v·v > 0).
    /// Timelike vectors point "inside" the light cone.
    fn is_timelike(&self) -> bool {
        self.is_positive_definite()
    }

    /// True if this element is spacelike (v·v < 0).
    /// Spacelike vectors point "outside" the light cone.
    fn is_spacelike(&self) -> bool {
        self.is_negative_definite()
    }

    /// True if this element is lightlike/null (v·v ≈ 0).
    /// Lightlike vectors point along the light cone.
    fn is_lightlike(&self) -> bool {
        self.is_null()
    }
}
```

### Wrappers for Indefinite Metrics

#### Proper<T> - Timelike and Normalized

The term "proper" comes from special relativity: proper time, proper velocity, proper acceleration. A "proper" vector is timelike and normalized.

```rust
/// A timelike element with unit magnitude.
///
/// "Proper" is the standard physics term for normalized timelike quantities:
/// - Proper time τ
/// - Proper velocity (4-velocity) with u·u = c² = 1 (natural units)
/// - Proper acceleration
pub struct Proper<T> {
    inner: T,
}

impl<T: IndefiniteNormed> Proper<T> {
    /// Creates a proper (normalized timelike) element.
    /// Returns None if the element is not timelike.
    pub fn try_new(inner: T) -> Option<Self> {
        let norm_sq = inner.bulk_norm_squared();
        if norm_sq > T::Scalar::epsilon() {
            // Timelike - can normalize
            let scale = T::Scalar::one() / norm_sq.sqrt();
            Some(Self { inner: inner.scale(scale) })
        } else {
            None  // Spacelike or null - can't be "proper"
        }
    }

    /// Returns the inner value, guaranteed to be timelike with unit norm.
    pub fn into_inner(self) -> T {
        self.inner
    }
}
```

#### Spacelike<T> - Spacelike and Normalized

```rust
/// A spacelike element with unit magnitude.
pub struct Spacelike<T> {
    inner: T,
}

impl<T: IndefiniteNormed> Spacelike<T> {
    /// Creates a normalized spacelike element.
    /// Returns None if the element is not spacelike.
    pub fn try_new(inner: T) -> Option<Self> {
        let norm_sq = inner.bulk_norm_squared();
        if norm_sq < -T::Scalar::epsilon() {
            // Spacelike - normalize by |v·v|
            let scale = T::Scalar::one() / (-norm_sq).sqrt();
            Some(Self { inner: inner.scale(scale) })
        } else {
            None  // Timelike or null
        }
    }
}
```

#### Null<T> - Constraint Wrapper for Lightlike

Like `Ideal<T>` for PGA, `Null<T>` is a **constraint** wrapper, not a normalization:

```rust
/// A lightlike/null element (v·v ≈ 0).
///
/// This is a constraint wrapper: it verifies the element is null,
/// but does NOT normalize (you can't normalize a null vector).
pub struct Null<T> {
    inner: T,
}

impl<T: IndefiniteNormed> Null<T> {
    /// Verifies the element is null. Does NOT normalize.
    pub fn try_new(inner: T) -> Option<Self> {
        if inner.is_null() {
            Some(Self { inner })
        } else {
            None  // Not null
        }
    }
}
```

### Type Aliases for Minkowski

```rust
// Spacetime Algebra (STA) in Cl(1,3,0) or Cl(3,1,0)
pub type ProperVector<T> = crate::wrappers::Proper<Vector<T>>;      // 4-velocity
pub type SpacelikeVector<T> = crate::wrappers::Spacelike<Vector<T>>; // Spatial direction
pub type NullVector<T> = crate::wrappers::Null<Vector<T>>;          // Photon direction

pub type ProperBivector<T> = crate::wrappers::Proper<Bivector<T>>;  // Electromagnetic field tensor?
```

### Sign Convention Handling

Different physics communities use different metric signatures:

| Convention | Signature | Timelike | Spacelike |
|------------|-----------|----------|-----------|
| Mostly-plus | (+,+,+,-) | v·v > 0 | v·v < 0 |
| Mostly-minus | (+,-,-,-) | v·v > 0 | v·v < 0 |
| West coast | (-,+,+,+) | v·v < 0 | v·v > 0 |

Rather than hardcoding conventions, we derive behavior from the signature:

```rust
impl<T: Float> CausalCharacter for Vector<T> {
    // The implementation checks which basis indices contribute positively
    // to the norm squared, based on the algebra's signature (p, q, r).
}
```

This allows the same `is_timelike()` method to work correctly regardless of which convention the user chose when defining their algebra.

### Combined Degenerate + Indefinite (Advanced)

Some algebras have both degenerate AND indefinite components (e.g., conformal space). These would need both:
- `DegenerateNormed` for bulk/weight decomposition
- `IndefiniteNormed` for causal character of the bulk

```rust
pub trait DegenerateIndefiniteNormed: DegenerateNormed + IndefiniteNormed {
    /// Bulk norm squared (may be negative for indefinite bulk).
    fn bulk_norm_squared(&self) -> Self::Scalar;

    /// Weight norm (always non-negative in current algebras).
    fn weight_norm(&self) -> Self::Scalar;

    /// Causal character of the bulk part.
    fn bulk_is_timelike(&self) -> bool {
        self.bulk_norm_squared() > Self::Scalar::epsilon()
    }

    fn bulk_is_spacelike(&self) -> bool {
        self.bulk_norm_squared() < -Self::Scalar::epsilon()
    }

    fn bulk_is_null(&self) -> bool {
        self.bulk_norm_squared().abs() < Self::Scalar::epsilon()
    }
}
```

### Testing Indefinite Signatures

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minkowski_causal_character() {
        // Cl(3,1,0) with metric (+,+,+,-)
        let timelike = Vector::new(0.0, 0.0, 0.0, 1.0);  // Pure time
        let spacelike = Vector::new(1.0, 0.0, 0.0, 0.0); // Pure space
        let lightlike = Vector::new(1.0, 0.0, 0.0, 1.0); // On light cone

        // v·v for timelike: 0 + 0 + 0 - 1 = -1 (wait, this depends on convention!)
        // Let's assume (+,+,+,-) means e₄² = -1
        assert!(timelike.is_timelike());  // Convention-dependent
        assert!(spacelike.is_spacelike());
        assert!(lightlike.is_lightlike());
    }

    #[test]
    fn proper_vector_normalization() {
        let v = Vector::new(0.0, 0.0, 0.0, 2.0);  // Timelike
        let proper = Proper::try_new(v).expect("should be timelike");

        // After normalization: |v·v| = 1
        assert!((proper.bulk_norm_squared().abs() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn null_vector_constraint() {
        let photon = Vector::new(1.0, 0.0, 0.0, 1.0);  // x = t, on light cone
        let null = Null::try_new(photon).expect("should be null");

        // No normalization applied - just constrained
        assert!(null.is_null());
    }

    #[test]
    fn spacelike_cannot_be_proper() {
        let spatial = Vector::new(1.0, 0.0, 0.0, 0.0);
        assert!(Proper::try_new(spatial).is_none());
    }
}
```

### Wrapper Summary Table

#### Euclidean Algebras Cl(n,0,0)

| Wrapper | Trait Required | Constraint | Normalization? | Use Case |
|---------|----------------|------------|----------------|----------|
| `Unit<T>` | `Normed` | ‖v‖ = 1 | Yes | Unit vectors, rotors |

#### Degenerate Algebras Cl(n,0,r) - PGA

| Wrapper | Trait Required | Constraint | Normalization? | Use Case |
|---------|----------------|------------|----------------|----------|
| `Bulk<T>` | `PositiveDefiniteBulk` | bulk = 1 | Yes | PGA versors |
| `Unitized<T>` | `DegenerateNormed` | weight = 1 | Yes | PGA finite elements |
| `Ideal<T>` | `DegenerateNormed` | weight = 0 | No | PGA at infinity |

#### Indefinite Algebras Cl(p,q,0) - Minkowski

| Wrapper | Trait Required | Constraint | Normalization? | Use Case |
|---------|----------------|------------|----------------|----------|
| `Proper<T>` | `IndefiniteNormed` | timelike, \|v·v\| = 1 | Yes | 4-velocity |
| `Spacelike<T>` | `IndefiniteNormed` | spacelike, \|v·v\| = 1 | Yes | Spatial directions |
| `Null<T>` | `IndefiniteNormed` | v·v = 0 | No | Photon worldlines |

#### Mixed Algebras Cl(p,q,r) - Conformal Minkowski

| Wrapper | Trait Required | Constraint | Normalization? | Use Case |
|---------|----------------|------------|----------------|----------|
| `ProperBulk<T>` | `BulkCausalCharacter` | bulk timelike, \|bulk²\| = 1 | Yes | Timelike + degenerate |
| `SpacelikeBulk<T>` | `BulkCausalCharacter` | bulk spacelike, \|bulk²\| = 1 | Yes | Spacelike + degenerate |
| `NullBulk<T>` | `BulkCausalCharacter` | bulk² = 0 | No | Null bulk (constraint) |
| `Unitized<T>` | `DegenerateNormed` | weight = 1 | Yes | Finite elements |
| `Ideal<T>` | `DegenerateNormed` | weight = 0 | No | At infinity |

## References

- [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
- [Rigid GA Wiki - Unitization](https://rigidgeometricalgebra.org/wiki/index.php?title=Unitization)
- [Wikipedia - Minkowski space](https://en.wikipedia.org/wiki/Minkowski_space) - Causal structure and conventions
- [Wikipedia - Four-velocity](https://en.wikipedia.org/wiki/Four-velocity) - "Proper" terminology in relativity
