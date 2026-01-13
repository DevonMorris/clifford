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

### Add Geometric Norm to DegenerateNormed Trait

The `DegenerateNormed` trait should also provide:

```rust
pub trait DegenerateNormed {
    type Scalar: Float;

    /// Bulk norm: magnitude of non-degenerate components.
    /// Returns √(self • self) where • is the dot product.
    fn bulk_norm(&self) -> Self::Scalar;

    /// Weight norm: magnitude of degenerate components.
    /// Returns √(self ⊗ self) where ⊗ is the antidot product.
    fn weight_norm(&self) -> Self::Scalar;

    /// Full geometric norm: bulk + weight as a dual number.
    /// Returns (bulk_norm, weight_norm) tuple.
    ///
    /// The geometric norm is homogeneous: ‖λu‖ = |λ|‖u‖
    fn geometric_norm(&self) -> (Self::Scalar, Self::Scalar) {
        (self.bulk_norm(), self.weight_norm())
    }

    /// Unitized geometric norm: bulk_norm / weight_norm.
    /// Returns None for ideal elements (weight = 0).
    ///
    /// Geometric interpretation:
    /// - Points: distance from origin
    /// - Lines: perpendicular distance from origin
    /// - Planes: perpendicular distance from origin
    /// - Motors: half the distance the origin moves
    fn unitized_geometric_norm(&self) -> Option<Self::Scalar> {
        let w = self.weight_norm();
        if w < Self::Scalar::epsilon() {
            None
        } else {
            Some(self.bulk_norm() / w)
        }
    }

    // ... existing methods: try_normalize, try_unitize, etc.
}
```

This provides the full geometric norm machinery from the Rigid GA literature.

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

1. `Unitized<T>` normalizes to weight_norm = 1 (finite elements)
2. `Ideal<T>` constrains weight_norm = 0 (elements at infinity)
3. Euclidean types get `Unit*` aliases
4. PGA versors get `Bulk*` and `Unitized*` aliases
5. PGA geometric entities get `Unitized*` (finite) and `Ideal*` (infinite) aliases
6. Invalid combinations are not generated
7. Documentation clearly explains the distinction

## Dependencies

- PRD-18.8 (Generate Normed trait impls)
- PRD-18.2 (Wrapper types)

## References

- [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
- [Rigid GA Wiki - Unitization](https://rigidgeometricalgebra.org/wiki/index.php?title=Unitization)
