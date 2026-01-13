# PRD-18.1: Normed Trait Hierarchy

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Add a `Normed` trait hierarchy to clifford that handles different algebra types (Euclidean, PGA, CGA, Minkowski)

## Overview

Different geometric algebras have different norm semantics:
- **Euclidean**: Standard positive-definite norm
- **PGA**: Degenerate metric with bulk norm and weight norm
- **CGA**: Null vectors, normalization via e∞ coefficient
- **Minkowski**: Indefinite metric (timelike/spacelike/lightlike)

This PRD adds a trait hierarchy to handle all these cases.

## Deliverables

### New File: `src/norm.rs`

#### Base `Normed` Trait

```rust
/// Trait for types that have a well-defined norm.
///
/// In geometric algebra, the norm of an element `u` is computed from
/// `u * rev(u)` (geometric product with reverse). The scalar part gives
/// the squared norm, which may be negative for indefinite metrics.
pub trait Normed {
    type Scalar: Float;

    /// Returns the squared norm of this element.
    ///
    /// **Note**: For indefinite metrics (Minkowski, etc.), this can be negative.
    fn norm_squared(&self) -> Self::Scalar;

    /// Returns the absolute norm of this element.
    ///
    /// Computed as `sqrt(|norm_squared()|)`. Always non-negative.
    fn norm(&self) -> Self::Scalar {
        self.norm_squared().abs().sqrt()
    }

    /// Attempts to normalize this element to unit norm.
    fn try_normalize(&self) -> Option<Self> where Self: Sized;

    /// Normalizes this element to unit norm.
    ///
    /// # Panics
    /// Panics if the norm is too small to normalize.
    fn normalize(&self) -> Self where Self: Sized {
        self.try_normalize().expect("cannot normalize zero element")
    }

    /// Scales this element by a scalar factor.
    fn scale(&self, factor: Self::Scalar) -> Self where Self: Sized;

    /// Returns true if this element can be normalized.
    fn is_normalizable(&self) -> bool where Self: Sized {
        self.norm() > Self::Scalar::epsilon()
    }
}
```

#### `DegenerateNormed` Trait (PGA)

```rust
/// Trait for types in algebras with degenerate metrics (e.g., PGA).
///
/// In PGA, the geometric norm has two components:
/// - **Bulk norm**: `‖u‖◉ = √(u • u)` - magnitude of non-degenerate part
/// - **Weight norm**: `‖u‖○ = √(u ⊙ u)` - magnitude via antidot product
///
/// Reference: [Rigid GA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
pub trait DegenerateNormed: Normed {
    /// Returns the squared bulk norm (from dot product).
    fn bulk_norm_squared(&self) -> Self::Scalar;

    /// Returns the bulk norm.
    fn bulk_norm(&self) -> Self::Scalar {
        self.bulk_norm_squared().abs().sqrt()
    }

    /// Returns the squared weight norm (from antidot product).
    fn weight_norm_squared(&self) -> Self::Scalar;

    /// Returns the weight norm.
    fn weight_norm(&self) -> Self::Scalar {
        self.weight_norm_squared().abs().sqrt()
    }

    /// Unitizes this element by dividing by its weight norm.
    ///
    /// Unitization is distinct from normalization:
    /// - **Normalize**: divide by bulk norm (makes `u • u = 1`)
    /// - **Unitize**: divide by weight norm (makes `u ⊙ u = 1`)
    fn try_unitize(&self) -> Option<Self> where Self: Sized;

    /// Unitizes this element.
    ///
    /// # Panics
    /// Panics if the weight norm is too small.
    fn unitize(&self) -> Self where Self: Sized {
        self.try_unitize().expect("cannot unitize element with zero weight")
    }
}
```

#### `IndefiniteNormed` Trait (Minkowski)

```rust
/// Trait for types in algebras with indefinite metrics (e.g., Minkowski).
///
/// In Minkowski spacetime (signature +---), vectors are classified by
/// the sign of their squared norm:
/// - **Timelike**: `v² > 0` (inside light cone)
/// - **Spacelike**: `v² < 0` (outside light cone)
/// - **Lightlike/Null**: `v² = 0` (on light cone)
///
/// Reference: [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
pub trait IndefiniteNormed: Normed {
    fn is_timelike(&self) -> bool {
        self.norm_squared() > Self::Scalar::zero()
    }

    fn is_spacelike(&self) -> bool {
        self.norm_squared() < Self::Scalar::zero()
    }

    fn is_lightlike(&self) -> bool {
        self.norm_squared().abs() < Self::Scalar::epsilon()
    }

    fn causal_character(&self) -> CausalCharacter {
        let ns = self.norm_squared();
        if ns.abs() < Self::Scalar::epsilon() {
            CausalCharacter::Lightlike
        } else if ns > Self::Scalar::zero() {
            CausalCharacter::Timelike
        } else {
            CausalCharacter::Spacelike
        }
    }
}

/// Classification of vectors in indefinite-metric spaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CausalCharacter {
    Timelike,
    Spacelike,
    Lightlike,
}
```

#### `ConformalNormed` Trait (CGA)

```rust
/// Trait for types in Conformal Geometric Algebra.
///
/// CGA represents points as null vectors (vectors with zero norm).
/// Normalization typically means setting the e∞ coefficient to 1.
///
/// Reference: [CGA docs](https://clifford.readthedocs.io/en/latest/tutorials/cga/index.html)
pub trait ConformalNormed: Normed {
    /// Returns the e∞ (infinity) coefficient of this element.
    fn einf_coefficient(&self) -> Self::Scalar;

    /// Returns true if this element represents a null vector.
    fn is_null(&self) -> bool {
        self.norm_squared().abs() < Self::Scalar::epsilon()
    }

    /// Normalizes by setting the e∞ coefficient to 1.
    fn try_normalize_cga(&self) -> Option<Self> where Self: Sized;
}
```

### Modified Files

- `src/lib.rs` - Add `pub mod norm;` and re-export traits
- `src/prelude.rs` - Include all norm traits

## Testing

```rust
#[test]
fn normed_trait_basic() {
    let v = Vector::new(3.0, 4.0, 0.0);
    assert!(relative_eq!(v.norm(), 5.0, epsilon = 1e-10, max_relative = 1e-10));

    let unit = v.normalize();
    assert!(relative_eq!(unit.norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
}

#[test]
fn degenerate_normed_pga() {
    let motor = Motor::from_rotation_z(0.5);
    assert!(relative_eq!(motor.bulk_norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
}

#[test]
fn indefinite_normed_classification() {
    // With +--- signature
    let timelike = FourVector::new(2.0, 1.0, 0.0, 0.0);  // 4 - 1 = 3 > 0
    let spacelike = FourVector::new(1.0, 2.0, 0.0, 0.0); // 1 - 4 = -3 < 0

    assert!(timelike.is_timelike());
    assert!(spacelike.is_spacelike());
}
```

## Success Criteria

1. `Normed` trait defined with default implementations
2. `DegenerateNormed` trait for PGA types
3. `IndefiniteNormed` trait for Minkowski types
4. `ConformalNormed` trait for CGA types
5. `CausalCharacter` enum defined
6. All traits exported in prelude
7. Documentation with mathematical references
