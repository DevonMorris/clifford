# PRD-18.2: Geometry-Specific Wrapper Types

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Add unambiguous wrapper types for normalized entities in different geometric algebras

## Overview

Each geometry type needs its own wrapper with a clear name describing what it guarantees:

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

## Deliverables

### New File: `src/wrappers.rs`

#### `Unit<T>` - Euclidean Norm Wrapper

```rust
/// A wrapper that guarantees the inner value has unit norm.
///
/// # Example
///
/// ```rust
/// use clifford::{Unit, specialized::euclidean::dim3::Vector};
///
/// let v = Vector::new(3.0, 4.0, 0.0);
/// let unit: Unit<Vector<f64>> = Unit::new_normalize(v);
/// assert!((unit.norm() - 1.0).abs() < 1e-10);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Unit<T> {
    inner: T,
}

impl<T: Normed> Unit<T> {
    /// Creates a `Unit<T>` by normalizing the input.
    pub fn try_new(inner: T) -> Option<Self> where T: Sized {
        inner.try_normalize().map(|normalized| Self { inner: normalized })
    }

    /// Creates a `Unit<T>` by normalizing the input.
    ///
    /// # Panics
    /// Panics if the input has zero norm.
    pub fn new_normalize(inner: T) -> Self where T: Sized {
        Self::try_new(inner).expect("cannot normalize zero element")
    }

    /// Creates a `Unit<T>` without checking normalization.
    ///
    /// # Safety
    /// The caller must ensure the inner value has unit norm.
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    pub fn into_inner(self) -> T { self.inner }
    pub fn as_inner(&self) -> &T { &self.inner }
}

impl<T> Deref for Unit<T> {
    type Target = T;
    fn deref(&self) -> &Self::Target { &self.inner }
}

impl<T> AsRef<T> for Unit<T> {
    fn as_ref(&self) -> &T { &self.inner }
}

impl<T> From<Unit<T>> for T {
    fn from(unit: Unit<T>) -> T { unit.inner }
}
```

#### `Bulk<T>` - PGA Versor Wrapper

```rust
/// A wrapper guaranteeing bulk norm = 1 (PGA rigid transformations).
///
/// In PGA, the bulk norm measures the non-degenerate part. For motors,
/// `bulk_norm() == 1` guarantees a proper rigid body transformation.
///
/// # Example
///
/// ```rust
/// use clifford::{Bulk, specialized::projective::dim3::Motor};
///
/// let motor = Motor::from_translation(1.0, 2.0, 3.0);
/// let rigid: Bulk<Motor<f64>> = Bulk::new_normalize(motor);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Bulk<T> {
    inner: T,
}

impl<T: DegenerateNormed> Bulk<T> {
    pub fn try_new(inner: T) -> Option<Self> where T: Sized {
        let n = inner.bulk_norm();
        if n < T::Scalar::epsilon() { return None; }
        Some(Self { inner: inner.scale(T::Scalar::one() / n) })
    }

    pub fn new_normalize(inner: T) -> Self where T: Sized {
        Self::try_new(inner).expect("cannot bulk-normalize zero element")
    }

    pub fn new_unchecked(inner: T) -> Self { Self { inner } }
    pub fn into_inner(self) -> T { self.inner }
    pub fn as_inner(&self) -> &T { &self.inner }
}

// Deref, AsRef, From implementations...
```

#### `Ideal<T>` - PGA Homogeneous Wrapper

```rust
/// A wrapper guaranteeing weight norm = 1 (standard homogeneous form).
///
/// In PGA, the weight norm measures the projective/ideal part. For points,
/// `weight_norm() == 1` means the point is in standard form with w = 1.
///
/// # Example
///
/// ```rust
/// use clifford::{Ideal, specialized::projective::dim3::Point};
///
/// let p = Point::new(2.0, 4.0, 6.0, 2.0);  // Homogeneous
/// let std: Ideal<Point<f64>> = Ideal::new_normalize(p);
/// // Now represents (1, 2, 3) with w = 1
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Ideal<T> {
    inner: T,
}

impl<T: DegenerateNormed> Ideal<T> {
    pub fn try_new(inner: T) -> Option<Self> where T: Sized {
        inner.try_unitize().map(|u| Self { inner: u })
    }

    pub fn new_normalize(inner: T) -> Self where T: Sized {
        Self::try_new(inner).expect("cannot weight-normalize ideal element")
    }

    pub fn new_unchecked(inner: T) -> Self { Self { inner } }
    pub fn into_inner(self) -> T { self.inner }
    pub fn as_inner(&self) -> &T { &self.inner }
}

// Deref, AsRef, From implementations...
```

#### `Proper<T>` - Minkowski Timelike Wrapper

```rust
/// A wrapper guaranteeing the vector is timelike (proper).
///
/// In Minkowski spacetime, timelike vectors have positive squared norm
/// (with +--- signature). This wrapper guarantees the element represents
/// a valid 4-velocity or similar timelike quantity.
///
/// # Example
///
/// ```rust
/// use clifford::{Proper, specialized::minkowski::dim4::FourVector};
///
/// let v = FourVector::new(2.0, 1.0, 0.0, 0.0);  // Timelike
/// let proper: Option<Proper<FourVector<f64>>> = Proper::try_new(v);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Proper<T> {
    inner: T,
}

impl<T: IndefiniteNormed> Proper<T> {
    pub fn try_new(inner: T) -> Option<Self> where T: Sized {
        if inner.is_timelike() {
            Some(Self { inner: inner.normalize() })
        } else {
            None
        }
    }

    pub fn new_unchecked(inner: T) -> Self { Self { inner } }
    pub fn into_inner(self) -> T { self.inner }
    pub fn as_inner(&self) -> &T { &self.inner }
}

// Deref, AsRef, From implementations...
```

### Additional Implementations

For each wrapper, implement:
- `Clone`, `Copy`, `Debug`, `PartialEq` (derive where T implements)
- `Serialize`/`Deserialize` with serde feature
- `Arbitrary` for property-based testing

### Modified Files

- `src/lib.rs` - Add `pub mod wrappers;` and re-export types
- `src/prelude.rs` - Include all wrapper types

## Testing

```rust
#[test]
fn unit_wrapper_guarantees_norm() {
    let v = Vector::new(3.0, 4.0, 0.0);
    let unit = Unit::new_normalize(v);
    assert!(relative_eq!(unit.norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
}

#[test]
fn bulk_wrapper_for_rigid_transforms() {
    let motor = Motor::new_unchecked(2.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0);
    let rigid = Bulk::new_normalize(motor);
    assert!(relative_eq!(rigid.bulk_norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
}

#[test]
fn ideal_wrapper_for_homogeneous_coords() {
    let p = Point::new(2.0, 4.0, 6.0, 2.0);
    let std = Ideal::new_normalize(p);
    assert!(relative_eq!(std.weight_norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
}

#[test]
fn proper_wrapper_rejects_spacelike() {
    let spacelike = FourVector::new(1.0, 2.0, 0.0, 0.0);  // 1 - 4 < 0
    assert!(Proper::try_new(spacelike).is_none());
}
```

## Success Criteria

1. `Unit<T>` wrapper for Euclidean types
2. `Bulk<T>` wrapper for PGA versors
3. `Ideal<T>` wrapper for PGA homogeneous types
4. `Proper<T>` wrapper for Minkowski timelike
5. All wrappers implement `Deref`, `AsRef`, `From`
6. Serde support with feature flag
7. Clear documentation with examples
