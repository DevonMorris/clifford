# PRD-6.4: Flat Points

**Status**: Pending
**Parent**: PRD-6 (Conformal Geometric Algebra)
**Depends On**: PRD-6.3 (Planes, Circles, Lines)
**Goal**: Implement flat points as the CGA analog of PGA points

## Reference

- https://conformalgeometricalgebra.org/wiki/index.php?title=Flat_point

## Background

### Flat Points vs Round Points

| Property | Round Point | Flat Point |
|----------|-------------|------------|
| Grade | 1 (vector) | 2 (bivector) |
| Representation | Null vector | Bivector with e∞ |
| Contains radius info | Yes (r=0) | No |
| PGA analog | - | Direct analog |

A **flat point** in CGA is a **grade-2 bivector** that represents a point without radius information. It's the precise analog of a point in Projective Geometric Algebra (PGA).

The relationship between flat and round points:
```
Flat Point = Round Point ∧ e∞
Round Point = Flat Point contracted with e₀
```

### Representation

A flat point has the form:
```
p = pₓ·e₁∞ + pᵧ·e₂∞ + pᵤ·e₃∞ + pᵥ·e₀∞
```

where each component contains the factor e∞ (point at infinity).

## Deliverables

### 1. Flat Point Type (`dim3/flat_point.rs`)

```rust
use crate::scalar::Float;
use super::Point;

/// A flat point in 3D CGA (grade-2 bivector).
///
/// Represents a point without radius/weight information, analogous to
/// a point in Projective Geometric Algebra (PGA).
///
/// Flat points are useful for:
/// - Line representation (dipole of flat points)
/// - Plane representation (tripole of flat points)
/// - When only position matters, not conformal weight
///
/// # Relationship to Round Points
///
/// - `FlatPoint = RoundPoint ∧ e∞`
/// - `RoundPoint = FlatPoint ⌋ e₀` (left contraction)
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, FlatPoint};
///
/// // Create a flat point at (1, 2, 3)
/// let fp = FlatPoint::new(1.0, 2.0, 3.0);
///
/// // Convert to round point
/// let rp: Point<f64> = fp.into();
/// assert!((rp.x() - 1.0).abs() < 1e-10);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Flat_point
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlatPoint<T: Float> {
    /// Coefficient of e₁∞ (x-coordinate times e∞)
    e1i: T,
    /// Coefficient of e₂∞ (y-coordinate times e∞)
    e2i: T,
    /// Coefficient of e₃∞ (z-coordinate times e∞)
    e3i: T,
    /// Coefficient of e₀∞ (weight times e∞)
    e0i: T,
}

impl<T: Float> FlatPoint<T> {
    /// Creates a flat point at Euclidean coordinates (x, y, z).
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        // FlatPoint = RoundPoint ∧ e∞
        // For unit weight: e1i = x, e2i = y, e3i = z, e0i = 1
        Self {
            e1i: x,
            e2i: y,
            e3i: z,
            e0i: T::one(),
        }
    }

    /// Creates a flat point from raw components (unchecked).
    #[inline]
    pub fn from_components_unchecked(e1i: T, e2i: T, e3i: T, e0i: T) -> Self {
        Self { e1i, e2i, e3i, e0i }
    }

    /// Creates a flat point from a round point.
    ///
    /// `FlatPoint = RoundPoint ∧ e∞`
    pub fn from_round(p: &Point<T>) -> Self {
        todo!("Derive from SymPy: derive_flat_from_round()")
    }

    /// The origin as a flat point.
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Extracts the x-coordinate.
    #[inline]
    pub fn x(&self) -> T {
        self.e1i / self.e0i
    }

    /// Extracts the y-coordinate.
    #[inline]
    pub fn y(&self) -> T {
        self.e2i / self.e0i
    }

    /// Extracts the z-coordinate.
    #[inline]
    pub fn z(&self) -> T {
        self.e3i / self.e0i
    }

    /// Returns the weight (e₀∞ component).
    #[inline]
    pub fn weight(&self) -> T {
        self.e0i
    }

    /// Normalizes the flat point so weight = 1.
    pub fn normalize(&self) -> Option<Self> {
        if self.e0i.abs() < T::epsilon() {
            None // Ideal point (at infinity)
        } else {
            Some(Self {
                e1i: self.e1i / self.e0i,
                e2i: self.e2i / self.e0i,
                e3i: self.e3i / self.e0i,
                e0i: T::one(),
            })
        }
    }

    /// Returns true if this is an ideal flat point (at infinity).
    #[inline]
    pub fn is_ideal(&self, epsilon: T) -> bool {
        self.e0i.abs() < epsilon
    }

    /// Converts to a round point.
    ///
    /// `RoundPoint = FlatPoint ⌋ e₀`
    pub fn to_round(&self) -> Point<T> {
        let (x, y, z) = (self.x(), self.y(), self.z());
        Point::new(x, y, z)
    }

    /// Component accessors.
    #[inline]
    pub fn e1i(&self) -> T { self.e1i }
    #[inline]
    pub fn e2i(&self) -> T { self.e2i }
    #[inline]
    pub fn e3i(&self) -> T { self.e3i }
    #[inline]
    pub fn e0i(&self) -> T { self.e0i }
}

impl<T: Float> From<Point<T>> for FlatPoint<T> {
    fn from(p: Point<T>) -> Self {
        FlatPoint::from_round(&p)
    }
}

impl<T: Float> From<FlatPoint<T>> for Point<T> {
    fn from(fp: FlatPoint<T>) -> Self {
        fp.to_round()
    }
}
```

### 2. Integration with Line Representation

Update `dim3/line.rs` to use flat points internally:

```rust
impl<T: Float> Line<T> {
    /// Creates a line as a dipole of flat points.
    ///
    /// This is an alternative representation that's compatible with PGA.
    pub fn from_flat_point_pair(p1: &FlatPoint<T>, p2: &FlatPoint<T>) -> Self {
        // L = p₁ ∧ p₂ (outer product of flat points)
        todo!("Derive from SymPy")
    }

    /// Extracts the line as a pair of flat points.
    pub fn to_flat_point_pair(&self) -> (FlatPoint<T>, FlatPoint<T>) {
        todo!("Derive from SymPy")
    }
}
```

## SymPy Derivations

Add to `derivations/src/clifford_derivations/cga.py`:

```python
@with_timeout(120)
def derive_flat_from_round():
    """Derive flat point from round point.

    FlatPoint = RoundPoint ∧ e∞

    The round point P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½r²·e∞
    becomes a flat point with components in e₁∞, e₂∞, e₃∞, e₀∞.
    """
    x, y, z = symbols('x y z')

    # Round point components
    e1 = x
    e2 = y
    e3 = z
    # e₀ = (e₋ - e₊)/2, e∞ = e₋ + e₊

    # Outer product with e∞ = e₋ + e₊
    # P ∧ e∞ has bivector components

    print("TODO: Implement flat from round derivation")


@with_timeout(120)
def derive_round_from_flat():
    """Derive round point from flat point.

    RoundPoint = FlatPoint ⌋ e₀

    Left contraction of flat point with origin.
    """
    print("TODO: Implement round from flat derivation")


@with_timeout(120)
def derive_line_from_flat_points():
    """Derive line as outer product of two flat points.

    L = p₁ ∧ p₂
    """
    print("TODO: Implement line from flat points derivation")
```

## Property Tests

```rust
proptest! {
    // ================================================================
    // Flat point roundtrip
    // ================================================================

    #[test]
    fn flat_point_roundtrip_coordinates(
        x in -100.0f64..100.0,
        y in -100.0f64..100.0,
        z in -100.0f64..100.0,
    ) {
        let fp = FlatPoint::new(x, y, z);
        prop_assert!(abs_diff_eq!(fp.x(), x, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(fp.y(), y, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(fp.z(), z, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Flat/Round conversion
    // ================================================================

    #[test]
    fn flat_round_roundtrip(
        x in -100.0f64..100.0,
        y in -100.0f64..100.0,
        z in -100.0f64..100.0,
    ) {
        let rp = Point::new(x, y, z);
        let fp = FlatPoint::from_round(&rp);
        let rp2 = fp.to_round();

        prop_assert!(abs_diff_eq!(rp.x(), rp2.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(rp.y(), rp2.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(rp.z(), rp2.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn round_flat_roundtrip(
        x in -100.0f64..100.0,
        y in -100.0f64..100.0,
        z in -100.0f64..100.0,
    ) {
        let fp = FlatPoint::new(x, y, z);
        let rp = fp.to_round();
        let fp2 = FlatPoint::from_round(&rp);

        prop_assert!(abs_diff_eq!(fp.x(), fp2.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(fp.y(), fp2.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(fp.z(), fp2.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Line from flat points
    // ================================================================

    #[test]
    fn line_from_flat_points_matches_round(
        x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
    ) {
        let rp1 = Point::new(x1, y1, z1);
        let rp2 = Point::new(x2, y2, z2);
        let fp1 = FlatPoint::new(x1, y1, z1);
        let fp2 = FlatPoint::new(x2, y2, z2);

        let line_round = Line::from_two_points(&rp1, &rp2);
        let line_flat = Line::from_flat_point_pair(&fp1, &fp2);

        // Both lines should be equivalent (up to scale)
        // Test by checking they contain the same points
        prop_assert!(line_flat.contains(&rp1, ABS_DIFF_EQ_EPS));
        prop_assert!(line_flat.contains(&rp2, ABS_DIFF_EQ_EPS));
    }
}
```

## Files to Create/Modify

### New Files
- `src/specialized/conformal/dim3/flat_point.rs`

### Modified Files
- `src/specialized/conformal/dim3/mod.rs` - Export FlatPoint
- `src/specialized/conformal/dim3/line.rs` - Add flat point methods
- `derivations/src/clifford_derivations/cga.py` - Add derivations

## Verification Checklist

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] Flat/Round point conversions roundtrip correctly
- [ ] Line from flat points matches line from round points
- [ ] Documentation explains relationship to PGA

## Dependencies

- PRD-6.3 (Planes, Circles, Lines) - must be complete

## Next Steps

After this PRD is complete, proceed to PRD-6.5 (Translator and Rotor Versors).
