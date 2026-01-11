# PRD-6.2: Round Points and Spheres

**Status**: Pending
**Parent**: PRD-6 (Conformal Geometric Algebra)
**Depends On**: PRD-6.1 (CGA Foundation)
**Goal**: Implement conformal points and spheres with the null constraint

## Reference

- https://conformalgeometricalgebra.org/wiki/index.php?title=Round_point
- https://conformalgeometricalgebra.org/wiki/index.php?title=Sphere

## Background

### Round Points (Conformal Points)

A round point in CGA is a **null vector** (grade-1 element that squares to zero). The conformal embedding maps a Euclidean point **p** = (x, y, z) to:

```
P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
```

The null constraint ensures `P · P = 0`, which is satisfied by this embedding.

### Spheres

A sphere in CGA is a **grade-4 quadvector**. The dual sphere representation:

```
S = c - ½r²·e∞
```

where `c` is the center point (as a conformal point) and `r` is the radius.

- **Real sphere**: r² > 0 (positive radius squared)
- **Point sphere**: r² = 0 (degenerate to a point)
- **Imaginary sphere**: r² < 0 (no real solutions)

## Deliverables

### 1. Module Structure

```
src/specialized/conformal/
├── mod.rs
└── dim3/
    ├── mod.rs
    ├── point.rs      # Round point type
    ├── sphere.rs     # Sphere type
    ├── arbitrary.rs  # Proptest support
    └── nalgebra.rs   # nalgebra conversions
```

### 2. Round Point Type (`dim3/point.rs`)

```rust
use crate::scalar::Float;

/// A round point in 3D CGA (null vector).
///
/// Represents a Euclidean point (x, y, z) embedded in the conformal model:
/// `P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞`
///
/// The null constraint `P · P = 0` is maintained by construction.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::Point;
///
/// // Create a point at (1, 2, 3)
/// let p = Point::new(1.0, 2.0, 3.0);
///
/// // Verify null constraint
/// assert!(p.is_null(1e-10));
///
/// // Distance between points
/// let q = Point::new(4.0, 6.0, 3.0);
/// let d = p.distance(&q);
/// assert!((d - 5.0).abs() < 1e-10);  // 3-4-5 triangle
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Round_point
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<T: Float> {
    /// Coefficient of e₁ (x-coordinate)
    e1: T,
    /// Coefficient of e₂ (y-coordinate)
    e2: T,
    /// Coefficient of e₃ (z-coordinate)
    e3: T,
    /// Coefficient of e₊ (part of null basis)
    ep: T,
    /// Coefficient of e₋ (part of null basis)
    em: T,
}

impl<T: Float> Point<T> {
    /// Creates a point at Euclidean coordinates (x, y, z).
    ///
    /// The conformal embedding ensures `P · P = 0`.
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        // P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
        // where e₀ = (e₋ - e₊)/2, e∞ = e₋ + e₊
        //
        // ep component: -1/2 + (x² + y² + z²)/2 = (x² + y² + z² - 1)/2
        // em component: 1/2 + (x² + y² + z²)/2 = (x² + y² + z² + 1)/2
        let r_sq = x * x + y * y + z * z;
        let half = T::one() / T::TWO;
        Self {
            e1: x,
            e2: y,
            e3: z,
            ep: (r_sq - T::one()) * half,
            em: (r_sq + T::one()) * half,
        }
    }

    /// Creates a point from raw conformal components (unchecked).
    ///
    /// Use this when you know the components satisfy the null constraint.
    #[inline]
    pub fn from_conformal_unchecked(e1: T, e2: T, e3: T, ep: T, em: T) -> Self {
        Self { e1, e2, e3, ep, em }
    }

    /// The origin point (0, 0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Extracts the Euclidean x-coordinate.
    ///
    /// For a properly normalized point, this returns the x-coordinate.
    #[inline]
    pub fn x(&self) -> T {
        // The e₀ component (weight) is (em - ep)/2 = 1 for normalized points
        // We divide by the weight to get Cartesian coordinates
        let weight = (self.em - self.ep) / T::TWO;
        self.e1 / weight
    }

    /// Extracts the Euclidean y-coordinate.
    #[inline]
    pub fn y(&self) -> T {
        let weight = (self.em - self.ep) / T::TWO;
        self.e2 / weight
    }

    /// Extracts the Euclidean z-coordinate.
    #[inline]
    pub fn z(&self) -> T {
        let weight = (self.em - self.ep) / T::TWO;
        self.e3 / weight
    }

    /// Returns the weight (e₀ component).
    ///
    /// For points created via `new()`, this is always 1.
    #[inline]
    pub fn weight(&self) -> T {
        (self.em - self.ep) / T::TWO
    }

    /// Checks if this is a valid null vector (P · P = 0).
    #[inline]
    pub fn is_null(&self, epsilon: T) -> bool {
        // P · P = e1² + e2² + e3² + ep² - em²
        //       = x² + y² + z² + ep² - em²
        let norm_sq = self.e1 * self.e1
            + self.e2 * self.e2
            + self.e3 * self.e3
            + self.ep * self.ep
            - self.em * self.em;
        norm_sq.abs() < epsilon
    }

    /// Euclidean distance to another point.
    ///
    /// Uses the inner product formula: d² = -2(P₁ · P₂)
    pub fn distance(&self, other: &Point<T>) -> T {
        // P₁ · P₂ = e1₁·e1₂ + e2₁·e2₂ + e3₁·e3₂ + ep₁·ep₂ - em₁·em₂
        let dot = self.e1 * other.e1
            + self.e2 * other.e2
            + self.e3 * other.e3
            + self.ep * other.ep
            - self.em * other.em;

        // d² = -2(P₁ · P₂) / (w₁ · w₂)
        // For unit weight points, d² = -2(P₁ · P₂)
        let w1 = self.weight();
        let w2 = other.weight();
        (T::from_i8(-2) * dot / (w1 * w2)).sqrt()
    }

    /// Accessor for e₁ component.
    #[inline]
    pub fn e1(&self) -> T { self.e1 }
    /// Accessor for e₂ component.
    #[inline]
    pub fn e2(&self) -> T { self.e2 }
    /// Accessor for e₃ component.
    #[inline]
    pub fn e3(&self) -> T { self.e3 }
    /// Accessor for e₊ component.
    #[inline]
    pub fn ep(&self) -> T { self.ep }
    /// Accessor for e₋ component.
    #[inline]
    pub fn em(&self) -> T { self.em }
}
```

### 3. Sphere Type (`dim3/sphere.rs`)

```rust
use crate::scalar::Float;
use super::Point;

/// A sphere in 3D CGA (grade-4 quadvector).
///
/// Represents a sphere with center (cx, cy, cz) and radius r.
///
/// In dual form: `S* = C - ½r²·e∞` where C is the center point.
///
/// # Real vs Imaginary
///
/// - `r² > 0`: Real sphere with positive radius
/// - `r² = 0`: Point sphere (degenerate)
/// - `r² < 0`: Imaginary sphere (no real points)
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Sphere};
///
/// // Unit sphere at origin
/// let sphere = Sphere::from_center_radius(0.0, 0.0, 0.0, 1.0);
///
/// // Point on sphere surface
/// let p = Point::new(1.0, 0.0, 0.0);
/// assert!(sphere.contains(&p, 1e-10));
///
/// // Point inside sphere
/// let q = Point::new(0.5, 0.0, 0.0);
/// assert!(sphere.signed_distance(&q) < 0.0);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Sphere
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere<T: Float> {
    /// Coefficient of e₁₂₃₊ (proportional to center x)
    e123p: T,
    /// Coefficient of e₁₂₃₋ (proportional to center y)
    e123m: T,
    /// Coefficient of e₁₂₊₋ (proportional to center z)
    e12pm: T,
    /// Coefficient of e₁₃₊₋ (related to radius and center)
    e13pm: T,
    /// Coefficient of e₂₃₊₋ (related to radius and center)
    e23pm: T,
}

impl<T: Float> Sphere<T> {
    /// Creates a sphere from center coordinates and radius.
    pub fn from_center_radius(cx: T, cy: T, cz: T, r: T) -> Self {
        // Implementation derived from SymPy
        // S* = C - ½r²·e∞
        todo!("Derive from SymPy: derive_sphere_from_center_radius()")
    }

    /// Creates a sphere from four points on its surface.
    ///
    /// The sphere is the outer product of four conformal points:
    /// `S = P₁ ∧ P₂ ∧ P₃ ∧ P₄`
    pub fn from_four_points(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>, p4: &Point<T>) -> Self {
        todo!("Derive from SymPy: derive_sphere_from_four_points()")
    }

    /// Extracts the center as a Euclidean point.
    pub fn center(&self) -> (T, T, T) {
        todo!("Derive from SymPy: derive_sphere_center_extraction()")
    }

    /// Returns the squared radius.
    ///
    /// - Positive: real sphere
    /// - Zero: point sphere
    /// - Negative: imaginary sphere
    pub fn radius_squared(&self) -> T {
        todo!("Derive from SymPy: derive_sphere_radius_squared()")
    }

    /// Returns the radius (panics if imaginary).
    pub fn radius(&self) -> T {
        let r_sq = self.radius_squared();
        assert!(r_sq >= T::zero(), "Cannot take radius of imaginary sphere");
        r_sq.sqrt()
    }

    /// Returns true if this is a real sphere (r² > 0).
    pub fn is_real(&self, epsilon: T) -> bool {
        self.radius_squared() > epsilon
    }

    /// Returns true if this is a point sphere (r² ≈ 0).
    pub fn is_point(&self, epsilon: T) -> bool {
        self.radius_squared().abs() < epsilon
    }

    /// Returns true if this is an imaginary sphere (r² < 0).
    pub fn is_imaginary(&self, epsilon: T) -> bool {
        self.radius_squared() < -epsilon
    }

    /// Returns true if the point lies on the sphere surface.
    ///
    /// Uses the inner product: `P · S = 0` iff P on S.
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        self.signed_distance(p).abs() < epsilon
    }

    /// Signed distance from point to sphere surface.
    ///
    /// - Negative: inside sphere
    /// - Zero: on surface
    /// - Positive: outside sphere
    pub fn signed_distance(&self, p: &Point<T>) -> T {
        todo!("Derive from SymPy: derive_sphere_point_distance()")
    }
}
```

### 4. SymPy Derivations

Add to `derivations/src/clifford_derivations/cga.py`:

```python
@with_timeout(120)
def derive_point_embedding():
    """Derive conformal point embedding.

    P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞

    Returns Rust code for Point::new().
    """
    x, y, z = symbols('x y z')

    # e₀ = (e₋ - e₊)/2, e∞ = e₋ + e₊
    # P = x·e₁ + y·e₂ + z·e₃ + (e₋ - e₊)/2 + ½(x² + y² + z²)·(e₋ + e₊)

    r_sq = x**2 + y**2 + z**2

    # ep coefficient: -1/2 + r²/2
    ep = simplify(-Rational(1, 2) + r_sq / 2)
    # em coefficient: 1/2 + r²/2
    em = simplify(Rational(1, 2) + r_sq / 2)

    print("Point embedding:")
    print(f"  e1 = {x}")
    print(f"  e2 = {y}")
    print(f"  e3 = {z}")
    print(f"  ep = {rust_code(ep)}")
    print(f"  em = {rust_code(em)}")

    return {'e1': x, 'e2': y, 'e3': z, 'ep': ep, 'em': em}


@with_timeout(120)
def derive_point_distance():
    """Derive distance formula between two conformal points.

    d² = -2(P₁ · P₂) for unit weight points.
    """
    x1, y1, z1, x2, y2, z2 = symbols('x1 y1 z1 x2 y2 z2')

    # Create conformal points
    r1_sq = x1**2 + y1**2 + z1**2
    r2_sq = x2**2 + y2**2 + z2**2

    ep1 = (r1_sq - 1) / 2
    em1 = (r1_sq + 1) / 2
    ep2 = (r2_sq - 1) / 2
    em2 = (r2_sq + 1) / 2

    # Inner product: P₁ · P₂ = e1₁·e1₂ + e2₁·e2₂ + e3₁·e3₂ + ep₁·ep₂ - em₁·em₂
    dot = x1*x2 + y1*y2 + z1*z2 + ep1*ep2 - em1*em2

    # d² = -2(P₁ · P₂)
    d_sq = simplify(-2 * dot)

    # Should equal (x1-x2)² + (y1-y2)² + (z1-z2)²
    expected = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2

    assert simplify(d_sq - expected) == 0, "Distance formula incorrect!"

    print("Distance squared formula verified:")
    print(f"  d² = {rust_code(d_sq)}")
    print(f"  Matches Euclidean: (x1-x2)² + (y1-y2)² + (z1-z2)²")

    return d_sq


@with_timeout(300)
def derive_sphere_from_center_radius():
    """Derive sphere representation from center and radius.

    S* = C - ½r²·e∞
    """
    cx, cy, cz, r = symbols('cx cy cz r')

    # TODO: Compute the full sphere components
    # This requires the grade-4 representation

    print("TODO: Implement sphere from center/radius derivation")


@with_timeout(300)
def derive_sphere_from_four_points():
    """Derive sphere from four points: S = P₁ ∧ P₂ ∧ P₃ ∧ P₄."""
    # This is computationally expensive
    print("TODO: Implement sphere from four points derivation")
```

### 5. nalgebra Conversions

```rust
#[cfg(any(
    feature = "nalgebra-0_32",
    feature = "nalgebra-0_33",
    feature = "nalgebra-0_34"
))]
mod nalgebra_impl {
    use super::*;

    cfg_if::cfg_if! {
        if #[cfg(feature = "nalgebra-0_34")] {
            use nalgebra_0_34 as na;
        } else if #[cfg(feature = "nalgebra-0_33")] {
            use nalgebra_0_33 as na;
        } else if #[cfg(feature = "nalgebra-0_32")] {
            use nalgebra_0_32 as na;
        }
    }

    impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T> {
        fn from(p: na::Point3<T>) -> Self {
            Point::new(p.x, p.y, p.z)
        }
    }

    impl<T: Float + na::Scalar> From<Point<T>> for na::Point3<T> {
        fn from(p: Point<T>) -> Self {
            na::Point3::new(p.x(), p.y(), p.z())
        }
    }

    impl<T: Float + na::Scalar> From<Sphere<T>> for (na::Point3<T>, T) {
        /// Converts sphere to (center, radius) tuple.
        fn from(s: Sphere<T>) -> Self {
            let (cx, cy, cz) = s.center();
            (na::Point3::new(cx, cy, cz), s.radius())
        }
    }
}
```

## Property Tests

```rust
proptest! {
    // ================================================================
    // Point null constraint
    // ================================================================

    #[test]
    fn point_is_null(
        x in -100.0f64..100.0,
        y in -100.0f64..100.0,
        z in -100.0f64..100.0,
    ) {
        let p = Point::new(x, y, z);
        prop_assert!(p.is_null(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn point_roundtrip_coordinates(
        x in -100.0f64..100.0,
        y in -100.0f64..100.0,
        z in -100.0f64..100.0,
    ) {
        let p = Point::new(x, y, z);
        prop_assert!(abs_diff_eq!(p.x(), x, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(p.y(), y, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(p.z(), z, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Distance
    // ================================================================

    #[test]
    fn distance_is_euclidean(
        x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
    ) {
        let p1 = Point::new(x1, y1, z1);
        let p2 = Point::new(x2, y2, z2);

        let cga_dist = p1.distance(&p2);
        let euclidean_dist = ((x2-x1).powi(2) + (y2-y1).powi(2) + (z2-z1).powi(2)).sqrt();

        prop_assert!(abs_diff_eq!(cga_dist, euclidean_dist, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn distance_matches_nalgebra(
        x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
    ) {
        let p1 = Point::new(x1, y1, z1);
        let p2 = Point::new(x2, y2, z2);

        let na_p1 = na::Point3::new(x1, y1, z1);
        let na_p2 = na::Point3::new(x2, y2, z2);

        let cga_dist = p1.distance(&p2);
        let na_dist = na::distance(&na_p1, &na_p2);

        prop_assert!(abs_diff_eq!(cga_dist, na_dist, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Sphere containment
    // ================================================================

    #[test]
    fn sphere_contains_surface_points(
        cx in -10.0f64..10.0, cy in -10.0f64..10.0, cz in -10.0f64..10.0,
        r in 0.1f64..10.0,
        theta in 0.0f64..std::f64::consts::TAU,
        phi in 0.0f64..std::f64::consts::PI,
    ) {
        let sphere = Sphere::from_center_radius(cx, cy, cz, r);
        let p = Point::new(
            cx + r * phi.sin() * theta.cos(),
            cy + r * phi.sin() * theta.sin(),
            cz + r * phi.cos(),
        );
        prop_assert!(sphere.contains(&p, ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn sphere_from_four_points_contains_them(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
        p3 in any::<Point<f64>>(),
        p4 in any::<Point<f64>>(),
    ) {
        let sphere = Sphere::from_four_points(&p1, &p2, &p3, &p4);
        prop_assert!(sphere.contains(&p1, ABS_DIFF_EQ_EPS));
        prop_assert!(sphere.contains(&p2, ABS_DIFF_EQ_EPS));
        prop_assert!(sphere.contains(&p3, ABS_DIFF_EQ_EPS));
        prop_assert!(sphere.contains(&p4, ABS_DIFF_EQ_EPS));
    }
}
```

## Files to Create/Modify

### New Files
- `src/specialized/conformal/mod.rs`
- `src/specialized/conformal/dim3/mod.rs`
- `src/specialized/conformal/dim3/point.rs`
- `src/specialized/conformal/dim3/sphere.rs`
- `src/specialized/conformal/dim3/arbitrary.rs`
- `src/specialized/conformal/dim3/nalgebra.rs`

### Modified Files
- `src/specialized/mod.rs` - Add `conformal` module
- `derivations/src/clifford_derivations/cga.py` - Add derivations

## Verification Checklist

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] Point null constraint verified by property tests
- [ ] Distance formula matches Euclidean/nalgebra
- [ ] Sphere containment tests pass
- [ ] nalgebra conversions work correctly

## Dependencies

- PRD-6.1 (CGA Foundation) - must be complete

## Next Steps

After this PRD is complete, proceed to PRD-6.3 (Planes, Circles, Lines).
