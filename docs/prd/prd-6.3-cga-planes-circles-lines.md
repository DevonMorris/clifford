# PRD-6.3: Planes, Circles, and Lines

**Status**: Pending
**Parent**: PRD-6 (Conformal Geometric Algebra)
**Depends On**: PRD-6.2 (Round Points and Spheres)
**Goal**: Implement planes, circles, and lines as CGA blades

## Reference

- https://conformalgeometricalgebra.org/wiki/index.php?title=Plane
- https://conformalgeometricalgebra.org/wiki/index.php?title=Circle
- https://conformalgeometricalgebra.org/wiki/index.php?title=Line
- https://conformalgeometricalgebra.org/wiki/index.php?title=Dipole

## Background

### Planes (Flat Spheres)

A plane in CGA is a **sphere that passes through infinity**. It's represented as a grade-4 quadvector where the component along e∞ is zero (or equivalently, the outer product with e∞ is zero).

```
Plane: ax + by + cz + d = 0
```

Constructed as the outer product of three non-collinear points.

### Circles

A circle in CGA is a **grade-3 trivector**. It can be constructed as:
- The intersection (meet) of two spheres
- The outer product of three points: `C = P₁ ∧ P₂ ∧ P₃`

### Lines (Flat Circles)

A line in CGA is a **circle that passes through infinity**. It's represented as a grade-3 trivector where the outer product with e∞ is zero.

Constructed as:
- The intersection of two planes
- The outer product of two points (with infinity): `L = P₁ ∧ P₂ ∧ e∞`

### Dipoles (Point Pairs)

A dipole is a **grade-2 bivector** representing a pair of points. It's the outer product of two conformal points: `D = P₁ ∧ P₂`.

## Deliverables

### 1. Plane Type (`dim3/plane.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere};

/// A plane in 3D CGA (flat sphere - grade-4 quadvector).
///
/// Represents the plane `ax + by + cz + d = 0` in implicit form.
/// A plane is a sphere that contains the point at infinity.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Plane};
///
/// // xy-plane (z = 0)
/// let plane = Plane::xy();
///
/// // Point on plane
/// let p = Point::new(1.0, 2.0, 0.0);
/// assert!(plane.contains(&p, 1e-10));
///
/// // Plane from three points
/// let p1 = Point::new(0.0, 0.0, 0.0);
/// let p2 = Point::new(1.0, 0.0, 0.0);
/// let p3 = Point::new(0.0, 1.0, 0.0);
/// let xy = Plane::from_three_points(&p1, &p2, &p3);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Plane
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane<T: Float> {
    // Grade-4 components (5 choose 4 = 5 components)
    // Basis: {e₁, e₂, e₃, e₊, e₋} where p=e₊ (positive) and m=e₋ (negative)
    e123p: T,  // Coefficient of e₁₂₃₊
    e123m: T,  // Coefficient of e₁₂₃₋
    e12pm: T,  // Coefficient of e₁₂₊₋
    e13pm: T,  // Coefficient of e₁₃₊₋
    e23pm: T,  // Coefficient of e₂₃₊₋
}

impl<T: Float> Plane<T> {
    /// Creates a plane from implicit equation ax + by + cz + d = 0.
    pub fn from_implicit(a: T, b: T, c: T, d: T) -> Self {
        todo!("Derive from SymPy")
    }

    /// Creates a plane from three non-collinear points.
    ///
    /// `Plane = P₁ ∧ P₂ ∧ P₃ ∧ e∞`
    pub fn from_three_points(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> Self {
        todo!("Derive from SymPy")
    }

    /// Creates a plane from a point and normal direction.
    pub fn from_point_normal(point: &Point<T>, normal: (T, T, T)) -> Self {
        todo!("Derive from SymPy")
    }

    /// The xy-plane (z = 0).
    #[inline]
    pub fn xy() -> Self {
        Self::from_implicit(T::zero(), T::zero(), T::one(), T::zero())
    }

    /// The xz-plane (y = 0).
    #[inline]
    pub fn xz() -> Self {
        Self::from_implicit(T::zero(), T::one(), T::zero(), T::zero())
    }

    /// The yz-plane (x = 0).
    #[inline]
    pub fn yz() -> Self {
        Self::from_implicit(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Returns the normal vector (a, b, c).
    pub fn normal(&self) -> (T, T, T) {
        todo!("Derive from SymPy")
    }

    /// Returns the signed distance from origin.
    pub fn distance_from_origin(&self) -> T {
        todo!("Derive from SymPy")
    }

    /// Signed distance from a point to the plane.
    pub fn signed_distance(&self, p: &Point<T>) -> T {
        todo!("Derive from SymPy")
    }

    /// Returns true if the point lies on the plane.
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        self.signed_distance(p).abs() < epsilon
    }

    /// Projects a point onto the plane.
    pub fn project(&self, p: &Point<T>) -> Point<T> {
        todo!("Derive from SymPy")
    }

    /// Meet of two planes: line of intersection.
    pub fn meet(&self, other: &Plane<T>) -> Line<T> {
        todo!("Derive from SymPy")
    }

    /// Meet of plane and sphere: circle of intersection.
    pub fn meet_sphere(&self, sphere: &Sphere<T>) -> Circle<T> {
        todo!("Derive from SymPy")
    }
}
```

### 2. Circle Type (`dim3/circle.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere, Plane};

/// A circle in 3D CGA (grade-3 trivector).
///
/// Represents the intersection of two spheres, or three points.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Circle, Sphere};
///
/// // Circle from three points
/// let p1 = Point::new(1.0, 0.0, 0.0);
/// let p2 = Point::new(0.0, 1.0, 0.0);
/// let p3 = Point::new(-1.0, 0.0, 0.0);
/// let circle = Circle::from_three_points(&p1, &p2, &p3);
///
/// // All points lie on circle
/// assert!(circle.contains(&p1, 1e-10));
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Circle
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle<T: Float> {
    // Grade-3 components (5 choose 3 = 10 components)
    // Basis: {e₁, e₂, e₃, e₊, e₋} where p=e₊ (positive) and m=e₋ (negative)
    e123: T,   // Coefficient of e₁₂₃
    e12p: T,   // Coefficient of e₁₂₊
    e12m: T,   // Coefficient of e₁₂₋
    e13p: T,   // Coefficient of e₁₃₊
    e13m: T,   // Coefficient of e₁₃₋
    e1pm: T,   // Coefficient of e₁₊₋
    e23p: T,   // Coefficient of e₂₃₊
    e23m: T,   // Coefficient of e₂₃₋
    e2pm: T,   // Coefficient of e₂₊₋
    e3pm: T,   // Coefficient of e₃₊₋
}

impl<T: Float> Circle<T> {
    /// Creates a circle from three points.
    ///
    /// `C = P₁ ∧ P₂ ∧ P₃`
    pub fn from_three_points(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> Self {
        todo!("Derive from SymPy")
    }

    /// Creates a circle as the intersection of two spheres.
    pub fn from_sphere_intersection(s1: &Sphere<T>, s2: &Sphere<T>) -> Self {
        // Circle = S₁ ∧ S₂ (outer product gives meet in dual space)
        todo!("Derive from SymPy")
    }

    /// Creates a circle from center, radius, and normal direction.
    pub fn from_center_radius_normal(
        center: &Point<T>,
        radius: T,
        normal: (T, T, T),
    ) -> Self {
        todo!("Derive from SymPy")
    }

    /// Extracts the center point.
    pub fn center(&self) -> Point<T> {
        todo!("Derive from SymPy")
    }

    /// Returns the radius.
    pub fn radius(&self) -> T {
        todo!("Derive from SymPy")
    }

    /// Returns the normal direction of the plane containing the circle.
    pub fn normal(&self) -> (T, T, T) {
        todo!("Derive from SymPy")
    }

    /// Returns true if this is a real circle (positive radius).
    pub fn is_real(&self, epsilon: T) -> bool {
        self.radius() > epsilon
    }

    /// Returns true if the point lies on the circle.
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        todo!("Derive from SymPy")
    }

    /// Returns the plane containing this circle.
    pub fn carrier_plane(&self) -> Plane<T> {
        todo!("Derive from SymPy")
    }
}
```

### 3. Line Type (`dim3/line.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Plane, Circle};

/// A line in 3D CGA (flat circle - grade-3 trivector through infinity).
///
/// Represents a line as a circle that passes through the point at infinity.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Line};
///
/// // Line through two points
/// let p1 = Point::new(0.0, 0.0, 0.0);
/// let p2 = Point::new(1.0, 1.0, 1.0);
/// let line = Line::from_two_points(&p1, &p2);
///
/// // Both points lie on line
/// assert!(line.contains(&p1, 1e-10));
/// assert!(line.contains(&p2, 1e-10));
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Line
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line<T: Float> {
    // Grade-3 components (only the "flat" subset that includes e∞)
    // This is similar to PGA line representation
    e123: T,
    e12i: T,  // e₁₂∞ component
    e13i: T,
    e23i: T,
    e1oi: T,  // e₁₀∞ component
    e2oi: T,
    e3oi: T,
}

impl<T: Float> Line<T> {
    /// Creates a line from two points.
    ///
    /// `L = P₁ ∧ P₂ ∧ e∞`
    pub fn from_two_points(p1: &Point<T>, p2: &Point<T>) -> Self {
        todo!("Derive from SymPy")
    }

    /// Creates a line as the intersection of two planes.
    pub fn from_plane_intersection(p1: &Plane<T>, p2: &Plane<T>) -> Self {
        todo!("Derive from SymPy")
    }

    /// Creates a line from a point and direction.
    pub fn from_point_direction(point: &Point<T>, direction: (T, T, T)) -> Self {
        todo!("Derive from SymPy")
    }

    /// The x-axis.
    pub fn x_axis() -> Self {
        Self::from_point_direction(&Point::origin(), (T::one(), T::zero(), T::zero()))
    }

    /// The y-axis.
    pub fn y_axis() -> Self {
        Self::from_point_direction(&Point::origin(), (T::zero(), T::one(), T::zero()))
    }

    /// The z-axis.
    pub fn z_axis() -> Self {
        Self::from_point_direction(&Point::origin(), (T::zero(), T::zero(), T::one()))
    }

    /// Returns the direction vector (unit direction).
    pub fn direction(&self) -> (T, T, T) {
        todo!("Derive from SymPy")
    }

    /// Returns a point on the line (closest to origin).
    pub fn point_on_line(&self) -> Point<T> {
        todo!("Derive from SymPy")
    }

    /// Returns true if the point lies on the line.
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        self.distance_to_point(p) < epsilon
    }

    /// Distance from a point to the line.
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        todo!("Derive from SymPy")
    }

    /// Closest point on line to given point.
    pub fn closest_point(&self, p: &Point<T>) -> Point<T> {
        todo!("Derive from SymPy")
    }

    /// Meet of line and plane: point of intersection.
    pub fn meet_plane(&self, plane: &Plane<T>) -> Point<T> {
        todo!("Derive from SymPy")
    }
}
```

### 4. Dipole Type (`dim3/dipole.rs`)

```rust
use crate::scalar::Float;
use super::Point;

/// A dipole (point pair) in 3D CGA (grade-2 bivector).
///
/// Represents two points as a single algebraic object.
/// `D = P₁ ∧ P₂`
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Dipole};
///
/// let p1 = Point::new(0.0, 0.0, 0.0);
/// let p2 = Point::new(1.0, 0.0, 0.0);
/// let dipole = Dipole::from_two_points(&p1, &p2);
///
/// // Extract the points
/// let (q1, q2) = dipole.extract_points();
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Dipole
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dipole<T: Float> {
    // Grade-2 components (5 choose 2 = 10 components)
    // Basis: {e₁, e₂, e₃, e₊, e₋} where p=e₊ (positive) and m=e₋ (negative)
    e12: T,   // Coefficient of e₁₂
    e13: T,   // Coefficient of e₁₃
    e1p: T,   // Coefficient of e₁₊
    e1m: T,   // Coefficient of e₁₋
    e23: T,   // Coefficient of e₂₃
    e2p: T,   // Coefficient of e₂₊
    e2m: T,   // Coefficient of e₂₋
    e3p: T,   // Coefficient of e₃₊
    e3m: T,   // Coefficient of e₃₋
    epm: T,   // Coefficient of e₊₋
}

impl<T: Float> Dipole<T> {
    /// Creates a dipole from two points.
    ///
    /// `D = P₁ ∧ P₂`
    pub fn from_two_points(p1: &Point<T>, p2: &Point<T>) -> Self {
        todo!("Derive from SymPy")
    }

    /// Extracts the two points from the dipole.
    ///
    /// Returns `None` if the dipole is degenerate (same point twice).
    pub fn extract_points(&self) -> Option<(Point<T>, Point<T>)> {
        todo!("Derive from SymPy - solve quadratic")
    }

    /// Distance between the two points.
    pub fn separation(&self) -> T {
        todo!("Derive from SymPy")
    }

    /// Midpoint of the two points.
    pub fn midpoint(&self) -> Point<T> {
        todo!("Derive from SymPy")
    }

    /// Returns true if this dipole is real (two distinct points).
    pub fn is_real(&self, epsilon: T) -> bool {
        self.separation() > epsilon
    }
}
```

## SymPy Derivations

Add to `derivations/src/clifford_derivations/cga.py`:

```python
@with_timeout(300)
def derive_plane_from_implicit():
    """Derive plane representation from ax + by + cz + d = 0."""
    a, b, c, d = symbols('a b c d')
    # Plane as grade-4 element
    print("TODO: Implement plane from implicit derivation")


@with_timeout(300)
def derive_circle_from_three_points():
    """Derive circle as C = P₁ ∧ P₂ ∧ P₃."""
    # This is an outer product of three conformal points
    print("TODO: Implement circle from three points derivation")


@with_timeout(300)
def derive_line_from_two_points():
    """Derive line as L = P₁ ∧ P₂ ∧ e∞."""
    print("TODO: Implement line from two points derivation")


@with_timeout(300)
def derive_dipole_extract_points():
    """Derive extraction of two points from dipole bivector."""
    # This requires solving a quadratic
    print("TODO: Implement dipole point extraction")
```

## Property Tests

```rust
proptest! {
    // ================================================================
    // Plane tests
    // ================================================================

    #[test]
    fn plane_from_three_points_contains_them(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
        p3 in any::<Point<f64>>(),
    ) {
        // Skip if collinear
        let plane = Plane::from_three_points(&p1, &p2, &p3);
        prop_assert!(plane.contains(&p1, ABS_DIFF_EQ_EPS));
        prop_assert!(plane.contains(&p2, ABS_DIFF_EQ_EPS));
        prop_assert!(plane.contains(&p3, ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn plane_distance_matches_euclidean(
        a in -1.0f64..1.0, b in -1.0f64..1.0, c in -1.0f64..1.0,
        d in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        // Normalize (a, b, c)
        let len = (a*a + b*b + c*c).sqrt();
        if len < ABS_DIFF_EQ_EPS {
            return Ok(());
        }
        let (a, b, c) = (a/len, b/len, c/len);
        let d = d / len;

        let plane = Plane::from_implicit(a, b, c, d);
        let p = Point::new(px, py, pz);

        let cga_dist = plane.signed_distance(&p);
        let euclidean_dist = a*px + b*py + c*pz + d;

        prop_assert!(abs_diff_eq!(cga_dist, euclidean_dist, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Circle tests
    // ================================================================

    #[test]
    fn circle_from_three_points_contains_them(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
        p3 in any::<Point<f64>>(),
    ) {
        let circle = Circle::from_three_points(&p1, &p2, &p3);
        prop_assert!(circle.contains(&p1, ABS_DIFF_EQ_EPS));
        prop_assert!(circle.contains(&p2, ABS_DIFF_EQ_EPS));
        prop_assert!(circle.contains(&p3, ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Line tests
    // ================================================================

    #[test]
    fn line_from_two_points_contains_them(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
    ) {
        let line = Line::from_two_points(&p1, &p2);
        prop_assert!(line.contains(&p1, ABS_DIFF_EQ_EPS));
        prop_assert!(line.contains(&p2, ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_distance_is_perpendicular(
        x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
        dx in -1.0f64..1.0, dy in -1.0f64..1.0, dz in -1.0f64..1.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let p1 = Point::new(x1, y1, z1);
        let len = (dx*dx + dy*dy + dz*dz).sqrt();
        if len < ABS_DIFF_EQ_EPS {
            return Ok(());
        }
        let p2 = Point::new(x1 + dx, y1 + dy, z1 + dz);
        let line = Line::from_two_points(&p1, &p2);

        let p = Point::new(px, py, pz);
        let closest = line.closest_point(&p);

        // Vector from closest to p should be perpendicular to line direction
        let to_p = (p.x() - closest.x(), p.y() - closest.y(), p.z() - closest.z());
        let dir = line.direction();
        let dot = to_p.0 * dir.0 + to_p.1 * dir.1 + to_p.2 * dir.2;

        prop_assert!(abs_diff_eq!(dot, 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Dipole tests
    // ================================================================

    #[test]
    fn dipole_roundtrip(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
    ) {
        let dipole = Dipole::from_two_points(&p1, &p2);
        if let Some((q1, q2)) = dipole.extract_points() {
            // Points should match (possibly swapped)
            let d1 = p1.distance(&q1) + p2.distance(&q2);
            let d2 = p1.distance(&q2) + p2.distance(&q1);
            prop_assert!(d1 < ABS_DIFF_EQ_EPS || d2 < ABS_DIFF_EQ_EPS);
        }
    }
}
```

## Files to Create/Modify

### New Files
- `src/specialized/conformal/dim3/plane.rs`
- `src/specialized/conformal/dim3/circle.rs`
- `src/specialized/conformal/dim3/line.rs`
- `src/specialized/conformal/dim3/dipole.rs`

### Modified Files
- `src/specialized/conformal/dim3/mod.rs` - Export new types
- `derivations/src/clifford_derivations/cga.py` - Add derivations

## Verification Checklist

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] Plane containment tests pass
- [ ] Circle from three points works correctly
- [ ] Line distance calculations match expected values
- [ ] Dipole point extraction roundtrips

## Dependencies

- PRD-6.2 (Round Points and Spheres) - must be complete

## Next Steps

After this PRD is complete, proceed to PRD-6.4 (Flat Points).
