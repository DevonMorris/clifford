# PRD-50: 2D Conformal Geometric Algebra Cl(3,1,0)

**Status**: Draft
**Goal**: Implement 2D Conformal Geometric Algebra via code generation

## Summary

Implement 2D Conformal Geometric Algebra (CGA) as Cl(3,1,0) via the code generation pipeline. CGA embeds 2D Euclidean space into a 4D algebra, enabling elegant representation of circles, lines, points, and conformal transformations in the plane.

## Motivation

2D Conformal Geometric Algebra is valuable for:
- **Circle geometry**: Circles and their intersections as first-class objects
- **Conformal transformations**: Rotations, translations, dilations, inversions in 2D
- **Educational purposes**: Simpler than 3D CGA, ideal for learning conformal algebra
- **2D graphics**: SVG-style geometry, font rendering, collision detection
- **Mechanical CAD**: 2D constraint solving with circles and lines

The algebra Cl(3,1,0) extends 2D Euclidean space with two additional basis vectors to represent the origin and point at infinity, analogous to how Cl(4,1,0) extends 3D.

## Algebra Structure

**Signature**: Cl(3,1,0)
- p = 3 (three basis vectors squaring to +1: e₁, e₂, e₃)
- q = 1 (one basis vector squaring to -1: e₄)
- r = 0 (non-degenerate)

**Dimension**: 2⁴ = 16 basis elements

| Grade | Count | Description |
|-------|-------|-------------|
| 0 | 1 | Scalar |
| 1 | 4 | Vectors (RoundPoint) |
| 2 | 6 | Bivectors (PointPair, FlatPoint) |
| 3 | 4 | Trivectors (Circle, Line) |
| 4 | 1 | Pseudoscalar |

## Basis Convention

Following the pattern established by conformal3:

- **e₁, e₂**: Euclidean basis vectors (spatial directions in 2D plane)
- **e₃**: Origin representation (e₃ = ½(e₋ - e₊))
- **e₄**: Point at infinity (e₄ = e₋ + e₊)

Where e₋² = -1 and e₊² = +1 are the conformal basis vectors.

### Alternative Representation

Some references use (e₊, e₋) directly:
- e₊² = +1 (positive signature)
- e₋² = -1 (negative signature)
- Origin: o = ½(e₋ - e₊)
- Infinity: n = e₋ + e₊

## Involution

**Primary involution**: `reverse`

The algebra is non-degenerate (r=0), so we use the reverse operation for computing norms. Since the signature is indefinite (3 positive, 1 negative), the norm can be positive, negative, or zero.

## Types

### Grade-1: RoundPoint (4 components)
A round point encodes both position and curvature in 2D:
```
p = pₓe₁ + pᵧe₂ + p_w e₃ + p_u e₄
```
- **Fields**: `x`, `y`, `w`, `u`
- **Geometric meaning**: 2D point with radius information
- **Null condition**: `x² + y² - 2wu = 0` for actual Euclidean points

### Grade-2: PointPair (6 components)
A point pair (dipole) represents two points or an oriented line segment:
```
d = m e₁₂ + vₓe₁₃ + vᵧe₂₃ + pₓe₁₄ + pᵧe₂₄ + p_w e₃₄
```
- **Fields**: `m`, `vx`, `vy`, `px`, `py`, `pw`
- **m**: Moment component (e₁₂ blade)
- **v**: Velocity components (e₁₃, e₂₃ blades)
- **p**: Position components (e₁₄, e₂₄, e₃₄ blades)
- **Subtype FlatPoint** (3 components): When m=0 and v=0 (blades e₁₄, e₂₄, e₃₄)

### Grade-3: Circle (4 components)
A circle in 2D space:
```
c = gₓe₂₃₄ + gᵧe₁₃₄ + g_w e₁₂₄ + g_u e₁₂₃
```
- **Fields**: `gx`, `gy`, `gw`, `gu`
- **Geometric meaning**: Circle with center and radius encoded
- **Subtype Line** (3 components): When gᵤ=0 (circle through infinity, i.e., a line)

### Grade-4: Pseudoscalar (1 component)
The oriented 4-volume element:
- **Fields**: `xywu`

### Even Subalgebra: Motor (8 components)
Grades 0, 2, 4 for conformal transformations:
- Grade 0: 1 scalar component
- Grade 2: 6 bivector components
- Grade 4: 1 pseudoscalar component
- **Description**: Represents general conformal motions in 2D

## Sparse Types

### FlatPoint (Grade 2, 3 blades)
Point pair with m=0 and v=0:
```
fp = pₓe₁₄ + pᵧe₂₄ + p_w e₃₄
```
- **Blades**: `e14`, `e24`, `e34`
- **Fields**: `px`, `py`, `pw`
- **Geometric meaning**: Euclidean point embedded in CGA

### Line (Grade 3, 3 blades)
Circle through infinity (gᵤ=0):
```
l = gₓe₂₃₄ + gᵧe₁₃₄ + g_w e₁₂₄
```
- **Blades**: `e234`, `e134`, `e124`
- **Fields**: `nx`, `ny`, `d` (normal and distance from origin)
- **Geometric meaning**: Directed line in 2D plane

## Key Operations

### Point Construction
Embed a Euclidean point (x, y) into CGA:
```
P = x·e₁ + y·e₂ + ½(x² + y²)·e₄ + e₃
```
Or using the normalized form:
```
P = o + x·e₁ + y·e₂ + ½|p|²·n
```
where o is the origin and n is the point at infinity.

### Circle from Three Points
Given three points P₁, P₂, P₃, the unique circle through them:
```
C = P₁ ∧ P₂ ∧ P₃
```

### Point Pair from Two Points
Given two points P₁, P₂:
```
PP = P₁ ∧ P₂
```

### Circle-Circle Intersection
Two circles C₁ and C₂ intersect at a point pair:
```
PP = C₁ ∨ C₂  (meet / antiwedge)
```

### Reflection through Circle
Reflect point P through circle C:
```
P' = C · P · C⁻¹  (sandwich product)
```

## Norm Trait

Implements `IndefiniteNormed` trait:
```rust
pub trait IndefiniteNormed {
    fn norm_squared(&self) -> f64;
    fn is_timelike(&self) -> bool { self.norm_squared() < 0.0 }
    fn is_spacelike(&self) -> bool { self.norm_squared() > 0.0 }
    fn is_null(&self) -> bool { self.norm_squared().abs() < f64::EPSILON }
}
```

For CGA, null vectors represent actual geometric points (lying on the null cone).

## Implementation Plan

1. Create `algebras/conformal2.toml` with algebra specification
2. Create module structure at `src/specialized/conformal/dim2/`
3. Update `build.rs` to include conformal2 algebra
4. Export types from `src/specialized/conformal/mod.rs`
5. Run `cargo build` to generate code
6. Add extension methods for point construction, circle extraction, etc.

## Algebra TOML Specification

```toml
# 2D Conformal Geometric Algebra
# Cl(3,1,0) - Embeds 2D Euclidean space into 4D conformal space

[algebra]
complete = false
name = "conformal2"
module_path = "conformal/dim2"
description = "2D Conformal Geometric Algebra Cl(3,1,0)"

[signature]
positive = ["e1", "e2", "e3"]
negative = ["e4"]
zero = []

[norm]
primary_involution = "reverse"

[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"

[types.RoundPoint]
grades = [1]
description = "Round point (null vectors represent actual 2D points)"
fields = ["x", "y", "w", "u"]

[types.PointPair]
grades = [2]
description = "Point pair (dipole, two points or oriented segment)"
fields = ["m", "vx", "vy", "px", "py", "pw"]

[types.Circle]
grades = [3]
description = "Circle in 2D space"
fields = ["gx", "gy", "gw", "gu"]

[types.Pseudoscalar]
grades = [4]
description = "Pseudoscalar (oriented 4-volume)"
fields = ["xywu"]

[types.Motor]
grades = [0, 2, 4]
description = "Motor (even subalgebra for conformal transformations)"
fields = ["s", "m", "vx", "vy", "px", "py", "pw", "ps"]

# Sparse types

[types.FlatPoint]
grades = [2]
fields = ["px", "py", "pw"]
blades = ["e14", "e24", "e34"]
description = "Flat point (m=0, v=0)"

[types.Line]
grades = [3]
fields = ["nx", "ny", "d"]
blades = ["e234", "e134", "e124"]
description = "Line (circle through infinity, gu=0)"
```

## File Structure

```
src/specialized/conformal/
  mod.rs              # Module root, exports both dim2 and dim3
  dim2/               # 2D Conformal GA
    mod.rs            # Documentation and re-exports
    extensions.rs     # Factory methods, coordinate extraction
    generated/        # Auto-generated code
      types.rs
      traits.rs
      conversions.rs
  dim3/               # 3D Conformal GA (existing)
    ...
```

## Extension Methods

Key methods for `extensions.rs`:

```rust
impl RoundPoint<T> {
    /// Create a point from Euclidean coordinates
    pub fn from_euclidean(x: T, y: T) -> Self;

    /// Extract Euclidean coordinates (if normalizable)
    pub fn to_euclidean(&self) -> Option<(T, T)>;

    /// Check if this is an actual point (null vector)
    pub fn is_euclidean_point(&self) -> bool;
}

impl Circle<T> {
    /// Create from center and radius
    pub fn from_center_radius(cx: T, cy: T, r: T) -> Self;

    /// Extract center and radius (if possible)
    pub fn center_radius(&self) -> Option<((T, T), T)>;

    /// Check if this is a line (circle through infinity)
    pub fn is_line(&self) -> bool;
}

impl Line<T> {
    /// Create from normal vector and distance
    pub fn from_normal_distance(nx: T, ny: T, d: T) -> Self;

    /// Create from two Euclidean points
    pub fn from_two_points(x1: T, y1: T, x2: T, y2: T) -> Self;
}

impl PointPair<T> {
    /// Extract the two points (if possible)
    pub fn to_points(&self) -> Option<((T, T), (T, T))>;
}
```

## Testing

Property-based tests to verify:

1. **Null cone**: Points constructed from Euclidean coordinates lie on null cone
2. **Circle through points**: `P₁ ∧ P₂ ∧ P₃` produces circle through all three points
3. **Incidence**: Points on a circle satisfy the incidence relation
4. **Circle-circle intersection**: Meet produces correct intersection points
5. **Reflection**: Reflection through a circle preserves angles (conformal)
6. **Motor composition**: Motors compose correctly for sequential transformations
7. **Duality**: Verify dual relationships between types

```rust
#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn point_is_null(x in -100.0f64..100.0, y in -100.0f64..100.0) {
            let p = RoundPoint::from_euclidean(x, y);
            prop_assert!(p.is_null(), "Euclidean point should be null vector");
        }

        #[test]
        fn circle_through_three_points(
            p1 in any_point(),
            p2 in any_point(),
            p3 in any_point()
        ) {
            let circle = p1.wedge(&p2).wedge(&p3);
            // Each point should be incident with the circle
            prop_assert!(is_incident(&p1, &circle));
            prop_assert!(is_incident(&p2, &circle));
            prop_assert!(is_incident(&p3, &circle));
        }
    }
}
```

## Success Criteria

- [ ] `algebras/conformal2.toml` defines all types correctly
- [ ] Code generation produces valid Rust code
- [ ] All basic operations (wedge, antiwedge, geometric product) work
- [ ] Point embedding preserves Euclidean coordinates
- [ ] Circle-from-three-points creates correct circles
- [ ] Circle intersection produces correct point pairs
- [ ] Motors represent conformal transformations correctly
- [ ] All tests pass with property-based verification

## References

- [Conformal Geometric Algebra Wiki](https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page)
- Dorst, Fontijne, Mann - "Geometric Algebra for Computer Science" (Chapter 13)
- Hestenes, Sobczyk - "Clifford Algebra to Geometric Calculus"
- [Ganja.js CGA2D Implementation](https://enkimute.github.io/ganja.js/)
