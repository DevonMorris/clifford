# PRD-22: Geometric Projections and Distance Operations

## Status: Proposed

## Problem

The `extensions.rs` file for 3D PGA is missing several geometric projection and distance operations documented in the [RGA wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Projections). While some operations exist (e.g., `Plane::project_point`), the implementations are ad-hoc rather than using the canonical PGA formulas, and several operations are missing entirely.

## Background

### Canonical Projection Formula

The RGA wiki defines projections using the **weight expansion** pattern:

```
project(a onto b) = b ∨ (a ∧ b⋆)
```

Where:
- `⋆` = complement (Hodge dual)
- `∧` = exterior product (join)
- `∨` = regressive product (meet)

This formula works uniformly for projecting any geometric object onto any other.

### Available Generated Products

The codegen already generates all necessary products:
- `complement_*` - Hodge complement for all types
- `exterior_*_*` - Exterior products between types
- `regressive_*_*` - Regressive products between types

## Current State

### Already Implemented (but some need refactoring)

| Method | Location | Notes |
|--------|----------|-------|
| `Point::distance(Point)` | extensions.rs:198 | Uses Cartesian coords, OK |
| `Point::distance_squared(Point)` | extensions.rs:203 | Uses Cartesian coords, OK |
| `Line::angle(Line)` | extensions.rs:405 | Manual formula |
| `Line::distance(Line)` | extensions.rs:413 | Manual formula |
| `Line::distance_to_point(Point)` | extensions.rs:467 | Manual formula |
| `Line::closest_point(Point)` | extensions.rs:508 | Manual formula |
| `Plane::distance_from_origin()` | extensions.rs:576 | Manual formula |
| `Plane::angle(Plane)` | extensions.rs:639 | Manual formula |
| `Plane::angle_to_line(Line)` | extensions.rs:647 | Manual formula |
| `Plane::project_point(Point)` | extensions.rs:665 | Manual formula |
| `Plane::project_line(Line)` | extensions.rs:677 | Manual formula |

### Missing Operations

| Operation | Formula (RGA wiki) | Description |
|-----------|-------------------|-------------|
| `Point::project_onto_plane(Plane)` | `g ∨ (p ∧ g⋆)` | Closest point on plane to point |
| `Point::project_onto_line(Line)` | `l ∨ (p ∧ l⋆)` | Closest point on line to point |
| `Point::distance_to_plane(Plane)` | `(p · g) / ‖g‖` | Signed distance from point to plane |
| `Point::distance_to_line(Line)` | `‖l_v × p_xyz + p_w·l_m‖ / ‖l‖` | Distance from point to line |
| `Line::angle_to_plane(Plane)` | `‖g_xyz × l_v‖ / (‖g‖·‖l‖)` | Angle between line and plane |
| `Line::project_onto_plane(Plane)` | `g ∨ (l ∧ g⋆)` | Project line onto plane |

## Proposed Changes

### Phase 1: Add Missing Point Operations

Add to `Point<T>`:

```rust
/// Project this point onto a plane.
///
/// Returns the point on the plane closest to this point.
/// Uses the canonical formula: `g ∨ (p ∧ g⋆)`
///
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Projections>
pub fn project_onto_plane(&self, plane: &Plane<T>) -> Point<T> {
    // g⋆ = complement_plane(g) -> Point
    // p ∧ g⋆ = exterior_point_point -> Line
    // g ∨ Line = regressive_plane_line -> Point
    let plane_complement = products::complement_plane(plane);
    let wedge = products::exterior_point_point(self, &plane_complement);
    products::regressive_plane_line(plane, &wedge)
}

/// Project this point onto a line.
///
/// Returns the point on the line closest to this point.
/// Uses the canonical formula: `l ∨ (p ∧ l⋆)`
///
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Projections>
pub fn project_onto_line(&self, line: &Line<T>) -> Point<T> {
    // l⋆ = complement_line(l) -> Line (self-dual in 4D)
    // p ∧ l⋆ = exterior_point_line -> Plane
    // l ∨ Plane = regressive_line_plane -> Point
    let line_complement = products::complement_line(line);
    let wedge = products::exterior_point_line(self, &line_complement);
    products::regressive_line_plane(line, &wedge)
}

/// Signed distance from this point to a plane.
///
/// Positive if point is on the side the plane normal points to.
///
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Euclidean_distance>
pub fn distance_to_plane(&self, plane: &Plane<T>) -> T {
    // (p · g) / ‖g‖
    let dot = products::dot_point_plane(self, plane);
    dot / plane.bulk_norm()
}

/// Distance from this point to a line.
///
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Euclidean_distance>
pub fn distance_to_line(&self, line: &Line<T>) -> T {
    // ‖l_v × p_xyz + p_w·l_m‖ / ‖l‖
    line.distance_to_point(self)
}
```

### Phase 2: Add Missing Line Operations

Add to `Line<T>`:

```rust
/// Angle between this line and a plane (in radians).
///
/// Returns the angle between the line and its projection onto the plane.
/// This is the complement of the angle between the line and the plane normal.
///
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Euclidean_angle>
pub fn angle_to_plane(&self, plane: &Plane<T>) -> T {
    // cos(φ) = ‖g_xyz × l_v‖ / (‖g‖·‖l‖)
    let n = plane.normal();
    let d = self.direction();
    let cross = n.cross(d);
    let cross_norm = cross.norm();
    let denom = plane.bulk_norm() * self.bulk_norm();
    if denom < T::epsilon() {
        T::zero()
    } else {
        (cross_norm / denom).acos()
    }
}

/// Project this line onto a plane.
///
/// Returns the line on the plane that is the orthogonal projection.
/// Uses the canonical formula: `g ∨ (l ∧ g⋆)`
///
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Projections>
pub fn project_onto_plane(&self, plane: &Plane<T>) -> Line<T> {
    // g⋆ = complement_plane(g) -> Point
    // l ∧ g⋆ = exterior_line_point -> Plane
    // g ∨ Plane = regressive_plane_plane -> Line
    let plane_complement = products::complement_plane(plane);
    let wedge = products::exterior_line_point(self, &plane_complement);
    products::regressive_plane_plane(plane, &wedge)
}
```

### Phase 3: Refactor Existing Operations (Optional)

Consider refactoring existing implementations to use the canonical formulas for consistency and to leverage the verified generated products. This is lower priority since the current implementations are tested and working.

## Implementation Notes

### Using Generated Products

All projection operations should use the generated products from `generated/products.rs`:
- `products::complement_*` for Hodge complement
- `products::exterior_*_*` for exterior products
- `products::regressive_*_*` for regressive products
- `products::dot_*_*` for inner products (distances)

### Norm Calculations

Distance and angle formulas require bulk norms:
- `Point`: Already has `bulk_norm()` via `DegenerateNormed` trait
- `Line`: Already has `bulk_norm()` via `DegenerateNormed` trait
- `Plane`: Already has `bulk_norm()` via `DegenerateNormed` trait

### Edge Cases

Handle degenerate cases:
- Zero-norm objects (ideal points, lines at infinity)
- Parallel lines/planes (zero cross product)
- Points on the geometry being projected onto

## Testing

### Property-Based Tests

```rust
proptest! {
    /// Projection of a point onto a plane lies on the plane.
    #[test]
    fn point_projection_onto_plane_lies_on_plane(
        p in any::<UnitizedPoint<f64>>(),
        g in any::<UnitizedPlane<f64>>(),
    ) {
        let proj = p.project_onto_plane(&g);
        // proj should satisfy the plane equation
        let residual = products::dot_point_plane(&proj, &g);
        prop_assert!(relative_eq!(residual, 0.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }

    /// Projection is the closest point (distance is minimized).
    #[test]
    fn point_projection_minimizes_distance(
        p in any::<UnitizedPoint<f64>>(),
        g in any::<UnitizedPlane<f64>>(),
    ) {
        let proj = p.project_onto_plane(&g);
        let dist_to_proj = p.distance(&proj);
        let dist_to_plane = p.distance_to_plane(&g).abs();
        prop_assert!(relative_eq!(dist_to_proj, dist_to_plane, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }

    /// Projection onto a line lies on the line.
    #[test]
    fn point_projection_onto_line_lies_on_line(
        p in any::<UnitizedPoint<f64>>(),
        l in any::<UnitizedLine<f64>>(),
    ) {
        let proj = p.project_onto_line(&l);
        // The distance from proj to line should be zero
        let dist = l.distance_to_point(&proj);
        prop_assert!(relative_eq!(dist, 0.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }
}
```

## Dependencies

- Generated products: `complement_*`, `exterior_*`, `regressive_*`, `dot_*`
- `DegenerateNormed` trait for bulk norms

## Acceptance Criteria

1. [ ] `Point::project_onto_plane()` implemented using canonical formula
2. [ ] `Point::project_onto_line()` implemented using canonical formula
3. [ ] `Point::distance_to_plane()` implemented
4. [ ] `Point::distance_to_line()` delegates to `Line::distance_to_point()`
5. [ ] `Line::angle_to_plane()` implemented
6. [ ] `Line::project_onto_plane()` implemented using canonical formula
7. [ ] Property-based tests verify geometric properties
8. [ ] Documentation links to RGA wiki
9. [ ] All implementations use generated products (no manual formulas)

## Future Considerations

- 2D PGA projections (similar patterns, simpler formulas)
- Antiprojections (if defined in RGA)
- Orthogonal complements and rejections
- Conformal GA projections (future algebra)

## References

- [Projections - RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Projections)
- [Euclidean Distance - RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Euclidean_distance)
- [Euclidean Angle - RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Euclidean_angle)
