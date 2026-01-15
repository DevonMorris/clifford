# PRD-30: Projections and Antiprojections

**Status**: Draft
**Goal**: Implement projection and antiprojection operations for all geometric algebras via code generation

## Background

In Geometric Algebra, projections and antiprojections are fundamental operations that allow computing the closest geometric relationship between objects. These operations are essential for:

- Finding the closest point on a line/plane to a given point
- Computing perpendicular distances
- Projecting geometry onto other geometry
- Many practical robotics and graphics applications

These operations should work across **all algebras** (Euclidean, Projective, Conformal, etc.) where the prerequisite operations are defined.

Reference: [RGA Wiki - Projections](https://rigidgeometricalgebra.org/wiki/index.php?title=Projections)

## Mathematical Definitions

### Projection

The **orthogonal projection** of object **a** onto object **b** is defined as:

```
proj_b(a) = b ∨ (a ∧ b☆)
```

Where:
- `∨` is the antiwedge (regressive) product
- `∧` is the wedge (exterior) product
- `b☆` is the weight dual of **b**

This can be interpreted as the "weight expansion of **a** into **b**".

### Antiprojection

The **orthogonal antiprojection** of object **a** onto object **b** is defined as:

```
antiproj_b(a) = b ∧ (a ∨ b☆)
```

Where:
- `∧` is the wedge (exterior) product
- `∨` is the antiwedge (regressive) product
- `b☆` is the weight dual of **b**

This can be interpreted as the "weight contraction of **a** with **b**".

### Geometric Interpretations

| Operation | Geometric Meaning |
|-----------|-------------------|
| `proj_line(point)` | Closest point on the line to the given point |
| `proj_plane(point)` | Closest point on the plane to the given point |
| `proj_plane(line)` | The line's projection onto the plane |
| `antiproj_point(plane)` | Plane through the point perpendicular to the original plane's normal |
| `antiproj_line(point)` | Plane containing the line and passing through the point |

### Key Properties

1. **Idempotence**: Projecting onto the same geometry twice gives the same result
2. **Orthogonality**: The difference between an object and its projection is orthogonal to the target
3. **Grade preservation**: Projection of a point onto a line/plane yields a point

## Prerequisites

The following operations must already be implemented (they are):
- Wedge product (`∧`) - `Wedge` trait
- Antiwedge product (`∨`) - `Antiwedge` trait
- Weight dual (`☆`) - `WeightDual` trait

## API Design

### Trait Definitions

```rust
/// Projection of self onto target geometry.
///
/// The projection finds the part of `self` that lies on `target`.
/// For a point projected onto a line or plane, this gives the closest
/// point on that geometry.
///
/// Formula: `target ∨ (self ∧ target☆)`
pub trait Project<Target> {
    /// The output type of the projection.
    type Output;

    /// Projects `self` onto `target`.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::ops::Project;
    /// use clifford::specialized::projective::dim3::{Point, Line};
    ///
    /// let point = Point::from_cartesian(1.0, 2.0, 3.0);
    /// let line = Line::x_axis();
    /// let closest = point.project(&line);  // Closest point on x-axis
    /// ```
    fn project(&self, target: &Target) -> Self::Output;
}

/// Antiprojection of self onto target geometry.
///
/// The antiprojection finds the geometry through `self` that is
/// perpendicular to `target`.
///
/// Formula: `target ∧ (self ∨ target☆)`
pub trait Antiproject<Target> {
    /// The output type of the antiprojection.
    type Output;

    /// Antiprojection of `self` onto `target`.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use clifford::ops::Antiproject;
    /// use clifford::specialized::projective::dim3::{Point, Plane};
    ///
    /// let point = Point::from_cartesian(1.0, 2.0, 3.0);
    /// let plane = Plane::xy();
    /// let perp_plane = point.antiproject(&plane);  // Plane through point, perpendicular to xy
    /// ```
    fn antiproject(&self, target: &Target) -> Self::Output;
}
```

### Type Combinations for 3D PGA

| Self | Target | `project()` Output | Geometric Meaning |
|------|--------|-------------------|-------------------|
| Point | Line | Point | Closest point on line |
| Point | Plane | Point | Closest point on plane |
| Line | Plane | Line | Line projected onto plane |

| Self | Target | `antiproject()` Output | Geometric Meaning |
|------|--------|------------------------|-------------------|
| Point | Line | Plane | Plane through point containing line |
| Point | Plane | Line | Line through point perpendicular to plane |
| Line | Point | Plane | Plane containing line and point |

### Type Combinations for 2D PGA

| Self | Target | `project()` Output | Geometric Meaning |
|------|--------|-------------------|-------------------|
| Point | Line | Point | Closest point on line |

| Self | Target | `antiproject()` Output | Geometric Meaning |
|------|--------|------------------------|-------------------|
| Point | Line | Line | Line through point perpendicular to target |

## Implementation Plan

Since Project and Antiproject are composite operations built from existing products (Wedge, Antiwedge, WeightDual), we implement them via **code generation** to work across all algebras.

### Step 1: Add Traits to `ops.rs`

Add `Project` and `Antiproject` traits to `src/ops.rs` with full documentation:

```rust
/// Projection of self onto target geometry.
///
/// Formula: `target ∨ (self ∧ target☆)`
pub trait Project<Target> {
    type Output;
    fn project(&self, target: &Target) -> Self::Output;
}

/// Antiprojection of self onto target geometry.
///
/// Formula: `target ∧ (self ∨ target☆)`
pub trait Antiproject<Target> {
    type Output;
    fn antiproject(&self, target: &Target) -> Self::Output;
}
```

### Step 2: Add Product Types to Codegen Discovery

In `crates/clifford-codegen/src/discovery/products.rs`:
- Add `Project` and `Antiproject` to `ProductType` enum
- Implement grade inference for both operations

### Step 3: Add to ProductsSpec

In `crates/clifford-codegen/src/spec/ir.rs`:
- Add `project: Vec<ProductEntry>` and `antiproject: Vec<ProductEntry>` to `ProductsSpec`

### Step 4: Add Symbolic Product Kinds

In `crates/clifford-codegen/src/symbolic/product.rs`:
- Add `Project` and `Antiproject` to `ProductKind` enum
- Implement symbolic computation as composite of existing operations

### Step 5: Add Trait Generation

In `crates/clifford-codegen/src/codegen/traits.rs`:
- Add `generate_project_trait()` and `generate_antiproject_trait()` methods
- Generate implementations for all valid type pairs

### Step 6: Update Parser

In `crates/clifford-codegen/src/spec/parser.rs`:
- Add project/antiproject inference in `infer_products_from_types()`

### Step 7: Regenerate All Algebras

Run build to regenerate all algebra code with new Project/Antiproject implementations.

### Step 8: Add Tests

Property-based tests to verify:
1. Projection is idempotent: `proj(proj(a, b), b) == proj(a, b)`
2. Projection lies on target geometry
3. Antiprojection contains original geometry

### Step 9: Update Documentation

- Add examples to trait documentation
- Document which type combinations are valid for each algebra

## Dependencies

The implementation requires these existing traits:
- `Wedge` - exterior product (already implemented)
- `Antiwedge` - regressive product (already implemented)
- `WeightDual` - weight dual operation (already implemented)

## Testing Strategy

### Unit Tests

```rust
#[test]
fn point_projection_onto_line() {
    let p = Point::from_cartesian(1.0, 1.0, 0.0);
    let line = Line::x_axis();
    let projected = p.project(&line);

    // Projected point should be on x-axis
    assert!(abs_diff_eq!(projected.y(), 0.0, epsilon = 1e-10));
    assert!(abs_diff_eq!(projected.z(), 0.0, epsilon = 1e-10));
    // x-coordinate preserved
    assert!(abs_diff_eq!(projected.x(), 1.0, epsilon = 1e-10));
}

#[test]
fn point_projection_onto_plane() {
    let p = Point::from_cartesian(1.0, 2.0, 3.0);
    let plane = Plane::xy();  // z = 0 plane
    let projected = p.project(&plane);

    // Projected point should be on xy-plane
    assert!(abs_diff_eq!(projected.z(), 0.0, epsilon = 1e-10));
    // x, y coordinates preserved
    assert!(abs_diff_eq!(projected.x(), 1.0, epsilon = 1e-10));
    assert!(abs_diff_eq!(projected.y(), 2.0, epsilon = 1e-10));
}
```

### Property-Based Tests

```rust
proptest! {
    #[test]
    fn projection_is_idempotent(
        p in any::<Point<f64>>(),
        line in any::<Line<f64>>(),
    ) {
        let proj1 = p.project(&line);
        let proj2 = proj1.project(&line);
        prop_assert!(relative_eq!(proj1, proj2, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn projection_is_closest_point(
        p in any::<Point<f64>>(),
        line in any::<Line<f64>>(),
    ) {
        let projected = p.project(&line);
        // The projected point should be on the line
        // (distance from projected point to line should be ~0)
        let dist = line.distance_to_point(&projected);
        prop_assert!(dist < RELATIVE_EQ_EPS);
    }
}
```

## Success Criteria

1. `Project` and `Antiproject` traits defined in `ops.rs`
2. Implementations for all relevant 3D PGA type combinations
3. Implementations for all relevant 2D PGA type combinations
4. Comprehensive test coverage including property-based tests
5. Documentation with geometric explanations and examples
6. All verification passes: `cargo fmt && cargo clippy && cargo nextest run && cargo doc --no-deps`

## Future Enhancements

1. **Codegen integration**: Automatically generate projection implementations for all valid type pairs
2. **Distance functions**: Use projections to compute distances between geometry
3. **Reflection operations**: Implement reflections using projections
4. **Rejection operation**: `a - proj_b(a)` gives the component of `a` perpendicular to `b`

## References

- [RGA Wiki - Projections](https://rigidgeometricalgebra.org/wiki/index.php?title=Projections)
- [RGA Wiki - Weight Expansion](https://rigidgeometricalgebra.org/wiki/index.php?title=Weight_expansion)
- [RGA Wiki - Weight Contraction](https://rigidgeometricalgebra.org/wiki/index.php?title=Weight_contraction)
- [RGA Wiki - Weight Dual](https://rigidgeometricalgebra.org/wiki/index.php?title=Weight_dual)
