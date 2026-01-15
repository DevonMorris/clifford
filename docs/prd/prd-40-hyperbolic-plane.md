# PRD-40: 2D Hyperbolic Projective Geometry via Code Generation

**Status**: Complete

## Overview

Implement 2D hyperbolic projective geometry as Cl(2,1,0) via the code generation pipeline.

## Algebra Details

**Signature**: Cl(2,1,0)
- p = 2 (two basis vectors squaring to +1: e₁, e₂)
- q = 1 (one basis vector squaring to -1: e₃)
- r = 0 (no null/degenerate basis vectors)

**Dimension**: 2³ = 8 elements

**Basis Elements**:
| Grade | Blades | Squares | Geometric Meaning |
|-------|--------|---------|-------------------|
| 0 | 1 | +1 | Scalar |
| 1 | e₁, e₂, e₃ | +1, +1, -1 | Points (homogeneous coords) |
| 2 | e₁₂, e₁₃, e₂₃ | -1, +1, +1 | Lines |
| 3 | e₁₂₃ | +1 | Pseudoscalar |

**Norm**: Uses **reverse** involution (non-degenerate but indefinite signature).
Implements `IndefiniteNormed` trait since norm² can be positive, negative, or zero.

## Geometric Interpretation

The hyperbolic plane (Lobachevsky plane) is a non-Euclidean geometry with constant negative curvature. It can be modeled using the hyperboloid model in Minkowski space.

Key properties:
- **Parallel postulate fails**: Through a point not on a line, infinitely many parallels exist
- **Angles of triangles**: Sum is always less than 180°
- **Geodesics**: Appear as hyperbolas in the Poincaré disk model
- **Ideal points**: Points "at infinity" on the boundary

### Causal Structure

Elements can have different causal characters based on their norm:
- **Timelike**: norm² < 0 (ordinary points inside the hyperboloid)
- **Null/Lightlike**: norm² = 0 (ideal points on the absolute)
- **Spacelike**: norm² > 0 (ultra-ideal points beyond infinity)

## Type Mapping

| Type | Grades | Fields | Description |
|------|--------|--------|-------------|
| Scalar | [0] | s | Scalar values |
| Point | [1] | x, y, t | Points (homogeneous: timelike=ordinary, null=ideal) |
| Line | [2] | xy, xt, yt | Geodesic lines |
| Pseudoscalar | [3] | xyt | Volume element |
| Rotor | [0, 2] | s, xy, xt, yt | Hyperbolic rotations and boosts |

## Applications

- **Non-Euclidean geometry**: Study of hyperbolic space
- **Special relativity**: Rapidity and Lorentz boosts
- **Hyperbolic tessellations**: Escher-like tilings
- **Complex analysis**: Poincaré disk, upper half-plane
- **Geodesy**: Surfaces with negative curvature

## Implementation

### Files to Create/Modify

1. `algebras/hyperbolic2.toml` - Algebra specification
2. `src/specialized/hyperbolic/dim2/mod.rs` - 2D hyperbolic plane module
3. `src/specialized/hyperbolic/mod.rs` - Update to export dim2
4. `build.rs` - Add to ALGEBRAS array

### TOML Specification

```toml
[algebra]
name = "hyperbolic2"
description = "2D Hyperbolic Projective Geometry"

[signature]
positive = ["e1", "e2"]
negative = ["e3"]
zero = []

[norm]
primary_involution = "reverse"

[types.Scalar]
grades = [0]
fields = ["s"]

[types.Point]
grades = [1]
fields = ["x", "y", "t"]

[types.Line]
grades = [2]
fields = ["xy", "xt", "yt"]

[types.Pseudoscalar]
grades = [3]
fields = ["xyt"]

[types.Rotor]
grades = [0, 2]
fields = ["s", "xy", "xt", "yt"]
```

## Relationship to Other Algebras

- **Hyperbolic numbers** (`hyperbolic/`): Cl(1,0,0) - split-complex numbers
- **Hyperbolic plane** (`hyperbolic/dim2/`): Cl(2,1,0) - this algebra
- **Elliptic plane** (`elliptic/dim2/`): Cl(3,0,0) - positive curvature analog
- **Minkowski plane** (`minkowski/dim2/`): Cl(1,1,0) - spacetime

The hyperbolic plane Cl(2,1,0) has one more positive dimension than Minkowski Cl(1,1,0),
giving it the structure needed for 2D hyperbolic geometry.

## Verification

1. `cargo build` - Codegen and compilation succeed
2. `cargo nextest run` - All tests pass
3. `cargo clippy` - No warnings
4. `cargo doc --no-deps` - Documentation builds

## Test Plan

- Verify indefinite norm properties (positive, negative, zero cases)
- Verify point-line incidence relations
- Verify rotor transformations preserve the metric
- Property tests for algebraic identities
