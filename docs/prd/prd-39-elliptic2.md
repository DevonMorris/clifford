# PRD-39: 2D Elliptic Projective Geometry via Code Generation

## Overview

Implement 2D elliptic projective geometry as Cl(3,0,0) via the code generation pipeline.

## Algebra Details

**Signature**: Cl(3,0,0)
- p = 3 (three basis vectors squaring to +1)
- q = 0 (no basis vectors squaring to -1)
- r = 0 (no null/degenerate basis vectors)

**Dimension**: 2³ = 8 elements

**Basis Elements**:
| Grade | Blades | Count | Geometric Meaning |
|-------|--------|-------|-------------------|
| 0 | 1 | 1 | Scalar |
| 1 | e₁, e₂, e₃ | 3 | Points (directions on sphere) |
| 2 | e₁₂, e₁₃, e₂₃ | 3 | Lines (great circles) |
| 3 | e₁₂₃ | 1 | Pseudoscalar |

**Norm**: Uses **reverse** involution (non-degenerate with all positive signature).

## Geometric Interpretation

The elliptic plane is the projective plane RP² equipped with an elliptic metric. It can be visualized as the unit sphere S² with antipodal points identified.

Key properties:
- **No parallel lines**: Any two lines intersect at exactly one point
- **Points**: Represented as unit vectors in R³ (antipodal pairs on sphere)
- **Lines**: Great circles on the sphere, represented as bivectors
- **Distance**: Angle between vectors (geodesic distance on sphere)
- **Angle**: Angle between bivector planes

## Type Mapping

| Type | Grades | Fields | Description |
|------|--------|--------|-------------|
| Scalar | [0] | s | Scalar values |
| Point | [1] | x, y, z | Points on elliptic plane |
| Line | [2] | xy, xz, yz | Lines (great circles) |
| Pseudoscalar | [3] | xyz | Volume element |
| Rotor | [0, 2] | s, xy, xz, yz | Rotations |

Note: For a general multivector with all grades, use `crate::algebra::Multivector<T, Euclidean3>`.

## Applications

- **Spherical geometry**: Navigation on Earth, astronomy
- **Computer graphics**: Omnidirectional cameras, environment mapping
- **Robotics**: Orientation representation without gimbal lock
- **Projective geometry**: Study of incidence relations

## Implementation

### Files to Create/Modify

1. `algebras/elliptic2.toml` - Algebra specification
2. `src/specialized/elliptic/mod.rs` - Module root
3. `src/specialized/elliptic/dim2/mod.rs` - 2D elliptic module
4. `src/specialized/mod.rs` - Add elliptic export
5. `build.rs` - Add to ALGEBRAS array

### TOML Specification

```toml
[algebra]
name = "elliptic2"
description = "2D Elliptic Projective Geometry"

[signature]
positive = ["e1", "e2", "e3"]
negative = []
zero = []

[norm]
primary_involution = "reverse"

[types.Scalar]
grades = [0]
fields = ["s"]

[types.Point]
grades = [1]
fields = ["x", "y", "z"]

[types.Line]
grades = [2]
fields = ["xy", "xz", "yz"]

[types.Pseudoscalar]
grades = [3]
fields = ["xyz"]

[types.Rotor]
grades = [0, 2]
fields = ["s", "xy", "xz", "yz"]
```

## Relationship to Other Algebras

- **Same signature as Euclidean3**: Cl(3,0,0) is algebraically identical to 3D Euclidean GA
- **Different interpretation**: Euclidean3 treats grade-1 as spatial vectors; Elliptic2 treats them as projective points
- **Type names reflect geometry**: `Point`, `Line` instead of `Vector`, `Bivector`

## Verification

1. `cargo build` - Codegen and compilation succeed
2. `cargo nextest run` - All tests pass
3. `cargo clippy` - No warnings
4. `cargo doc --no-deps` - Documentation builds

## Test Plan

- Verify point-line incidence: `point.wedge(line) = 0` when point lies on line
- Verify line intersection: `line1.wedge(line2)` gives intersection point
- Verify rotor transformations preserve distances
- Property tests for algebraic identities
