# PRD-42: 3D Conformal Geometric Algebra Cl(4,1,0)

**Status**: Draft
**Goal**: Implement 3D Conformal Geometric Algebra via code generation

## Summary

Implement 3D Conformal Geometric Algebra (CGA) as Cl(4,1,0) via the code generation pipeline. CGA embeds 3D Euclidean space into a 5D algebra, enabling elegant representation of spheres, circles, planes, and conformal transformations.

## Motivation

Conformal Geometric Algebra is powerful for:
- **Round geometry**: Spheres, circles, and their intersections
- **Conformal transformations**: Rotations, translations, dilations, inversions
- **Computer graphics**: Ray tracing, collision detection
- **Robotics**: Inverse kinematics, motion planning
- **Physics**: Electromagnetism, quantum mechanics

The algebra Cl(4,1,0) extends 3D Euclidean space with two additional basis vectors to represent the origin and point at infinity.

## Algebra Structure

**Signature**: Cl(4,1,0)
- p = 4 (four basis vectors squaring to +1: e₁, e₂, e₃, e₄)
- q = 1 (one basis vector squaring to -1: e₅)
- r = 0 (non-degenerate)

**Dimension**: 2⁵ = 32 basis elements

| Grade | Count | Description |
|-------|-------|-------------|
| 0 | 1 | Scalar |
| 1 | 5 | Vectors (RoundPoint) |
| 2 | 10 | Bivectors (Dipole, FlatPoint) |
| 3 | 10 | Trivectors (Circle, Line) |
| 4 | 5 | Quadrivectors (Sphere, Plane) |
| 5 | 1 | Pseudoscalar |

## Basis Convention

Following the Conformal GA wiki (conformalgeometricalgebra.org):

- **e₁, e₂, e₃**: Euclidean basis vectors (spatial directions)
- **e₄**: Origin representation (e₄ = ½(e₋ - e₊))
- **e₅**: Point at infinity (e₅ = e₋ + e₊)

Where e₋² = -1 and e₊² = +1 are the conformal basis vectors.

## Involution

**Primary involution**: `reverse`

The algebra is non-degenerate (r=0), so we use the reverse operation for computing norms. Since the signature is indefinite (4 positive, 1 negative), the norm can be positive, negative, or zero.

## Types (from conformalgeometricalgebra.org)

### Grade-1: RoundPoint (5 components)
A round point encodes both position and curvature:
```
a = aₓe₁ + aᵧe₂ + a_z e₃ + a_w e₄ + a_u e₅
```
- **Fields**: `x`, `y`, `z`, `w`, `u`
- **Geometric meaning**: 3D point with radius information

### Grade-2: Dipole (10 components)
A dipole represents a point pair or oriented line segment:
```
d = vₓe₁₄ + vᵧe₂₄ + v_z e₃₄ + mₓe₂₃ + mᵧe₃₁ + m_z e₁₂ + pₓe₁₅ + pᵧe₂₅ + p_z e₃₅ + p_w e₄₅
```
- **Fields**: `vx`, `vy`, `vz`, `mx`, `my`, `mz`, `px`, `py`, `pz`, `pw`
- **Subtype FlatPoint** (4 components): When v=0 and m=0

### Grade-3: Circle (10 components)
A circle in 3D space:
```
c = cᵍₓe₄₂₃ + cᵍᵧe₄₃₁ + cᵍ_z e₄₁₂ + cᵍ_w e₃₂₁ + cᵛₓe₄₁₅ + cᵛᵧe₄₂₅ + cᵛ_z e₄₃₅ + cᵐₓe₂₃₅ + cᵐᵧe₃₁₅ + cᵐ_z e₁₂₅
```
- **Fields**: `gx`, `gy`, `gz`, `gw`, `vx`, `vy`, `vz`, `mx`, `my`, `mz`
- **Subtype Line** (6 components): When g=0 (circle through infinity)

### Grade-4: Sphere (5 components)
A sphere in 3D space:
```
s = s_u e₁₂₃₄ + sₓe₄₂₃₅ + sᵧe₄₃₁₅ + s_z e₄₁₂₅ + s_w e₃₂₁₅
```
- **Fields**: `u`, `x`, `y`, `z`, `w`
- **Subtype Plane** (4 components): When u=0 (sphere through infinity)

### Grade-5: Pseudoscalar (1 component)
The oriented 5-volume element:
- **Fields**: `xyzwu`

### Even Subalgebra: Motor (16 components)
Grades 0, 2, 4 for conformal transformations:
- **Fields**: scalar + 10 bivector + 5 quadrivector components

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

1. Create `algebras/conformal3.toml` with algebra specification
2. Add `Cl4_1_0` signature type to `src/signature/euclidean.rs`
3. Create module structure at `src/specialized/conformal/dim3/`
4. Update `build.rs` to include conformal3 algebra
5. Export types from `src/specialized/conformal/mod.rs`
6. Run `cargo build` to generate code

## File Structure

```
src/specialized/conformal/
  mod.rs              # Module root
  dim3/               # 3D Conformal GA
    mod.rs            # Documentation and re-exports
    generated/        # Auto-generated code
      types.rs
      traits.rs
      conversions.rs
```

## Testing

- Verify null vectors have zero norm (actual CGA points)
- Test sphere-sphere intersection produces circle
- Test plane-sphere intersection produces circle
- Verify conformal transformations preserve null cone

## References

- [Conformal Geometric Algebra Wiki](https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page)
- [Geometric Algebra for Computer Science](https://geometricalgebra.org/)
- Dorst, Fontijne, Mann - "Geometric Algebra for Computer Science"
