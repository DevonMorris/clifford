# PRD-41: 3D Minkowski Spacetime Cl(3,1,0)

**Status**: Complete

## Summary

Implement 3+1 Minkowski spacetime algebra Cl(3,1,0) via code generation. This is the geometric algebra of special relativity with 3 spatial dimensions and 1 temporal dimension.

## Motivation

Minkowski spacetime is fundamental to special relativity and provides a natural framework for:
- Lorentz transformations and boosts
- Relativistic mechanics
- Electromagnetic field representations
- Relativistic particle physics

The algebra Cl(3,1,0) (also known as the spacetime algebra or STA) encodes the metric signature of special relativity where spatial vectors square to +1 and the temporal vector squares to -1.

## Algebra Structure

**Signature**: Cl(3,1,0) - 3 positive, 1 negative, 0 zero

| Grade | Blades | Squares | Count | Description |
|-------|--------|---------|-------|-------------|
| 0 | 1 | +1 | 1 | Scalar |
| 1 | e₁, e₂, e₃, e₄ | +1, +1, +1, -1 | 4 | Spacetime vectors |
| 2 | e₁₂, e₁₃, e₁₄, e₂₃, e₂₄, e₃₄ | -1, -1, +1, -1, +1, +1 | 6 | Bivectors |
| 3 | e₁₂₃, e₁₂₄, e₁₃₄, e₂₃₄ | -1, +1, +1, +1 | 4 | Trivectors |
| 4 | e₁₂₃₄ | -1 | 1 | Pseudoscalar |

**Total**: 16 basis elements (2⁴)

## Causal Structure

The indefinite metric creates a causal structure for vectors:

| Condition | Type | Physical Interpretation |
|-----------|------|------------------------|
| v² < 0 | Timelike | Massive particles, observers |
| v² = 0 | Null/Lightlike | Light rays, photons |
| v² > 0 | Spacelike | Spatial separations |

Note: Convention follows (-,+,+,+) signature where e₄² = -1 is the temporal direction.

## Involution

**Primary involution**: `reverse`

The algebra is non-degenerate (r=0), so we use the reverse operation for computing norms:
```
reverse(a + be₁ + ce₁₂ + de₁₂₃ + ee₁₂₃₄) = a + be₁ - ce₁₂ - de₁₂₃ + ee₁₂₃₄
```

The norm is computed as: `norm² = x * reverse(x)`

This produces an **indefinite norm** that can be positive, negative, or zero.

## Types

### Grade-1: Vector (Spacetime Vector)
- **Fields**: `x`, `y`, `z`, `t` (spatial + temporal components)
- **Interpretation**: Events, spacetime displacements, 4-velocities
- **Causal character**: Timelike (t² > x² + y² + z²), null, or spacelike

### Grade-2: Bivector
- **Fields**: `xy`, `xz`, `xt`, `yz`, `yt`, `zt`
- **Interpretation**: Spacetime planes, electromagnetic field tensor
- **Structure**:
  - Spatial bivectors (xy, xz, yz): Magnetic field components
  - Spacetime bivectors (xt, yt, zt): Electric field components

### Grade-3: Trivector
- **Fields**: `xyz`, `xyt`, `xzt`, `yzt`
- **Interpretation**: Spacetime volumes, current densities

### Grade-4: Pseudoscalar
- **Fields**: `xyzt`
- **Interpretation**: Oriented spacetime volume element

### Even Subalgebra: Eventor (grades 0, 2, 4)
- **Fields**: `s`, `xy`, `xz`, `xt`, `yz`, `yt`, `zt`, `xyzt`
- **Interpretation**: Lorentz transformations (rotations + boosts)

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

## Implementation Plan

1. Create `algebras/minkowski3.toml` with algebra specification
2. Add `Cl3_1_0` signature type to `src/signature/euclidean.rs`
3. Create module structure at `src/specialized/minkowski/dim3/`
4. Update `build.rs` to include minkowski3 algebra
5. Export types from `src/specialized/minkowski/mod.rs`
6. Run `cargo build` to generate code

## File Structure

```
src/specialized/minkowski/
  mod.rs              # Module root (update to include dim3)
  dim2/               # Existing 2D Minkowski
  dim3/               # New 3D Minkowski
    mod.rs            # Documentation and re-exports
    generated/        # Auto-generated code
      types.rs
      traits.rs
```

## Testing

- Verify Lorentz boost composition
- Test electromagnetic field duality (F ↔ *F)
- Property tests for norm signs matching causal character
- Verify metric signature: e₁² = e₂² = e₃² = +1, e₄² = -1

## References

- [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
- [Geometric Algebra for Physicists](https://www.cambridge.org/core/books/geometric-algebra-for-physicists/FB8D3ACB76AB3AB10BA7F27505925091)
- Hestenes, D. "Space-Time Algebra" (1966)
