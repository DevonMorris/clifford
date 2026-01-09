# PRD-6: Conformal GA (CGA) & Polish

**Status**: Pending
**Goal**: Conformal GA support and final documentation polish

## Background

Conformal Geometric Algebra embeds Euclidean space into a higher-dimensional
space where circles, spheres, and planes are represented uniformly as blades.
Transformations including rotations, translations, and dilations are all
represented as versors.

CGA adds two extra basis vectors: e+ (squares to +1) and e- (squares to -1).
Often written as e∞ (point at infinity) and e₀ (origin).

## Deliverables

### 1. CGA Signatures (`src/signature/cga.rs`)

```rust
/// 2D Conformal GA: Cl(3,1,0)
/// Embeds 2D Euclidean space
#[derive(Clone, Copy, Debug, Default)]
pub struct CGA2D;

impl Signature for CGA2D {
    const P: usize = 3;
    const Q: usize = 1;
    const R: usize = 0;
    const DIM: usize = 4;
    const NUM_BLADES: usize = 16;

    fn metric(i: usize) -> i8 {
        if i < 3 { 1 } else { -1 }
    }
}

/// 3D Conformal GA: Cl(4,1,0)
/// Embeds 3D Euclidean space
#[derive(Clone, Copy, Debug, Default)]
pub struct CGA3D;

impl Signature for CGA3D {
    const P: usize = 4;
    const Q: usize = 1;
    const R: usize = 0;
    const DIM: usize = 5;
    const NUM_BLADES: usize = 32;

    fn metric(i: usize) -> i8 {
        if i < 4 { 1 } else { -1 }
    }
}
```

### 2. CGA Utilities (`src/specialized/cga3d/`)

```rust
/// Null basis vectors
impl<T: Float> CGA3D {
    /// Origin point: e₀ = (e- - e+) / 2
    pub fn origin<T: Float>() -> Multivector<T, CGA3D>;

    /// Point at infinity: e∞ = e- + e+
    pub fn infinity<T: Float>() -> Multivector<T, CGA3D>;
}

/// Conformal point from Euclidean coordinates
pub fn point<T: Float>(x: T, y: T, z: T) -> Multivector<T, CGA3D> {
    // P = x*e1 + y*e2 + z*e3 + (x² + y² + z²)/2 * e∞ + e₀
}

/// Extract Euclidean coordinates from conformal point
pub fn to_euclidean<T: Float>(p: &Multivector<T, CGA3D>) -> Option<[T; 3]>;

/// Sphere from center and radius
pub fn sphere<T: Float>(cx: T, cy: T, cz: T, r: T) -> Multivector<T, CGA3D>;

/// Plane from normal and distance
pub fn plane<T: Float>(nx: T, ny: T, nz: T, d: T) -> Multivector<T, CGA3D>;

/// Circle (intersection of two spheres)
pub fn circle_from_spheres<T: Float>(
    s1: &Multivector<T, CGA3D>,
    s2: &Multivector<T, CGA3D>,
) -> Multivector<T, CGA3D>;
```

### 3. CGA Transformations

```rust
/// Translator in CGA
pub fn translator<T: Float>(dx: T, dy: T, dz: T) -> Multivector<T, CGA3D> {
    // T = 1 - (d/2) * e∞ where d = dx*e1 + dy*e2 + dz*e3
}

/// Dilator (uniform scaling) in CGA
pub fn dilator<T: Float>(scale: T) -> Multivector<T, CGA3D>;

/// General CGA versor transformation
pub fn transform<T: Float>(
    versor: &Multivector<T, CGA3D>,
    element: &Multivector<T, CGA3D>,
) -> Multivector<T, CGA3D> {
    versor.sandwich(element)
}
```

### 4. Documentation Polish

- [ ] All public items have comprehensive rustdoc
- [ ] Mathematical notation explained (∧, ·, ⟨⟩, ~, etc.)
- [ ] Geometric intuition for every operation
- [ ] Doc examples compile and run
- [ ] Module-level documentation with concepts
- [ ] Links to learning resources

### 5. Optional: SIMD (`src/simd/`)

Feature-gated SIMD support using `portable_simd` (nightly):

```rust
#[cfg(feature = "simd")]
pub mod simd {
    use std::simd::*;

    /// SIMD-accelerated f32x4 operations
    pub fn geometric_product_4x4(a: &[f32; 16], b: &[f32; 16]) -> [f32; 16];
}
```

## Files to Create

- `src/signature/cga.rs`
- `src/specialized/cga3d/mod.rs`
- `src/specialized/cga3d/embed.rs`
- `src/specialized/cga3d/objects.rs`
- `src/simd/mod.rs` (optional)

## Testing (proptest)

```rust
proptest! {
    #[test]
    fn conformal_point_roundtrip(x in -100.0..100.0f64, y in -100.0..100.0, z in -100.0..100.0) {
        let p = point(x, y, z);
        let [rx, ry, rz] = to_euclidean(&p).unwrap();
        prop_assert!((x - rx).abs() < 1e-10);
        prop_assert!((y - ry).abs() < 1e-10);
        prop_assert!((z - rz).abs() < 1e-10);
    }

    #[test]
    fn translation_versor(
        dx in -10.0..10.0f64,
        dy in -10.0..10.0,
        dz in -10.0..10.0,
        x in -100.0..100.0,
        y in -100.0..100.0,
        z in -100.0..100.0,
    ) {
        let t = translator(dx, dy, dz);
        let p = point(x, y, z);
        let tp = transform(&t, &p);
        let [rx, ry, rz] = to_euclidean(&tp).unwrap();
        prop_assert!((rx - (x + dx)).abs() < 1e-10);
        prop_assert!((ry - (y + dy)).abs() < 1e-10);
        prop_assert!((rz - (z + dz)).abs() < 1e-10);
    }
}
```

## Final Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - all tests pass
- [ ] `cargo clippy` - no warnings
- [ ] `cargo doc` - builds without warnings
- [ ] All doc examples run
- [ ] README updated with usage examples
- [ ] CLAUDE.md status updated

## Release Checklist

- [ ] Version bump in Cargo.toml
- [ ] CHANGELOG updated
- [ ] All PRDs marked complete
- [ ] `CARGO_REGISTRY_TOKEN` set
- [ ] Tag and push: `git tag v0.1.0 && git push --tags`
