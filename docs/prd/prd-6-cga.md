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

### 2. CGA Utilities (`src/specialized/ceuclidean::dim3/`)

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

### 4. nalgebra Interoperability

When `nalgebra-0_33` or `nalgebra-0_34` feature is enabled, CGA types provide
conversions to/from nalgebra types:

#### Point Conversions

```rust
// CGA conformal point <-> nalgebra Point3 (via Euclidean extraction)
impl<T: Float + na::Scalar> From<na::Point3<T>> for ConformalPoint<T> {
    fn from(p: na::Point3<T>) -> Self {
        point(p.x, p.y, p.z)
    }
}

impl<T: Float + na::Scalar> TryFrom<ConformalPoint<T>> for na::Point3<T> {
    type Error = NalgebraConversionError;

    fn try_from(cp: ConformalPoint<T>) -> Result<Self, Self::Error> {
        let [x, y, z] = to_euclidean(&cp).ok_or(NalgebraConversionError::InvalidConformalPoint)?;
        Ok(na::Point3::new(x, y, z))
    }
}
```

#### Sphere Conversions

```rust
// Sphere <-> center point + radius
impl<T: Float + na::Scalar> From<Sphere<T>> for (na::Point3<T>, T) {
    fn from(s: Sphere<T>) -> Self {
        let (center, radius) = s.center_radius();
        (na::Point3::new(center[0], center[1], center[2]), radius)
    }
}

impl<T: Float + na::Scalar> From<(na::Point3<T>, T)> for Sphere<T> {
    fn from((center, radius): (na::Point3<T>, T)) -> Self {
        sphere(center.x, center.y, center.z, radius)
    }
}
```

#### Plane Conversions

```rust
// CGA Plane <-> nalgebra normal + distance
impl<T: Float + na::Scalar> From<CGAPlane<T>> for (na::Unit<na::Vector3<T>>, T) {
    fn from(p: CGAPlane<T>) -> Self {
        let (normal, dist) = p.normal_distance();
        (na::Unit::new_normalize(na::Vector3::new(normal[0], normal[1], normal[2])), dist)
    }
}
```

#### Transformation Conversions

```rust
// CGA versors can produce nalgebra transformations
impl<T: Float + na::RealField> CGATranslator<T> {
    /// Extract as nalgebra Translation3
    pub fn to_nalgebra(&self) -> na::Translation3<T> {
        let [dx, dy, dz] = self.displacement();
        na::Translation3::new(dx, dy, dz)
    }
}

impl<T: Float + na::RealField> CGADilator<T> {
    /// Extract as nalgebra uniform scale factor
    pub fn scale_factor(&self) -> T { ... }
}
```

#### Error Type Extension

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NalgebraConversionError {
    InvalidConformalPoint,  // Not a valid conformal point embedding
    InvalidSphere,          // Imaginary radius or degenerate
    InvalidPlane,           // Zero normal vector
}
```

### 5. Documentation Polish

- [ ] All public items have comprehensive rustdoc
- [ ] Mathematical notation explained (∧, ·, ⟨⟩, ~, etc.)
- [ ] Geometric intuition for every operation
- [ ] Doc examples compile and run
- [ ] Module-level documentation with concepts
- [ ] Links to learning resources

### 6. Optional: SIMD (`src/simd/`)

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
- `src/specialized/ceuclidean::dim3/mod.rs`
- `src/specialized/ceuclidean::dim3/embed.rs`
- `src/specialized/ceuclidean::dim3/objects.rs`
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

## Benchmarks

CGA operations should be benchmarked against:
1. Generic `Multivector` implementation
2. nalgebra equivalents where applicable
3. Traditional approaches (e.g., sphere intersection via algebraic methods)

### Benchmark File (`benches/cga.rs`)

```rust
//! Benchmarks for CGA operations.
//!
//! Run with: `cargo bench --bench cga`

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use clifford::specialized::cga3d::{point, sphere, plane, translator, dilator, transform, to_euclidean};

// Point embedding/extraction
fn bench_point_embed(c: &mut Criterion) {
    c.bench_function("cga3d/point_embed", |bencher| {
        bencher.iter(|| point(black_box(1.0), black_box(2.0), black_box(3.0)))
    });
}

fn bench_point_extract(c: &mut Criterion) {
    let p = point(1.0, 2.0, 3.0);

    c.bench_function("cga3d/point_extract", |bencher| {
        bencher.iter(|| to_euclidean(&black_box(p)))
    });
}

// Transformation versors
fn bench_translator_create(c: &mut Criterion) {
    c.bench_function("cga3d/translator_create", |bencher| {
        bencher.iter(|| translator(black_box(1.0), black_box(2.0), black_box(3.0)))
    });
}

fn bench_translator_apply(c: &mut Criterion) {
    let t = translator(1.0, 2.0, 3.0);
    let p = point(0.0, 0.0, 0.0);

    c.bench_function("cga3d/translator_apply", |bencher| {
        bencher.iter(|| transform(&black_box(t), &black_box(p)))
    });
}

fn bench_dilator_apply(c: &mut Criterion) {
    let d = dilator(2.0);
    let p = point(1.0, 1.0, 1.0);

    c.bench_function("cga3d/dilator_apply", |bencher| {
        bencher.iter(|| transform(&black_box(d), &black_box(p)))
    });
}

// Geometric object operations
fn bench_sphere_create(c: &mut Criterion) {
    c.bench_function("cga3d/sphere_create", |bencher| {
        bencher.iter(|| sphere(black_box(0.0), black_box(0.0), black_box(0.0), black_box(1.0)))
    });
}

fn bench_sphere_intersect(c: &mut Criterion) {
    let s1 = sphere(0.0, 0.0, 0.0, 2.0);
    let s2 = sphere(1.0, 0.0, 0.0, 2.0);

    c.bench_function("cga3d/sphere_intersect_circle", |bencher| {
        bencher.iter(|| black_box(s1).meet(&black_box(s2)))
    });
}

fn bench_plane_sphere_intersect(c: &mut Criterion) {
    let plane = plane(0.0, 0.0, 1.0, 0.0); // z = 0
    let s = sphere(0.0, 0.0, 0.0, 1.0);

    c.bench_function("cga3d/plane_sphere_intersect", |bencher| {
        bencher.iter(|| black_box(plane).meet(&black_box(s)))
    });
}

criterion_group!(
    benches,
    bench_point_embed,
    bench_point_extract,
    bench_translator_create,
    bench_translator_apply,
    bench_dilator_apply,
    bench_sphere_create,
    bench_sphere_intersect,
    bench_plane_sphere_intersect,
);
criterion_main!(benches);
```

### nalgebra Comparison (`benches/cga_nalgebra.rs`)

Compare CGA operations with traditional nalgebra approaches:

```rust
#![cfg(any(feature = "nalgebra-0_33", feature = "nalgebra-0_34"))]

use clifford::specialized::cga3d::*;

// Translation: CGA translator vs nalgebra Translation3
fn bench_translate_cga(c: &mut Criterion) {
    let t = translator(1.0, 2.0, 3.0);
    let p = point(0.0, 0.0, 0.0);

    c.bench_function("comparison/translate/cga", |bencher| {
        bencher.iter(|| {
            let tp = transform(&black_box(t), &black_box(p));
            to_euclidean(&tp)
        })
    });
}

fn bench_translate_nalgebra(c: &mut Criterion) {
    let t = na::Translation3::new(1.0, 2.0, 3.0);
    let p = na::Point3::new(0.0, 0.0, 0.0);

    c.bench_function("comparison/translate/nalgebra", |bencher| {
        bencher.iter(|| black_box(t) * black_box(p))
    });
}

// Sphere-sphere intersection: CGA meet vs algebraic
fn bench_sphere_intersect_cga(c: &mut Criterion) {
    let s1 = sphere(0.0, 0.0, 0.0, 2.0);
    let s2 = sphere(1.0, 0.0, 0.0, 2.0);

    c.bench_function("comparison/sphere_intersect/cga", |bencher| {
        bencher.iter(|| black_box(s1).meet(&black_box(s2)))
    });
}

fn bench_sphere_intersect_algebraic(c: &mut Criterion) {
    // Traditional approach: solve quadratic system
    let c1 = na::Point3::new(0.0, 0.0, 0.0);
    let r1 = 2.0f64;
    let c2 = na::Point3::new(1.0, 0.0, 0.0);
    let r2 = 2.0f64;

    c.bench_function("comparison/sphere_intersect/algebraic", |bencher| {
        bencher.iter(|| {
            // Compute intersection circle center and radius
            let d = (black_box(c2) - black_box(c1)).norm();
            let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
            let h = (r1 * r1 - a * a).sqrt();
            let center = c1 + (c2 - c1).normalize() * a;
            (center, h) // circle center and radius
        })
    });
}

// Point-on-sphere test: CGA inner product vs distance check
fn bench_point_on_sphere_cga(c: &mut Criterion) {
    let s = sphere(0.0, 0.0, 0.0, 1.0);
    let p = point(1.0, 0.0, 0.0);

    c.bench_function("comparison/point_on_sphere/cga", |bencher| {
        bencher.iter(|| {
            // In CGA, point lies on sphere if inner product is zero
            black_box(p).inner(&black_box(s)).norm() < 1e-10
        })
    });
}

fn bench_point_on_sphere_algebraic(c: &mut Criterion) {
    let center = na::Point3::new(0.0, 0.0, 0.0);
    let radius = 1.0f64;
    let p = na::Point3::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/point_on_sphere/algebraic", |bencher| {
        bencher.iter(|| {
            ((black_box(p) - black_box(center)).norm() - radius).abs() < 1e-10
        })
    });
}
```

### Expected Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| Point embed | < 5 ns | Simple arithmetic |
| Point extract | < 5 ns | Division + extraction |
| Translator create | < 5 ns | Vector + scalar construction |
| Translator apply | < 20 ns | Sandwich product |
| Dilator apply | < 20 ns | Sandwich product |
| Sphere intersect | < 15 ns | Meet operation |
| CGA vs algebraic | Comparable | CGA may have overhead but more uniform |

### Key Insights

CGA benchmarks should reveal:
1. **Embedding overhead**: Cost of moving between Euclidean and conformal representations
2. **Uniformity benefit**: Single meet/join vs case-by-case geometric algorithms
3. **Composition advantage**: Chaining CGA versors vs composing affine transforms

### Benchmark README Section

Add to `benches/README.md`:

```markdown
## CGA Benchmarks

Operations on specialized CGA types.

### Point Operations
![cga3d_point_embed](reports/cga3d_point_embed_pdf.svg)
![cga3d_point_extract](reports/cga3d_point_extract_pdf.svg)

### Transformation Versors
![cga3d_translator_apply](reports/cga3d_translator_apply_pdf.svg)
![cga3d_dilator_apply](reports/cga3d_dilator_apply_pdf.svg)

### Geometric Intersections
![cga3d_sphere_intersect](reports/cga3d_sphere_intersect_circle_pdf.svg)

### Comparison with Traditional Methods
| Operation | CGA | Traditional | Notes |
|-----------|-----|-------------|-------|
| Translate point | ~X ns | ~Y ns | CGA has embed/extract overhead |
| Sphere intersect | ~X ns | ~Y ns | CGA is uniform, traditional is specialized |
| Point-on-sphere | ~X ns | ~Y ns | CGA uses inner product |
```

## Final Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - all tests pass
- [ ] `cargo clippy` - no warnings
- [ ] `cargo doc` - builds without warnings
- [ ] All doc examples run
- [ ] README updated with usage examples
- [ ] CLAUDE.md status updated
- [ ] nalgebra conversions tested with feature flags
- [ ] Benchmarks run and documented

## Release Checklist

- [ ] Version bump in Cargo.toml
- [ ] CHANGELOG updated
- [ ] All PRDs marked complete
- [ ] `CARGO_REGISTRY_TOKEN` set
- [ ] Tag and push: `git tag v0.1.0 && git push --tags`
