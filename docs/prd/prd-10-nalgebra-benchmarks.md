# PRD-10: nalgebra vs clifford Benchmarks

**Status**: Pending
**Goal**: Comprehensive benchmarks comparing clifford specialized types with nalgebra equivalents

## Motivation

With nalgebra interoperability implemented (PRD-8), users need to understand the
performance characteristics when choosing between:

1. Pure clifford operations
2. Pure nalgebra operations
3. Converting between libraries for specific operations

This PRD establishes benchmarks to measure:
- Raw operation speed comparisons
- Conversion overhead costs
- Realistic workflow scenarios

## Benchmark Categories

### 1. Vector Operations

Compare `dim2::Vector`/`dim3::Vector` vs `na::Vector2`/`na::Vector3`:

| Operation | clifford | nalgebra | Notes |
|-----------|----------|----------|-------|
| Addition | `a + b` | `a + b` | Same semantics |
| Dot product | `a.dot(b)` | `a.dot(&b)` | Same semantics |
| Cross product (3D) | `a.cross(b)` | `a.cross(&b)` | Same semantics |
| Normalization | `a.normalized()` | `a.normalize()` | Same semantics |
| Scalar multiply | `a * s` | `a * s` | Same semantics |

```rust
fn bench_vec3_dot_clifford(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 2.0, 3.0);
    let b = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("comparison/vec3_dot/clifford", |bencher| {
        bencher.iter(|| black_box(a).dot(black_box(b)))
    });
}

fn bench_vec3_dot_nalgebra(c: &mut Criterion) {
    let a = na::Vector3::new(1.0, 2.0, 3.0);
    let b = na::Vector3::new(4.0, 5.0, 6.0);

    c.bench_function("comparison/vec3_dot/nalgebra", |bencher| {
        bencher.iter(|| black_box(a).dot(&black_box(b)))
    });
}
```

### 2. Rotation Operations

Compare `dim3::Rotor` vs `na::UnitQuaternion` and `na::Rotation3`:

| Operation | clifford | nalgebra (quat) | nalgebra (rot) |
|-----------|----------|-----------------|----------------|
| Create from angle/axis | `Rotor::from_angle_plane()` | `UnitQuaternion::from_axis_angle()` | `Rotation3::from_axis_angle()` |
| Rotate vector | `r.rotate(v)` | `q * v` | `r * v` |
| Compose rotations | `r1.compose(r2)` | `q1 * q2` | `r1 * r2` |
| Inverse | `r.inverse()` | `q.inverse()` | `r.inverse()` |
| Slerp | `r1.slerp(r2, t)` | `q1.slerp(&q2, t)` | N/A (via quat) |

```rust
fn bench_rotation_clifford(c: &mut Criterion) {
    let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());
    let v = dim3::Vector::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rotate_vec3/clifford", |bencher| {
        bencher.iter(|| black_box(rotor).rotate(black_box(v)))
    });
}

fn bench_rotation_quaternion(c: &mut Criterion) {
    let q = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let v = na::Vector3::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rotate_vec3/nalgebra_quat", |bencher| {
        bencher.iter(|| black_box(q) * black_box(v))
    });
}

fn bench_rotation_matrix(c: &mut Criterion) {
    let r = na::Rotation3::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let v = na::Vector3::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rotate_vec3/nalgebra_rot", |bencher| {
        bencher.iter(|| black_box(r) * black_box(v))
    });
}
```

### 3. 2D Rotation Operations

Compare `dim2::Rotor` vs `na::Rotation2` and `na::UnitComplex`:

| Operation | clifford | nalgebra (rot) | nalgebra (complex) |
|-----------|----------|----------------|-------------------|
| Create from angle | `Rotor::from_angle()` | `Rotation2::new()` | `UnitComplex::new()` |
| Rotate vector | `r.rotate(v)` | `r * v` | `c * v` |
| Compose | `r1.compose(r2)` | `r1 * r2` | `c1 * c2` |

### 4. Conversion Overhead

Measure the cost of converting between clifford and nalgebra types:

```rust
fn bench_vec3_to_nalgebra(c: &mut Criterion) {
    let v = dim3::Vector::new(1.0, 2.0, 3.0);

    c.bench_function("conversion/vec3_to_nalgebra", |bencher| {
        bencher.iter(|| {
            let na_v: na::Vector3<f64> = black_box(v).into();
            black_box(na_v)
        })
    });
}

fn bench_vec3_from_nalgebra(c: &mut Criterion) {
    let v = na::Vector3::new(1.0, 2.0, 3.0);

    c.bench_function("conversion/vec3_from_nalgebra", |bencher| {
        bencher.iter(|| {
            let cliff_v: dim3::Vector<f64> = black_box(v).into();
            black_box(cliff_v)
        })
    });
}

fn bench_rotor_to_quaternion(c: &mut Criterion) {
    let r = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());

    c.bench_function("conversion/rotor_to_quaternion", |bencher| {
        bencher.iter(|| {
            let q: na::UnitQuaternion<f64> = black_box(r).into();
            black_box(q)
        })
    });
}

fn bench_quaternion_to_rotor(c: &mut Criterion) {
    let q = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);

    c.bench_function("conversion/quaternion_to_rotor", |bencher| {
        bencher.iter(|| {
            let r: dim3::Rotor<f64> = black_box(q).into();
            black_box(r)
        })
    });
}
```

### 5. Realistic Workflows

Benchmark common usage patterns that involve multiple operations:

```rust
/// Scenario: Transform multiple points with a rotation
fn bench_batch_rotation_clifford(c: &mut Criterion) {
    let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());
    let points: Vec<dim3::Vector<f64>> = (0..100)
        .map(|i| dim3::Vector::new(i as f64, 0.0, 0.0))
        .collect();

    c.bench_function("workflow/batch_rotate_100/clifford", |bencher| {
        bencher.iter(|| {
            points.iter().map(|p| rotor.rotate(*p)).collect::<Vec<_>>()
        })
    });
}

fn bench_batch_rotation_nalgebra(c: &mut Criterion) {
    let q = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let points: Vec<na::Vector3<f64>> = (0..100)
        .map(|i| na::Vector3::new(i as f64, 0.0, 0.0))
        .collect();

    c.bench_function("workflow/batch_rotate_100/nalgebra", |bencher| {
        bencher.iter(|| {
            points.iter().map(|p| q * p).collect::<Vec<_>>()
        })
    });
}

/// Scenario: Mixed workflow - receive nalgebra data, process with GA, return nalgebra
fn bench_mixed_workflow(c: &mut Criterion) {
    let na_points: Vec<na::Vector3<f64>> = (0..10)
        .map(|i| na::Vector3::new(i as f64, 0.0, 0.0))
        .collect();
    let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());

    c.bench_function("workflow/mixed_convert_rotate_convert", |bencher| {
        bencher.iter(|| {
            na_points.iter().map(|p| {
                // Convert from nalgebra
                let cliff_p: dim3::Vector<f64> = (*p).into();
                // Process with GA
                let rotated = rotor.rotate(cliff_p);
                // Convert back to nalgebra
                let result: na::Vector3<f64> = rotated.into();
                result
            }).collect::<Vec<_>>()
        })
    });
}
```

### 6. GA-Specific Operations (No nalgebra Equivalent)

Benchmark operations unique to GA that nalgebra cannot express directly:

```rust
/// Wedge product (exterior algebra) - no direct nalgebra equivalent
fn bench_wedge_product(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 0.0, 0.0);
    let b = dim3::Vector::new(0.0, 1.0, 0.0);

    c.bench_function("ga_only/vec3_wedge", |bencher| {
        bencher.iter(|| black_box(a).wedge(black_box(b)))
    });
}

/// Rotor composition is more elegant than quaternion multiplication
fn bench_rotor_chain(c: &mut Criterion) {
    let r1 = dim3::Rotor::from_angle_plane(0.1, dim3::Bivector::unit_xy());
    let r2 = dim3::Rotor::from_angle_plane(0.2, dim3::Bivector::unit_xz());
    let r3 = dim3::Rotor::from_angle_plane(0.3, dim3::Bivector::unit_yz());

    c.bench_function("ga_only/rotor_chain_3", |bencher| {
        bencher.iter(|| {
            black_box(r1).compose(black_box(r2)).compose(black_box(r3))
        })
    });
}

/// Reflection through a plane (natural in GA)
fn bench_reflection(c: &mut Criterion) {
    let plane_normal = dim3::Vector::new(0.0, 1.0, 0.0); // y-plane
    let v = dim3::Vector::new(1.0, 2.0, 3.0);

    c.bench_function("ga_only/vec3_reflect", |bencher| {
        bencher.iter(|| {
            // v' = -n * v * n (sandwich with grade-1 element)
            let reflected = -black_box(plane_normal).sandwich(black_box(v));
            black_box(reflected)
        })
    });
}
```

## Deliverables

### 1. Benchmark File (`benches/nalgebra_comparison.rs`)

```rust
//! Benchmarks comparing clifford operations with nalgebra equivalents.
//!
//! Run with: `cargo bench --bench nalgebra_comparison --features nalgebra-0_34`

#![cfg(any(feature = "nalgebra-0_33", feature = "nalgebra-0_34"))]
#![allow(missing_docs)]

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use clifford::specialized::euclidean::{dim2, dim3};

// ... benchmark functions ...

criterion_group!(
    benches,
    // Vector comparisons
    bench_vec3_dot_clifford,
    bench_vec3_dot_nalgebra,
    bench_vec3_cross_clifford,
    bench_vec3_cross_nalgebra,
    // Rotation comparisons
    bench_rotation_clifford,
    bench_rotation_quaternion,
    bench_rotation_matrix,
    // Conversion overhead
    bench_vec3_to_nalgebra,
    bench_vec3_from_nalgebra,
    bench_rotor_to_quaternion,
    bench_quaternion_to_rotor,
    // Workflows
    bench_batch_rotation_clifford,
    bench_batch_rotation_nalgebra,
    bench_mixed_workflow,
    // GA-only
    bench_wedge_product,
    bench_rotor_chain,
);
criterion_main!(benches);
```

### 2. Cargo.toml Update

```toml
[[bench]]
name = "nalgebra_comparison"
harness = false
required-features = ["nalgebra-0_34"]
```

### 3. README Updates (`benches/README.md`)

Add new section:

```markdown
## nalgebra Comparison Benchmarks

Compare clifford's specialized types against nalgebra equivalents.

Run with nalgebra-0_34:
```bash
cargo bench --bench nalgebra_comparison --features nalgebra-0_34
```

### Vector Operations

| Operation | clifford | nalgebra | Difference |
|-----------|----------|----------|------------|
| Vec3 dot | ~X ns | ~Y ns | ... |
| Vec3 cross | ~X ns | ~Y ns | ... |
| Vec3 add | ~X ns | ~Y ns | ... |

### Rotation Operations

| Operation | clifford (Rotor) | nalgebra (Quat) | nalgebra (Rot) |
|-----------|------------------|-----------------|----------------|
| Rotate vec3 | ~X ns | ~Y ns | ~Z ns |
| Compose | ~X ns | ~Y ns | ~Z ns |
| Slerp | ~X ns | ~Y ns | N/A |

### Conversion Overhead

| Conversion | Time |
|------------|------|
| Vec3 -> Vector3 | ~X ns |
| Vector3 -> Vec3 | ~X ns |
| Rotor -> Quaternion | ~X ns |
| Quaternion -> Rotor | ~X ns |

### Key Findings

- **Vector operations**: [summary of findings]
- **Rotations**: [summary of findings]
- **Conversions**: [summary of findings]
- **Recommendation**: [when to use which]
```

## Expected Results

Based on implementation characteristics, we expect:

1. **Vector operations**: Near-identical performance (both use simple arithmetic)
2. **Rotation (rotor vs quaternion)**: Similar performance (same underlying math)
3. **Rotation (rotor vs rotation matrix)**: Matrix may be faster for single rotations, rotor faster for composition
4. **Conversions**: Very low overhead (~1-2ns for vectors, ~5-10ns for rotors)
5. **GA-specific ops**: Unique to clifford, no comparison needed

## CI Integration

Add benchmark job to CI (optional, for regression tracking):

```yaml
bench-nalgebra:
  name: Benchmark nalgebra comparison
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@stable
    - uses: Swatinem/rust-cache@v2
    - run: cargo bench --bench nalgebra_comparison --features nalgebra-0_34 -- --noplot
```

## Verification

- [ ] `cargo bench --bench nalgebra_comparison --features nalgebra-0_33` runs
- [ ] `cargo bench --bench nalgebra_comparison --features nalgebra-0_34` runs
- [ ] Results documented in `benches/README.md`
- [ ] SVG reports captured in `benches/reports/`
- [ ] Key findings summarized

## Future Considerations

- **SIMD comparison**: When SIMD is implemented (PRD-9), compare against nalgebra's SIMD paths
- **PGA/CGA benchmarks**: When implemented, compare Motor vs Isometry, conformal point ops
- **Batch operations**: Compare vectorized batch processing approaches
