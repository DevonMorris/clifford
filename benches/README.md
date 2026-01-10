# Benchmarks

Criterion benchmarks comparing generic `Multivector` operations vs specialized 2D/3D types.

## Performance Summary

The specialized types provide significant speedups over the generic implementation:

| Operation | Generic | Specialized 3D | Specialized 2D | Speedup |
|-----------|---------|----------------|----------------|---------|
| Vector dot | ~59 ns | ~8.5 ns | ~0.85 ns | **7-70x** |
| Vector wedge | ~61 ns | ~8.5 ns | ~0.88 ns | **7-70x** |
| Vector add | ~1.7 ns | ~1.7 ns | ~0.83 ns | ~1-2x |
| Rotor rotation | ~50 ns | ~5.5 ns | ~1.5 ns | **9-33x** |
| Rotor compose | - | ~10.7 ns | ~0.98 ns | - |

## Generic Multivector Benchmarks

Operations on the sparse `Multivector<T, S>` type.

### Vector Addition (~1.7 ns)
![generic/vector_add](reports/generic/vector_add.svg)

### Vector Dot Product (~59 ns)
![generic/vector_dot](reports/generic/vector_dot.svg)

### Vector Wedge Product (~61 ns)
![generic/vector_wedge](reports/generic/vector_wedge.svg)

### Geometric Product (~41 ns)
![generic/vector_geometric](reports/generic/vector_geometric.svg)

### Rotor Sandwich (~50 ns)
![generic/rotor_sandwich](reports/generic/rotor_sandwich.svg)

### Full Multivector Geometric (~282 ns)
![generic/full_geometric](reports/generic/full_geometric.svg)

## Specialized 2D Benchmarks

Operations on fixed-size `dim2::Vector`, `dim2::Rotor` types.

### Vector Dot (~0.85 ns)
![euclidean/dim2/vector_dot](reports/euclidean/dim2/vector_dot.svg)

### Vector Wedge (~0.88 ns)
![euclidean/dim2/vector_wedge](reports/euclidean/dim2/vector_wedge.svg)

### Vector Add (~0.83 ns)
![euclidean/dim2/vector_add](reports/euclidean/dim2/vector_add.svg)

### Rotor Rotate (~1.5 ns)
![euclidean/dim2/rotor_rotate](reports/euclidean/dim2/rotor_rotate.svg)

### Rotor Compose (~0.98 ns)
![euclidean/dim2/rotor_compose](reports/euclidean/dim2/rotor_compose.svg)

### Rotor Slerp (~35 ns)
![euclidean/dim2/rotor_slerp](reports/euclidean/dim2/rotor_slerp.svg)

## Specialized 3D Benchmarks

Operations on fixed-size `dim3::Vector`, `dim3::Bivector`, `dim3::Rotor` types.

### Vector Dot (~8.5 ns)
![euclidean/dim3/vector_dot](reports/euclidean/dim3/vector_dot.svg)

### Vector Wedge (~8.5 ns)
![euclidean/dim3/vector_wedge](reports/euclidean/dim3/vector_wedge.svg)

### Vector Cross (~8.8 ns)
![euclidean/dim3/vector_cross](reports/euclidean/dim3/vector_cross.svg)

### Vector Add (~1.7 ns)
![euclidean/dim3/vector_add](reports/euclidean/dim3/vector_add.svg)

### Rotor Rotate (~5.5 ns)
![euclidean/dim3/rotor_rotate](reports/euclidean/dim3/rotor_rotate.svg)

### Rotor Compose (~10.7 ns)
![euclidean/dim3/rotor_compose](reports/euclidean/dim3/rotor_compose.svg)

### Rotor Slerp (~36 ns)
![euclidean/dim3/rotor_slerp](reports/euclidean/dim3/rotor_slerp.svg)

### Rotor from Vectors (~23 ns)
![euclidean/dim3/rotor_from_vectors](reports/euclidean/dim3/rotor_from_vectors.svg)

## nalgebra Comparison Benchmarks

Head-to-head comparison of clifford specialized types vs nalgebra equivalents.

### Performance Summary

| Operation | clifford | nalgebra | Winner |
|-----------|----------|----------|--------|
| **3D Rotation** | 4.4 ns | 10.5 ns (quat) / 7.4 ns (rot) | **clifford 2.4x faster** |
| **Batch rotate 100** | 464 ns | 1094 ns | **clifford 2.4x faster** |
| **Vec3 dot** | 6.9 ns | 7.5 ns | ~equal |
| **Vec3 cross** | 7.1 ns | 7.0 ns | ~equal |
| **Vec3 normalize** | 3.5 ns | 3.4 ns | ~equal |
| **Rotation compose** | 8.6 ns | 8.1 ns | ~equal |
| **Rotation slerp** | 26.8 ns | 18.6 ns | nalgebra 1.4x faster |
| **2D rotation** | 1.4 ns | 1.1 ns | nalgebra slightly faster |
| **Vec2 dot** | 0.6 ns | 0.7 ns | ~equal |

### Key Findings

1. **Rotation is clifford's strength**: The rotor-based rotation is 2.4x faster than quaternion rotation
2. **Vector operations are equivalent**: Both libraries optimize these to similar performance
3. **Slerp favors nalgebra**: nalgebra's quaternion slerp is more optimized
4. **Conversion overhead is minimal**: ~1 ns for vectors, ~6 ns for rotor↔quaternion

### 3D Vector Operations

#### Vec3 Dot - clifford (~6.9 ns)
![comparison_vec3_dot_clifford](reports/comparison/comparison_vec3_dot_clifford.svg)

#### Vec3 Dot - nalgebra (~7.5 ns)
![comparison_vec3_dot_nalgebra](reports/comparison/comparison_vec3_dot_nalgebra.svg)

#### Vec3 Cross - clifford (~7.1 ns)
![comparison_vec3_cross_clifford](reports/comparison/comparison_vec3_cross_clifford.svg)

#### Vec3 Cross - nalgebra (~7.0 ns)
![comparison_vec3_cross_nalgebra](reports/comparison/comparison_vec3_cross_nalgebra.svg)

#### Vec3 Normalize - clifford (~3.5 ns)
![comparison_vec3_normalize_clifford](reports/comparison/comparison_vec3_normalize_clifford.svg)

#### Vec3 Normalize - nalgebra (~3.4 ns)
![comparison_vec3_normalize_nalgebra](reports/comparison/comparison_vec3_normalize_nalgebra.svg)

### 3D Rotation Operations

#### Rotate Vec3 - clifford Rotor (~4.4 ns)
![comparison_rotate_vec3_clifford](reports/comparison/comparison_rotate_vec3_clifford.svg)

#### Rotate Vec3 - nalgebra Quaternion (~10.5 ns)
![comparison_rotate_vec3_nalgebra_quat](reports/comparison/comparison_rotate_vec3_nalgebra_quat.svg)

#### Rotate Vec3 - nalgebra Rotation3 (~7.4 ns)
![comparison_rotate_vec3_nalgebra_rot](reports/comparison/comparison_rotate_vec3_nalgebra_rot.svg)

#### Rotation Compose - clifford (~8.6 ns)
![comparison_rotation_compose_clifford](reports/comparison/comparison_rotation_compose_clifford.svg)

#### Rotation Compose - nalgebra (~8.1 ns)
![comparison_rotation_compose_nalgebra](reports/comparison/comparison_rotation_compose_nalgebra.svg)

#### Rotation Slerp - clifford (~26.8 ns)
![comparison_rotation_slerp_clifford](reports/comparison/comparison_rotation_slerp_clifford.svg)

#### Rotation Slerp - nalgebra (~18.6 ns)
![comparison_rotation_slerp_nalgebra](reports/comparison/comparison_rotation_slerp_nalgebra.svg)

### 2D Operations

#### Vec2 Dot - clifford (~0.6 ns)
![comparison_vec2_dot_clifford](reports/comparison/comparison_vec2_dot_clifford.svg)

#### Vec2 Dot - nalgebra (~0.7 ns)
![comparison_vec2_dot_nalgebra](reports/comparison/comparison_vec2_dot_nalgebra.svg)

#### Rotate Vec2 - clifford (~1.4 ns)
![comparison_rotate_vec2_clifford](reports/comparison/comparison_rotate_vec2_clifford.svg)

#### Rotate Vec2 - nalgebra (~1.1 ns)
![comparison_rotate_vec2_nalgebra](reports/comparison/comparison_rotate_vec2_nalgebra.svg)

### Conversion Overhead

| Conversion | Time |
|------------|------|
| Vec3 → nalgebra | ~1.1 ns |
| nalgebra → Vec3 | ~1.1 ns |
| Rotor → Quaternion | ~6.3 ns |
| Quaternion → Rotor | ~1.5 ns |

#### Vec3 to nalgebra (~1.1 ns)
![conversion_vec3_to_nalgebra](reports/comparison/conversion_vec3_to_nalgebra.svg)

#### Vec3 from nalgebra (~1.1 ns)
![conversion_vec3_from_nalgebra](reports/comparison/conversion_vec3_from_nalgebra.svg)

#### Rotor to Quaternion (~6.3 ns)
![conversion_rotor_to_quaternion](reports/comparison/conversion_rotor_to_quaternion.svg)

#### Quaternion to Rotor (~1.5 ns)
![conversion_quaternion_to_rotor](reports/comparison/conversion_quaternion_to_rotor.svg)

### Workflow Benchmarks

#### Batch Rotate 100 Points - clifford (~464 ns)
![workflow_batch_rotate_100_clifford](reports/comparison/workflow_batch_rotate_100_clifford.svg)

#### Batch Rotate 100 Points - nalgebra (~1094 ns)
![workflow_batch_rotate_100_nalgebra](reports/comparison/workflow_batch_rotate_100_nalgebra.svg)

#### Mixed Workflow (convert → rotate → convert) (~53 ns for 10 points)
![workflow_mixed_convert_rotate_convert](reports/comparison/workflow_mixed_convert_rotate_convert.svg)

### GA-Specific Operations (no nalgebra equivalent)

These operations demonstrate geometric algebra capabilities unique to clifford.

#### Vec3 Wedge Product (~6.6 ns)
![ga_only_vec3_wedge](reports/comparison/ga_only_vec3_wedge.svg)

#### Vec3 Geometric Product (~7.6 ns)
![ga_only_vec3_geometric](reports/comparison/ga_only_vec3_geometric.svg)

#### Rotor Chain (3 rotations) (~17.6 ns)
![ga_only_rotor_chain_3](reports/comparison/ga_only_rotor_chain_3.svg)

## Running Benchmarks

```bash
# Run all benchmarks (excluding nalgebra comparison)
cargo bench

# Run only generic benchmarks
cargo bench --bench generic

# Run only specialized benchmarks
cargo bench --bench specialized

# Run nalgebra comparison benchmarks (requires feature flag)
cargo bench --bench nalgebra_comparison --features nalgebra-0_34

# Run specific benchmark by name
cargo bench -- "rotor"
```

After running, detailed HTML reports are generated at `target/criterion/report/index.html`.

## Regenerating Report Images

To update the SVG plots in this directory after running benchmarks:

```bash
# Copy generic benchmarks
for svg in target/criterion/generic_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#generic_}
  cp "$svg" "benches/reports/generic/${name}.svg"
done

# Copy euclidean/dim2 benchmarks
for svg in target/criterion/euclidean_dim2_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#euclidean_dim2_}
  cp "$svg" "benches/reports/euclidean/dim2/${name}.svg"
done

# Copy euclidean/dim3 benchmarks
for svg in target/criterion/euclidean_dim3_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#euclidean_dim3_}
  cp "$svg" "benches/reports/euclidean/dim3/${name}.svg"
done

# Copy nalgebra comparison benchmarks
for dir in target/criterion/comparison_* target/criterion/conversion_* \
           target/criterion/workflow_* target/criterion/ga_only_*; do
  if [ -d "$dir" ]; then
    bench=$(basename "$dir")
    cp "$dir/report/pdf.svg" "benches/reports/comparison/${bench}.svg"
  fi
done
```

## Adding New Benchmarks

1. Add benchmark function to appropriate file (`generic.rs` or `specialized.rs`):

```rust
fn bench_new_operation(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 2.0, 3.0);

    c.bench_function("euclidean/dim3/new_operation", |bencher| {
        bencher.iter(|| black_box(a).some_operation())
    });
}
```

2. Add the function to `criterion_group!` at the bottom of the file.

3. Run benchmarks: `cargo bench`

4. Copy the new SVG and update this README.
