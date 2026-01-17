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

## Projective Geometric Algebra (PGA) Benchmarks

Operations on specialized PGA types for rigid body transformations.

### 2D PGA Performance Summary

| Operation | Time | Description |
|-----------|------|-------------|
| Point join | ~8 ns | Join two points to get a line |
| Line meet | ~2 ns | Meet two lines to get intersection |
| Motor transform point | ~3 ns | Apply rigid transform to point |
| Motor compose | ~10 ns | Compose two motors |
| Motor inverse | ~3 ns | Compute motor inverse |

### 3D PGA Performance Summary

| Operation | Time | Description |
|-----------|------|-------------|
| Motor transform point | ~16 ns | Apply rigid transform to point |
| Motor transform line | ~24 ns | Apply rigid transform to line |
| Motor compose | ~15 ns | Compose two motors |
| Motor inverse | ~11 ns | Compute motor inverse |
| Line meet plane | ~10 ns | Intersection point |
| Line join point | ~10 ns | Plane through line and point |
| Motor commutator | ~20 ns | [M₁, M₂] = M₁M₂ - M₂M₁ |

### 2D PGA Benchmarks

#### Point Join (~8 ns)
![projective/dim2/point_join](reports/projective/dim2/point_join.svg)

#### Line Meet (~2 ns)
![projective/dim2/line_meet](reports/projective/dim2/line_meet.svg)

#### Motor Transform Point (~3 ns)
![projective/dim2/motor_transform_point](reports/projective/dim2/motor_transform_point.svg)

#### Motor Transform Line (~3 ns)
![projective/dim2/motor_transform_line](reports/projective/dim2/motor_transform_line.svg)

#### Motor Compose (~10 ns)
![projective/dim2/motor_compose](reports/projective/dim2/motor_compose.svg)

#### Motor Inverse (~3 ns)
![projective/dim2/motor_inverse](reports/projective/dim2/motor_inverse.svg)

#### Point Distance (~3 ns)
![projective/dim2/point_distance](reports/projective/dim2/point_distance.svg)

#### Line Distance to Point (~7 ns)
![projective/dim2/line_distance_to_point](reports/projective/dim2/line_distance_to_point.svg)

#### Batch Transform 100 Points (~330 ns)
![projective/dim2/batch_transform_100_points](reports/projective/dim2/batch_transform_100_points.svg)

### 3D PGA Benchmarks

#### Motor Transform Point (~16 ns)
![projective/dim3/motor_transform_point](reports/projective/dim3/motor_transform_point.svg)

#### Motor Transform Line (~24 ns)
![projective/dim3/motor_transform_line](reports/projective/dim3/motor_transform_line.svg)

#### Motor Compose (~15 ns)
![projective/dim3/motor_compose](reports/projective/dim3/motor_compose.svg)

#### Motor Inverse (~11 ns)
![projective/dim3/motor_inverse](reports/projective/dim3/motor_inverse.svg)

#### Motor Commutator (~20 ns)
![projective/dim3/motor_commutator](reports/projective/dim3/motor_commutator.svg)

#### Line Meet Plane (~10 ns)
![projective/dim3/line_meet_plane](reports/projective/dim3/line_meet_plane.svg)

#### Line Join Point (~10 ns)
![projective/dim3/line_join_point](reports/projective/dim3/line_join_point.svg)

#### Line Distance to Point (~14 ns)
![projective/dim3/line_distance_to_point](reports/projective/dim3/line_distance_to_point.svg)

#### Line Closest Point (~15 ns)
![projective/dim3/line_closest_point](reports/projective/dim3/line_closest_point.svg)

#### Point Distance (~4 ns)
![projective/dim3/point_distance](reports/projective/dim3/point_distance.svg)

#### Point Dot Product (~2 ns)
![projective/dim3/point_dot](reports/projective/dim3/point_dot.svg)

#### Line Dot Product (~12 ns)
![projective/dim3/line_dot](reports/projective/dim3/line_dot.svg)

#### Line Distance (~21 ns)
![projective/dim3/line_distance](reports/projective/dim3/line_distance.svg)

#### Point Left Contract Line (~9 ns)
![projective/dim3/point_left_contract_line](reports/projective/dim3/point_left_contract_line.svg)

#### Point Left Contract Plane (~9 ns)
![projective/dim3/point_left_contract_plane](reports/projective/dim3/point_left_contract_plane.svg)

#### Line Left Contract Plane (~10 ns)
![projective/dim3/line_left_contract_plane](reports/projective/dim3/line_left_contract_plane.svg)

#### Batch Transform 100 Points (~1.6 µs)
![projective/dim3/batch_transform_100_points](reports/projective/dim3/batch_transform_100_points.svg)

## PGA vs nalgebra Comparison

Head-to-head comparison of PGA Motor operations vs nalgebra Isometry.

### Performance Summary

| Operation | clifford PGA | nalgebra | Notes |
|-----------|-------------|----------|-------|
| **2D transform point** | 3.3 ns | 9.1 ns | **PGA 3x faster** |
| **2D compose** | 10.4 ns | 8.9 ns | nalgebra ~1.2x faster |
| **2D inverse** | 2.8 ns | 1.4 ns | nalgebra 2x faster |
| **3D transform point** | 16 ns | 15 ns | Comparable |
| **3D compose** | 15.4 ns | 15.6 ns | Comparable |
| **3D inverse** | 11.5 ns | 13.2 ns | **PGA 15% faster** |
| **3D batch 100** | 1.6 µs | 1.5 µs | Comparable |

### Key Findings

1. **2D transform is PGA's strength**: 3x faster than nalgebra Isometry2
2. **3D operations are comparable**: Similar performance to quaternion-based approach
3. **PGA inverse is faster in 3D**: Motor inverse is simpler than Isometry inverse
4. **PGA provides unique operations**: meet, join, contraction have no nalgebra equivalent

### 2D Comparisons

#### 2D Transform Point - clifford (~3.3 ns)
![pga_comparison_2d_transform_point_clifford](reports/pga_comparison/pga_comparison_2d_transform_point_clifford.svg)

#### 2D Transform Point - nalgebra (~9.1 ns)
![pga_comparison_2d_transform_point_nalgebra](reports/pga_comparison/pga_comparison_2d_transform_point_nalgebra.svg)

#### 2D Compose - clifford (~10.4 ns)
![pga_comparison_2d_compose_clifford](reports/pga_comparison/pga_comparison_2d_compose_clifford.svg)

#### 2D Compose - nalgebra (~8.9 ns)
![pga_comparison_2d_compose_nalgebra](reports/pga_comparison/pga_comparison_2d_compose_nalgebra.svg)

#### 2D Inverse - clifford (~2.8 ns)
![pga_comparison_2d_inverse_clifford](reports/pga_comparison/pga_comparison_2d_inverse_clifford.svg)

#### 2D Inverse - nalgebra (~1.4 ns)
![pga_comparison_2d_inverse_nalgebra](reports/pga_comparison/pga_comparison_2d_inverse_nalgebra.svg)

#### 2D Batch 100 - clifford (~380 ns)
![pga_comparison_2d_batch_100_clifford](reports/pga_comparison/pga_comparison_2d_batch_100_clifford.svg)

#### 2D Batch 100 - nalgebra (~920 ns)
![pga_comparison_2d_batch_100_nalgebra](reports/pga_comparison/pga_comparison_2d_batch_100_nalgebra.svg)

### 3D Comparisons

#### 3D Transform Point - clifford (~16 ns)
![pga_comparison_3d_transform_point_clifford](reports/pga_comparison/pga_comparison_3d_transform_point_clifford.svg)

#### 3D Transform Point - nalgebra (~15 ns)
![pga_comparison_3d_transform_point_nalgebra](reports/pga_comparison/pga_comparison_3d_transform_point_nalgebra.svg)

#### 3D Compose - clifford (~15.4 ns)
![pga_comparison_3d_compose_clifford](reports/pga_comparison/pga_comparison_3d_compose_clifford.svg)

#### 3D Compose - nalgebra (~15.6 ns)
![pga_comparison_3d_compose_nalgebra](reports/pga_comparison/pga_comparison_3d_compose_nalgebra.svg)

#### 3D Inverse - clifford (~11.5 ns)
![pga_comparison_3d_inverse_clifford](reports/pga_comparison/pga_comparison_3d_inverse_clifford.svg)

#### 3D Inverse - nalgebra (~13.2 ns)
![pga_comparison_3d_inverse_nalgebra](reports/pga_comparison/pga_comparison_3d_inverse_nalgebra.svg)

#### 3D Batch 100 - clifford (~1.6 µs)
![pga_comparison_3d_batch_100_clifford](reports/pga_comparison/pga_comparison_3d_batch_100_clifford.svg)

#### 3D Batch 100 - nalgebra (~1.5 µs)
![pga_comparison_3d_batch_100_nalgebra](reports/pga_comparison/pga_comparison_3d_batch_100_nalgebra.svg)

#### 3D Point Distance - clifford (~4 ns)
![pga_comparison_3d_point_distance_clifford](reports/pga_comparison/pga_comparison_3d_point_distance_clifford.svg)

#### 3D Point Distance - nalgebra (~3 ns)
![pga_comparison_3d_point_distance_nalgebra](reports/pga_comparison/pga_comparison_3d_point_distance_nalgebra.svg)

### PGA Conversions

| Conversion | Time |
|------------|------|
| Motor → Isometry3 | ~6 ns |
| Isometry3 → Motor | ~19 ns |
| Point → nalgebra | ~2.6 ns |
| nalgebra → Point | ~2.2 ns |

#### Motor to Isometry3 (~6 ns)
![pga_conversion_motor_to_isometry3](reports/pga_comparison/pga_conversion_motor_to_isometry3.svg)

#### Isometry3 to Motor (~19 ns)
![pga_conversion_isometry3_to_motor](reports/pga_comparison/pga_conversion_isometry3_to_motor.svg)

#### Point to nalgebra (~2.6 ns)
![pga_conversion_point3_to_nalgebra](reports/pga_comparison/pga_conversion_point3_to_nalgebra.svg)

#### Point from nalgebra (~2.2 ns)
![pga_conversion_point3_from_nalgebra](reports/pga_comparison/pga_conversion_point3_from_nalgebra.svg)

### PGA-Only Operations (no nalgebra equivalent)

These operations demonstrate PGA capabilities unique to geometric algebra.

#### Line Meet Plane (~10 ns)
![pga_only_3d_line_meet_plane](reports/pga_comparison/pga_only_3d_line_meet_plane.svg)

#### Line Join Point (~10 ns)
![pga_only_3d_line_join_point](reports/pga_comparison/pga_only_3d_line_join_point.svg)

#### Motor Commutator (~20 ns)
![pga_only_3d_motor_commutator](reports/pga_comparison/pga_only_3d_motor_commutator.svg)

#### Point Left Contract Line (~9 ns)
![pga_only_3d_point_left_contract_line](reports/pga_comparison/pga_only_3d_point_left_contract_line.svg)

## Euclidean GA vs nalgebra Comparison

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
# Run all benchmarks (excluding nalgebra comparisons)
cargo bench

# Run only generic benchmarks
cargo bench --bench generic

# Run only specialized (Euclidean) benchmarks
cargo bench --bench specialized

# Run only PGA benchmarks
cargo bench --bench projective

# Run Euclidean vs nalgebra comparison benchmarks
cargo bench --bench nalgebra_comparison --no-default-features --features "proptest-support,nalgebra-0_34"

# Run PGA vs nalgebra comparison benchmarks
cargo bench --bench pga_nalgebra_comparison --no-default-features --features "proptest-support,nalgebra-0_34"

# Run specific benchmark by name
cargo bench -- "motor"
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

# Copy Euclidean nalgebra comparison benchmarks
for dir in target/criterion/comparison_* target/criterion/conversion_* \
           target/criterion/workflow_* target/criterion/ga_only_*; do
  if [ -d "$dir" ]; then
    bench=$(basename "$dir")
    cp "$dir/report/pdf.svg" "benches/reports/comparison/${bench}.svg"
  fi
done

# Copy projective/dim2 benchmarks
for svg in target/criterion/projective_dim2_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#projective_dim2_}
  cp "$svg" "benches/reports/projective/dim2/${name}.svg"
done

# Copy projective/dim3 benchmarks
for svg in target/criterion/projective_dim3_*/report/pdf.svg; do
  bench=$(basename $(dirname $(dirname "$svg")))
  name=${bench#projective_dim3_}
  cp "$svg" "benches/reports/projective/dim3/${name}.svg"
done

# Copy PGA nalgebra comparison benchmarks
for dir in target/criterion/pga_comparison_* target/criterion/pga_conversion_* \
           target/criterion/pga_only_*; do
  if [ -d "$dir" ]; then
    bench=$(basename "$dir")
    cp "$dir/report/pdf.svg" "benches/reports/pga_comparison/${bench}.svg"
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
