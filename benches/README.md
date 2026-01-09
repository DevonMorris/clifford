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

## Running Benchmarks

```bash
# Run all benchmarks
cargo bench

# Run only generic benchmarks
cargo bench --bench generic

# Run only specialized benchmarks
cargo bench --bench specialized

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
