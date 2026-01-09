# Benchmarks

Criterion benchmarks for Clifford's `Multivector` operations on `Euclidean3` (3D space).

## Results Summary

| Benchmark | Time | Description |
|-----------|------|-------------|
| `multivector_add` | ~1.6 ns | Addition of two 3D multivectors |
| `vector_geometric_generic` | ~40 ns | Geometric product of two 3D vectors |
| `rotor_sandwich_generic` | ~46 ns | Rotor rotation via sandwich product |
| `vector_dot_generic` | ~58 ns | Inner product of two 3D vectors |
| `vector_wedge_generic` | ~61 ns | Outer (wedge) product of two 3D vectors |
| `full_multivector_geometric` | ~278 ns | Geometric product with all 8 components |

## Timing Distributions

### Multivector Addition (~1.6 ns)
![multivector_add](reports/multivector_add_pdf.svg)

### Geometric Product - Vectors (~40 ns)
![vector_geometric_generic](reports/vector_geometric_generic_pdf.svg)

### Rotor Sandwich Product (~46 ns)
![rotor_sandwich_generic](reports/rotor_sandwich_generic_pdf.svg)

### Inner Product (~58 ns)
![vector_dot_generic](reports/vector_dot_generic_pdf.svg)

### Outer/Wedge Product (~61 ns)
![vector_wedge_generic](reports/vector_wedge_generic_pdf.svg)

### Full Multivector Geometric Product (~278 ns)
![full_multivector_geometric](reports/full_multivector_geometric_pdf.svg)

## Running Benchmarks

```bash
cargo bench
```

After running, detailed HTML reports are generated at `target/criterion/report/index.html`.

## Regenerating Report Images

To update the SVG plots in this directory after running benchmarks:

```bash
for bench in vector_dot_generic vector_wedge_generic vector_geometric_generic \
             rotor_sandwich_generic multivector_add full_multivector_geometric; do
  cp "target/criterion/$bench/report/pdf.svg" "benches/reports/${bench}_pdf.svg"
done
```

## Adding New Benchmarks

When adding new operations, add benchmarks to `multivector.rs`:

```rust
fn bench_new_operation(c: &mut Criterion) {
    let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);

    c.bench_function("new_operation", |bencher| {
        bencher.iter(|| black_box(&a).some_operation())
    });
}
```

Add the function to `criterion_group!` at the bottom of the file, run benchmarks, and copy the new SVG.
