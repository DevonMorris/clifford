# Benchmarks

This directory contains criterion benchmarks for Clifford's `Multivector` operations.

## Running Benchmarks

```bash
cargo bench
```

## Viewing HTML Reports

After running benchmarks, HTML reports are generated in `target/criterion/`.

**Main report index:**
```
target/criterion/report/index.html
```

**Individual benchmark reports:**
```
target/criterion/<benchmark_name>/report/index.html
```

Open these files in a browser to view detailed performance analysis, including:
- Timing distributions
- Iteration time plots
- Comparison with previous runs
- Statistical analysis

## Current Benchmarks

| Benchmark | Description |
|-----------|-------------|
| `vector_dot_generic` | Inner product of two 3D vectors |
| `vector_wedge_generic` | Outer (wedge) product of two 3D vectors |
| `vector_geometric_generic` | Geometric product of two 3D vectors |
| `rotor_sandwich_generic` | Rotor rotation via sandwich product |
| `multivector_add` | Addition of two 3D multivectors |
| `full_multivector_geometric` | Geometric product with all 8 components non-zero |

## Adding New Benchmarks

When adding new operations to the library, add corresponding benchmarks to `multivector.rs`:

```rust
fn bench_new_operation(c: &mut Criterion) {
    // Setup test data outside the benchmark loop
    let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);

    c.bench_function("new_operation", |bencher| {
        bencher.iter(|| {
            // Only the operation being measured goes here
            black_box(&a).some_operation()
        })
    });
}
```

Then add the function to the `criterion_group!` macro at the bottom of the file.
