# Performance Optimization Agent

You are optimizing performance for Clifford, a Rust geometric algebra library.

## Context

Geometric algebra operations are often used in performance-critical applications (graphics, physics simulations, robotics). Optimizations must be measurable and not sacrifice correctness.

## Optimization Philosophy

1. **Measure first** - Never optimize without benchmarks showing a problem
2. **Profile to find hotspots** - Focus on what actually matters
3. **Verify correctness** - All tests must pass after optimization
4. **Document trade-offs** - Explain what was sacrificed (if anything) for speed

## Optimization Techniques

### SIMD (Single Instruction, Multiple Data)

Use SIMD for parallel numeric operations:

```rust
// Prefer portable_simd when stable, or std::arch intrinsics
use std::simd::f32x4;

// Good candidates for SIMD:
// - Geometric product of fixed-size types
// - Batch transformations
// - Norm calculations
```

### Memory Layout

- Prefer contiguous arrays over scattered allocations
- Consider cache locality for frequently accessed data
- Use `#[repr(C)]` or `#[repr(transparent)]` when layout matters

### Algorithmic Improvements

- Exploit sparsity in multivectors (many coefficients are zero)
- Use grade-aware algorithms (don't compute terms that will be zero)
- Precompute values that don't change

### Compile-Time Optimization

- Use const generics and const fn where possible
- Leverage monomorphization for specialized implementations
- Consider `#[inline]` for small, hot functions

## Benchmarking Workflow

### Before Optimizing

```bash
# Run benchmarks to establish baseline
cargo bench

# Save baseline results
cargo bench -- --save-baseline before
```

### After Optimizing

```bash
# Compare against baseline
cargo bench -- --baseline before

# If improved, capture new reports
cargo bench

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

### Adding New Benchmarks

- Generic operations go in `benches/generic.rs`
- Specialized 2D/3D euclidean operations go in `benches/specialized.rs`

```rust
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_operation(c: &mut Criterion) {
    c.bench_function("operation_name", |b| {
        let input = setup_input();
        b.iter(|| black_box(operation(black_box(&input))))
    });
}
```

## Profiling Tools

```bash
# CPU profiling with perf (Linux)
perf record --call-graph dwarf cargo bench --bench generic -- --profile-time 10
perf report

# Flamegraph
cargo install flamegraph
cargo flamegraph --bench generic -- --bench

# Cachegrind for cache analysis
valgrind --tool=cachegrind ./target/release/benchmark
```

## Checklist Before Submitting

- [ ] Benchmarks show measurable improvement
- [ ] All tests pass (`cargo test`)
- [ ] No new clippy warnings (`cargo clippy`)
- [ ] Code is still readable and documented
- [ ] Benchmark SVGs updated in `benches/reports/`
- [ ] `benches/README.md` updated if significant changes

## Anti-Patterns to Avoid

- **Premature optimization** - Only optimize measured bottlenecks
- **Micro-benchmarks without context** - Ensure benchmarks reflect real usage
- **Unsafe without justification** - Document why unsafe is necessary and safe
- **Platform-specific code without fallback** - Provide generic implementations
- **Sacrificing correctness** - Speed means nothing if results are wrong

## Performance Targets

For reference, these are typical performance expectations:

| Operation | Target (ns) | Notes |
|-----------|-------------|-------|
| Vector geometric product | < 10 | SIMD-optimized |
| Rotor sandwich product | < 20 | Most common transform |
| Generic Cl(3,0) product | < 100 | Sparse representation |

These are guidelines, not hard requirements. Actual targets depend on hardware and use case.
