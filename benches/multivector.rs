//! Benchmarks for Multivector operations.
//!
//! Run with: `cargo bench`

#![allow(missing_docs)]

use std::hint::black_box;

use clifford::algebra::Multivector;
use clifford::signature::Euclidean3;
use criterion::{Criterion, criterion_group, criterion_main};

fn bench_vector_dot(c: &mut Criterion) {
    let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    let b: Multivector<f64, Euclidean3> = Multivector::vector(&[4.0, 5.0, 6.0]);

    c.bench_function("vector_dot_generic", |bencher| {
        bencher.iter(|| black_box(&a).inner(black_box(&b)).scalar_part())
    });
}

fn bench_vector_wedge(c: &mut Criterion) {
    let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    let b: Multivector<f64, Euclidean3> = Multivector::vector(&[4.0, 5.0, 6.0]);

    c.bench_function("vector_wedge_generic", |bencher| {
        bencher.iter(|| black_box(&a).outer(black_box(&b)))
    });
}

fn bench_geometric_product(c: &mut Criterion) {
    let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    let b: Multivector<f64, Euclidean3> = Multivector::vector(&[4.0, 5.0, 6.0]);

    c.bench_function("vector_geometric_generic", |bencher| {
        bencher.iter(|| black_box(&a) * black_box(&b))
    });
}

fn bench_rotor_rotation(c: &mut Criterion) {
    // Create a rotor: cos(θ/2) + sin(θ/2) * e₁₂
    let angle = std::f64::consts::PI / 4.0; // 45 degrees
    let half = angle / 2.0;
    let mut rotor: Multivector<f64, Euclidean3> = Multivector::scalar(half.cos());
    let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    let e12 = &e1 * &e2;
    rotor = &rotor + &(&e12 * half.sin());

    let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 0.0, 0.0]);

    c.bench_function("rotor_sandwich_generic", |bencher| {
        bencher.iter(|| black_box(&rotor).sandwich(black_box(&v)))
    });
}

fn bench_multivector_add(c: &mut Criterion) {
    let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    let b: Multivector<f64, Euclidean3> = Multivector::vector(&[4.0, 5.0, 6.0]);

    c.bench_function("multivector_add", |bencher| {
        bencher.iter(|| black_box(&a) + black_box(&b))
    });
}

fn bench_full_multivector_product(c: &mut Criterion) {
    // Full multivector with all 8 components non-zero
    let mut a: Multivector<f64, Euclidean3> = Multivector::zero();
    let mut b: Multivector<f64, Euclidean3> = Multivector::zero();
    for i in 0..8 {
        a.set(clifford::basis::Blade::from_index(i), (i + 1) as f64);
        b.set(clifford::basis::Blade::from_index(i), (8 - i) as f64);
    }

    c.bench_function("full_multivector_geometric", |bencher| {
        bencher.iter(|| black_box(&a) * black_box(&b))
    });
}

criterion_group!(
    benches,
    bench_vector_dot,
    bench_vector_wedge,
    bench_geometric_product,
    bench_rotor_rotation,
    bench_multivector_add,
    bench_full_multivector_product,
);
criterion_main!(benches);
