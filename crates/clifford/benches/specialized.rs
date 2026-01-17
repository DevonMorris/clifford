//! Benchmarks for specialized 2D/3D geometric algebra types.
//!
//! Run with: `cargo bench --bench specialized`

#![allow(missing_docs, clippy::missing_docs_in_private_items)]

use std::f64::consts::FRAC_PI_4;
use std::hint::black_box;

use clifford::ops::{Transform, Wedge};
use clifford::specialized::euclidean::{dim2, dim3};
use criterion::{Criterion, criterion_group, criterion_main};

// ============================================================================
// 2D Benchmarks
// ============================================================================

fn bench_vec2_dot(c: &mut Criterion) {
    let a: dim2::Vector<f64> = dim2::Vector::new(1.0, 2.0);
    let b: dim2::Vector<f64> = dim2::Vector::new(3.0, 4.0);

    c.bench_function("euclidean/dim2/vector_dot", |bencher| {
        bencher.iter(|| black_box(a).dot(black_box(b)))
    });
}

fn bench_vec2_wedge(c: &mut Criterion) {
    let a: dim2::Vector<f64> = dim2::Vector::new(1.0, 2.0);
    let b: dim2::Vector<f64> = dim2::Vector::new(3.0, 4.0);

    c.bench_function("euclidean/dim2/vector_wedge", |bencher| {
        bencher.iter(|| black_box(a).wedge(&black_box(b)))
    });
}

fn bench_vec2_add(c: &mut Criterion) {
    let a: dim2::Vector<f64> = dim2::Vector::new(1.0, 2.0);
    let b: dim2::Vector<f64> = dim2::Vector::new(3.0, 4.0);

    c.bench_function("euclidean/dim2/vector_add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });
}

fn bench_rotor2_rotation(c: &mut Criterion) {
    let rotor = dim2::Rotor::<f64>::from_angle(FRAC_PI_4);
    let v = dim2::Vector::new(1.0, 0.0);

    c.bench_function("euclidean/dim2/rotor_rotate", |bencher| {
        bencher.iter(|| black_box(rotor).transform(&black_box(v)))
    });
}

fn bench_rotor2_compose(c: &mut Criterion) {
    let r1 = dim2::Rotor::<f64>::from_angle(FRAC_PI_4);
    let r2 = dim2::Rotor::<f64>::from_angle(FRAC_PI_4 / 2.0);

    c.bench_function("euclidean/dim2/rotor_compose", |bencher| {
        bencher.iter(|| black_box(r1).compose(black_box(r2)))
    });
}

fn bench_rotor2_slerp(c: &mut Criterion) {
    let r1 = dim2::Rotor::<f64>::from_angle(0.0);
    let r2 = dim2::Rotor::<f64>::from_angle(FRAC_PI_4);

    c.bench_function("euclidean/dim2/rotor_slerp", |bencher| {
        bencher.iter(|| black_box(r1).slerp(black_box(r2), black_box(0.5)))
    });
}

// ============================================================================
// 3D Benchmarks
// ============================================================================

fn bench_vec3_dot(c: &mut Criterion) {
    let a: dim3::Vector<f64> = dim3::Vector::new(1.0, 2.0, 3.0);
    let b: dim3::Vector<f64> = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("euclidean/dim3/vector_dot", |bencher| {
        bencher.iter(|| black_box(a).dot(black_box(b)))
    });
}

fn bench_vec3_wedge(c: &mut Criterion) {
    let a: dim3::Vector<f64> = dim3::Vector::new(1.0, 2.0, 3.0);
    let b: dim3::Vector<f64> = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("euclidean/dim3/vector_wedge", |bencher| {
        bencher.iter(|| black_box(a).wedge(&black_box(b)))
    });
}

fn bench_vec3_cross(c: &mut Criterion) {
    let a: dim3::Vector<f64> = dim3::Vector::new(1.0, 2.0, 3.0);
    let b: dim3::Vector<f64> = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("euclidean/dim3/vector_cross", |bencher| {
        bencher.iter(|| black_box(a).cross(black_box(b)))
    });
}

fn bench_vec3_add(c: &mut Criterion) {
    let a: dim3::Vector<f64> = dim3::Vector::new(1.0, 2.0, 3.0);
    let b: dim3::Vector<f64> = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("euclidean/dim3/vector_add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });
}

fn bench_rotor3_rotation(c: &mut Criterion) {
    let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_rz());
    let v = dim3::Vector::new(1.0, 0.0, 0.0);

    c.bench_function("euclidean/dim3/rotor_rotate", |bencher| {
        bencher.iter(|| black_box(rotor).transform(&black_box(v)))
    });
}

fn bench_rotor3_compose(c: &mut Criterion) {
    let r1 = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_rz());
    let r2 = dim3::Rotor::from_angle_plane(FRAC_PI_4 / 2.0, dim3::Bivector::unit_ry());

    c.bench_function("euclidean/dim3/rotor_compose", |bencher| {
        bencher.iter(|| black_box(r1).compose(black_box(r2)))
    });
}

fn bench_rotor3_slerp(c: &mut Criterion) {
    let r1 = dim3::Rotor::from_angle_plane(0.0, dim3::Bivector::unit_rz());
    let r2 = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_rz());

    c.bench_function("euclidean/dim3/rotor_slerp", |bencher| {
        bencher.iter(|| black_box(r1).slerp(black_box(r2), black_box(0.5)))
    });
}

fn bench_rotor3_from_vectors(c: &mut Criterion) {
    let a = dim3::Vector::<f64>::unit_x();
    let b = dim3::Vector::new(0.707, 0.707, 0.0).normalized();

    c.bench_function("euclidean/dim3/rotor_from_vectors", |bencher| {
        bencher.iter(|| dim3::Rotor::from_vectors(black_box(a), black_box(b)))
    });
}

criterion_group!(
    benches,
    // 2D benchmarks
    bench_vec2_dot,
    bench_vec2_wedge,
    bench_vec2_add,
    bench_rotor2_rotation,
    bench_rotor2_compose,
    bench_rotor2_slerp,
    // 3D benchmarks
    bench_vec3_dot,
    bench_vec3_wedge,
    bench_vec3_cross,
    bench_vec3_add,
    bench_rotor3_rotation,
    bench_rotor3_compose,
    bench_rotor3_slerp,
    bench_rotor3_from_vectors,
);
criterion_main!(benches);
