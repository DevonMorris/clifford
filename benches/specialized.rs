//! Benchmarks for specialized 2D/3D geometric algebra types.
//!
//! Run with: `cargo bench --bench specialized`

#![allow(missing_docs)]

use std::f64::consts::FRAC_PI_4;
use std::hint::black_box;

use clifford::specialized::ga2d::{Rotor2, Vec2};
use clifford::specialized::ga3d::{Bivec3, Rotor3, Vec3};
use criterion::{Criterion, criterion_group, criterion_main};

// ============================================================================
// 2D Benchmarks
// ============================================================================

fn bench_vec2_dot(c: &mut Criterion) {
    let a: Vec2<f64> = Vec2::new(1.0, 2.0);
    let b: Vec2<f64> = Vec2::new(3.0, 4.0);

    c.bench_function("ga2d/vector_dot", |bencher| {
        bencher.iter(|| black_box(a).dot(black_box(b)))
    });
}

fn bench_vec2_wedge(c: &mut Criterion) {
    let a: Vec2<f64> = Vec2::new(1.0, 2.0);
    let b: Vec2<f64> = Vec2::new(3.0, 4.0);

    c.bench_function("ga2d/vector_wedge", |bencher| {
        bencher.iter(|| black_box(a).wedge(black_box(b)))
    });
}

fn bench_vec2_add(c: &mut Criterion) {
    let a: Vec2<f64> = Vec2::new(1.0, 2.0);
    let b: Vec2<f64> = Vec2::new(3.0, 4.0);

    c.bench_function("ga2d/vector_add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });
}

fn bench_rotor2_rotation(c: &mut Criterion) {
    let rotor = Rotor2::<f64>::from_angle(FRAC_PI_4);
    let v = Vec2::new(1.0, 0.0);

    c.bench_function("ga2d/rotor_rotate", |bencher| {
        bencher.iter(|| black_box(rotor).rotate(black_box(v)))
    });
}

fn bench_rotor2_compose(c: &mut Criterion) {
    let r1 = Rotor2::<f64>::from_angle(FRAC_PI_4);
    let r2 = Rotor2::<f64>::from_angle(FRAC_PI_4 / 2.0);

    c.bench_function("ga2d/rotor_compose", |bencher| {
        bencher.iter(|| black_box(r1).compose(black_box(r2)))
    });
}

fn bench_rotor2_slerp(c: &mut Criterion) {
    let r1 = Rotor2::<f64>::identity();
    let r2 = Rotor2::<f64>::from_angle(FRAC_PI_4);

    c.bench_function("ga2d/rotor_slerp", |bencher| {
        bencher.iter(|| black_box(r1).slerp(black_box(r2), black_box(0.5)))
    });
}

// ============================================================================
// 3D Benchmarks
// ============================================================================

fn bench_vec3_dot(c: &mut Criterion) {
    let a: Vec3<f64> = Vec3::new(1.0, 2.0, 3.0);
    let b: Vec3<f64> = Vec3::new(4.0, 5.0, 6.0);

    c.bench_function("ga3d/vector_dot", |bencher| {
        bencher.iter(|| black_box(a).dot(black_box(b)))
    });
}

fn bench_vec3_wedge(c: &mut Criterion) {
    let a: Vec3<f64> = Vec3::new(1.0, 2.0, 3.0);
    let b: Vec3<f64> = Vec3::new(4.0, 5.0, 6.0);

    c.bench_function("ga3d/vector_wedge", |bencher| {
        bencher.iter(|| black_box(a).wedge(black_box(b)))
    });
}

fn bench_vec3_cross(c: &mut Criterion) {
    let a: Vec3<f64> = Vec3::new(1.0, 2.0, 3.0);
    let b: Vec3<f64> = Vec3::new(4.0, 5.0, 6.0);

    c.bench_function("ga3d/vector_cross", |bencher| {
        bencher.iter(|| black_box(a).cross(black_box(b)))
    });
}

fn bench_vec3_add(c: &mut Criterion) {
    let a: Vec3<f64> = Vec3::new(1.0, 2.0, 3.0);
    let b: Vec3<f64> = Vec3::new(4.0, 5.0, 6.0);

    c.bench_function("ga3d/vector_add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });
}

fn bench_rotor3_rotation(c: &mut Criterion) {
    let rotor = Rotor3::from_angle_plane(FRAC_PI_4, Bivec3::unit_xy());
    let v = Vec3::new(1.0, 0.0, 0.0);

    c.bench_function("ga3d/rotor_rotate", |bencher| {
        bencher.iter(|| black_box(rotor).rotate(black_box(v)))
    });
}

fn bench_rotor3_compose(c: &mut Criterion) {
    let r1 = Rotor3::from_angle_plane(FRAC_PI_4, Bivec3::unit_xy());
    let r2 = Rotor3::from_angle_plane(FRAC_PI_4 / 2.0, Bivec3::unit_xz());

    c.bench_function("ga3d/rotor_compose", |bencher| {
        bencher.iter(|| black_box(r1).compose(black_box(r2)))
    });
}

fn bench_rotor3_slerp(c: &mut Criterion) {
    let r1 = Rotor3::<f64>::identity();
    let r2 = Rotor3::from_angle_plane(FRAC_PI_4, Bivec3::unit_xy());

    c.bench_function("ga3d/rotor_slerp", |bencher| {
        bencher.iter(|| black_box(r1).slerp(black_box(r2), black_box(0.5)))
    });
}

fn bench_rotor3_from_vectors(c: &mut Criterion) {
    let a = Vec3::<f64>::unit_x();
    let b = Vec3::new(0.707, 0.707, 0.0).normalized();

    c.bench_function("ga3d/rotor_from_vectors", |bencher| {
        bencher.iter(|| Rotor3::from_vectors(black_box(a), black_box(b)))
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
