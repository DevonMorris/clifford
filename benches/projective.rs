//! Benchmarks for Projective Geometric Algebra (PGA) types.
//!
//! Run with:
//! ```bash
//! cargo bench --bench projective
//! ```
//!
//! Run with nalgebra comparison:
//! ```bash
//! cargo bench --bench projective --features nalgebra-0_34
//! ```

#![allow(missing_docs)]

use std::f64::consts::FRAC_PI_4;
use std::hint::black_box;

use clifford::specialized::projective::{dim2, dim3};
use criterion::{Criterion, criterion_group, criterion_main};

// ============================================================================
// 2D PGA Benchmarks
// ============================================================================

fn bench_pga2_point_join(c: &mut Criterion) {
    let p1 = dim2::Point::from_cartesian(0.0, 0.0);
    let p2 = dim2::Point::from_cartesian(1.0, 1.0);

    c.bench_function("projective/dim2/point_join", |bencher| {
        bencher.iter(|| black_box(p1).join(&black_box(p2)))
    });
}

fn bench_pga2_line_meet(c: &mut Criterion) {
    let l1: dim2::Line<f64> = dim2::Line::x_axis();
    let l2: dim2::Line<f64> = dim2::Line::y_axis();

    c.bench_function("projective/dim2/line_meet", |bencher| {
        bencher.iter(|| black_box(l1).meet(&black_box(l2)))
    });
}

fn bench_pga2_motor_transform_point(c: &mut Criterion) {
    let motor = dim2::Motor::from_rotation(FRAC_PI_4);
    let p = dim2::Point::from_cartesian(3.0, 4.0);

    c.bench_function("projective/dim2/motor_transform_point", |bencher| {
        bencher.iter(|| black_box(motor).transform_point(&black_box(p)))
    });
}

fn bench_pga2_motor_transform_line(c: &mut Criterion) {
    let motor = dim2::Motor::from_rotation(FRAC_PI_4);
    let line: dim2::Line<f64> = dim2::Line::from_implicit(1.0, 1.0, 0.0);

    c.bench_function("projective/dim2/motor_transform_line", |bencher| {
        bencher.iter(|| black_box(motor).transform_line(&black_box(line)))
    });
}

fn bench_pga2_motor_inverse(c: &mut Criterion) {
    let motor = dim2::Motor::from_rotation(FRAC_PI_4);

    c.bench_function("projective/dim2/motor_inverse", |bencher| {
        bencher.iter(|| black_box(motor).inverse())
    });
}

fn bench_pga2_point_distance(c: &mut Criterion) {
    let p1 = dim2::Point::from_cartesian(0.0, 0.0);
    let p2 = dim2::Point::from_cartesian(3.0, 4.0);

    c.bench_function("projective/dim2/point_distance", |bencher| {
        bencher.iter(|| black_box(p1).distance(&black_box(p2)))
    });
}

fn bench_pga2_line_distance_to_point(c: &mut Criterion) {
    let line: dim2::Line<f64> = dim2::Line::x_axis();
    let p = dim2::Point::from_cartesian(3.0, 4.0);

    c.bench_function("projective/dim2/line_distance_to_point", |bencher| {
        bencher.iter(|| black_box(line).distance_to_point(&black_box(p)))
    });
}

// ============================================================================
// 3D PGA Benchmarks
// ============================================================================

fn bench_pga3_motor_transform_point(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4);
    let p = dim3::Point::from_cartesian(1.0, 0.0, 0.0);

    c.bench_function("projective/dim3/motor_transform_point", |bencher| {
        bencher.iter(|| black_box(motor).transform_point(&black_box(p)))
    });
}

fn bench_pga3_motor_transform_line(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4);
    let line: dim3::Line<f64> = dim3::Line::x_axis();

    c.bench_function("projective/dim3/motor_transform_line", |bencher| {
        bencher.iter(|| black_box(motor).transform_line(&black_box(line)))
    });
}

fn bench_pga3_motor_inverse(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4);

    c.bench_function("projective/dim3/motor_inverse", |bencher| {
        bencher.iter(|| black_box(motor).inverse())
    });
}

fn bench_pga3_line_meet_plane(c: &mut Criterion) {
    let line: dim3::Line<f64> = dim3::Line::z_axis();
    let plane: dim3::Plane<f64> = dim3::Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);

    c.bench_function("projective/dim3/line_meet_plane", |bencher| {
        bencher.iter(|| black_box(line).meet(&black_box(plane)))
    });
}

fn bench_pga3_line_join_point(c: &mut Criterion) {
    let line: dim3::Line<f64> = dim3::Line::z_axis();
    let p = dim3::Point::from_cartesian(1.0, 0.0, 0.0);

    c.bench_function("projective/dim3/line_join_point", |bencher| {
        bencher.iter(|| black_box(line).join_point(&black_box(p)))
    });
}

fn bench_pga3_line_distance_to_point(c: &mut Criterion) {
    let line: dim3::Line<f64> = dim3::Line::z_axis();
    let p = dim3::Point::from_cartesian(3.0, 4.0, 5.0);

    c.bench_function("projective/dim3/line_distance_to_point", |bencher| {
        bencher.iter(|| black_box(line).distance_to_point(&black_box(p)))
    });
}

fn bench_pga3_line_closest_point(c: &mut Criterion) {
    let line: dim3::Line<f64> = dim3::Line::z_axis();
    let p = dim3::Point::from_cartesian(3.0, 4.0, 5.0);

    c.bench_function("projective/dim3/line_closest_point", |bencher| {
        bencher.iter(|| black_box(line).closest_point(&black_box(p)))
    });
}

fn bench_pga3_point_distance(c: &mut Criterion) {
    let p1 = dim3::Point::from_cartesian(0.0, 0.0, 0.0);
    let p2 = dim3::Point::from_cartesian(1.0, 2.0, 3.0);

    c.bench_function("projective/dim3/point_distance", |bencher| {
        bencher.iter(|| black_box(p1).distance(&black_box(p2)))
    });
}

fn bench_pga3_point_dot(c: &mut Criterion) {
    let p1 = dim3::Point::from_cartesian(1.0, 2.0, 3.0);
    let p2 = dim3::Point::from_cartesian(4.0, 5.0, 6.0);

    c.bench_function("projective/dim3/point_dot", |bencher| {
        bencher.iter(|| black_box(p1).dot(&black_box(p2)))
    });
}

fn bench_pga3_line_dot(c: &mut Criterion) {
    let l1: dim3::Line<f64> = dim3::Line::x_axis();
    let l2: dim3::Line<f64> = dim3::Line::y_axis();

    c.bench_function("projective/dim3/line_dot", |bencher| {
        bencher.iter(|| black_box(l1).dot(&black_box(l2)))
    });
}

fn bench_pga3_line_distance(c: &mut Criterion) {
    let l1: dim3::Line<f64> = dim3::Line::z_axis();
    let l2 = dim3::Line::from_point_and_direction(
        &dim3::Point::from_cartesian(1.0, 0.0, 0.0),
        &clifford::specialized::euclidean::dim3::Vector::new(0.0, 1.0, 0.0),
    );

    c.bench_function("projective/dim3/line_distance", |bencher| {
        bencher.iter(|| black_box(l1).distance(&black_box(l2)))
    });
}

fn bench_pga3_point_left_contract_line(c: &mut Criterion) {
    let p = dim3::Point::from_cartesian(1.0, 0.0, 0.0);
    let line: dim3::Line<f64> = dim3::Line::z_axis();

    c.bench_function("projective/dim3/point_left_contract_line", |bencher| {
        bencher.iter(|| black_box(p).left_contract_line(&black_box(line)))
    });
}

fn bench_pga3_motor_transform_plane(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4);
    let plane: dim3::Plane<f64> = dim3::Plane::from_normal_and_distance(1.0, 0.0, 0.0, 5.0);

    c.bench_function("projective/dim3/motor_transform_plane", |bencher| {
        bencher.iter(|| black_box(motor).transform_plane(&black_box(plane)))
    });
}

// ============================================================================
// Batch operations
// ============================================================================

fn bench_pga3_batch_transform_points(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4);
    let points: Vec<dim3::Point<f64>> = (0..100)
        .map(|i| dim3::Point::from_cartesian(i as f64, 0.0, 0.0))
        .collect();

    c.bench_function("projective/dim3/batch_transform_100_points", |bencher| {
        bencher.iter(|| {
            points
                .iter()
                .map(|p| black_box(motor).transform_point(p))
                .collect::<Vec<_>>()
        })
    });
}

fn bench_pga2_batch_transform_points(c: &mut Criterion) {
    let motor = dim2::Motor::from_rotation(FRAC_PI_4);
    let points: Vec<dim2::Point<f64>> = (0..100)
        .map(|i| dim2::Point::from_cartesian(i as f64, 0.0))
        .collect();

    c.bench_function("projective/dim2/batch_transform_100_points", |bencher| {
        bencher.iter(|| {
            points
                .iter()
                .map(|p| black_box(motor).transform_point(p))
                .collect::<Vec<_>>()
        })
    });
}

// ============================================================================
// Criterion Groups
// ============================================================================

criterion_group!(
    pga2_benches,
    bench_pga2_point_join,
    bench_pga2_line_meet,
    bench_pga2_motor_transform_point,
    bench_pga2_motor_transform_line,
    bench_pga2_motor_inverse,
    bench_pga2_point_distance,
    bench_pga2_line_distance_to_point,
    bench_pga2_batch_transform_points,
);

criterion_group!(
    pga3_benches,
    bench_pga3_motor_transform_point,
    bench_pga3_motor_transform_line,
    bench_pga3_motor_transform_plane,
    bench_pga3_motor_inverse,
    bench_pga3_line_meet_plane,
    bench_pga3_line_join_point,
    bench_pga3_line_distance_to_point,
    bench_pga3_line_closest_point,
    bench_pga3_point_distance,
    bench_pga3_point_dot,
    bench_pga3_line_dot,
    bench_pga3_line_distance,
    bench_pga3_point_left_contract_line,
    bench_pga3_batch_transform_points,
);

criterion_main!(pga2_benches, pga3_benches);
