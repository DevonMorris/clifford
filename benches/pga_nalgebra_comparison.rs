//! Benchmarks comparing PGA operations with nalgebra equivalents.
//!
//! Run with nalgebra 0.34:
//! ```bash
//! cargo bench --bench pga_nalgebra_comparison --features nalgebra-0_34
//! ```
//!
//! Run with nalgebra 0.33:
//! ```bash
//! cargo bench --bench pga_nalgebra_comparison --features nalgebra-0_33
//! ```

#![allow(missing_docs)]

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::FRAC_PI_4;
use std::hint::black_box;

#[cfg(feature = "nalgebra-0_32")]
use nalgebra_0_32 as na;
#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use clifford::specialized::projective::dim3;

// ============================================================================
// 3D Motor vs Isometry3 - Transform Point
// ============================================================================

fn bench_pga3_motor_transform_clifford(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4)
        .compose(&dim3::Motor::from_translation(1.0, 2.0, 3.0));
    let p = dim3::Point::new(1.0, 0.0, 0.0);

    c.bench_function("pga_comparison/3d_transform_point/clifford", |bencher| {
        bencher.iter(|| black_box(motor).transform_point(&black_box(p)))
    });
}

fn bench_pga3_isometry_transform_nalgebra(c: &mut Criterion) {
    let rotation = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let translation = na::Translation3::new(1.0, 2.0, 3.0);
    let iso = na::Isometry3::from_parts(translation, rotation);
    let p = na::Point3::new(1.0, 0.0, 0.0);

    c.bench_function("pga_comparison/3d_transform_point/nalgebra", |bencher| {
        bencher.iter(|| black_box(iso) * black_box(p))
    });
}

// ============================================================================
// 3D Motor vs Isometry3 - Compose
// ============================================================================

fn bench_pga3_motor_compose_clifford(c: &mut Criterion) {
    let m1 = dim3::Motor::from_rotation_z(FRAC_PI_4);
    let m2 = dim3::Motor::from_translation(1.0, 2.0, 3.0);

    c.bench_function("pga_comparison/3d_compose/clifford", |bencher| {
        bencher.iter(|| black_box(m1).compose(&black_box(m2)))
    });
}

fn bench_pga3_isometry_compose_nalgebra(c: &mut Criterion) {
    let rotation = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let iso1 = na::Isometry3::from_parts(na::Translation3::identity(), rotation);
    let iso2 = na::Isometry3::translation(1.0, 2.0, 3.0);

    c.bench_function("pga_comparison/3d_compose/nalgebra", |bencher| {
        bencher.iter(|| black_box(iso1) * black_box(iso2))
    });
}

// ============================================================================
// 3D Motor vs Isometry3 - Inverse
// ============================================================================

fn bench_pga3_motor_inverse_clifford(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4)
        .compose(&dim3::Motor::from_translation(1.0, 2.0, 3.0));

    c.bench_function("pga_comparison/3d_inverse/clifford", |bencher| {
        bencher.iter(|| black_box(motor).inverse())
    });
}

fn bench_pga3_isometry_inverse_nalgebra(c: &mut Criterion) {
    let rotation = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let translation = na::Translation3::new(1.0, 2.0, 3.0);
    let iso = na::Isometry3::from_parts(translation, rotation);

    c.bench_function("pga_comparison/3d_inverse/nalgebra", |bencher| {
        bencher.iter(|| black_box(iso).inverse())
    });
}

// ============================================================================
// Batch Transform - 100 points
// ============================================================================

fn bench_pga3_batch_clifford(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4)
        .compose(&dim3::Motor::from_translation(1.0, 2.0, 3.0));
    let points: Vec<dim3::Point<f64>> = (0..100)
        .map(|i| dim3::Point::new(i as f64, 0.0, 0.0))
        .collect();

    c.bench_function("pga_comparison/3d_batch_100/clifford", |bencher| {
        bencher.iter(|| {
            points
                .iter()
                .map(|p| black_box(motor).transform_point(p))
                .collect::<Vec<_>>()
        })
    });
}

fn bench_pga3_batch_nalgebra(c: &mut Criterion) {
    let rotation = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let translation = na::Translation3::new(1.0, 2.0, 3.0);
    let iso = na::Isometry3::from_parts(translation, rotation);
    let points: Vec<na::Point3<f64>> = (0..100)
        .map(|i| na::Point3::new(i as f64, 0.0, 0.0))
        .collect();

    c.bench_function("pga_comparison/3d_batch_100/nalgebra", |bencher| {
        bencher.iter(|| {
            points
                .iter()
                .map(|p| black_box(iso) * p)
                .collect::<Vec<_>>()
        })
    });
}

// ============================================================================
// Point Distance
// ============================================================================

fn bench_pga3_point_distance_clifford(c: &mut Criterion) {
    let p1 = dim3::Point::new(0.0, 0.0, 0.0);
    let p2 = dim3::Point::new(1.0, 2.0, 3.0);

    c.bench_function("pga_comparison/3d_point_distance/clifford", |bencher| {
        bencher.iter(|| black_box(p1).distance(&black_box(p2)))
    });
}

fn bench_pga3_point_distance_nalgebra(c: &mut Criterion) {
    let p1 = na::Point3::new(0.0, 0.0, 0.0);
    let p2 = na::Point3::new(1.0, 2.0, 3.0);

    c.bench_function("pga_comparison/3d_point_distance/nalgebra", |bencher| {
        bencher.iter(|| na::distance(&black_box(p1), &black_box(p2)))
    });
}

// ============================================================================
// Conversions
// ============================================================================

fn bench_pga3_motor_to_isometry(c: &mut Criterion) {
    let motor = dim3::Motor::from_rotation_z(FRAC_PI_4)
        .compose(&dim3::Motor::from_translation(1.0, 2.0, 3.0));

    c.bench_function("pga_conversion/motor_to_isometry3", |bencher| {
        bencher.iter(|| {
            let iso: na::Isometry3<f64> = black_box(motor).into();
            black_box(iso)
        })
    });
}

fn bench_pga3_isometry_to_motor(c: &mut Criterion) {
    let rotation = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let translation = na::Translation3::new(1.0, 2.0, 3.0);
    let iso = na::Isometry3::from_parts(translation, rotation);

    c.bench_function("pga_conversion/isometry3_to_motor", |bencher| {
        bencher.iter(|| {
            let motor: dim3::Motor<f64> = black_box(iso).into();
            black_box(motor)
        })
    });
}

fn bench_pga3_point_to_nalgebra(c: &mut Criterion) {
    let p = dim3::Point::new(1.0, 2.0, 3.0);

    c.bench_function("pga_conversion/point3_to_nalgebra", |bencher| {
        bencher.iter(|| {
            let na_p: na::Point3<f64> = black_box(p).try_into().unwrap();
            black_box(na_p)
        })
    });
}

fn bench_pga3_point_from_nalgebra(c: &mut Criterion) {
    let p = na::Point3::new(1.0, 2.0, 3.0);

    c.bench_function("pga_conversion/point3_from_nalgebra", |bencher| {
        bencher.iter(|| {
            let pga_p: dim3::Point<f64> = black_box(p).into();
            black_box(pga_p)
        })
    });
}

// ============================================================================
// PGA-specific operations (no nalgebra equivalent)
// ============================================================================

fn bench_pga3_line_meet_plane(c: &mut Criterion) {
    let line: dim3::Line<f64> = dim3::Line::z_axis();
    let plane: dim3::Plane<f64> = dim3::Plane::from_normal_and_distance(0.0, 0.0, 1.0, -5.0);

    c.bench_function("pga_only/3d_line_meet_plane", |bencher| {
        bencher.iter(|| black_box(line).meet(&black_box(plane)))
    });
}

fn bench_pga3_line_join_point(c: &mut Criterion) {
    let line: dim3::Line<f64> = dim3::Line::z_axis();
    let p = dim3::Point::new(1.0, 0.0, 0.0);

    c.bench_function("pga_only/3d_line_join_point", |bencher| {
        bencher.iter(|| black_box(line).join_point(&black_box(p)))
    });
}

fn bench_pga3_motor_commutator(c: &mut Criterion) {
    let m1 = dim3::Motor::from_rotation_x(0.1);
    let m2 = dim3::Motor::from_rotation_z(0.2);

    c.bench_function("pga_only/3d_motor_commutator", |bencher| {
        bencher.iter(|| black_box(m1).commutator(&black_box(m2)))
    });
}

fn bench_pga3_left_contraction(c: &mut Criterion) {
    let p = dim3::Point::new(1.0, 0.0, 0.0);
    let line: dim3::Line<f64> = dim3::Line::z_axis();

    c.bench_function("pga_only/3d_point_left_contract_line", |bencher| {
        bencher.iter(|| black_box(p).left_contract_line(&black_box(line)))
    });
}

// ============================================================================
// Criterion Groups
// ============================================================================

criterion_group!(
    pga3_comparison,
    bench_pga3_motor_transform_clifford,
    bench_pga3_isometry_transform_nalgebra,
    bench_pga3_motor_compose_clifford,
    bench_pga3_isometry_compose_nalgebra,
    bench_pga3_motor_inverse_clifford,
    bench_pga3_isometry_inverse_nalgebra,
    bench_pga3_batch_clifford,
    bench_pga3_batch_nalgebra,
    bench_pga3_point_distance_clifford,
    bench_pga3_point_distance_nalgebra,
);

criterion_group!(
    pga_conversions,
    bench_pga3_motor_to_isometry,
    bench_pga3_isometry_to_motor,
    bench_pga3_point_to_nalgebra,
    bench_pga3_point_from_nalgebra,
);

criterion_group!(
    pga_specific,
    bench_pga3_line_meet_plane,
    bench_pga3_line_join_point,
    bench_pga3_motor_commutator,
    bench_pga3_left_contraction,
);

criterion_main!(pga3_comparison, pga_conversions, pga_specific);
