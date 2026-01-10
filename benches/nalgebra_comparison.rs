//! Benchmarks comparing clifford operations with nalgebra equivalents.
//!
//! Run with nalgebra 0.34:
//! ```bash
//! cargo bench --bench nalgebra_comparison --features nalgebra-0_34
//! ```
//!
//! Run with nalgebra 0.33:
//! ```bash
//! cargo bench --bench nalgebra_comparison --features nalgebra-0_33
//! ```

#![allow(missing_docs)]

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::FRAC_PI_4;
use std::hint::black_box;

#[cfg(feature = "nalgebra-0_33")]
use nalgebra_0_33 as na;
#[cfg(feature = "nalgebra-0_34")]
use nalgebra_0_34 as na;

use clifford::specialized::euclidean::{dim2, dim3};

// ============================================================================
// 3D Vector Operations - Comparison
// ============================================================================

fn bench_vec3_dot_clifford(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 2.0, 3.0);
    let b = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("comparison/vec3_dot/clifford", |bencher| {
        bencher.iter(|| black_box(a).dot(black_box(b)))
    });
}

fn bench_vec3_dot_nalgebra(c: &mut Criterion) {
    let a = na::Vector3::new(1.0, 2.0, 3.0);
    let b = na::Vector3::new(4.0, 5.0, 6.0);

    c.bench_function("comparison/vec3_dot/nalgebra", |bencher| {
        bencher.iter(|| black_box(a).dot(&black_box(b)))
    });
}

fn bench_vec3_cross_clifford(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 2.0, 3.0);
    let b = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("comparison/vec3_cross/clifford", |bencher| {
        bencher.iter(|| black_box(a).cross(black_box(b)))
    });
}

fn bench_vec3_cross_nalgebra(c: &mut Criterion) {
    let a = na::Vector3::new(1.0, 2.0, 3.0);
    let b = na::Vector3::new(4.0, 5.0, 6.0);

    c.bench_function("comparison/vec3_cross/nalgebra", |bencher| {
        bencher.iter(|| black_box(a).cross(&black_box(b)))
    });
}

fn bench_vec3_normalize_clifford(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 2.0, 3.0);

    c.bench_function("comparison/vec3_normalize/clifford", |bencher| {
        bencher.iter(|| black_box(a).normalized())
    });
}

fn bench_vec3_normalize_nalgebra(c: &mut Criterion) {
    let a = na::Vector3::new(1.0, 2.0, 3.0);

    c.bench_function("comparison/vec3_normalize/nalgebra", |bencher| {
        bencher.iter(|| black_box(a).normalize())
    });
}

// ============================================================================
// 3D Rotation Operations - Comparison
// ============================================================================

fn bench_rotate_vec3_clifford(c: &mut Criterion) {
    let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());
    let v = dim3::Vector::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rotate_vec3/clifford", |bencher| {
        bencher.iter(|| black_box(rotor).rotate(black_box(v)))
    });
}

fn bench_rotate_vec3_nalgebra_quat(c: &mut Criterion) {
    let q = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let v = na::Vector3::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rotate_vec3/nalgebra_quat", |bencher| {
        bencher.iter(|| black_box(q) * black_box(v))
    });
}

fn bench_rotate_vec3_nalgebra_rot(c: &mut Criterion) {
    let r = na::Rotation3::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let v = na::Vector3::new(1.0, 0.0, 0.0);

    c.bench_function("comparison/rotate_vec3/nalgebra_rot", |bencher| {
        bencher.iter(|| black_box(r) * black_box(v))
    });
}

fn bench_rotation_compose_clifford(c: &mut Criterion) {
    let r1 = dim3::Rotor::from_angle_plane(0.1, dim3::Bivector::unit_xy());
    let r2 = dim3::Rotor::from_angle_plane(0.2, dim3::Bivector::unit_xz());

    c.bench_function("comparison/rotation_compose/clifford", |bencher| {
        bencher.iter(|| black_box(r1).compose(black_box(r2)))
    });
}

fn bench_rotation_compose_nalgebra(c: &mut Criterion) {
    let q1 = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), 0.1);
    let q2 = na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), 0.2);

    c.bench_function("comparison/rotation_compose/nalgebra", |bencher| {
        bencher.iter(|| black_box(q1) * black_box(q2))
    });
}

fn bench_rotation_slerp_clifford(c: &mut Criterion) {
    let r1 = dim3::Rotor::from_angle_plane(0.0, dim3::Bivector::unit_xy());
    let r2 = dim3::Rotor::from_angle_plane(1.0, dim3::Bivector::unit_xy());

    c.bench_function("comparison/rotation_slerp/clifford", |bencher| {
        bencher.iter(|| black_box(r1).slerp(black_box(r2), 0.5))
    });
}

fn bench_rotation_slerp_nalgebra(c: &mut Criterion) {
    let q1 = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), 0.0);
    let q2 = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), 1.0);

    c.bench_function("comparison/rotation_slerp/nalgebra", |bencher| {
        bencher.iter(|| black_box(q1).slerp(&black_box(q2), 0.5))
    });
}

// ============================================================================
// 2D Operations - Comparison
// ============================================================================

fn bench_vec2_dot_clifford(c: &mut Criterion) {
    let a = dim2::Vector::new(1.0, 2.0);
    let b = dim2::Vector::new(3.0, 4.0);

    c.bench_function("comparison/vec2_dot/clifford", |bencher| {
        bencher.iter(|| black_box(a).dot(black_box(b)))
    });
}

fn bench_vec2_dot_nalgebra(c: &mut Criterion) {
    let a = na::Vector2::new(1.0, 2.0);
    let b = na::Vector2::new(3.0, 4.0);

    c.bench_function("comparison/vec2_dot/nalgebra", |bencher| {
        bencher.iter(|| black_box(a).dot(&black_box(b)))
    });
}

fn bench_rotate_vec2_clifford(c: &mut Criterion) {
    let rotor = dim2::Rotor::from_angle(FRAC_PI_4);
    let v = dim2::Vector::new(1.0, 0.0);

    c.bench_function("comparison/rotate_vec2/clifford", |bencher| {
        bencher.iter(|| black_box(rotor).rotate(black_box(v)))
    });
}

fn bench_rotate_vec2_nalgebra(c: &mut Criterion) {
    let r = na::Rotation2::new(FRAC_PI_4);
    let v = na::Vector2::new(1.0, 0.0);

    c.bench_function("comparison/rotate_vec2/nalgebra", |bencher| {
        bencher.iter(|| black_box(r) * black_box(v))
    });
}

// ============================================================================
// Conversion Overhead
// ============================================================================

fn bench_vec3_to_nalgebra(c: &mut Criterion) {
    let v = dim3::Vector::new(1.0, 2.0, 3.0);

    c.bench_function("conversion/vec3_to_nalgebra", |bencher| {
        bencher.iter(|| {
            let na_v: na::Vector3<f64> = black_box(v).into();
            black_box(na_v)
        })
    });
}

fn bench_vec3_from_nalgebra(c: &mut Criterion) {
    let v = na::Vector3::new(1.0, 2.0, 3.0);

    c.bench_function("conversion/vec3_from_nalgebra", |bencher| {
        bencher.iter(|| {
            let cliff_v: dim3::Vector<f64> = black_box(v).into();
            black_box(cliff_v)
        })
    });
}

fn bench_rotor_to_quaternion(c: &mut Criterion) {
    let r = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());

    c.bench_function("conversion/rotor_to_quaternion", |bencher| {
        bencher.iter(|| {
            let q: na::UnitQuaternion<f64> = black_box(r).into();
            black_box(q)
        })
    });
}

fn bench_quaternion_to_rotor(c: &mut Criterion) {
    let q = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);

    c.bench_function("conversion/quaternion_to_rotor", |bencher| {
        bencher.iter(|| {
            let r: dim3::Rotor<f64> = black_box(q).into();
            black_box(r)
        })
    });
}

// ============================================================================
// Realistic Workflows
// ============================================================================

fn bench_batch_rotate_clifford(c: &mut Criterion) {
    let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());
    let points: Vec<dim3::Vector<f64>> = (0..100)
        .map(|i| dim3::Vector::new(i as f64, 0.0, 0.0))
        .collect();

    c.bench_function("workflow/batch_rotate_100/clifford", |bencher| {
        bencher.iter(|| {
            points
                .iter()
                .map(|p| black_box(rotor).rotate(*p))
                .collect::<Vec<_>>()
        })
    });
}

fn bench_batch_rotate_nalgebra(c: &mut Criterion) {
    let q = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), FRAC_PI_4);
    let points: Vec<na::Vector3<f64>> = (0..100)
        .map(|i| na::Vector3::new(i as f64, 0.0, 0.0))
        .collect();

    c.bench_function("workflow/batch_rotate_100/nalgebra", |bencher| {
        bencher.iter(|| points.iter().map(|p| black_box(q) * p).collect::<Vec<_>>())
    });
}

fn bench_mixed_workflow(c: &mut Criterion) {
    let na_points: Vec<na::Vector3<f64>> = (0..10)
        .map(|i| na::Vector3::new(i as f64, 0.0, 0.0))
        .collect();
    let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_4, dim3::Bivector::unit_xy());

    c.bench_function("workflow/mixed_convert_rotate_convert", |bencher| {
        bencher.iter(|| {
            na_points
                .iter()
                .map(|p| {
                    // Convert from nalgebra
                    let cliff_p: dim3::Vector<f64> = (*p).into();
                    // Process with GA
                    let rotated = black_box(rotor).rotate(cliff_p);
                    // Convert back to nalgebra
                    let result: na::Vector3<f64> = rotated.into();
                    result
                })
                .collect::<Vec<_>>()
        })
    });
}

// ============================================================================
// GA-Specific Operations (no nalgebra equivalent)
// ============================================================================

fn bench_vec3_wedge(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 0.0, 0.0);
    let b = dim3::Vector::new(0.0, 1.0, 0.0);

    c.bench_function("ga_only/vec3_wedge", |bencher| {
        bencher.iter(|| black_box(a).wedge(black_box(b)))
    });
}

fn bench_rotor_chain(c: &mut Criterion) {
    let r1 = dim3::Rotor::from_angle_plane(0.1, dim3::Bivector::unit_xy());
    let r2 = dim3::Rotor::from_angle_plane(0.2, dim3::Bivector::unit_xz());
    let r3 = dim3::Rotor::from_angle_plane(0.3, dim3::Bivector::unit_yz());

    c.bench_function("ga_only/rotor_chain_3", |bencher| {
        bencher.iter(|| black_box(r1).compose(black_box(r2)).compose(black_box(r3)))
    });
}

fn bench_vec3_geometric(c: &mut Criterion) {
    let a = dim3::Vector::new(1.0, 2.0, 3.0);
    let b = dim3::Vector::new(4.0, 5.0, 6.0);

    c.bench_function("ga_only/vec3_geometric", |bencher| {
        bencher.iter(|| black_box(a).geometric(black_box(b)))
    });
}

// ============================================================================
// Criterion Groups
// ============================================================================

criterion_group!(
    vec3_comparison,
    bench_vec3_dot_clifford,
    bench_vec3_dot_nalgebra,
    bench_vec3_cross_clifford,
    bench_vec3_cross_nalgebra,
    bench_vec3_normalize_clifford,
    bench_vec3_normalize_nalgebra,
);

criterion_group!(
    rotation_comparison,
    bench_rotate_vec3_clifford,
    bench_rotate_vec3_nalgebra_quat,
    bench_rotate_vec3_nalgebra_rot,
    bench_rotation_compose_clifford,
    bench_rotation_compose_nalgebra,
    bench_rotation_slerp_clifford,
    bench_rotation_slerp_nalgebra,
);

criterion_group!(
    vec2_comparison,
    bench_vec2_dot_clifford,
    bench_vec2_dot_nalgebra,
    bench_rotate_vec2_clifford,
    bench_rotate_vec2_nalgebra,
);

criterion_group!(
    conversion_overhead,
    bench_vec3_to_nalgebra,
    bench_vec3_from_nalgebra,
    bench_rotor_to_quaternion,
    bench_quaternion_to_rotor,
);

criterion_group!(
    workflows,
    bench_batch_rotate_clifford,
    bench_batch_rotate_nalgebra,
    bench_mixed_workflow,
);

criterion_group!(
    ga_specific,
    bench_vec3_wedge,
    bench_rotor_chain,
    bench_vec3_geometric,
);

criterion_main!(
    vec3_comparison,
    rotation_comparison,
    vec2_comparison,
    conversion_overhead,
    workflows,
    ga_specific,
);
