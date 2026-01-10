//! Rerun visualization example.
//!
//! This example demonstrates how to visualize clifford geometric algebra
//! types using the Rerun SDK.
//!
//! # Running
//!
//! ```bash
//! cargo run --example rerun_visualization --features rerun-0_28
//! ```
//!
//! This will spawn a Rerun viewer window and display the visualizations.

use std::f32::consts::{FRAC_PI_4, PI};

use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
use clifford::specialized::projective::dim3::{Motor, Point};
use clifford::specialized::visualization::{rerun, AsPosition};
use tracing::info;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let rec = rerun::RecordingStreamBuilder::new("clifford_demo").spawn()?;

    // =========================================================================
    // Euclidean 3D Examples
    // =========================================================================

    // Log a point using AsPosition wrapper
    let point = Vector::new(1.0_f32, 2.0, 3.0);
    rec.log("euclidean/point", &rerun::Points3D::new([AsPosition(point)]))?;

    // Log direction arrows
    let directions = [
        Vector::new(1.0_f32, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 0.0, 1.0),
    ];
    rec.log(
        "euclidean/axes",
        &rerun::Arrows3D::from_vectors(directions),
    )?;

    // Log a rotation as a transform
    let rotor = Rotor::from_angle_plane(FRAC_PI_4, Bivector::unit_xy());
    rec.log("euclidean/rotation", &rerun::Transform3D::from(rotor))?;

    // Visualize rotation by showing rotated vectors
    let original = Vector::new(1.0_f32, 0.0, 0.0);
    let mut rotated_points = vec![AsPosition(original)];

    for i in 1..=8 {
        let angle = (i as f32) * FRAC_PI_4;
        let r = Rotor::from_angle_plane(angle, Bivector::unit_xy());
        let rotated = r.rotate(original);
        rotated_points.push(AsPosition(rotated));
    }
    rec.log(
        "euclidean/rotation_path",
        &rerun::Points3D::new(rotated_points),
    )?;

    // =========================================================================
    // Projective 3D (PGA) Examples
    // =========================================================================

    // Log PGA points
    let pga_points = [
        Point::new(0.0_f32, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    rec.log("pga/square", &rerun::Points3D::new(pga_points))?;

    // Log a motor (rigid transformation)
    let motor = Motor::from_rotation_z(FRAC_PI_4).compose(&Motor::from_translation(2.0, 0.0, 0.0));
    rec.log("pga/motor_transform", &rerun::Transform3D::from(motor))?;

    // Visualize motor action on a point
    let start = Point::new(0.0_f32, 0.0, 0.0);
    let mut trajectory = vec![start];

    for i in 1..=16 {
        let t = (i as f32) / 16.0;
        let angle = t * 2.0 * PI;
        let m = Motor::from_rotation_z(angle).compose(&Motor::from_translation(t * 3.0, 0.0, t));
        let transformed = m.transform_point(&start);
        trajectory.push(transformed);
    }
    rec.log(
        "pga/spiral_trajectory",
        &rerun::Points3D::new(trajectory.iter().copied()),
    )?;

    info!("Logged visualization data to Rerun viewer");

    Ok(())
}
