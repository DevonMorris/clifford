//! Animated Rerun visualization example.
//!
//! This example demonstrates animating geometric algebra transformations
//! using the Rerun SDK with time-varying data.
//!
//! # Running
//!
//! ```bash
//! cargo run --example rerun_animation --features rerun-0_28
//! ```
//!
//! This will spawn a Rerun viewer showing:
//! - A vector being rotated by a rotor
//! - A cube being transformed by a motor (rotation + translation)
//! - A coordinate frame following a helical path

use std::f32::consts::TAU;

use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
use clifford::specialized::projective::dim3::{Motor, Point};
use clifford::specialized::visualization::{AsPosition, rerun};
use tracing::info;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let rec = rerun::RecordingStreamBuilder::new("clifford_animation").spawn()?;

    // Number of frames
    let num_frames = 200;

    for frame in 0..num_frames {
        let t = frame as f32 / num_frames as f32;

        // Set the timeline for this frame
        rec.set_time_sequence("frame", frame as i64);

        // =====================================================================
        // Animation 1: Rotating vector (rotor action)
        // =====================================================================
        {
            let angle = t * TAU; // Full rotation
            let rotor = Rotor::from_angle_plane(angle, Bivector::unit_xy());

            // Original vector along x-axis
            let original = Vector::new(2.0_f32, 0.0, 0.0);
            let rotated = rotor.rotate(original);

            // Log the rotating arrow
            rec.log(
                "rotor/rotating_vector",
                &rerun::Arrows3D::from_vectors([rotated])
                    .with_origins([[0.0, 0.0, 0.0]])
                    .with_colors([[255, 100, 100]]),
            )?;

            // Log a trail of the tip position
            rec.log(
                "rotor/trail",
                &rerun::Points3D::new([AsPosition(rotated)])
                    .with_colors([[255, 100, 100, 128]])
                    .with_radii([0.02]),
            )?;
        }

        // =====================================================================
        // Animation 2: Cube transformed by motor
        // =====================================================================
        {
            // Rotating and translating motor
            let rotation = Motor::from_rotation_z(t * TAU);
            let translation = Motor::from_translation(
                3.0 * (t * TAU).cos(),
                3.0 * (t * TAU).sin(),
                (t * TAU * 2.0).sin(),
            );
            let motor = rotation.compose(&translation);

            // Define a cube centered at origin
            let half = 0.5_f32;
            let cube_vertices = [
                Point::from_cartesian(-half, -half, -half),
                Point::from_cartesian(half, -half, -half),
                Point::from_cartesian(half, half, -half),
                Point::from_cartesian(-half, half, -half),
                Point::from_cartesian(-half, -half, half),
                Point::from_cartesian(half, -half, half),
                Point::from_cartesian(half, half, half),
                Point::from_cartesian(-half, half, half),
            ];

            // Transform all vertices
            let transformed: Vec<Point<f32>> = cube_vertices
                .iter()
                .map(|p| motor.transform_point(p))
                .collect();

            // Log transformed vertices
            rec.log(
                "motor/cube_vertices",
                &rerun::Points3D::new(transformed.iter().copied())
                    .with_colors([[100, 255, 100]])
                    .with_radii([0.08]),
            )?;

            // Log edges as line strips
            let edges = [
                // Bottom face
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 0),
                // Top face
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 4),
                // Vertical edges
                (0, 4),
                (1, 5),
                (2, 6),
                (3, 7),
            ];

            for (i, (a, b)) in edges.iter().enumerate() {
                let p1 = transformed[*a];
                let p2 = transformed[*b];
                rec.log(
                    format!("motor/cube_edges/{}", i),
                    &rerun::LineStrips3D::new([[
                        [p1.x(), p1.y(), p1.z()],
                        [p2.x(), p2.y(), p2.z()],
                    ]])
                    .with_colors([[100, 255, 100]]),
                )?;
            }

            // Log the motor as a transform
            rec.log("motor/transform", &rerun::Transform3D::from(motor))?;
        }

        // =====================================================================
        // Animation 3: Helical motion (screw motion via composed motors)
        // =====================================================================
        {
            // Create a helical path: rotation around z combined with z translation
            let angle = t * TAU * 2.0; // Two full rotations
            let z_offset = t * 4.0; // Rise 4 units

            let motor = Motor::from_rotation_z(angle)
                .compose(&Motor::from_translation(-5.0, 0.0, z_offset));

            // Show coordinate frame following the helix
            let origin = Point::origin();
            let tx = motor.transform_point(&Point::from_cartesian(0.5, 0.0, 0.0));
            let ty = motor.transform_point(&Point::from_cartesian(0.0, 0.5, 0.0));
            let tz = motor.transform_point(&Point::from_cartesian(0.0, 0.0, 0.5));
            let o = motor.transform_point(&origin);

            rec.log(
                "helix/frame_x",
                &rerun::Arrows3D::from_vectors([[tx.x() - o.x(), tx.y() - o.y(), tx.z() - o.z()]])
                    .with_origins([[o.x(), o.y(), o.z()]])
                    .with_colors([[255, 0, 0]]),
            )?;
            rec.log(
                "helix/frame_y",
                &rerun::Arrows3D::from_vectors([[ty.x() - o.x(), ty.y() - o.y(), ty.z() - o.z()]])
                    .with_origins([[o.x(), o.y(), o.z()]])
                    .with_colors([[0, 255, 0]]),
            )?;
            rec.log(
                "helix/frame_z",
                &rerun::Arrows3D::from_vectors([[tz.x() - o.x(), tz.y() - o.y(), tz.z() - o.z()]])
                    .with_origins([[o.x(), o.y(), o.z()]])
                    .with_colors([[0, 0, 255]]),
            )?;

            // Trail point
            rec.log(
                "helix/trail",
                &rerun::Points3D::new([[o.x(), o.y(), o.z()]])
                    .with_colors([[100, 100, 255, 128]])
                    .with_radii([0.03]),
            )?;
        }
    }

    info!("Animation complete! Use the timeline scrubber in Rerun to replay.");

    Ok(())
}
