//! Bivector and Hodge dual visualization.
//!
//! This example demonstrates the relationship between bivectors (oriented areas)
//! and their Hodge dual vectors in 3D Euclidean space.
//!
//! # Running
//!
//! ```bash
//! cargo run --example rerun_bivector --features rerun-0_28
//! ```
//!
//! # Geometric Algebra Concepts
//!
//! In 3D, a bivector represents an oriented plane segment (like an oriented area).
//! The Hodge dual of a bivector is a vector perpendicular to that plane.
//!
//! The bivector u ∧ v has:
//! - Magnitude = area of parallelogram spanned by u and v
//! - Zero when u and v are parallel (no area!)
//! - Maximum when u and v are perpendicular
//!
//! The Hodge dual vector is perpendicular to the plane containing u and v.

use std::f32::consts::TAU;

use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
use clifford::specialized::visualization::rerun;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let rec = rerun::RecordingStreamBuilder::new("clifford_bivector").spawn()?;

    let num_frames = 300;

    for frame in 0..num_frames {
        let t = frame as f32 / num_frames as f32;
        rec.set_time_sequence("frame", frame as i64);

        // Fixed vector u along x-axis
        let u = Vector::new(2.0_f32, 0.0, 0.0);

        // Rotating vector v - starts along y, rotates in xy-plane then tilts up
        let angle_xy = t * TAU * 2.0; // Two full rotations in xy
        let angle_tilt = (t * TAU).sin() * 0.5; // Gentle tilt up and down

        let rotor_xy = Rotor::from_angle_plane(angle_xy, Bivector::unit_xy());
        let rotor_tilt = Rotor::from_angle_plane(angle_tilt, Bivector::unit_xz());

        let base_v = Vector::new(0.0, 1.5_f32, 0.0);
        let v = rotor_tilt.rotate(rotor_xy.rotate(base_v));

        // The bivector is u ∧ v (wedge product)
        let bivector = u.wedge(v);

        // The Hodge dual is a vector perpendicular to the plane
        let dual = bivector.dual();

        // =====================================================================
        // Draw the parallelogram (bivector visualization)
        // =====================================================================

        // Four corners of the parallelogram centered at origin
        let p0 = Vector::new(0.0, 0.0, 0.0); // origin
        let p1 = u; // origin + u
        let p2 = u + v; // origin + u + v
        let p3 = v; // origin + v

        // Draw the parallelogram edges
        rec.log(
            "bivector/edges",
            &rerun::LineStrips3D::new([[
                [p0.x(), p0.y(), p0.z()],
                [p1.x(), p1.y(), p1.z()],
                [p2.x(), p2.y(), p2.z()],
                [p3.x(), p3.y(), p3.z()],
                [p0.x(), p0.y(), p0.z()], // close the loop
            ]])
            .with_colors([[255, 200, 100]])
            .with_radii([0.02]),
        )?;

        // Draw the filled parallelogram as two triangles
        // Triangle 1: p0, p1, p2
        // Triangle 2: p0, p2, p3
        rec.log(
            "bivector/surface",
            &rerun::Mesh3D::new([
                [p0.x(), p0.y(), p0.z()],
                [p1.x(), p1.y(), p1.z()],
                [p2.x(), p2.y(), p2.z()],
                [p3.x(), p3.y(), p3.z()],
            ])
            .with_triangle_indices([[0_u32, 1, 2], [0, 2, 3]])
            .with_vertex_colors([
                [255, 150, 50, 180],
                [255, 150, 50, 180],
                [255, 150, 50, 180],
                [255, 150, 50, 180],
            ]),
        )?;

        // Draw the spanning vectors u and v
        // u is fixed (red), v is rotating (green)
        rec.log(
            "bivector/vector_u_fixed",
            &rerun::Arrows3D::from_vectors([[u.x(), u.y(), u.z()]])
                .with_origins([[0.0, 0.0, 0.0]])
                .with_colors([[255, 80, 80]]),
        )?;
        rec.log(
            "bivector/vector_v_rotating",
            &rerun::Arrows3D::from_vectors([[v.x(), v.y(), v.z()]])
                .with_origins([[0.0, 0.0, 0.0]])
                .with_colors([[80, 255, 80]]),
        )?;

        // =====================================================================
        // Draw the Hodge dual vector
        // =====================================================================

        // Scale the dual for visibility (normalize and scale)
        let dual_norm = dual.norm();
        let dual_scaled = if dual_norm > 1e-6 {
            let scale = 2.0 / dual_norm;
            Vector::new(dual.x() * scale, dual.y() * scale, dual.z() * scale)
        } else {
            dual
        };

        // Draw from center of parallelogram
        let center = (u + v) * 0.5;
        rec.log(
            "dual/vector",
            &rerun::Arrows3D::from_vectors([[dual_scaled.x(), dual_scaled.y(), dual_scaled.z()]])
                .with_origins([[center.x(), center.y(), center.z()]])
                .with_colors([[100, 150, 255]]),
        )?;

        // Also draw in opposite direction (bivectors are oriented, dual can point either way)
        rec.log(
            "dual/vector_neg",
            &rerun::Arrows3D::from_vectors([[
                -dual_scaled.x(),
                -dual_scaled.y(),
                -dual_scaled.z(),
            ]])
            .with_origins([[center.x(), center.y(), center.z()]])
            .with_colors([[100, 150, 255, 100]]),
        )?;

        // =====================================================================
        // Add labels with bivector components
        // =====================================================================

        rec.log(
            "info/bivector_xy",
            &rerun::TextLog::new(format!("e₁₂: {:.2}", bivector.xy())),
        )?;
        rec.log(
            "info/bivector_xz",
            &rerun::TextLog::new(format!("e₁₃: {:.2}", bivector.xz())),
        )?;
        rec.log(
            "info/bivector_yz",
            &rerun::TextLog::new(format!("e₂₃: {:.2}", bivector.yz())),
        )?;
        rec.log(
            "info/dual",
            &rerun::TextLog::new(format!(
                "dual: ({:.2}, {:.2}, {:.2})",
                dual.x(),
                dual.y(),
                dual.z()
            )),
        )?;

        // =====================================================================
        // Reference coordinate axes
        // =====================================================================

        if frame == 0 {
            rec.log_static(
                "axes/x",
                &rerun::Arrows3D::from_vectors([[3.0_f32, 0.0, 0.0]])
                    .with_origins([[0.0, 0.0, 0.0]])
                    .with_colors([[255, 0, 0, 80]]),
            )?;
            rec.log_static(
                "axes/y",
                &rerun::Arrows3D::from_vectors([[0.0_f32, 3.0, 0.0]])
                    .with_origins([[0.0, 0.0, 0.0]])
                    .with_colors([[0, 255, 0, 80]]),
            )?;
            rec.log_static(
                "axes/z",
                &rerun::Arrows3D::from_vectors([[0.0_f32, 0.0, 3.0]])
                    .with_origins([[0.0, 0.0, 0.0]])
                    .with_colors([[0, 0, 255, 80]]),
            )?;
        }
    }

    println!("Bivector visualization complete!");
    println!();
    println!("Watch how the bivector changes as v rotates around fixed u:");
    println!("- Red arrow (u): fixed vector along x-axis");
    println!("- Green arrow (v): rotating vector");
    println!("- Orange parallelogram: bivector u ∧ v (oriented area)");
    println!("- Blue arrow: Hodge dual (perpendicular to plane)");
    println!();
    println!("Notice: when u and v align, the bivector collapses to zero!");

    Ok(())
}
