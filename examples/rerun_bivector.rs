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
//! - Zero when u ∥ v (parallel = no area!)
//! - Maximum when u ⊥ v (perpendicular)
//! - Sign flips when v crosses u (orientation reverses)
//!
//! The Hodge dual vector is perpendicular to the xy-plane (points along z).

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

        // Rotating vector v in the xy-plane
        let angle = t * TAU; // One full rotation
        let rotor = Rotor::from_angle_plane(angle, Bivector::unit_xy());

        let base_v = Vector::new(0.0, 1.5_f32, 0.0);
        let v = rotor.rotate(base_v);

        // The bivector is u ∧ v (wedge product)
        let bivector = u.wedge(v);

        // The Hodge dual is a vector perpendicular to the plane
        let dual = bivector.dual();

        // Determine bivector sign/orientation for coloring
        // When the bivector flips sign, the orientation reverses
        // We use the z-component of the dual (or xy component of bivector)
        let positive_orientation = dual.z() >= 0.0;
        let (color_primary, color_secondary) = if positive_orientation {
            ([255_u8, 150, 50, 180], [255_u8, 100, 100]) // Orange/red for positive
        } else {
            ([100_u8, 150, 255, 180], [100_u8, 100, 255]) // Blue/purple for negative
        };

        // =====================================================================
        // Draw the parallelogram (bivector visualization)
        // =====================================================================

        // Four corners of the parallelogram centered at origin
        let p0 = Vector::new(0.0, 0.0, 0.0); // origin
        let p1 = u; // origin + u
        let p2 = u + v; // origin + u + v
        let p3 = v; // origin + v

        // Draw the parallelogram edges (color changes with orientation)
        rec.log(
            "bivector/edges",
            &rerun::LineStrips3D::new([[
                [p0.x(), p0.y(), p0.z()],
                [p1.x(), p1.y(), p1.z()],
                [p2.x(), p2.y(), p2.z()],
                [p3.x(), p3.y(), p3.z()],
                [p0.x(), p0.y(), p0.z()], // close the loop
            ]])
            .with_colors([color_secondary])
            .with_radii([0.02]),
        )?;

        // Draw the filled parallelogram as two triangles
        // Triangle 1: p0, p1, p2
        // Triangle 2: p0, p2, p3
        // Color indicates orientation (orange = positive, blue = negative)
        rec.log(
            "bivector/surface",
            &rerun::Mesh3D::new([
                [p0.x(), p0.y(), p0.z()],
                [p1.x(), p1.y(), p1.z()],
                [p2.x(), p2.y(), p2.z()],
                [p3.x(), p3.y(), p3.z()],
            ])
            .with_triangle_indices([[0_u32, 1, 2], [0, 2, 3]])
            .with_vertex_colors([color_primary, color_primary, color_primary, color_primary]),
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
    println!("- Orange surface: positive orientation (u ∧ v > 0)");
    println!("- Blue surface: negative orientation (u ∧ v < 0)");
    println!("- Blue arrow: Hodge dual (perpendicular to plane)");
    println!();
    println!("Notice:");
    println!("  - When u and v align, the bivector collapses to zero!");
    println!("  - When v crosses u, the orientation (color) flips!");

    Ok(())
}
