//! 2D Robot Arm - Motor Kinematics Visualization
//!
//! This demo demonstrates forward kinematics of a 2-link robot arm using
//! motor composition in 2D Projective Geometric Algebra.
//!
//! Run with: `cargo run -p clifford-viz --example projective2_robot --release`
//!
//! ## Mathematical Background
//!
//! In 2D PGA, each joint transformation is represented as a motor:
//! - Rotation motors: `M = cos(θ/2) + sin(θ/2)·e₁₂₃`
//! - Translation motors: `M = 1 + (d/2)·direction`
//!
//! The total transformation is the composition of all motors:
//! `M_total = M_base · T_link1 · M_elbow · T_link2`
//!
//! The end effector position is found by transforming the origin:
//! `P_end = M_total⁻¹ · Origin · M_total`

use clifford::ops::{Transform, Versor};
use clifford::specialized::projective::dim2::{Flector, Motor, Point};
use clifford_viz::common::prelude::*;
use egui_plot::{Plot, Points};
use std::f64::consts::PI;

/// State for the 2D robot arm demonstration.
struct Robot2DDemo {
    /// Base joint angle (radians).
    joint1_angle: f32,
    /// Elbow joint angle (radians).
    joint2_angle: f32,
    /// Length of the first link.
    link1_length: f32,
    /// Length of the second link.
    link2_length: f32,
    /// Animation state for joint 1.
    animation1: Animation,
    /// Animation state for joint 2.
    animation2: Animation,
    /// Whether to animate joint 1.
    animate_joint1: bool,
    /// Whether to animate joint 2.
    animate_joint2: bool,
    /// Whether to show the workspace boundary.
    show_workspace: bool,
    /// Whether to show intermediate frames.
    show_frames: bool,
    /// Whether to show the grid.
    show_grid: bool,
    /// Trail of end effector positions.
    trail: Vec<(f64, f64)>,
    /// Whether to show the trail.
    show_trail: bool,
    /// Maximum trail length.
    max_trail_length: usize,
}

impl Default for Robot2DDemo {
    fn default() -> Self {
        Self {
            joint1_angle: 0.5,
            joint2_angle: -0.8,
            link1_length: 2.0,
            link2_length: 1.5,
            animation1: Animation::with_duration(4.0),
            animation2: Animation::with_duration(3.0),
            animate_joint1: false,
            animate_joint2: false,
            show_workspace: true,
            show_frames: true,
            show_grid: true,
            trail: Vec::new(),
            show_trail: true,
            max_trail_length: 200,
        }
    }
}

impl Robot2DDemo {
    /// Compute forward kinematics using motor composition.
    ///
    /// Returns (base_pos, elbow_pos, end_effector_pos).
    ///
    /// In 2D PGA, Motor × Motor = Flector, so we track types carefully.
    fn forward_kinematics(&self) -> (Point<f64>, Point<f64>, Point<f64>) {
        let j1 = f64::from(self.joint1_angle);
        let j2 = f64::from(self.joint2_angle);
        let l1 = f64::from(self.link1_length);
        let l2 = f64::from(self.link2_length);

        // Base is at origin
        let base = Point::origin();

        // Motor for base rotation
        let m1 = Motor::from_rotation(j1);

        // Motor for translation along link1 (in local frame, along x-axis)
        let t1 = Motor::from_translation(l1, 0.0);

        // Motor for elbow rotation
        let m2 = Motor::from_rotation(j2);

        // Motor for translation along link2 (in local frame, along x-axis)
        let t2 = Motor::from_translation(l2, 0.0);

        // Compose motors: M_total = M1 · T1 · M2 · T2
        // Note: In 2D PGA, Motor × Motor = Flector, Flector × Flector = Flector

        // Position after link1 (elbow position)
        // m1 · t1 → Flector
        let m_to_elbow: Flector<f64> = m1.compose(&t1);
        let elbow = m_to_elbow.transform(&Point::origin());

        // Position after link2 (end effector)
        // m2 · t2 → Flector
        // Flector · Flector → Flector
        let m2_t2: Flector<f64> = m2.compose(&t2);
        let m_total: Flector<f64> = m_to_elbow.compose(&m2_t2);
        let end_effector = m_total.transform(&Point::origin());

        (base, elbow, end_effector)
    }

    /// Get the transformation Flectors for frame visualization.
    ///
    /// Returns a list of transformations at each stage of the kinematic chain.
    /// Note: Due to 2D PGA algebra, some compositions result in Flectors.
    fn get_frame_transforms(&self) -> Vec<FrameTransform> {
        let j1 = f64::from(self.joint1_angle);
        let j2 = f64::from(self.joint2_angle);
        let l1 = f64::from(self.link1_length);

        let m1 = Motor::from_rotation(j1);
        let t1 = Motor::from_translation(l1, 0.0);
        let m2 = Motor::from_rotation(j2);

        // Build transformation chain
        let at_elbow: Flector<f64> = m1.compose(&t1);
        let after_elbow: Motor<f64> = at_elbow.compose(&m2);

        vec![
            FrameTransform::Identity,
            FrameTransform::Motor(m1),
            FrameTransform::Flector(at_elbow),
            FrameTransform::Motor(after_elbow),
        ]
    }
}

/// A transformation that can be either a Motor or Flector.
///
/// In 2D PGA, Motor × Motor = Flector, so kinematic chains
/// alternate between these types.
enum FrameTransform {
    /// Identity transformation.
    Identity,
    /// Motor transformation (odd subalgebra).
    Motor(Motor<f64>),
    /// Flector transformation (even subalgebra).
    Flector(Flector<f64>),
}

impl FrameTransform {
    /// Transform a point using this transformation.
    fn transform_point(&self, p: &Point<f64>) -> Point<f64> {
        match self {
            FrameTransform::Identity => *p,
            FrameTransform::Motor(m) => m.transform(p),
            FrameTransform::Flector(f) => f.transform(p),
        }
    }
}

impl Robot2DDemo {
    /// Compute the workspace boundary (reachable area).
    fn workspace_boundary(&self, segments: usize) -> Vec<[f64; 2]> {
        let l1 = f64::from(self.link1_length);
        let l2 = f64::from(self.link2_length);

        let r_max = l1 + l2; // Maximum reach
        let r_min = (l1 - l2).abs(); // Minimum reach

        // Outer boundary
        let mut points = Vec::new();
        for i in 0..=segments {
            let angle = 2.0 * PI * i as f64 / segments as f64;
            points.push([r_max * angle.cos(), r_max * angle.sin()]);
        }

        // Inner boundary (if links have different lengths)
        if r_min > 0.01 {
            // Add a gap
            points.push([f64::NAN, f64::NAN]);
            for i in 0..=segments {
                let angle = 2.0 * PI * i as f64 / segments as f64;
                points.push([r_min * angle.cos(), r_min * angle.sin()]);
            }
        }

        points
    }

    /// Update the end effector trail.
    fn update_trail(&mut self, end_effector: &Point<f64>) {
        if let Some((x, y)) = end_effector.to_cartesian() {
            self.trail.push((x, y));
            if self.trail.len() > self.max_trail_length {
                self.trail.remove(0);
            }
        }
    }
}

impl VisualizationApp for Robot2DDemo {
    fn name(&self) -> &'static str {
        "2D Robot Arm - Motor Kinematics"
    }

    fn update(&mut self, dt: f32) {
        // Update animations
        self.animation1.update(dt);
        self.animation2.update(dt);

        if self.animate_joint1 && self.animation1.playing {
            // Oscillate between -π/2 and π/2
            let t = self.animation1.progress();
            self.joint1_angle = (PI as f32) * (t * 2.0 - 1.0) * 0.5;
        }

        if self.animate_joint2 && self.animation2.playing {
            // Oscillate between -π and π
            let t = self.animation2.progress();
            self.joint2_angle = (PI as f32) * (t * 2.0 - 1.0) * 0.75;
        }

        // Update trail if animating
        if (self.animate_joint1 && self.animation1.playing)
            || (self.animate_joint2 && self.animation2.playing)
        {
            let (_, _, end_effector) = self.forward_kinematics();
            self.update_trail(&end_effector);
        }
    }

    fn render(&self, ui: &mut egui::Ui) {
        let (base, elbow, end_effector) = self.forward_kinematics();
        let l1 = f64::from(self.link1_length);
        let l2 = f64::from(self.link2_length);
        let bounds = (l1 + l2 + 1.0).max(5.0);

        Plot::new("robot2d_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .show(ui, |plot_ui| {
                // Draw grid
                if self.show_grid {
                    for line in grid_2d(bounds, 1.0) {
                        plot_ui.line(line);
                    }
                    for axis in axes_2d(bounds) {
                        plot_ui.line(axis);
                    }
                }

                // Draw workspace boundary
                if self.show_workspace {
                    let boundary = self.workspace_boundary(64);
                    // Filter out NaN segments
                    let mut current_segment: Vec<[f64; 2]> = Vec::new();
                    for point in &boundary {
                        if point[0].is_nan() {
                            // End current segment, start new one
                            if current_segment.len() >= 2 {
                                plot_ui.line(
                                    egui_plot::Line::new(egui_plot::PlotPoints::new(
                                        current_segment.clone(),
                                    ))
                                    .color(with_alpha(palette::GRID, 60))
                                    .width(1.0),
                                );
                            }
                            current_segment.clear();
                        } else {
                            current_segment.push(*point);
                        }
                    }
                    // Draw final segment
                    if current_segment.len() >= 2 {
                        plot_ui.line(
                            egui_plot::Line::new(egui_plot::PlotPoints::new(current_segment))
                                .color(with_alpha(palette::GRID, 60))
                                .width(1.0)
                                .name("Workspace"),
                        );
                    }
                }

                // Draw trail
                if self.show_trail && !self.trail.is_empty() {
                    let trail_points: Vec<[f64; 2]> =
                        self.trail.iter().map(|(x, y)| [*x, *y]).collect();
                    plot_ui.line(
                        egui_plot::Line::new(egui_plot::PlotPoints::new(trail_points))
                            .color(with_alpha(palette::PLANE, 100))
                            .width(1.5)
                            .name("Trail"),
                    );
                }

                // Draw coordinate frames at each joint
                if self.show_frames {
                    let transforms = self.get_frame_transforms();
                    let frame_size = 0.4;

                    for (i, transform) in transforms.iter().enumerate() {
                        // Get frame origin
                        let origin = transform.transform_point(&Point::origin());
                        if let Some((ox, oy)) = origin.to_cartesian() {
                            // X-axis of frame (transform unit x vector direction)
                            let x_dir = transform.transform_point(&Point::from_cartesian(1.0, 0.0));
                            let y_dir = transform.transform_point(&Point::from_cartesian(0.0, 1.0));

                            if let (Some((xx, xy)), Some((yx, yy))) =
                                (x_dir.to_cartesian(), y_dir.to_cartesian())
                            {
                                let dx = (xx - ox) * frame_size;
                                let dy = (xy - oy) * frame_size;
                                let dyx = (yx - ox) * frame_size;
                                let dyy = (yy - oy) * frame_size;

                                // X-axis (red)
                                for arrow in
                                    arrow_2d(ox, oy, dx, dy, with_alpha(palette::X_AXIS, 150))
                                {
                                    plot_ui.line(arrow);
                                }

                                // Y-axis (green)
                                for arrow in
                                    arrow_2d(ox, oy, dyx, dyy, with_alpha(palette::Y_AXIS, 150))
                                {
                                    plot_ui.line(arrow);
                                }

                                // Label
                                plot_ui.points(
                                    Points::new(vec![[ox, oy]])
                                        .color(with_alpha(palette::GRID, 100))
                                        .radius(3.0)
                                        .name(format!("Frame {}", i)),
                                );
                            }
                        }
                    }
                }

                // Get joint positions
                let (bx, by) = base.to_cartesian().unwrap_or((0.0, 0.0));
                let (ex, ey) = elbow.to_cartesian().unwrap_or((0.0, 0.0));
                let (px, py) = end_effector.to_cartesian().unwrap_or((0.0, 0.0));

                // Draw links
                // Link 1: base to elbow
                plot_ui
                    .line(labeled_line_segment(bx, by, ex, ey, palette::LINE, "Link 1").width(4.0));

                // Link 2: elbow to end effector
                plot_ui.line(
                    labeled_line_segment(ex, ey, px, py, palette::MOTOR, "Link 2").width(4.0),
                );

                // Draw joints
                // Base joint (fixed to ground)
                plot_ui.line(line_segment(
                    bx - 0.3,
                    by - 0.1,
                    bx + 0.3,
                    by - 0.1,
                    palette::GRID,
                ));
                plot_ui.line(line_segment(
                    bx - 0.2,
                    by - 0.2,
                    bx + 0.2,
                    by - 0.2,
                    palette::GRID,
                ));
                plot_ui.points(
                    Points::new(vec![[bx, by]])
                        .color(palette::ROTOR)
                        .radius(10.0)
                        .filled(true)
                        .name("Base"),
                );

                // Elbow joint
                plot_ui.points(
                    Points::new(vec![[ex, ey]])
                        .color(palette::ROTOR)
                        .radius(8.0)
                        .filled(true)
                        .name("Elbow"),
                );

                // End effector
                plot_ui.points(
                    Points::new(vec![[px, py]])
                        .color(palette::PLANE)
                        .radius(10.0)
                        .filled(true)
                        .name("End Effector"),
                );

                // Draw joint angle arcs
                let arc_radius = 0.5;

                // Base joint arc
                if self.joint1_angle.abs() > 0.05 {
                    for arc in arc_with_arrow(
                        bx,
                        by,
                        arc_radius,
                        0.0,
                        f64::from(self.joint1_angle),
                        with_alpha(palette::ROTOR, 180),
                    ) {
                        plot_ui.line(arc);
                    }
                }

                // Elbow joint arc (relative to link1 direction)
                if self.joint2_angle.abs() > 0.05 {
                    let link1_angle = f64::from(self.joint1_angle);
                    for arc in arc_with_arrow(
                        ex,
                        ey,
                        arc_radius * 0.8,
                        link1_angle,
                        link1_angle + f64::from(self.joint2_angle),
                        with_alpha(palette::MOTOR, 180),
                    ) {
                        plot_ui.line(arc);
                    }
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // Joint controls
        section_separator(ui, Some("Joint Angles"));

        // Joint 1
        angle_slider_range(ui, "Joint 1 (base)", &mut self.joint1_angle, -180.0, 180.0);
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.animate_joint1, "Animate");
            if self.animate_joint1
                && ui
                    .button(if self.animation1.playing {
                        "⏸"
                    } else {
                        "▶"
                    })
                    .clicked()
            {
                self.animation1.playing = !self.animation1.playing;
            }
        });

        ui.add_space(4.0);

        // Joint 2
        angle_slider_range(ui, "Joint 2 (elbow)", &mut self.joint2_angle, -180.0, 180.0);
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.animate_joint2, "Animate");
            if self.animate_joint2
                && ui
                    .button(if self.animation2.playing {
                        "⏸"
                    } else {
                        "▶"
                    })
                    .clicked()
            {
                self.animation2.playing = !self.animation2.playing;
            }
        });

        // Link lengths
        section_separator(ui, Some("Link Lengths"));

        ui.horizontal(|ui| {
            ui.label("Link 1:");
            ui.add(egui::Slider::new(&mut self.link1_length, 0.5..=4.0).fixed_decimals(1));
        });

        ui.horizontal(|ui| {
            ui.label("Link 2:");
            ui.add(egui::Slider::new(&mut self.link2_length, 0.5..=3.0).fixed_decimals(1));
        });

        // End effector position
        section_separator(ui, Some("End Effector"));

        let (_, _, end_effector) = self.forward_kinematics();
        if let Some((x, y)) = end_effector.to_cartesian() {
            point2_display(ui, "Position", x as f32, y as f32);

            let dist = (x * x + y * y).sqrt();
            value_display(ui, "Distance", dist as f32, 2);
        }

        // Motor composition display
        section_separator(ui, Some("Motor Composition"));

        let j1 = f64::from(self.joint1_angle);
        let j2 = f64::from(self.joint2_angle);
        let l1 = f64::from(self.link1_length);
        let l2 = f64::from(self.link2_length);

        let m1 = Motor::from_rotation(j1);
        let t1 = Motor::from_translation(l1, 0.0);
        let m2 = Motor::from_rotation(j2);
        let t2 = Motor::from_translation(l2, 0.0);
        // Motor × Motor = Flector, Flector × Flector = Flector
        let m_total: Flector<f64> = m1.compose(&t1).compose(&m2.compose(&t2));

        info_box(
            ui,
            &format!(
                "M_total = R₁ · T₁ · R₂ · T₂\n\nR₁ = rot({:.1}°)\nT₁ = trans({:.1}, 0)\nR₂ = rot({:.1}°)\nT₂ = trans({:.1}, 0)\n\nNote: Motor × Motor = Flector",
                self.joint1_angle.to_degrees(),
                self.link1_length,
                self.joint2_angle.to_degrees(),
                self.link2_length
            ),
        );

        // Flector has grades 0 and 2 (even subalgebra)
        ga_value_display(
            ui,
            "F_total",
            &[
                ("1", m_total.s() as f32),
                ("d", m_total.d() as f32),
                ("nx", m_total.nx() as f32),
                ("ny", m_total.ny() as f32),
            ],
        );

        // Display options
        section_separator(ui, Some("Display Options"));
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_workspace, "Show workspace boundary");
        ui.checkbox(&mut self.show_frames, "Show coordinate frames");
        ui.checkbox(&mut self.show_trail, "Show end effector trail");

        if self.show_trail {
            if ui.button("Clear Trail").clicked() {
                self.trail.clear();
            }
        }

        // Preset poses
        section_separator(ui, Some("Preset Poses"));

        ui.horizontal(|ui| {
            if ui.button("Home").clicked() {
                self.joint1_angle = 0.0;
                self.joint2_angle = 0.0;
                self.trail.clear();
            }
            if ui.button("Reach Up").clicked() {
                self.joint1_angle = PI as f32 / 2.0;
                self.joint2_angle = 0.0;
                self.trail.clear();
            }
            if ui.button("Folded").clicked() {
                self.joint1_angle = PI as f32 / 4.0;
                self.joint2_angle = -PI as f32 / 2.0;
                self.trail.clear();
            }
        });
    }

    fn info(&self, ui: &mut egui::Ui) {
        let (_, _, end_effector) = self.forward_kinematics();
        let l1 = f64::from(self.link1_length);
        let l2 = f64::from(self.link2_length);

        ui.horizontal(|ui| {
            ui.label(format!(
                "θ₁={:.1}° θ₂={:.1}°",
                self.joint1_angle.to_degrees(),
                self.joint2_angle.to_degrees()
            ));
            ui.separator();
            if let Some((x, y)) = end_effector.to_cartesian() {
                ui.label(format!("End: ({:.2}, {:.2})", x, y));
            }
            ui.separator();
            ui.label(format!("Reach: {:.1} - {:.1}", (l1 - l2).abs(), l1 + l2));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(ROBOT2D_EDUCATION)
    }
}

/// Educational content for the 2D robot arm visualization.
const ROBOT2D_EDUCATION: EducationalContent = EducationalContent {
    title: "Robot Kinematics with Motors",

    overview: "\
This demo shows how 2D Projective GA motors elegantly solve forward kinematics \
for a 2-link robot arm.

Each joint transformation (rotation or translation) is a motor. The total \
transformation from base to end effector is simply the COMPOSITION of all \
motors along the kinematic chain.

This is cleaner than traditional approaches using rotation matrices and \
homogeneous transformations because:
• Motors compose naturally: M_total = M₁ · M₂ · ... · Mₙ
• No singularities (unlike Euler angles)
• Rotation and translation unified in one algebraic element",

    math_background: "\
MOTOR KINEMATICS

A 2-link robot arm has 4 transformations:
1. R₁ - Base rotation (joint 1)
2. T₁ - Translation along link 1
3. R₂ - Elbow rotation (joint 2)
4. T₂ - Translation along link 2

Each is a motor in 2D PGA (grades 1 and 3):
    Rotation: R = cos(θ/2)·e₁₂₃ - sin(θ/2)·e₃
    Translation: T = e₁₂₃ + (d/2)·direction

The total motor is their composition:
    M_total = R₁ · T₁ · R₂ · T₂

Note: Composition order matters! We apply transformations right-to-left \
(like matrix multiplication).

FORWARD KINEMATICS

The end effector position is found by transforming the origin:
    P_end = M_total⁻¹ · Origin · M_total

This antisandwich product propagates the point through all the \
transformations in the chain.",

    how_to_use: "\
• JOINT SLIDERS: Adjust θ₁ (base) and θ₂ (elbow) angles
• ANIMATE: Check 'Animate' and press play to see continuous motion
• LINK LENGTHS: Adjust arm segment lengths
• WORKSPACE: Gray circles show reachable area
• FRAMES: Coordinate frames at each joint show the local orientation
• TRAIL: See the path traced by the end effector
• PRESETS: Quick poses for common configurations",

    key_concepts: "\
• Motors encode rigid transformations (rotation + translation)
• Motor composition gives the total transformation
• Antisandwich product transforms points: P' = M⁻¹PM
• Forward kinematics: compose motors base → end effector
• Workspace = all reachable end effector positions
• Joint limits would constrain valid angle ranges (not shown)",

    resources: &[
        (
            "Robot Kinematics with GA",
            "https://rigidgeometricalgebra.org/wiki/",
        ),
        (
            "Steven De Keninck - Kinematics",
            "https://bivector.net/PGA4CS.pdf",
        ),
        (
            "Geometric Algebra for Computer Science",
            "https://geometricalgebra.org/",
        ),
    ],
};

fn main() -> eframe::Result<()> {
    run_app::<Robot2DDemo>()
}
