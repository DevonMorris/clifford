//! 2D Robot Arm Visualization with PGA Motors
//!
//! This demo demonstrates forward kinematics of a 2-link robot arm using
//! 2D Projective Geometric Algebra (PGA). Each joint rotation is represented
//! by a Motor, and positions are computed by sequential transformation.
//!
//! Run with: `cargo run -p clifford-viz --example projective2_robot --release`
//!
//! ## Mathematical Background
//!
//! In 2D PGA, Motors encode rigid transformations:
//! - **Rotation** around origin: `M = cos(θ/2) + sin(θ/2)·e₃`
//! - **Translation**: `M = 1 + (dx·e₁ + dy·e₂)/2`
//!
//! For a robot arm, we use motors to rotate link endpoints around joints,
//! computing positions through the kinematic chain.

use std::collections::VecDeque;
use std::f32::consts::{PI, TAU};

use clifford::ops::Transform;
use clifford::specialized::projective::dim2::{Motor, Point};
use clifford_viz::common::prelude::*;
use egui_plot::{Plot, Points};

/// Fixed viewport bounds for stable viewing experience.
const VIEWPORT_BOUNDS: f64 = 5.0;

/// Maximum number of trail points to keep.
const MAX_TRAIL_LENGTH: usize = 200;

/// Animation mode for the robot arm.
#[derive(Clone, Copy, PartialEq, Default)]
enum AnimationMode {
    /// Manual control via sliders.
    #[default]
    Manual,
    /// Trace a circle at the workspace edge.
    CircleTrace,
    /// Figure-8 pattern using Lissajous curves.
    Figure8,
}

/// Demo state for 2D Robot Arm visualization.
struct RobotArmDemo {
    /// Length of first link (base to elbow).
    link1_length: f64,
    /// Length of second link (elbow to end effector).
    link2_length: f64,
    /// Base joint angle (radians).
    theta1: f32,
    /// Elbow joint angle (radians).
    theta2: f32,
    /// Animation state.
    animation: Animation,
    /// Current animation mode.
    animation_mode: AnimationMode,
    /// Whether to show the coordinate grid.
    show_grid: bool,
    /// Whether to show workspace boundary circles.
    show_workspace: bool,
    /// Whether to show joint angle arcs.
    show_angles: bool,
    /// Whether to show end effector trail.
    show_trail: bool,
    /// Trail of end effector positions.
    trail: VecDeque<(f64, f64)>,
}

impl Default for RobotArmDemo {
    fn default() -> Self {
        Self {
            link1_length: 2.0,
            link2_length: 1.5,
            theta1: 0.5,  // ~29 degrees
            theta2: -0.3, // ~-17 degrees
            animation: Animation::with_duration(4.0),
            animation_mode: AnimationMode::Manual,
            show_grid: true,
            show_workspace: true,
            show_angles: true,
            show_trail: false,
            trail: VecDeque::new(),
        }
    }
}

impl RobotArmDemo {
    /// Computes all joint positions using forward kinematics.
    ///
    /// Returns (base, elbow, end_effector) as Points.
    ///
    /// # Forward Kinematics Algorithm
    ///
    /// 1. Base is at the origin
    /// 2. Elbow = rotate (link1_length, 0) around origin by θ₁
    /// 3. End effector = elbow + rotate (link2_length, 0) by (θ₁ + θ₂)
    fn forward_kinematics(&self) -> (Point<f64>, Point<f64>, Point<f64>) {
        let base = Point::origin();

        // Elbow: rotate link1 endpoint around origin by theta1
        let link1_local = Point::from_cartesian(self.link1_length, 0.0);
        let r1 = Motor::from_rotation(f64::from(self.theta1));
        let elbow = r1.transform(&link1_local);

        // End effector: rotate link2 by total angle, then translate to elbow
        let link2_local = Point::from_cartesian(self.link2_length, 0.0);
        let r_total = Motor::from_rotation(f64::from(self.theta1 + self.theta2));
        let link2_rotated = r_total.transform(&link2_local);

        // Add elbow position to get end effector in world coordinates
        let (ex, ey) = elbow.to_cartesian().unwrap();
        let (dx, dy) = link2_rotated.to_cartesian().unwrap();
        let end_effector = Point::from_cartesian(ex + dx, ey + dy);

        (base, elbow, end_effector)
    }

    /// Renders the robot arm links.
    fn render_arm(
        &self,
        plot_ui: &mut egui_plot::PlotUi,
        base: (f64, f64),
        elbow: (f64, f64),
        end: (f64, f64),
    ) {
        // Link 1: Base to Elbow (thicker, primary color)
        plot_ui.line(
            labeled_line_segment(base.0, base.1, elbow.0, elbow.1, palette::LINE, "Link 1")
                .width(4.0),
        );

        // Link 2: Elbow to End Effector (slightly thinner, secondary color)
        plot_ui.line(
            labeled_line_segment(elbow.0, elbow.1, end.0, end.1, palette::MOTOR, "Link 2")
                .width(3.5),
        );
    }

    /// Renders the joint markers and end effector.
    fn render_joints(
        &self,
        plot_ui: &mut egui_plot::PlotUi,
        base: (f64, f64),
        elbow: (f64, f64),
        end: (f64, f64),
    ) {
        // Base joint (fixed pivot) - larger
        plot_ui.line(circle_2d(base.0, base.1, 0.15, palette::POINT, 32));
        plot_ui.points(
            Points::new(vec![[base.0, base.1]])
                .color(palette::POINT)
                .radius(8.0)
                .filled(true)
                .name("Base Joint"),
        );

        // Elbow joint
        plot_ui.line(circle_2d(elbow.0, elbow.1, 0.12, palette::ROTOR, 32));
        plot_ui.points(
            Points::new(vec![[elbow.0, elbow.1]])
                .color(palette::ROTOR)
                .radius(6.0)
                .filled(true)
                .name("Elbow Joint"),
        );

        // End effector (highlighted)
        plot_ui.points(
            Points::new(vec![[end.0, end.1]])
                .color(palette::SELECTED)
                .radius(10.0)
                .filled(true)
                .name("End Effector"),
        );
    }

    /// Renders the workspace boundary (reachable area).
    fn render_workspace(&self, plot_ui: &mut egui_plot::PlotUi) {
        let r_max = self.link1_length + self.link2_length; // Fully extended
        let r_min = (self.link1_length - self.link2_length).abs(); // Fully folded

        // Outer reachable boundary
        plot_ui.line(
            circle_2d(0.0, 0.0, r_max, with_alpha(palette::CIRCLE, 80), 64).name("Outer workspace"),
        );

        // Inner unreachable boundary (if links differ in length)
        if r_min > 0.01 {
            plot_ui.line(
                circle_2d(0.0, 0.0, r_min, with_alpha(palette::LINE_SECONDARY, 60), 64)
                    .name("Inner boundary"),
            );
        }
    }

    /// Renders joint angle arcs.
    fn render_angle_arcs(
        &self,
        plot_ui: &mut egui_plot::PlotUi,
        base: (f64, f64),
        elbow: (f64, f64),
    ) {
        // Base joint angle arc (from x-axis)
        if self.theta1.abs() > 0.05 {
            let arc_radius = 0.5;
            for arc in arc_with_arrow(
                base.0,
                base.1,
                arc_radius,
                0.0,
                f64::from(self.theta1),
                with_alpha(palette::ROTOR, 180),
            ) {
                plot_ui.line(arc);
            }
        }

        // Elbow joint angle arc (relative to link1 direction)
        if self.theta2.abs() > 0.05 {
            let arc_radius = 0.4;
            for arc in arc_with_arrow(
                elbow.0,
                elbow.1,
                arc_radius,
                f64::from(self.theta1), // Start from link1 direction
                f64::from(self.theta1 + self.theta2), // End at link2 direction
                with_alpha(palette::MOTOR, 180),
            ) {
                plot_ui.line(arc);
            }
        }
    }

    /// Renders the end effector trail.
    fn render_trail(&self, plot_ui: &mut egui_plot::PlotUi) {
        if self.trail.len() < 2 {
            return;
        }

        // Draw trail as connected line segments with fading opacity
        let points: Vec<[f64; 2]> = self.trail.iter().map(|(x, y)| [*x, *y]).collect();

        plot_ui.line(
            egui_plot::Line::new(egui_plot::PlotPoints::new(points))
                .color(with_alpha(palette::SELECTED, 100))
                .width(1.5)
                .name("Trail"),
        );
    }

    /// Updates animation and trail.
    fn update_animation(&mut self) {
        if self.animation_mode != AnimationMode::Manual && self.animation.playing {
            let t = self.animation.progress();

            match self.animation_mode {
                AnimationMode::CircleTrace => {
                    // Trace a circle - rotate base, keep elbow straight
                    self.theta1 = t * TAU;
                    self.theta2 = 0.0;
                }
                AnimationMode::Figure8 => {
                    // Figure-8 using Lissajous-like motion
                    self.theta1 = (t * TAU).sin() * PI / 2.0;
                    self.theta2 = (t * TAU * 2.0).sin() * PI / 3.0;
                }
                AnimationMode::Manual => {}
            }
        }
    }

    /// Updates the end effector trail.
    fn update_trail(&mut self) {
        if self.show_trail {
            let (_, _, end) = self.forward_kinematics();
            if let Some((x, y)) = end.to_cartesian() {
                self.trail.push_back((x, y));
                while self.trail.len() > MAX_TRAIL_LENGTH {
                    self.trail.pop_front();
                }
            }
        } else {
            self.trail.clear();
        }
    }
}

impl VisualizationApp for RobotArmDemo {
    fn name(&self) -> &'static str {
        "2D Robot Arm - Motor Kinematics"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        self.update_animation();
        self.update_trail();
    }

    fn render(&mut self, ui: &mut egui::Ui) {
        let (_base_pt, elbow_pt, end_pt) = self.forward_kinematics();

        // Convert to Cartesian for rendering
        let base = (0.0, 0.0);
        let elbow = elbow_pt.to_cartesian().unwrap();
        let end = end_pt.to_cartesian().unwrap();

        Plot::new("robot_arm_plot")
            .data_aspect(1.0)
            .show_axes(false)
            .show_grid(false)
            .allow_zoom(true)
            .allow_drag(true)
            .show(ui, |plot_ui| {
                // Grid (using fixed bounds for stability)
                if self.show_grid {
                    for line in grid_2d(VIEWPORT_BOUNDS, 1.0) {
                        plot_ui.line(line);
                    }
                    for axis in axes_2d(VIEWPORT_BOUNDS) {
                        plot_ui.line(axis);
                    }
                }

                // Workspace boundary
                if self.show_workspace {
                    self.render_workspace(plot_ui);
                }

                // End effector trail
                if self.show_trail {
                    self.render_trail(plot_ui);
                }

                // Robot arm links
                self.render_arm(plot_ui, base, elbow, end);

                // Joint angle arcs
                if self.show_angles {
                    self.render_angle_arcs(plot_ui, base, elbow);
                }

                // Joints and end effector
                self.render_joints(plot_ui, base, elbow, end);

                // Base mount (ground symbol)
                let mount_width = 0.4;
                plot_ui.line(line_segment(
                    -mount_width,
                    -0.15,
                    mount_width,
                    -0.15,
                    palette::GRID_MAJOR,
                ));
                // Ground hatch marks
                for i in 0..5 {
                    let x = -mount_width + (i as f64 + 0.5) * mount_width / 2.5;
                    plot_ui.line(line_segment(x, -0.15, x - 0.1, -0.25, palette::GRID_MAJOR));
                }
            });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Joint Angles ===
        group_header(ui, "Joint Angles");

        angle_slider_range(ui, "Base \u{03b8}\u{2081}", &mut self.theta1, -180.0, 180.0);
        angle_slider_range(
            ui,
            "Elbow \u{03b8}\u{2082}",
            &mut self.theta2,
            -180.0,
            180.0,
        );

        if ui.button("Reset Angles").clicked() {
            self.theta1 = 0.5;
            self.theta2 = -0.3;
        }

        // === Link Configuration ===
        group_header(ui, "Link Configuration");

        ui.horizontal(|ui| {
            ui.label("Link 1:");
            ui.add(
                egui::DragValue::new(&mut self.link1_length)
                    .speed(0.05)
                    .range(0.5..=4.0)
                    .suffix(" units"),
            );
        });

        ui.horizontal(|ui| {
            ui.label("Link 2:");
            ui.add(
                egui::DragValue::new(&mut self.link2_length)
                    .speed(0.05)
                    .range(0.5..=4.0)
                    .suffix(" units"),
            );
        });

        // === Animation ===
        group_header(ui, "Animation");

        ui.horizontal(|ui| {
            if ui
                .selectable_label(self.animation_mode == AnimationMode::Manual, "Manual")
                .clicked()
            {
                self.animation_mode = AnimationMode::Manual;
                self.animation.playing = false;
            }
            if ui
                .selectable_label(self.animation_mode == AnimationMode::CircleTrace, "Circle")
                .clicked()
            {
                self.animation_mode = AnimationMode::CircleTrace;
                self.animation.playing = true;
            }
            if ui
                .selectable_label(self.animation_mode == AnimationMode::Figure8, "Figure-8")
                .clicked()
            {
                self.animation_mode = AnimationMode::Figure8;
                self.animation.playing = true;
            }
        });

        if self.animation_mode != AnimationMode::Manual {
            animation_controls(ui, &mut self.animation);
            progress_slider(ui, &mut self.animation);
        }

        // === Display Options ===
        group_header(ui, "Display");

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_grid, "Grid");
            ui.checkbox(&mut self.show_workspace, "Workspace");
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_angles, "Angles");
            ui.checkbox(&mut self.show_trail, "Trail");
        });

        if self.show_trail && ui.button("Clear Trail").clicked() {
            self.trail.clear();
        }

        // === End Effector Position ===
        group_header(ui, "End Effector");

        let (_, _, end_pt) = self.forward_kinematics();
        if let Some((x, y)) = end_pt.to_cartesian() {
            ui.horizontal(|ui| {
                ui.label("Position:");
                ui.monospace(format!("({:.3}, {:.3})", x, y));
            });

            let reach = (x * x + y * y).sqrt();
            ui.horizontal(|ui| {
                ui.label("Reach:");
                ui.monospace(format!("{:.3}", reach));
            });
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        let (_, elbow, end) = self.forward_kinematics();
        let elbow_coords = elbow.to_cartesian().unwrap();
        let end_coords = end.to_cartesian().unwrap();

        ui.horizontal(|ui| {
            ui.label(format!("End: ({:.2}, {:.2})", end_coords.0, end_coords.1));
            ui.separator();
            ui.label(format!(
                "Elbow: ({:.2}, {:.2})",
                elbow_coords.0, elbow_coords.1
            ));
            ui.separator();
            ui.label(format!(
                "\u{03b8}\u{2081}={:.1}\u{00b0} \u{03b8}\u{2082}={:.1}\u{00b0}",
                self.theta1.to_degrees(),
                self.theta2.to_degrees()
            ));
            ui.separator();
            let reach = (end_coords.0 * end_coords.0 + end_coords.1 * end_coords.1).sqrt();
            ui.label(format!("Reach: {:.2}", reach));
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(ROBOT_ARM_EDUCATION)
    }
}

/// Educational content for the 2D Robot Arm visualization.
const ROBOT_ARM_EDUCATION: EducationalContent = EducationalContent {
    title: "2D Robot Arm with PGA Motors",

    overview: "\
This visualization demonstrates forward kinematics of a 2-link robot arm \
using 2D Projective Geometric Algebra (PGA). Each joint rotation is \
represented by a Motor, and the end effector position is computed by \
sequentially transforming points through the kinematic chain.

Key insight: In 2D PGA, Motors encode rotation elegantly. We use \
Motor::from_rotation(\u{03b8}) to create a rotation around the origin, \
then apply it via the antisandwich product (transform).",

    math_background: "\
FORWARD KINEMATICS computes positions from joint angles:

Base joint at origin: P\u{2080} = (0, 0)

Elbow position:
    L\u{2081} = (link1_length, 0)  // Link 1 in local coords
    M\u{2081} = Motor::from_rotation(\u{03b8}\u{2081})
    Elbow = M\u{2081}.transform(L\u{2081})

End effector position:
    L\u{2082} = (link2_length, 0)  // Link 2 in local coords
    M\u{2082} = Motor::from_rotation(\u{03b8}\u{2081} + \u{03b8}\u{2082})  // Total rotation
    End = Elbow + M\u{2082}.transform(L\u{2082})

MOTOR FORMULA for rotation around origin:
    M = cos(\u{03b8}/2) - sin(\u{03b8}/2)\u{00b7}e\u{2083}

The antisandwich product M\u{207b}\u{00b9}PM transforms point P.

WORKSPACE is the annular region between:
    r_max = L\u{2081} + L\u{2082}  (fully extended)
    r_min = |L\u{2081} - L\u{2082}|  (fully folded)

Any point in this region is reachable by some (\u{03b8}\u{2081}, \u{03b8}\u{2082}).",

    how_to_use: "\
\u{2022} JOINT SLIDERS: Drag \u{03b8}\u{2081} and \u{03b8}\u{2082} to control the arm
\u{2022} LINK LENGTHS: Adjust to change arm proportions and workspace
\u{2022} ANIMATION: Select Circle or Figure-8 to see motion patterns
\u{2022} WORKSPACE: Enable to see the reachable area boundaries
\u{2022} TRAIL: Enable to visualize the end effector path over time
\u{2022} ANGLES: Show joint angle arcs for visual feedback

Try setting both angles to 0\u{00b0} to see the arm fully extended along the x-axis, \
or set \u{03b8}\u{2082} to \u{00b1}180\u{00b0} to fold the arm back on itself.",

    key_concepts: "\
\u{2022} Motors encode rotation via antisandwich product
\u{2022} Sequential transformation follows kinematic chain
\u{2022} Joint angles are relative: \u{03b8}\u{2082} is relative to link 1
\u{2022} Workspace = annular region of reachable positions
\u{2022} Forward kinematics: (\u{03b8}\u{2081}, \u{03b8}\u{2082}) \u{2192} (x, y)
\u{2022} Inverse kinematics: (x, y) \u{2192} (\u{03b8}\u{2081}, \u{03b8}\u{2082}) (not covered)",

    resources: &[
        (
            "Rigid Geometric Algebra Wiki - Motors",
            "https://rigidgeometricalgebra.org/wiki/index.php?title=Motor",
        ),
        (
            "Forward Kinematics",
            "https://en.wikipedia.org/wiki/Forward_kinematics",
        ),
        (
            "Look, Ma, No Matrices!",
            "https://enkimute.github.io/LookMaNoMatrices/",
        ),
    ],
};

fn main() -> eframe::Result<()> {
    run_app::<RobotArmDemo>()
}
