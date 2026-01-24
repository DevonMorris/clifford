//! 3D Robot Arm Visualization with PGA Motors.
//!
//! This demo demonstrates forward kinematics of a 3-joint robot arm using
//! 3D Projective Geometric Algebra (PGA). Each joint rotation is represented
//! by a Motor, and positions are computed by sequential motor composition.
//!
//! ## Robot Structure
//!
//! - Joint 1 (Base): Rotation around Z-axis (yaw)
//! - Joint 2 (Shoulder): Rotation around Y-axis (pitch)
//! - Joint 3 (Elbow): Rotation around Y-axis (pitch)
//!
//! ## Mathematical Background
//!
//! Forward kinematics uses motor composition:
//! ```text
//! M_total = M_base * T_1 * M_shoulder * T_2 * M_elbow * T_3
//! End_effector = M_total * Origin
//! ```

use crate::common::prelude::*;
use clifford::ops::Transform;
use clifford::specialized::projective::dim3::{Motor, Point};
use three_d::*;

/// Demo for 3D Robot Arm using native three-d rendering.
pub struct Projective3RobotDemo {
    /// 3D camera for rendering.
    camera: Option<Camera>,
    /// Orbit control for camera interaction.
    control: Option<OrbitControl>,
    /// World coordinate axes.
    world_axes: Option<Axes>,

    /// Base joint mesh (sphere).
    base_mesh: Option<Gm<Mesh, PhysicalMaterial>>,
    /// Link 1 mesh (thin cylinder).
    link1_mesh: Option<Gm<Mesh, PhysicalMaterial>>,
    /// Link 2 mesh (thin cylinder).
    link2_mesh: Option<Gm<Mesh, PhysicalMaterial>>,
    /// Link 3 mesh (thin cylinder).
    link3_mesh: Option<Gm<Mesh, PhysicalMaterial>>,
    /// End effector mesh (sphere).
    end_effector_mesh: Option<Gm<Mesh, PhysicalMaterial>>,

    /// Joint 1 (shoulder) mesh.
    joint1_mesh: Option<Gm<Mesh, PhysicalMaterial>>,
    /// Joint 2 (elbow) mesh.
    joint2_mesh: Option<Gm<Mesh, PhysicalMaterial>>,
    /// Joint 3 (wrist) mesh.
    joint3_mesh: Option<Gm<Mesh, PhysicalMaterial>>,

    /// Ambient lighting.
    ambient_light: Option<AmbientLight>,
    /// Key directional light.
    key_light: Option<DirectionalLight>,
    /// Fill light.
    fill_light: Option<DirectionalLight>,

    /// Joint 1 angle (base rotation around Z).
    theta1: f32,
    /// Joint 2 angle (shoulder rotation around Y).
    theta2: f32,
    /// Joint 3 angle (elbow rotation around Y).
    theta3: f32,

    /// Link 1 length.
    link1_length: f32,
    /// Link 2 length.
    link2_length: f32,
    /// Link 3 length.
    link3_length: f32,

    /// Animation controller.
    animation: Animation,
    /// Whether to show world coordinate axes.
    show_world_axes: bool,
    /// Whether to show joint spheres.
    show_joints: bool,
}

impl Default for Projective3RobotDemo {
    fn default() -> Self {
        Self {
            camera: None,
            control: None,
            world_axes: None,
            base_mesh: None,
            link1_mesh: None,
            link2_mesh: None,
            link3_mesh: None,
            end_effector_mesh: None,
            joint1_mesh: None,
            joint2_mesh: None,
            joint3_mesh: None,
            ambient_light: None,
            key_light: None,
            fill_light: None,
            theta1: 0.3,
            theta2: 0.8,  // ~45 degrees
            theta3: -0.7, // ~-40 degrees (elbow bent)
            link1_length: 1.0,
            link2_length: 1.5,
            link3_length: 1.0,
            animation: Animation::with_duration(4.0),
            show_world_axes: false,
            show_joints: true,
        }
    }
}

impl Projective3RobotDemo {
    /// Computes forward kinematics returning joint positions.
    ///
    /// Uses cumulative rotations to find link directions in world frame,
    /// then computes positions by translating along those directions.
    ///
    /// Key insight: ideal points (directions) are only affected by rotation,
    /// not translation. So we track cumulative rotation separately.
    fn forward_kinematics(&self) -> (Point<f64>, Point<f64>, Point<f64>, Point<f64>) {
        // Direction vectors as ideal points (w=0)
        let x_dir = Point::ideal(1.0, 0.0, 0.0);
        let z_dir = Point::ideal(0.0, 0.0, 1.0);

        // Joint rotations (around axes through origin)
        let r1 = Motor::from_rotation_z(f64::from(self.theta1));
        let r2 = Motor::from_rotation_y(f64::from(self.theta2));
        let r3 = Motor::from_rotation_y(f64::from(self.theta3));

        // Cumulative rotations for each frame
        let rot_1 = r1;
        let rot_12 = rot_1 * r2;
        let rot_123 = rot_12 * r3;

        // Link 1: along Z in frame 1 (Z rotation doesn't change Z direction)
        let link1_dir = rot_1.transform(&z_dir);
        let l1 = f64::from(self.link1_length);

        // Shoulder = base + L1 * link1_direction
        let shoulder =
            Point::from_cartesian(l1 * link1_dir.x(), l1 * link1_dir.y(), l1 * link1_dir.z());

        // Link 2: along X in frame 2, rotated to world frame
        let link2_dir = rot_12.transform(&x_dir);
        let l2 = f64::from(self.link2_length);

        // Elbow = shoulder + L2 * link2_direction
        let elbow = Point::from_cartesian(
            shoulder.cartesian_x() + l2 * link2_dir.x(),
            shoulder.cartesian_y() + l2 * link2_dir.y(),
            shoulder.cartesian_z() + l2 * link2_dir.z(),
        );

        // Link 3: along X in frame 3, rotated to world frame
        let link3_dir = rot_123.transform(&x_dir);
        let l3 = f64::from(self.link3_length);

        // End effector = elbow + L3 * link3_direction
        let end_effector = Point::from_cartesian(
            elbow.cartesian_x() + l3 * link3_dir.x(),
            elbow.cartesian_y() + l3 * link3_dir.y(),
            elbow.cartesian_z() + l3 * link3_dir.z(),
        );

        (shoulder, elbow, end_effector, end_effector)
    }

    /// Gets the motor for a specific joint configuration.
    fn get_joint_motor(&self, joint_index: usize) -> Motor<f64> {
        let m1 = Motor::from_rotation_z(f64::from(self.theta1));
        let t1 = Motor::from_translation(0.0, 0.0, f64::from(self.link1_length));
        let m2 = Motor::from_rotation_y(f64::from(self.theta2));
        let t2 = Motor::from_translation(f64::from(self.link2_length), 0.0, 0.0);
        let m3 = Motor::from_rotation_y(f64::from(self.theta3));
        let t3 = Motor::from_translation(f64::from(self.link3_length), 0.0, 0.0);

        match joint_index {
            0 => Motor::identity(),
            1 => (m1 * t1).unitized(),
            2 => (m1 * t1 * m2 * t2).unitized(),
            3 => (m1 * t1 * m2 * t2 * m3 * t3).unitized(),
            _ => Motor::identity(),
        }
    }

    /// Build a rotation/translation transform for a pre-scaled arrow.
    ///
    /// The arrow mesh is X-aligned and pre-scaled to the link length.
    /// This function just rotates it to point from start toward end and translates to start.
    fn arrow_transform(start: Vec3, end: Vec3) -> Mat4 {
        let dir = end - start;
        let length = dir.magnitude();
        if length < 1e-6 {
            return Mat4::from_translation(start);
        }
        let dir_normalized = dir / length;

        // Rotate from X-axis to target direction
        let x_axis = vec3(1.0, 0.0, 0.0);
        let dot = x_axis.dot(dir_normalized);

        let rotation = if dot > 0.9999 {
            Mat4::identity()
        } else if dot < -0.9999 {
            Mat4::from_angle_z(Rad(std::f32::consts::PI))
        } else {
            let rot_axis = x_axis.cross(dir_normalized);
            let rot_axis_len = rot_axis.magnitude();
            if rot_axis_len < 1e-6 {
                Mat4::identity()
            } else {
                let rot_axis_normalized = rot_axis / rot_axis_len;
                let angle = dot.clamp(-1.0, 1.0).acos();
                Mat4::from_axis_angle(rot_axis_normalized, Rad(angle))
            }
        };

        // Translate to start point (arrow starts at origin)
        Mat4::from_translation(start) * rotation
    }
}

impl VisualizationApp for Projective3RobotDemo {
    fn name(&self) -> &'static str {
        "3D PGA Robot Arm"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            // Animate joints with visible amplitude
            let t = self.animation.angle(); // 0 to 2π over loop duration
            // Base swings ±60 degrees (±1.05 rad)
            self.theta1 = t.sin() * 1.05;
            // Shoulder oscillates between 30 and 80 degrees
            self.theta2 = 0.96 + (t * 2.0).cos() * 0.44;
            // Elbow oscillates between -90 and -30 degrees
            self.theta3 = -1.05 + (t * 1.5).sin() * 0.52;
        }
    }

    fn render(&mut self, _ui: &mut egui::Ui) {
        // 3D rendering is handled by render_3d()
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Joint Controls ===
        group_header(ui, "Joint Angles");

        angle_slider_range(ui, "Base (Z-rot)", &mut self.theta1, -180.0, 180.0);
        angle_slider_range(ui, "Shoulder (Y-rot)", &mut self.theta2, -90.0, 90.0);
        angle_slider_range(ui, "Elbow (Y-rot)", &mut self.theta3, -135.0, 135.0);

        // === Link Lengths ===
        section_separator(ui, Some("Link Lengths"));
        ui.horizontal(|ui| {
            ui.label("Link 1:");
            ui.add(
                egui::DragValue::new(&mut self.link1_length)
                    .speed(0.1)
                    .range(0.5..=3.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Link 2:");
            ui.add(
                egui::DragValue::new(&mut self.link2_length)
                    .speed(0.1)
                    .range(0.5..=3.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Link 3:");
            ui.add(
                egui::DragValue::new(&mut self.link3_length)
                    .speed(0.1)
                    .range(0.5..=3.0),
            );
        });

        // === Animation ===
        section_separator(ui, Some("Animation"));
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // === Joint Positions ===
        section_separator(ui, Some("Joint Positions"));
        let (shoulder, elbow, end_effector, _) = self.forward_kinematics();
        if let Some((x, y, z)) = shoulder.to_cartesian() {
            ui.label(format!("Shoulder: ({:.2}, {:.2}, {:.2})", x, y, z));
        }
        if let Some((x, y, z)) = elbow.to_cartesian() {
            ui.label(format!("Elbow: ({:.2}, {:.2}, {:.2})", x, y, z));
        }
        if let Some((x, y, z)) = end_effector.to_cartesian() {
            ui.label(format!("End: ({:.2}, {:.2}, {:.2})", x, y, z));
        }

        // === Motor Components ===
        section_separator(ui, Some("Forward Kinematics"));
        let m_end = self.get_joint_motor(3);
        ga_value_display(
            ui,
            "M_end",
            &[("s", m_end.s() as f32), ("ps", m_end.ps() as f32)],
        );

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_world_axes, "World Axes");
            ui.checkbox(&mut self.show_joints, "Joint Spheres");
        });

        // === Camera Info ===
        section_separator(ui, Some("Camera"));
        ui.label(
            egui::RichText::new("Drag to orbit, scroll to zoom")
                .small()
                .weak(),
        );
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.label("3-DOF robot arm with motor-based forward kinematics");
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(PROJECTIVE3_ROBOT_EDUCATION)
    }
}

impl VisualizationApp3D for Projective3RobotDemo {
    fn init_3d(&mut self, context: &Context) {
        // Create perspective camera
        self.camera = Some(Camera::new_perspective(
            Viewport::new_at_origo(1, 1),
            vec3(5.0, 4.0, 6.0),
            vec3(0.0, 1.0, 0.0),
            vec3(0.0, 1.0, 0.0),
            degrees(45.0),
            0.1,
            100.0,
        ));

        // Orbit control
        self.control = Some(OrbitControl::new(vec3(0.0, 1.0, 0.0), 1.0, 30.0));

        // World axes
        self.world_axes = Some(Axes::new(context, 0.02, 2.5));

        // Create arrow meshes for links (pre-scaled like Axes does)
        let arrow_radius = 0.02;

        let mut arrow1 = CpuMesh::arrow(0.9, 0.6, 16);
        arrow1
            .transform(Mat4::from_nonuniform_scale(
                self.link1_length,
                arrow_radius,
                arrow_radius,
            ))
            .unwrap();

        let mut arrow2 = CpuMesh::arrow(0.9, 0.6, 16);
        arrow2
            .transform(Mat4::from_nonuniform_scale(
                self.link2_length,
                arrow_radius,
                arrow_radius,
            ))
            .unwrap();

        let mut arrow3 = CpuMesh::arrow(0.9, 0.6, 16);
        arrow3
            .transform(Mat4::from_nonuniform_scale(
                self.link3_length,
                arrow_radius,
                arrow_radius,
            ))
            .unwrap();

        let sphere = CpuMesh::sphere(16);

        // Base - small sphere at origin
        self.base_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(150, 150, 170, 255),
                    roughness: 0.6,
                    metallic: 0.2,
                    ..Default::default()
                },
            ),
        ));

        // Link meshes - pre-scaled arrows
        self.link1_mesh = Some(Gm::new(
            Mesh::new(context, &arrow1),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(100, 180, 255, 255), // Light blue
                    roughness: 0.5,
                    metallic: 0.1,
                    ..Default::default()
                },
            ),
        ));
        self.link2_mesh = Some(Gm::new(
            Mesh::new(context, &arrow2),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(100, 255, 180, 255), // Cyan
                    roughness: 0.5,
                    metallic: 0.1,
                    ..Default::default()
                },
            ),
        ));
        self.link3_mesh = Some(Gm::new(
            Mesh::new(context, &arrow3),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(255, 200, 100, 255), // Yellow-orange
                    roughness: 0.5,
                    metallic: 0.1,
                    ..Default::default()
                },
            ),
        ));

        // End effector (green sphere)
        self.end_effector_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(50, 255, 100, 255),
                    roughness: 0.4,
                    metallic: 0.3,
                    ..Default::default()
                },
            ),
        ));

        // Joint spheres (orange)
        self.joint1_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(255, 150, 50, 255),
                    roughness: 0.4,
                    metallic: 0.2,
                    ..Default::default()
                },
            ),
        ));
        self.joint2_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(255, 150, 50, 255),
                    roughness: 0.4,
                    metallic: 0.2,
                    ..Default::default()
                },
            ),
        ));
        self.joint3_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(255, 150, 50, 255),
                    roughness: 0.4,
                    metallic: 0.2,
                    ..Default::default()
                },
            ),
        ));

        // Lighting
        self.ambient_light = Some(AmbientLight::new(context, 0.3, Srgba::WHITE));
        self.key_light = Some(DirectionalLight::new(
            context,
            1.0,
            Srgba::WHITE,
            vec3(-1.0, -1.0, -1.0).normalize(),
        ));
        self.fill_light = Some(DirectionalLight::new(
            context,
            0.4,
            Srgba::new(200, 200, 255, 255),
            vec3(1.0, -0.5, 0.5).normalize(),
        ));
    }

    fn render_3d(&mut self, frame: &mut FrameInput) {
        // Compute joint positions
        let (shoulder, elbow, end_effector, _) = self.forward_kinematics();

        let base_pos = vec3(0.0, 0.0, 0.0);
        let shoulder_pos = shoulder
            .to_cartesian()
            .map(|(x, y, z)| vec3(x as f32, y as f32, z as f32))
            .unwrap_or(vec3(0.0, 0.0, 1.0));
        let elbow_pos = elbow
            .to_cartesian()
            .map(|(x, y, z)| vec3(x as f32, y as f32, z as f32))
            .unwrap_or(vec3(1.0, 0.0, 1.0));
        let end_pos = end_effector
            .to_cartesian()
            .map(|(x, y, z)| vec3(x as f32, y as f32, z as f32))
            .unwrap_or(vec3(2.0, 0.0, 1.0));

        let camera = self.camera.as_mut().unwrap();
        let control = self.control.as_mut().unwrap();

        camera.set_viewport(frame.viewport);
        control.handle_events(camera, &mut frame.events);

        // Update base (small sphere at origin)
        if let Some(mesh) = &mut self.base_mesh {
            mesh.set_transformation(Mat4::from_scale(0.08));
        }

        // Update link meshes (pre-scaled arrows, just rotate and translate)
        if let Some(mesh) = &mut self.link1_mesh {
            mesh.set_transformation(Self::arrow_transform(base_pos, shoulder_pos));
        }

        if let Some(mesh) = &mut self.link2_mesh {
            mesh.set_transformation(Self::arrow_transform(shoulder_pos, elbow_pos));
        }

        if let Some(mesh) = &mut self.link3_mesh {
            mesh.set_transformation(Self::arrow_transform(elbow_pos, end_pos));
        }

        // Update end effector
        if let Some(mesh) = &mut self.end_effector_mesh {
            mesh.set_transformation(Mat4::from_translation(end_pos) * Mat4::from_scale(0.1));
        }

        // Update joint spheres
        let joint_scale = 0.06;
        if let Some(mesh) = &mut self.joint1_mesh {
            mesh.set_transformation(
                Mat4::from_translation(shoulder_pos) * Mat4::from_scale(joint_scale),
            );
        }
        if let Some(mesh) = &mut self.joint2_mesh {
            mesh.set_transformation(
                Mat4::from_translation(elbow_pos) * Mat4::from_scale(joint_scale),
            );
        }
        if let Some(mesh) = &mut self.joint3_mesh {
            mesh.set_transformation(
                Mat4::from_translation(end_pos) * Mat4::from_scale(joint_scale),
            );
        }

        // Collect objects to render
        let mut objects: Vec<&dyn Object> = Vec::new();

        if self.show_world_axes {
            if let Some(axes) = &self.world_axes {
                objects.push(axes);
            }
        }

        if let Some(mesh) = &self.base_mesh {
            objects.push(mesh);
        }
        if let Some(mesh) = &self.link1_mesh {
            objects.push(mesh);
        }
        if let Some(mesh) = &self.link2_mesh {
            objects.push(mesh);
        }
        if let Some(mesh) = &self.link3_mesh {
            objects.push(mesh);
        }
        if let Some(mesh) = &self.end_effector_mesh {
            objects.push(mesh);
        }

        if self.show_joints {
            if let Some(mesh) = &self.joint1_mesh {
                objects.push(mesh);
            }
            if let Some(mesh) = &self.joint2_mesh {
                objects.push(mesh);
            }
            if let Some(mesh) = &self.joint3_mesh {
                objects.push(mesh);
            }
        }

        // Collect lights
        let lights: Vec<&dyn Light> = vec![
            self.ambient_light.as_ref().unwrap(),
            self.key_light.as_ref().unwrap(),
            self.fill_light.as_ref().unwrap(),
        ];

        // Render
        frame.screen().render(camera, objects, &lights);
    }
}

/// Educational content for the 3D Robot Arm visualization.
const PROJECTIVE3_ROBOT_EDUCATION: EducationalContent = EducationalContent {
    title: "3D Robot Arm with PGA Motors",

    overview: "\
This demo shows forward kinematics of a 3-DOF robot arm using PGA motors. \
Each joint is a rotation motor, and link lengths are encoded as translation motors. \
The end effector position is computed by composing all motors in the kinematic chain.

Motors provide a compact, numerically stable representation for robot kinematics \
that avoids gimbal lock and matrix drift issues.",

    math_background: "\
Forward kinematics computes end effector pose from joint angles:

    M_end = M_1 * T_1 * M_2 * T_2 * M_3 * T_3

Where:
  - M_i = rotation motor for joint i
  - T_i = translation motor for link i

Joint rotation motors:
  - M_1 = rotation around Z (base yaw)
  - M_2 = rotation around local Y (shoulder pitch)
  - M_3 = rotation around local Y (elbow pitch)

Position computation:
  - P_shoulder = (M_1 * T_1) * Origin
  - P_elbow = (M_1 * T_1 * M_2 * T_2) * Origin
  - P_end = (M_1 * T_1 * M_2 * T_2 * M_3 * T_3) * Origin",

    how_to_use: "\
- Adjust joint angles with sliders
- Modify link lengths to change arm geometry
- Click Play to animate the arm
- Observe end effector position update
- Drag to orbit camera, scroll to zoom",

    key_concepts: "\
- Motors compose multiplicatively: M_total = M_2 * M_1
- Rotation motors: M = cos(t/2) + sin(t/2)*axis
- Translation motors: T = 1 + d/2 where d is displacement bivector
- Unitization prevents numerical drift
- No gimbal lock (unlike Euler angles)
- Interpolation via motor slerp is natural",

    resources: &[
        (
            "RGA Wiki - Motors",
            "https://rigidgeometricalgebra.org/wiki/index.php?title=Motor",
        ),
        (
            "PGA for Computer Scientists",
            "https://bivector.net/PGA4CS.pdf",
        ),
        (
            "Geometric Algebra for Robotics",
            "https://geometricalgebra.org/",
        ),
    ],
};
