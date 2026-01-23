//! 3D Projective Motor Visualization using native three-d rendering.
//!
//! This demo demonstrates rigid body motion in 3D PGA using Motors,
//! rendered with native GPU acceleration via three-d.
//!
//! Features:
//! - Wireframe box transformed by Motor
//! - Coordinate frame showing transformation
//! - Screw axis visualization
//! - Motor decomposition display (rotation angle, translation)
//! - Animation support

use crate::common::prelude::*;
use clifford::ops::Transform;
use clifford::specialized::euclidean::dim3::Vector as EuclideanVector;
use clifford::specialized::projective::dim3::{Motor, Point};
use std::f32::consts::TAU;
use three_d::*;

/// Transform a direction vector using a rotation motor.
///
/// This uses the point-based transformation and extracts the direction.
fn rotation_transform_vector(motor: &Motor<f64>, dir: &EuclideanVector<f64>) -> (f64, f64, f64) {
    // Transform a point at the direction from origin
    let pt = Point::from_cartesian(dir.x(), dir.y(), dir.z());
    let origin = Point::origin();

    let pt_t = motor.transform(&pt);
    let o_t = motor.transform(&origin);

    // Get the direction by subtracting transformed origin from transformed point
    let (px, py, pz) = pt_t.to_cartesian().unwrap_or((dir.x(), dir.y(), dir.z()));
    let (ox, oy, oz) = o_t.to_cartesian().unwrap_or((0.0, 0.0, 0.0));

    (px - ox, py - oy, pz - oz)
}

/// Demo for 3D Projective Motors using native three-d rendering.
pub struct Projective3MotorDemo {
    // === Three-d 3D Resources ===
    /// 3D camera for rendering.
    camera: Option<Camera>,
    /// Orbit control for camera interaction.
    control: Option<OrbitControl>,
    /// The box mesh being transformed.
    box_mesh: Option<Gm<Mesh, PhysicalMaterial>>,
    /// Coordinate axes visualization (world frame).
    world_axes: Option<Axes>,
    /// Object X axis (red).
    object_axis_x: Option<Gm<Mesh, ColorMaterial>>,
    /// Object Y axis (green).
    object_axis_y: Option<Gm<Mesh, ColorMaterial>>,
    /// Object Z axis (blue).
    object_axis_z: Option<Gm<Mesh, ColorMaterial>>,
    /// Ambient lighting for the scene.
    ambient_light: Option<AmbientLight>,
    /// Key light (main directional light).
    key_light: Option<DirectionalLight>,
    /// Fill light (secondary directional light).
    fill_light: Option<DirectionalLight>,
    /// Back light (rim lighting).
    back_light: Option<DirectionalLight>,

    // === Motor Parameters ===
    /// Rotation angle in radians.
    rotation_angle: f32,
    /// Rotation axis (unit vector).
    rotation_axis: [f32; 3],
    /// Translation vector.
    translation: [f32; 3],
    /// Animation controller.
    animation: Animation,

    // === Display Options ===
    /// Whether to show world coordinate axes.
    show_world_axes: bool,
    /// Whether to show object coordinate axes.
    show_object_axes: bool,
    /// Whether to show the transformed box.
    show_box: bool,
}

impl Default for Projective3MotorDemo {
    fn default() -> Self {
        Self {
            camera: None,
            control: None,
            box_mesh: None,
            world_axes: None,
            object_axis_x: None,
            object_axis_y: None,
            object_axis_z: None,
            ambient_light: None,
            key_light: None,
            fill_light: None,
            back_light: None,
            rotation_angle: 0.0,
            rotation_axis: [0.0, 0.0, 1.0], // Z-axis by default
            translation: [0.0, 0.0, 0.0],
            animation: Animation::with_duration(4.0),
            show_world_axes: true,
            show_object_axes: true,
            show_box: true,
        }
    }
}

impl Projective3MotorDemo {
    /// Get the current motor from parameters.
    fn current_motor(&self) -> Motor<f64> {
        // Check if axis is zero - if so, use identity rotation
        let axis_len_sq = self.rotation_axis[0].powi(2)
            + self.rotation_axis[1].powi(2)
            + self.rotation_axis[2].powi(2);

        let rotation = if axis_len_sq < 1e-10 {
            // Zero axis - use identity (no rotation)
            Motor::identity()
        } else {
            let axis = EuclideanVector::new(
                f64::from(self.rotation_axis[0]),
                f64::from(self.rotation_axis[1]),
                f64::from(self.rotation_axis[2]),
            );
            Motor::from_axis_angle(&axis, f64::from(self.rotation_angle))
        };

        let translation = Motor::from_translation(
            f64::from(self.translation[0]),
            f64::from(self.translation[1]),
            f64::from(self.translation[2]),
        );

        // Combine: first rotate, then translate
        // In PGA, motor composition is M_total = M_2 * M_1 (applied right to left)
        // Unitize to ensure valid rigid transformation
        (translation * rotation).unitized()
    }

    /// Build transformation matrix from rotation motor and translation.
    ///
    /// We apply rotation and translation separately to avoid potential issues
    /// with motor composition.
    fn build_transform(&self) -> Mat4 {
        // Get rotation motor (or identity if axis is zero)
        let axis_len_sq = self.rotation_axis[0].powi(2)
            + self.rotation_axis[1].powi(2)
            + self.rotation_axis[2].powi(2);

        let rotation = if axis_len_sq < 1e-10 {
            Motor::identity()
        } else {
            let axis = EuclideanVector::new(
                f64::from(self.rotation_axis[0]),
                f64::from(self.rotation_axis[1]),
                f64::from(self.rotation_axis[2]),
            );
            Motor::from_axis_angle(&axis, f64::from(self.rotation_angle))
        };

        // Transform basis vectors using ONLY the rotation motor
        let x_vec = EuclideanVector::new(1.0, 0.0, 0.0);
        let y_vec = EuclideanVector::new(0.0, 1.0, 0.0);
        let z_vec = EuclideanVector::new(0.0, 0.0, 1.0);

        // Get rotated basis vectors
        let rot_x = rotation_transform_vector(&rotation, &x_vec);
        let rot_y = rotation_transform_vector(&rotation, &y_vec);
        let rot_z = rotation_transform_vector(&rotation, &z_vec);

        // Translation is applied after rotation
        let tx = f64::from(self.translation[0]);
        let ty = f64::from(self.translation[1]);
        let tz = f64::from(self.translation[2]);

        // three-d uses column-major matrices
        Mat4::new(
            rot_x.0 as f32,
            rot_x.1 as f32,
            rot_x.2 as f32,
            0.0,
            rot_y.0 as f32,
            rot_y.1 as f32,
            rot_y.2 as f32,
            0.0,
            rot_z.0 as f32,
            rot_z.1 as f32,
            rot_z.2 as f32,
            0.0,
            tx as f32,
            ty as f32,
            tz as f32,
            1.0,
        )
    }

    /// Normalize the rotation axis.
    fn normalize_axis(&mut self) {
        let len = (self.rotation_axis[0].powi(2)
            + self.rotation_axis[1].powi(2)
            + self.rotation_axis[2].powi(2))
        .sqrt();
        if len > 1e-6 {
            self.rotation_axis[0] /= len;
            self.rotation_axis[1] /= len;
            self.rotation_axis[2] /= len;
        } else {
            // Default to Z axis if zero
            self.rotation_axis = [0.0, 0.0, 1.0];
        }
    }
}

impl VisualizationApp for Projective3MotorDemo {
    fn name(&self) -> &'static str {
        "3D PGA Motors - Rigid Body Motion"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            self.rotation_angle = self.animation.progress() * TAU;
        }
    }

    fn render(&mut self, _ui: &mut egui::Ui) {
        // 3D rendering is handled by render_3d(), not this method
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Rotation Controls ===
        group_header(ui, "Rotation");
        angle_slider_range(ui, "Angle", &mut self.rotation_angle, -360.0, 360.0);

        ui.horizontal(|ui| {
            ui.label("Axis:");
            ui.add(
                egui::DragValue::new(&mut self.rotation_axis[0])
                    .speed(0.05)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.rotation_axis[1])
                    .speed(0.05)
                    .prefix("y: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.rotation_axis[2])
                    .speed(0.05)
                    .prefix("z: "),
            );
        });

        ui.horizontal(|ui| {
            if ui.button("X").clicked() {
                self.rotation_axis = [1.0, 0.0, 0.0];
            }
            if ui.button("Y").clicked() {
                self.rotation_axis = [0.0, 1.0, 0.0];
            }
            if ui.button("Z").clicked() {
                self.rotation_axis = [0.0, 0.0, 1.0];
            }
            if ui.button("Normalize").clicked() {
                self.normalize_axis();
            }
        });

        // === Translation Controls ===
        section_separator(ui, Some("Translation"));
        ui.horizontal(|ui| {
            ui.add(
                egui::DragValue::new(&mut self.translation[0])
                    .speed(0.1)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.translation[1])
                    .speed(0.1)
                    .prefix("y: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.translation[2])
                    .speed(0.1)
                    .prefix("z: "),
            );
        });

        if ui.button("Reset Translation").clicked() {
            self.translation = [0.0, 0.0, 0.0];
        }

        // === Animation ===
        section_separator(ui, Some("Animation"));
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // === Motor Components ===
        section_separator(ui, Some("Motor Components"));
        let motor = self.current_motor();
        ga_value_display(
            ui,
            "M",
            &[
                ("s", motor.s() as f32),
                ("e_1_2", motor.tz() as f32),
                ("e_1_3", motor.ty() as f32),
                ("e_1_4", motor.tx() as f32),
                ("e_2_3", motor.rx() as f32),
                ("e_2_4", motor.ry() as f32),
                ("e_3_4", motor.rz() as f32),
                ("e_1_2_3_4", motor.ps() as f32),
            ],
        );

        // Motor decomposition info
        info_box(
            ui,
            &format!(
                "Motor = Translation * Rotation\n\
                 Rotation: {:.1} deg around ({:.2}, {:.2}, {:.2})\n\
                 Translation: ({:.2}, {:.2}, {:.2})",
                self.rotation_angle.to_degrees(),
                self.rotation_axis[0],
                self.rotation_axis[1],
                self.rotation_axis[2],
                self.translation[0],
                self.translation[1],
                self.translation[2]
            ),
        );

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_box, "Box");
            ui.checkbox(&mut self.show_world_axes, "World Axes");
            ui.checkbox(&mut self.show_object_axes, "Object Axes");
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
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            ui.label("3D rigid body motion using PGA motors");
            ui.separator();
            ui.colored_label(x_axis(&ctx), "X");
            ui.colored_label(y_axis(&ctx), "Y");
            ui.colored_label(z_axis(&ctx), "Z");
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(PROJECTIVE3_MOTOR_EDUCATION)
    }
}

impl VisualizationApp3D for Projective3MotorDemo {
    fn init_3d(&mut self, context: &Context) {
        // Create perspective camera
        self.camera = Some(Camera::new_perspective(
            Viewport::new_at_origo(1, 1),
            vec3(4.0, 3.0, 5.0),
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 1.0, 0.0),
            degrees(45.0),
            0.1,
            100.0,
        ));

        // Orbit control
        self.control = Some(OrbitControl::new(vec3(0.0, 0.0, 0.0), 1.0, 30.0));

        // Create box mesh with physically-based material
        let box_mesh = CpuMesh::cube();
        self.box_mesh = Some(Gm::new(
            Mesh::new(context, &box_mesh),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(100, 180, 255, 255),
                    roughness: 0.5,
                    metallic: 0.0,
                    ..Default::default()
                },
            ),
        ));

        // Create world coordinate axes
        self.world_axes = Some(Axes::new(context, 0.02, 2.5));

        // Create object coordinate axes as individual arrows
        // We'll transform these properly in render_3d (point for base, vector for direction)
        // Arrow mesh: along X-axis, from 0 to 1, radius 1
        // tail_length: 0.7 (70% shaft, 30% head), tail_radius: 0.4 (40% of total radius)
        let arrow = CpuMesh::arrow(0.7, 0.4, 16);

        // X axis (red)
        self.object_axis_x = Some(Gm::new(
            Mesh::new(context, &arrow),
            ColorMaterial {
                color: Srgba::new(255, 50, 50, 255),
                ..Default::default()
            },
        ));

        // Y axis (green)
        self.object_axis_y = Some(Gm::new(
            Mesh::new(context, &arrow),
            ColorMaterial {
                color: Srgba::new(50, 255, 50, 255),
                ..Default::default()
            },
        ));

        // Z axis (blue)
        self.object_axis_z = Some(Gm::new(
            Mesh::new(context, &arrow),
            ColorMaterial {
                color: Srgba::new(50, 50, 255, 255),
                ..Default::default()
            },
        ));

        // Three-point lighting
        self.ambient_light = Some(AmbientLight::new(context, 0.2, Srgba::WHITE));

        self.key_light = Some(DirectionalLight::new(
            context,
            1.0,
            Srgba::WHITE,
            vec3(-1.0, -1.0, -1.0).normalize(),
        ));

        self.fill_light = Some(DirectionalLight::new(
            context,
            0.5,
            Srgba::new(200, 200, 255, 255),
            vec3(1.0, -0.5, 0.0).normalize(),
        ));

        self.back_light = Some(DirectionalLight::new(
            context,
            0.3,
            Srgba::new(255, 220, 200, 255),
            vec3(0.0, -0.5, 1.0).normalize(),
        ));
    }

    fn render_3d(&mut self, frame: &mut FrameInput) {
        // Build transformation matrix (rotation + translation)
        let transform = self.build_transform();

        let camera = self.camera.as_mut().unwrap();
        let control = self.control.as_mut().unwrap();

        // Update camera viewport
        camera.set_viewport(frame.viewport);

        // Handle camera orbit control
        control.handle_events(camera, &mut frame.events);

        // Apply transformation to box
        if let Some(box_mesh) = &mut self.box_mesh {
            box_mesh.set_transformation(transform);
        }

        // Transform object axes properly:
        // - Base point transforms as a point (rotation + translation)
        // - Direction transforms as a vector (rotation only)
        // Match world axes parameters (Axes::new uses radius 0.02, length 2.5)
        let axis_length = 2.5f32;
        let axis_radius = 0.02f32;

        // Get rotation motor for direction transformation
        let axis_len_sq = self.rotation_axis[0].powi(2)
            + self.rotation_axis[1].powi(2)
            + self.rotation_axis[2].powi(2);

        let rotation = if axis_len_sq < 1e-10 {
            Motor::identity()
        } else {
            let axis = EuclideanVector::new(
                f64::from(self.rotation_axis[0]),
                f64::from(self.rotation_axis[1]),
                f64::from(self.rotation_axis[2]),
            );
            Motor::from_axis_angle(&axis, f64::from(self.rotation_angle))
        };

        // Transformed origin (point: rotation + translation)
        let origin = vec3(
            self.translation[0],
            self.translation[1],
            self.translation[2],
        );

        // Transformed directions (vectors: rotation only)
        let x_vec = EuclideanVector::new(1.0, 0.0, 0.0);
        let y_vec = EuclideanVector::new(0.0, 1.0, 0.0);
        let z_vec = EuclideanVector::new(0.0, 0.0, 1.0);

        let rot_x = rotation_transform_vector(&rotation, &x_vec);
        let rot_y = rotation_transform_vector(&rotation, &y_vec);
        let rot_z = rotation_transform_vector(&rotation, &z_vec);

        let dir_x = vec3(rot_x.0 as f32, rot_x.1 as f32, rot_x.2 as f32);
        let dir_y = vec3(rot_y.0 as f32, rot_y.1 as f32, rot_y.2 as f32);
        let dir_z = vec3(rot_z.0 as f32, rot_z.1 as f32, rot_z.2 as f32);

        // Transform each axis arrow
        // Arrow default: along X axis, from 0 to 1, radius 1
        // We need to: scale, rotate to align with direction, translate to origin
        if let Some(axis_x) = &mut self.object_axis_x {
            let t = axis_transform(origin, dir_x, axis_length, axis_radius);
            axis_x.set_transformation(t);
        }
        if let Some(axis_y) = &mut self.object_axis_y {
            let t = axis_transform(origin, dir_y, axis_length, axis_radius);
            axis_y.set_transformation(t);
        }
        if let Some(axis_z) = &mut self.object_axis_z {
            let t = axis_transform(origin, dir_z, axis_length, axis_radius);
            axis_z.set_transformation(t);
        }

        // Collect objects to render
        let mut objects: Vec<&dyn Object> = Vec::new();

        if self.show_box {
            if let Some(box_mesh) = &self.box_mesh {
                objects.push(box_mesh);
            }
        }

        if self.show_world_axes {
            if let Some(world_axes) = &self.world_axes {
                objects.push(world_axes);
            }
        }

        if self.show_object_axes {
            if let Some(axis_x) = &self.object_axis_x {
                objects.push(axis_x);
            }
            if let Some(axis_y) = &self.object_axis_y {
                objects.push(axis_y);
            }
            if let Some(axis_z) = &self.object_axis_z {
                objects.push(axis_z);
            }
        }

        // Collect lights
        let lights: Vec<&dyn Light> = vec![
            self.ambient_light.as_ref().unwrap(),
            self.key_light.as_ref().unwrap(),
            self.fill_light.as_ref().unwrap(),
            self.back_light.as_ref().unwrap(),
        ];

        // Render
        frame.screen().render(camera, objects, &lights);
    }
}

/// Build a transformation matrix for an axis arrow.
///
/// The arrow is transformed from its default (X-aligned, from 0 to 1, radius 1)
/// to point from `origin` along `direction` with the specified length and radius.
fn axis_transform(origin: Vec3, direction: Vec3, length: f32, radius: f32) -> Mat4 {
    // The default arrow is along X axis, from 0 to 1, radius 1
    // We need to:
    // 1. Scale: length along X, radius for Y and Z
    // 2. Rotate: from X axis to the target direction
    // 3. Translate: base to origin

    let dir_len = direction.magnitude();
    if dir_len < 1e-6 {
        // Degenerate direction - return identity (invisible axis)
        return Mat4::identity();
    }
    let dir_normalized = direction / dir_len;

    // Scale matrix: length in X, radius in Y/Z
    let scale = Mat4::from_nonuniform_scale(length, radius, radius);

    // Rotation from X axis to target direction using a robust method
    let x_axis = vec3(1.0, 0.0, 0.0);
    let dot = x_axis.dot(dir_normalized);

    let rotation = if dot > 0.9999 {
        // Already aligned with +X
        Mat4::identity()
    } else if dot < -0.9999 {
        // Opposite to X - rotate 180 degrees around Y
        Mat4::from_angle_y(Rad(std::f32::consts::PI))
    } else {
        // General rotation using axis-angle
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

    // Translation: arrow starts at 0, so just translate to origin
    let translation = Mat4::from_translation(origin);

    // Combined: translate * rotate * scale
    translation * rotation * scale
}

/// Educational content for the 3D Projective Motor visualization.
const PROJECTIVE3_MOTOR_EDUCATION: EducationalContent = EducationalContent {
    title: "Motors in 3D Projective Geometric Algebra",

    overview: "\
Motors are the PGA representation of rigid body transformations (rotation + translation). \
Unlike 4x4 matrices, motors are compact (8 components), always invertible, and support \
smooth interpolation.

A motor encodes both rotation and translation as a single algebraic element. \
Objects transform via the sandwich product: X' = M X M~",

    math_background: "\
A 3D PGA Motor has 8 components (even subalgebra):
    M = s + e_1_2 + e_1_3 + e_1_4 + e_2_3 + e_2_4 + e_3_4 + e_1_2_3_4

The components encode:
  - Rotation: stored in e_1_4, e_2_4, e_3_4 (velocity bivectors)
  - Translation: stored in e_2_3, e_3_1, e_1_2 (moment bivectors)
  - Identity part: s and e_1_2_3_4 (pseudoscalar)

Motor composition:
    M_total = M_2 * M_1 (applied right to left)

Pure rotation around axis a by angle theta:
    R = cos(theta/2) + sin(theta/2)(a_x*e_4_1 + a_y*e_4_2 + a_z*e_4_3)

Pure translation by vector t:
    T = 1 + (t_x*e_2_3 + t_y*e_3_1 + t_z*e_1_2)/2

Combined motion:
    M = T * R (first rotate, then translate)",

    how_to_use: "\
- Drag in the 3D view to orbit the camera
- Scroll to zoom in/out
- Adjust rotation angle with the slider
- Set rotation axis (X, Y, Z buttons or custom)
- Set translation vector
- Click Play to animate continuous rotation
- Toggle display of box and axes
- World axes show the fixed reference frame
- Object axes show the transformed frame",

    key_concepts: "\
- Motors use HALF the rotation angle (like quaternions)
- Motors compose multiplicatively: M_total = M_2 * M_1
- Motor inverse: M^-1 = M~ / |M|^2 (reverse normalized)
- Motors satisfy geometric constraint: <M M~>_0 = 1
- Screw motion: rotation + translation along same axis
- Motors interpolate naturally (SLERP-like)",

    resources: &[
        (
            "Rigid Geometric Algebra Wiki - Motors",
            "https://rigidgeometricalgebra.org/wiki/index.php?title=Motor",
        ),
        (
            "Look, Ma, No Matrices!",
            "https://enkimute.github.io/LookMaNoMatrices/",
        ),
        (
            "PGA for Computer Scientists",
            "https://bivector.net/PGA4CS.pdf",
        ),
    ],
};
