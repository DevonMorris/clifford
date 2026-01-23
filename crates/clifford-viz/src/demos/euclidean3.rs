//! 3D Euclidean rotor visualization using native three-d rendering.
//!
//! This demo demonstrates rotation in 3D Euclidean space using rotors from
//! Geometric Algebra, rendered with native GPU acceleration via three-d.
//!
//! Features:
//! - Wireframe cube rotating under rotor transformation
//! - Coordinate axes (X=red, Y=green, Z=blue)
//! - Camera orbit controls (drag to rotate view)
//! - Rotation plane selection (XY, XZ, YZ)
//! - Animation support

use crate::common::prelude::*;
use clifford::ops::Transform;
use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
use std::f32::consts::TAU;
use three_d::*;

/// Which plane to rotate in.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum RotationPlane {
    /// XY plane (rotation around Z axis)
    #[default]
    XY,
    /// XZ plane (rotation around Y axis)
    XZ,
    /// YZ plane (rotation around X axis)
    YZ,
}

impl RotationPlane {
    /// Convert to a unit bivector representing this plane.
    fn to_bivector(self) -> Bivector<f64> {
        match self {
            Self::XY => Bivector::unit_rz(), // rz = e12 = XY plane
            Self::XZ => Bivector::unit_ry(), // ry = e13 = XZ plane
            Self::YZ => Bivector::unit_rx(), // rx = e23 = YZ plane
        }
    }

    /// Human-readable name.
    fn name(self) -> &'static str {
        match self {
            Self::XY => "XY (around Z)",
            Self::XZ => "XZ (around Y)",
            Self::YZ => "YZ (around X)",
        }
    }
}

/// Demo for 3D Euclidean rotors using native three-d rendering.
pub struct Euclidean3Demo {
    /// 3D camera for rendering.
    camera: Option<Camera>,
    /// Orbit control for camera interaction.
    control: Option<OrbitControl>,
    /// The cube mesh being rotated.
    cube: Option<Gm<Mesh, PhysicalMaterial>>,
    /// Coordinate axes visualization.
    axes: Option<Axes>,
    /// Ambient lighting for the scene.
    ambient_light: Option<AmbientLight>,
    /// Key light (main directional light).
    key_light: Option<DirectionalLight>,
    /// Fill light (secondary directional light).
    fill_light: Option<DirectionalLight>,
    /// Back light (rim lighting).
    back_light: Option<DirectionalLight>,

    /// Current rotation angle in radians.
    rotation_angle: f32,
    /// The plane of rotation (determines rotation axis).
    rotation_plane: RotationPlane,
    /// Animation controller.
    animation: Animation,

    /// Whether to show coordinate axes.
    show_axes: bool,
    /// Whether to show the rotating cube.
    show_cube: bool,
}

impl Default for Euclidean3Demo {
    fn default() -> Self {
        Self {
            camera: None,
            control: None,
            cube: None,
            axes: None,
            ambient_light: None,
            key_light: None,
            fill_light: None,
            back_light: None,
            rotation_angle: 0.0,
            rotation_plane: RotationPlane::XY,
            animation: Animation::with_duration(4.0),
            show_axes: true,
            show_cube: true,
        }
    }
}

impl Euclidean3Demo {
    /// Get the current rotor from angle and plane.
    fn current_rotor(&self) -> Rotor<f64> {
        Rotor::from_angle_plane(
            f64::from(self.rotation_angle),
            self.rotation_plane.to_bivector(),
        )
    }

    /// Convert a rotor to a three-d Mat4 transformation matrix.
    fn rotor_to_mat4(rotor: &Rotor<f64>) -> Mat4 {
        // Transform basis vectors to get rotation matrix columns
        let x = rotor.transform(&Vector::unit_x());
        let y = rotor.transform(&Vector::unit_y());
        let z = rotor.transform(&Vector::unit_z());

        // three-d uses column-major matrices
        Mat4::new(
            x.x() as f32,
            x.y() as f32,
            x.z() as f32,
            0.0,
            y.x() as f32,
            y.y() as f32,
            y.z() as f32,
            0.0,
            z.x() as f32,
            z.y() as f32,
            z.z() as f32,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        )
    }
}

impl VisualizationApp for Euclidean3Demo {
    fn name(&self) -> &'static str {
        "3D Euclidean Rotors"
    }

    fn update(&mut self, dt: f32) {
        self.animation.update(dt);
        if self.animation.playing {
            self.rotation_angle = self.animation.progress() * TAU;
        }
    }

    fn render(&mut self, _ui: &mut egui::Ui) {
        // 3D rendering is handled by render_3d(), not this method
        // This is called but we don't need to do anything
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Rotation Controls ===
        group_header(ui, "Rotation");
        angle_slider_range(ui, "Angle", &mut self.rotation_angle, -360.0, 360.0);

        ui.horizontal(|ui| {
            ui.label("Plane:");
            egui::ComboBox::from_id_salt("rotation_plane")
                .selected_text(self.rotation_plane.name())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.rotation_plane,
                        RotationPlane::XY,
                        "XY (around Z)",
                    );
                    ui.selectable_value(
                        &mut self.rotation_plane,
                        RotationPlane::XZ,
                        "XZ (around Y)",
                    );
                    ui.selectable_value(
                        &mut self.rotation_plane,
                        RotationPlane::YZ,
                        "YZ (around X)",
                    );
                });
        });

        // === Animation ===
        ui.add_space(spacing::XS);
        animation_controls(ui, &mut self.animation);
        progress_slider(ui, &mut self.animation);

        // === Rotor Components ===
        section_separator(ui, Some("Rotor Components"));
        let rotor = self.current_rotor();
        ga_value_display(
            ui,
            "R",
            &[
                ("1", rotor.s() as f32),
                ("e_2_3", rotor.rx() as f32),
                ("e_1_3", rotor.ry() as f32),
                ("e_1_2", rotor.rz() as f32),
            ],
        );

        info_box(
            ui,
            &format!(
                "R = cos(theta/2) + sin(theta/2)B\n\
                 theta = {:.1} deg\n\
                 B = {} plane",
                self.rotation_angle.to_degrees(),
                match self.rotation_plane {
                    RotationPlane::XY => "e_1_2",
                    RotationPlane::XZ => "e_1_3",
                    RotationPlane::YZ => "e_2_3",
                }
            ),
        );

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_cube, "Cube");
            ui.checkbox(&mut self.show_axes, "Axes");
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
            ui.label("3D rotation using rotor (GA equivalent of quaternion)");
            ui.separator();
            ui.colored_label(x_axis(&ctx), "X");
            ui.colored_label(y_axis(&ctx), "Y");
            ui.colored_label(z_axis(&ctx), "Z");
        });
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(EUCLIDEAN3_EDUCATION)
    }
}

impl VisualizationApp3D for Euclidean3Demo {
    fn init_3d(&mut self, context: &Context) {
        // Create perspective camera
        self.camera = Some(Camera::new_perspective(
            Viewport::new_at_origo(1, 1),
            vec3(3.0, 2.0, 4.0),
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 1.0, 0.0),
            degrees(45.0),
            0.1,
            100.0,
        ));

        // Orbit control
        self.control = Some(OrbitControl::new(vec3(0.0, 0.0, 0.0), 1.0, 20.0));

        // Create cube mesh with physically-based material for proper shading
        let cube_mesh = CpuMesh::cube();
        self.cube = Some(Gm::new(
            Mesh::new(context, &cube_mesh),
            PhysicalMaterial::new_opaque(
                context,
                &CpuMaterial {
                    albedo: Srgba::new(100, 150, 255, 255),
                    roughness: 0.5,
                    metallic: 0.0,
                    ..Default::default()
                },
            ),
        ));

        // Create coordinate axes
        self.axes = Some(Axes::new(context, 0.02, 2.0));

        // Three-point lighting for good shading on all cube faces
        self.ambient_light = Some(AmbientLight::new(context, 0.2, Srgba::WHITE));

        // Key light: main light from upper-right-front
        self.key_light = Some(DirectionalLight::new(
            context,
            1.0,
            Srgba::WHITE,
            vec3(-1.0, -1.0, -1.0).normalize(),
        ));

        // Fill light: softer light from left side
        self.fill_light = Some(DirectionalLight::new(
            context,
            0.5,
            Srgba::new(200, 200, 255, 255),
            vec3(1.0, -0.5, 0.0).normalize(),
        ));

        // Back light: rim light from behind
        self.back_light = Some(DirectionalLight::new(
            context,
            0.3,
            Srgba::new(255, 220, 200, 255),
            vec3(0.0, -0.5, 1.0).normalize(),
        ));
    }

    fn render_3d(&mut self, frame: &mut FrameInput) {
        // Compute rotor transformation first (before mutable borrows)
        let rotor = self.current_rotor();
        let transform = Self::rotor_to_mat4(&rotor);

        let camera = self.camera.as_mut().unwrap();
        let control = self.control.as_mut().unwrap();

        // Update camera viewport
        camera.set_viewport(frame.viewport);

        // Handle camera orbit control
        control.handle_events(camera, &mut frame.events);

        if let Some(cube) = &mut self.cube {
            cube.set_transformation(transform);
        }

        // Collect objects to render
        let mut objects: Vec<&dyn Object> = Vec::new();

        if self.show_cube {
            if let Some(cube) = &self.cube {
                objects.push(cube);
            }
        }

        if self.show_axes {
            if let Some(axes) = &self.axes {
                objects.push(axes);
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

/// Educational content for the 3D Euclidean rotor visualization.
const EUCLIDEAN3_EDUCATION: EducationalContent = EducationalContent {
    title: "3D Rotors in Geometric Algebra",

    overview: "\
This visualization demonstrates 3D rotations using Geometric Algebra (GA) rotors. \
Rotors are the GA equivalent of quaternions, but with clearer geometric meaning.

The key insight is that rotations happen in PLANES, not around axes. \
A rotor encodes rotation in a specific plane (bivector) by a specific angle.",

    math_background: "\
A 3D rotor R encoding rotation by angle theta in plane B is:

    R = cos(theta/2) + sin(theta/2)B

where B is a unit bivector representing the rotation plane:
  - e_1_2 (XY plane) rotates around Z axis
  - e_1_3 (XZ plane) rotates around Y axis
  - e_2_3 (YZ plane) rotates around X axis

Vectors transform via the sandwich product:
    v' = R v R~
where R~ is the reverse (conjugate).

Rotor components relate to quaternions:
    R = s + rz*e_1_2 + ry*e_1_3 + rx*e_2_3
    q = w + xi + yj + zk

    s = w, rx = x, ry = y, rz = z",

    how_to_use: "\
- Drag in the 3D view to orbit the camera
- Scroll to zoom in/out
- Use the angle slider to rotate the cube
- Change the rotation plane with the dropdown
- Click Play to animate continuous rotation
- Toggle display of cube and axes",

    key_concepts: "\
- Rotors use HALF the rotation angle (like quaternions)
- Rotation happens in a PLANE (bivector), not around an axis
- The axis of rotation is perpendicular to the plane
- Rotor composition: R_total = R_2 * R_1
- Unit rotors preserve vector lengths
- In 3D Euclidean GA, rotors ARE quaternions",

    resources: &[
        (
            "Rigid Geometric Algebra Wiki",
            "https://rigidgeometricalgebra.org/wiki/",
        ),
        (
            "Look, Ma, No Matrices!",
            "https://enkimute.github.io/LookMaNoMatrices/",
        ),
    ],
};
