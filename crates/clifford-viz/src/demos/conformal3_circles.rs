//! Circle Operations - 3D Conformal GA Visualization
//!
//! This demo demonstrates circle operations in 3D CGA:
//! - Creating circles from 3 points via wedge product: C = P1 ^ P2 ^ P3
//! - Visualizing the circle's plane normal
//! - Detecting when a circle degenerates into a line (collinear points)
//!
//! In CGA, a circle is a grade-3 multivector (trivector) created by
//! the outer product of three points.

use crate::common::prelude::*;
use clifford::ops::Wedge;
use clifford::specialized::conformal::dim3::{Circle, RoundPoint};
use three_d::*;

/// Scene size for visualization.
const SCENE_SIZE: f32 = 5.0;

/// Epsilon for numerical comparisons.
const EPSILON: f64 = 1e-8;

/// Number of segments to approximate a circle.
const CIRCLE_SEGMENTS: usize = 32;

/// Demo for 3D CGA Circle operations.
pub struct Conformal3CirclesDemo {
    // === Three-d 3D Resources ===
    /// Graphics context.
    context: Option<Context>,
    /// 3D camera.
    camera: Option<Camera>,
    /// Orbit control.
    control: Option<OrbitControl>,
    /// Ambient light.
    ambient_light: Option<AmbientLight>,
    /// Key light.
    key_light: Option<DirectionalLight>,
    /// World axes.
    world_axes: Option<Axes>,

    /// Point 1 mesh.
    point1_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point 2 mesh.
    point2_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point 3 mesh.
    point3_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Circle segment meshes (cylinders arranged in a ring).
    circle_segments: Vec<Gm<Mesh, ColorMaterial>>,
    /// Normal arrow mesh.
    normal_mesh: Option<Gm<Mesh, ColorMaterial>>,

    // === Point Positions ===
    /// Point 1 position.
    point1: [f32; 3],
    /// Point 2 position.
    point2: [f32; 3],
    /// Point 3 position.
    point3: [f32; 3],

    // === Display Options ===
    /// Show world axes.
    show_axes: bool,
    /// Show the circle.
    show_circle: bool,
    /// Show the normal vector.
    show_normal: bool,
}

impl Default for Conformal3CirclesDemo {
    fn default() -> Self {
        Self {
            context: None,
            camera: None,
            control: None,
            ambient_light: None,
            key_light: None,
            world_axes: None,
            point1_mesh: None,
            point2_mesh: None,
            point3_mesh: None,
            circle_segments: Vec::new(),
            normal_mesh: None,
            point1: [1.5, 0.0, 0.0],
            point2: [-0.75, 1.3, 0.0],
            point3: [-0.75, -1.3, 0.0],
            show_axes: true,
            show_circle: true,
            show_normal: true,
        }
    }
}

impl Conformal3CirclesDemo {
    /// Create a CGA circle from three points.
    fn create_circle(&self) -> Circle<f64> {
        let p1 = RoundPoint::from_euclidean(
            f64::from(self.point1[0]),
            f64::from(self.point1[1]),
            f64::from(self.point1[2]),
        );
        let p2 = RoundPoint::from_euclidean(
            f64::from(self.point2[0]),
            f64::from(self.point2[1]),
            f64::from(self.point2[2]),
        );
        let p3 = RoundPoint::from_euclidean(
            f64::from(self.point3[0]),
            f64::from(self.point3[1]),
            f64::from(self.point3[2]),
        );

        // Circle = P1 ^ P2 ^ P3
        let pair = p1.wedge(&p2);
        pair.wedge(&p3)
    }

    /// Check if the circle is actually a line (points are collinear).
    fn is_line(&self, circle: &Circle<f64>) -> bool {
        circle.is_line(EPSILON)
    }

    /// Calculate circle center and radius from three points using classical geometry.
    fn circle_params(&self) -> Option<(Vector3<f32>, Vector3<f32>, f32)> {
        let p1 = vec3(self.point1[0], self.point1[1], self.point1[2]);
        let p2 = vec3(self.point2[0], self.point2[1], self.point2[2]);
        let p3 = vec3(self.point3[0], self.point3[1], self.point3[2]);

        // Vectors from p1 to p2 and p3
        let v12 = p2 - p1;
        let v13 = p3 - p1;

        // Normal to the plane
        let normal = v12.cross(v13);
        let normal_len = normal.magnitude();
        if normal_len < 1e-6 {
            return None; // Collinear points
        }
        let normal = normal / normal_len;

        // Find circumcenter using the formula
        // Center = P1 + ((|v13|^2 * (v12 . v12xv13) x v12) + (|v12|^2 * (v13 . v13xv12) x v13)) / (2 * |v12 x v13|^2)
        let v12_sq = v12.dot(v12);
        let v13_sq = v13.dot(v13);
        let v12_cross_v13 = v12.cross(v13);
        let denom = 2.0 * v12_cross_v13.dot(v12_cross_v13);
        if denom.abs() < 1e-12 {
            return None;
        }

        let center =
            p1 + (v13_sq * v12_cross_v13.cross(v12) + v12_sq * v13.cross(v12_cross_v13)) / denom;

        // Radius is distance from center to any point
        let radius = (center - p1).magnitude();

        Some((center, normal, radius))
    }

    /// Create a transformation matrix for a cylinder (for the normal arrow).
    fn cylinder_transform(start: Vector3<f32>, end: Vector3<f32>, radius: f32) -> Mat4 {
        let diff = end - start;
        let length = diff.magnitude();
        if length < 1e-6 {
            return Mat4::from_scale(0.0);
        }

        let dir = diff / length;
        let midpoint = (start + end) * 0.5;

        // Default cylinder is along Y axis with height 2
        let y_axis = vec3(0.0, 1.0, 0.0);
        let dot = y_axis.dot(dir);

        let rotation = if dot > 0.9999 {
            Mat4::identity()
        } else if dot < -0.9999 {
            Mat4::from_angle_x(Rad(std::f32::consts::PI))
        } else {
            let axis = y_axis.cross(dir).normalize();
            let angle = dot.acos();
            Mat4::from_axis_angle(axis, Rad(angle))
        };

        let translation = Mat4::from_translation(midpoint);
        let scale = Mat4::from_nonuniform_scale(radius, length / 2.0, radius);

        translation * rotation * scale
    }

    /// Generate points on a circle given center, normal, and radius.
    fn circle_points(
        center: Vector3<f32>,
        normal: Vector3<f32>,
        radius: f32,
        num_points: usize,
    ) -> Vec<Vector3<f32>> {
        // Find two perpendicular vectors in the circle plane
        let z_axis = vec3(0.0, 0.0, 1.0);
        let u = if normal.dot(z_axis).abs() < 0.9 {
            z_axis.cross(normal).normalize()
        } else {
            vec3(1.0, 0.0, 0.0).cross(normal).normalize()
        };
        let v = normal.cross(u).normalize();

        let mut points = Vec::with_capacity(num_points);
        for i in 0..num_points {
            let angle = 2.0 * std::f32::consts::PI * (i as f32) / (num_points as f32);
            let point = center + u * (radius * angle.cos()) + v * (radius * angle.sin());
            points.push(point);
        }
        points
    }
}

impl VisualizationApp for Conformal3CirclesDemo {
    fn name(&self) -> &'static str {
        "Conformal 3D - Circle Operations"
    }

    fn update(&mut self, _dt: f32) {
        // No animation
    }

    fn render(&mut self, _ui: &mut egui::Ui) {
        // 3D rendering in render_3d
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Point Controls ===
        group_header(ui, "Points");

        ui.horizontal(|ui| {
            ui.label("P1:");
            ui.add(
                egui::DragValue::new(&mut self.point1[0])
                    .speed(0.1)
                    .prefix("x:"),
            );
            ui.add(
                egui::DragValue::new(&mut self.point1[1])
                    .speed(0.1)
                    .prefix("y:"),
            );
            ui.add(
                egui::DragValue::new(&mut self.point1[2])
                    .speed(0.1)
                    .prefix("z:"),
            );
        });

        ui.horizontal(|ui| {
            ui.label("P2:");
            ui.add(
                egui::DragValue::new(&mut self.point2[0])
                    .speed(0.1)
                    .prefix("x:"),
            );
            ui.add(
                egui::DragValue::new(&mut self.point2[1])
                    .speed(0.1)
                    .prefix("y:"),
            );
            ui.add(
                egui::DragValue::new(&mut self.point2[2])
                    .speed(0.1)
                    .prefix("z:"),
            );
        });

        ui.horizontal(|ui| {
            ui.label("P3:");
            ui.add(
                egui::DragValue::new(&mut self.point3[0])
                    .speed(0.1)
                    .prefix("x:"),
            );
            ui.add(
                egui::DragValue::new(&mut self.point3[1])
                    .speed(0.1)
                    .prefix("y:"),
            );
            ui.add(
                egui::DragValue::new(&mut self.point3[2])
                    .speed(0.1)
                    .prefix("z:"),
            );
        });

        ui.add_space(12.0);

        // === Presets ===
        group_header(ui, "Presets");

        ui.horizontal(|ui| {
            if ui.button("XY Circle").clicked() {
                self.point1 = [1.5, 0.0, 0.0];
                self.point2 = [-0.75, 1.3, 0.0];
                self.point3 = [-0.75, -1.3, 0.0];
            }
            if ui.button("XZ Circle").clicked() {
                self.point1 = [1.5, 0.0, 0.0];
                self.point2 = [-0.75, 0.0, 1.3];
                self.point3 = [-0.75, 0.0, -1.3];
            }
        });

        ui.horizontal(|ui| {
            if ui.button("YZ Circle").clicked() {
                self.point1 = [0.0, 1.5, 0.0];
                self.point2 = [0.0, -0.75, 1.3];
                self.point3 = [0.0, -0.75, -1.3];
            }
            if ui.button("Tilted").clicked() {
                self.point1 = [1.0, 0.5, 0.3];
                self.point2 = [-0.5, 1.2, -0.4];
                self.point3 = [-0.3, -0.8, 0.9];
            }
        });

        if ui.button("Collinear (Line)").clicked() {
            self.point1 = [-2.0, 0.0, 0.0];
            self.point2 = [0.0, 0.0, 0.0];
            self.point3 = [2.0, 0.0, 0.0];
        }

        ui.add_space(12.0);

        // === Circle Info ===
        group_header(ui, "Circle Info");

        let circle = self.create_circle();
        let is_line = self.is_line(&circle);

        if is_line {
            info_box(ui, "Points are collinear - degenerates to LINE");
        } else if let Some((center, normal, radius)) = self.circle_params() {
            info_box(
                ui,
                &format!(
                    "Center: ({:.2}, {:.2}, {:.2})\n\
                     Radius: {:.2}\n\
                     Normal: ({:.2}, {:.2}, {:.2})",
                    center.x, center.y, center.z, radius, normal.x, normal.y, normal.z
                ),
            );
        }

        ui.add_space(12.0);

        // === Display Options ===
        group_header(ui, "Display");
        ui.checkbox(&mut self.show_axes, "Show Axes");
        ui.checkbox(&mut self.show_circle, "Show Circle");
        ui.checkbox(&mut self.show_normal, "Show Normal");
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.label("Create a circle from 3 points: C = P1 ^ P2 ^ P3");
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(CONFORMAL3_CIRCLES_EDUCATION)
    }
}

impl VisualizationApp3D for Conformal3CirclesDemo {
    fn init_3d(&mut self, context: &Context) {
        self.context = Some(context.clone());

        // Camera
        self.camera = Some(Camera::new_perspective(
            Viewport::new_at_origo(1, 1),
            vec3(5.0, 4.0, 4.0),
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 0.0, 1.0),
            degrees(45.0),
            0.1,
            100.0,
        ));

        // Controls
        self.control = Some(OrbitControl::new(
            vec3(0.0, 0.0, 0.0),
            1.0,
            SCENE_SIZE * 4.0,
        ));

        // Axes
        self.world_axes = Some(Axes::new(context, 0.02, 2.0));

        // Lights
        self.ambient_light = Some(AmbientLight::new(context, 0.4, Srgba::WHITE));
        self.key_light = Some(DirectionalLight::new(
            context,
            0.8,
            Srgba::WHITE,
            vec3(1.0, 0.5, 1.0).normalize(),
        ));

        // Meshes
        let sphere_cpu = CpuMesh::sphere(16);
        let cylinder_cpu = CpuMesh::cylinder(8);

        // Point meshes
        self.point1_mesh = Some(Gm::new(
            Mesh::new(context, &sphere_cpu),
            ColorMaterial {
                color: Srgba::new(255, 100, 100, 255),
                ..Default::default()
            },
        ));
        self.point2_mesh = Some(Gm::new(
            Mesh::new(context, &sphere_cpu),
            ColorMaterial {
                color: Srgba::new(100, 255, 100, 255),
                ..Default::default()
            },
        ));
        self.point3_mesh = Some(Gm::new(
            Mesh::new(context, &sphere_cpu),
            ColorMaterial {
                color: Srgba::new(100, 100, 255, 255),
                ..Default::default()
            },
        ));

        // Circle segments (cylinders arranged in a ring)
        for _ in 0..CIRCLE_SEGMENTS {
            self.circle_segments.push(Gm::new(
                Mesh::new(context, &cylinder_cpu),
                ColorMaterial {
                    color: Srgba::new(255, 200, 50, 255),
                    ..Default::default()
                },
            ));
        }

        // Normal arrow mesh (cylinder)
        self.normal_mesh = Some(Gm::new(
            Mesh::new(context, &cylinder_cpu),
            ColorMaterial {
                color: Srgba::new(200, 100, 255, 255),
                ..Default::default()
            },
        ));
    }

    fn render_3d(&mut self, frame: &mut FrameInput) {
        // Pre-compute circle parameters
        let circle_params = self.circle_params();
        let is_line = self.is_line(&self.create_circle());

        // Copy display options
        let show_axes = self.show_axes;
        let show_circle = self.show_circle;
        let show_normal = self.show_normal;
        let p1 = self.point1;
        let p2 = self.point2;
        let p3 = self.point3;

        // Update camera
        let camera = self.camera.as_mut().unwrap();
        let control = self.control.as_mut().unwrap();
        camera.set_viewport(frame.viewport);
        control.handle_events(camera, &mut frame.events);

        let hide = Mat4::from_scale(0.0);
        let point_size = 0.12;

        // Update point meshes
        if let Some(mesh) = &mut self.point1_mesh {
            mesh.set_transformation(
                Mat4::from_translation(vec3(p1[0], p1[1], p1[2])) * Mat4::from_scale(point_size),
            );
        }
        if let Some(mesh) = &mut self.point2_mesh {
            mesh.set_transformation(
                Mat4::from_translation(vec3(p2[0], p2[1], p2[2])) * Mat4::from_scale(point_size),
            );
        }
        if let Some(mesh) = &mut self.point3_mesh {
            mesh.set_transformation(
                Mat4::from_translation(vec3(p3[0], p3[1], p3[2])) * Mat4::from_scale(point_size),
            );
        }

        // Update circle segments
        if show_circle && !is_line {
            if let Some((center, normal, radius)) = circle_params {
                let points = Self::circle_points(center, normal, radius, CIRCLE_SEGMENTS);
                for (i, mesh) in self.circle_segments.iter_mut().enumerate() {
                    let start = points[i];
                    let end = points[(i + 1) % CIRCLE_SEGMENTS];
                    mesh.set_transformation(Self::cylinder_transform(start, end, 0.03));
                }
            } else {
                for mesh in &mut self.circle_segments {
                    mesh.set_transformation(hide);
                }
            }
        } else {
            for mesh in &mut self.circle_segments {
                mesh.set_transformation(hide);
            }
        }

        // Update normal arrow
        if let Some(mesh) = &mut self.normal_mesh {
            if show_normal && !is_line {
                if let Some((center, normal, _radius)) = circle_params {
                    let arrow_length = 1.5;
                    let arrow_end = center + normal * arrow_length;
                    mesh.set_transformation(Self::cylinder_transform(center, arrow_end, 0.03));
                } else {
                    mesh.set_transformation(hide);
                }
            } else {
                mesh.set_transformation(hide);
            }
        }

        // Collect objects
        let mut objects: Vec<&dyn Object> = Vec::new();

        if show_axes {
            if let Some(axes) = &self.world_axes {
                objects.push(axes);
            }
        }

        if let Some(mesh) = &self.point1_mesh {
            objects.push(mesh);
        }
        if let Some(mesh) = &self.point2_mesh {
            objects.push(mesh);
        }
        if let Some(mesh) = &self.point3_mesh {
            objects.push(mesh);
        }
        for mesh in &self.circle_segments {
            objects.push(mesh);
        }
        if let Some(mesh) = &self.normal_mesh {
            objects.push(mesh);
        }

        // Lights
        let lights: Vec<&dyn Light> = vec![
            self.ambient_light.as_ref().unwrap(),
            self.key_light.as_ref().unwrap(),
        ];

        // Render
        frame.screen().render(camera, objects, &lights);
    }
}

/// Educational content for circle operations.
const CONFORMAL3_CIRCLES_EDUCATION: EducationalContent = EducationalContent {
    title: "Circles in 3D Conformal Geometric Algebra",
    overview: "In CGA, a circle is represented as a grade-3 multivector (trivector) \
               created by the outer product of three points.",
    math_background: "Circle from three points:\n\
                      C = P1 ^ P2 ^ P3\n\n\
                      Where each point P is embedded as:\n\
                      P = x*e1 + y*e2 + z*e3 + o + 0.5*|p|^2*inf\n\n\
                      The circle lies in the plane containing P1, P2, P3.\n\
                      If points are collinear, the result is a line (degenerate circle).",
    how_to_use: "Drag the point coordinates to move P1, P2, P3.\n\
                 The circle passing through all three points is shown.\n\
                 Use presets for standard configurations.\n\
                 Try 'Collinear' to see line detection.",
    key_concepts: "- Wedge product creates higher-grade elements\n\
                   - Three points determine a unique circle\n\
                   - Collinear points give a line (infinite radius circle)\n\
                   - Circle normal is perpendicular to the containing plane",
    resources: &[
        ("CGA Wiki", "https://conformalgeometricalgebra.org"),
        ("CGA Primer", "https://bivector.net/CGA.html"),
    ],
};
