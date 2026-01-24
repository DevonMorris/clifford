//! 3D Projective Plücker Lines Visualization using native three-d rendering.
//!
//! This demo demonstrates Plücker line coordinates in 3D PGA,
//! rendered with native GPU acceleration via three-d.
//!
//! Features:
//! - Line creation via two points
//! - Plücker coordinate display (direction and moment)
//! - Line-line intersection/distance
//! - Line-plane intersection
//! - Interactive point manipulation

use crate::common::prelude::*;
use clifford::ops::{Antiwedge, Wedge};
use clifford::specialized::projective::dim3::{Line, Plane, Point};
use three_d::*;

/// Demo for 3D Plücker Lines using native three-d rendering.
pub struct Projective3LinesDemo {
    /// 3D camera for rendering.
    camera: Option<Camera>,
    /// Orbit control for camera interaction.
    control: Option<OrbitControl>,
    /// World coordinate axes.
    world_axes: Option<Axes>,

    /// Line 1 cylinder mesh.
    line1_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point A sphere mesh.
    point_a_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point B sphere mesh.
    point_b_mesh: Option<Gm<Mesh, ColorMaterial>>,

    /// Line 2 cylinder mesh.
    line2_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point C sphere mesh.
    point_c_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point D sphere mesh.
    point_d_mesh: Option<Gm<Mesh, ColorMaterial>>,

    /// Reference plane mesh.
    plane_mesh: Option<Gm<Mesh, ColorMaterial>>,

    /// Intersection point sphere mesh.
    intersection_mesh: Option<Gm<Mesh, ColorMaterial>>,

    /// Ambient lighting.
    ambient_light: Option<AmbientLight>,
    /// Key directional light.
    key_light: Option<DirectionalLight>,

    /// Point A coordinates for Line 1.
    point_a: [f32; 3],
    /// Point B coordinates for Line 1.
    point_b: [f32; 3],

    /// Point C coordinates for Line 2.
    point_c: [f32; 3],
    /// Point D coordinates for Line 2.
    point_d: [f32; 3],

    /// Reference plane normal vector.
    plane_normal: [f32; 3],
    /// Reference plane distance from origin.
    plane_distance: f32,

    /// Whether to show world coordinate axes.
    show_world_axes: bool,
    /// Whether to show Line 1.
    show_line1: bool,
    /// Whether to show Line 2.
    show_line2: bool,
    /// Whether to show the reference plane.
    show_plane: bool,
    /// Whether to show intersection points.
    show_intersection: bool,
}

impl Default for Projective3LinesDemo {
    fn default() -> Self {
        Self {
            camera: None,
            control: None,
            world_axes: None,
            line1_mesh: None,
            point_a_mesh: None,
            point_b_mesh: None,
            line2_mesh: None,
            point_c_mesh: None,
            point_d_mesh: None,
            plane_mesh: None,
            intersection_mesh: None,
            ambient_light: None,
            key_light: None,
            // Line 1: from (1, 0, 0) to (0, 1, 1)
            point_a: [1.0, 0.0, 0.0],
            point_b: [0.0, 1.0, 1.0],
            // Line 2: from (0, 0, 0) to (1, 1, 0)
            point_c: [0.0, 0.0, 0.0],
            point_d: [1.0, 1.0, 0.0],
            // XY plane (z = 0)
            plane_normal: [0.0, 0.0, 1.0],
            plane_distance: 0.0,
            show_world_axes: true,
            show_line1: true,
            show_line2: true,
            show_plane: true,
            show_intersection: true,
        }
    }
}

impl Projective3LinesDemo {
    /// Get Line 1 from points A and B.
    fn line1(&self) -> Line<f64> {
        let a = Point::from_cartesian(
            f64::from(self.point_a[0]),
            f64::from(self.point_a[1]),
            f64::from(self.point_a[2]),
        );
        let b = Point::from_cartesian(
            f64::from(self.point_b[0]),
            f64::from(self.point_b[1]),
            f64::from(self.point_b[2]),
        );
        a.wedge(&b)
    }

    /// Get Line 2 from points C and D.
    fn line2(&self) -> Line<f64> {
        let c = Point::from_cartesian(
            f64::from(self.point_c[0]),
            f64::from(self.point_c[1]),
            f64::from(self.point_c[2]),
        );
        let d = Point::from_cartesian(
            f64::from(self.point_d[0]),
            f64::from(self.point_d[1]),
            f64::from(self.point_d[2]),
        );
        c.wedge(&d)
    }

    /// Get the reference plane.
    fn plane(&self) -> Plane<f64> {
        Plane::from_normal_and_distance(
            f64::from(self.plane_normal[0]),
            f64::from(self.plane_normal[1]),
            f64::from(self.plane_normal[2]),
            f64::from(self.plane_distance),
        )
    }

    /// Compute line-plane intersection.
    fn line_plane_intersection(&self, line: &Line<f64>) -> Option<Point<f64>> {
        let plane = self.plane();
        let intersection: Point<f64> = line.antiwedge(&plane);

        // Check if intersection is finite (not at infinity)
        if intersection.is_finite(1e-6) {
            Some(intersection)
        } else {
            None // Line parallel to plane
        }
    }

    /// Build a cylinder mesh transformation for a line segment.
    fn line_transform(p1: Vec3, p2: Vec3, radius: f32) -> Mat4 {
        let dir = p2 - p1;
        let length = dir.magnitude();
        if length < 1e-6 {
            return Mat4::identity();
        }
        let dir_normalized = dir / length;

        // Scale: length in X (cylinder along X), radius in Y/Z
        let scale = Mat4::from_nonuniform_scale(length, radius, radius);

        // Rotate from X-axis to target direction
        let x_axis = vec3(1.0, 0.0, 0.0);
        let dot = x_axis.dot(dir_normalized);

        let rotation = if dot > 0.9999 {
            Mat4::identity()
        } else if dot < -0.9999 {
            Mat4::from_angle_y(Rad(std::f32::consts::PI))
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

        // Translate to start point (cylinder goes from 0 to length along X after scale)
        let translation = Mat4::from_translation(p1);

        translation * rotation * scale
    }
}

impl VisualizationApp for Projective3LinesDemo {
    fn name(&self) -> &'static str {
        "3D PGA Plucker Lines"
    }

    fn update(&mut self, _dt: f32) {
        // No animation for this demo
    }

    fn render(&mut self, _ui: &mut egui::Ui) {
        // 3D rendering is handled by render_3d()
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Line 1 Controls ===
        group_header(ui, "Line 1 (Yellow)");

        ui.horizontal(|ui| {
            ui.label("Point A:");
            ui.add(
                egui::DragValue::new(&mut self.point_a[0])
                    .speed(0.1)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_a[1])
                    .speed(0.1)
                    .prefix("y: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_a[2])
                    .speed(0.1)
                    .prefix("z: "),
            );
        });

        ui.horizontal(|ui| {
            ui.label("Point B:");
            ui.add(
                egui::DragValue::new(&mut self.point_b[0])
                    .speed(0.1)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_b[1])
                    .speed(0.1)
                    .prefix("y: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_b[2])
                    .speed(0.1)
                    .prefix("z: "),
            );
        });

        // Display Plücker coordinates for Line 1
        let line1 = self.line1();
        let dir1 = line1.direction();
        let mom1 = line1.moment();

        info_box(
            ui,
            &format!(
                "Direction: ({:.3}, {:.3}, {:.3})\n\
                 Moment:    ({:.3}, {:.3}, {:.3})\n\
                 d.m = {:.6} {}",
                dir1.x(),
                dir1.y(),
                dir1.z(),
                mom1.x(),
                mom1.y(),
                mom1.z(),
                line1.plucker_residual(),
                if line1.satisfies_plucker_condition(1e-6) {
                    "(valid)"
                } else {
                    "(INVALID)"
                }
            ),
        );

        // === Line 2 Controls ===
        section_separator(ui, Some("Line 2 (Cyan)"));

        ui.horizontal(|ui| {
            ui.label("Point C:");
            ui.add(
                egui::DragValue::new(&mut self.point_c[0])
                    .speed(0.1)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_c[1])
                    .speed(0.1)
                    .prefix("y: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_c[2])
                    .speed(0.1)
                    .prefix("z: "),
            );
        });

        ui.horizontal(|ui| {
            ui.label("Point D:");
            ui.add(
                egui::DragValue::new(&mut self.point_d[0])
                    .speed(0.1)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_d[1])
                    .speed(0.1)
                    .prefix("y: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.point_d[2])
                    .speed(0.1)
                    .prefix("z: "),
            );
        });

        // Display Plücker coordinates for Line 2
        let line2 = self.line2();
        let dir2 = line2.direction();
        let mom2 = line2.moment();

        info_box(
            ui,
            &format!(
                "Direction: ({:.3}, {:.3}, {:.3})\n\
                 Moment:    ({:.3}, {:.3}, {:.3})\n\
                 d.m = {:.6} {}",
                dir2.x(),
                dir2.y(),
                dir2.z(),
                mom2.x(),
                mom2.y(),
                mom2.z(),
                line2.plucker_residual(),
                if line2.satisfies_plucker_condition(1e-6) {
                    "(valid)"
                } else {
                    "(INVALID)"
                }
            ),
        );

        // === Line-Line Analysis ===
        section_separator(ui, Some("Line-Line Analysis"));

        let distance = line1.distance(&line2);
        let angle = line1.angle(&line2).to_degrees();
        let parallel = line1.is_parallel(&line2);

        info_box(
            ui,
            &format!(
                "Distance: {:.4}\n\
                 Angle: {:.2} deg\n\
                 Parallel: {}",
                distance,
                angle,
                if parallel { "Yes" } else { "No" }
            ),
        );

        // === Plane Controls ===
        section_separator(ui, Some("Reference Plane"));

        ui.horizontal(|ui| {
            ui.label("Normal:");
            ui.add(
                egui::DragValue::new(&mut self.plane_normal[0])
                    .speed(0.1)
                    .prefix("x: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.plane_normal[1])
                    .speed(0.1)
                    .prefix("y: "),
            );
            ui.add(
                egui::DragValue::new(&mut self.plane_normal[2])
                    .speed(0.1)
                    .prefix("z: "),
            );
        });

        ui.horizontal(|ui| {
            ui.label("Distance:");
            ui.add(egui::DragValue::new(&mut self.plane_distance).speed(0.1));
        });

        ui.horizontal(|ui| {
            if ui.button("XY").clicked() {
                self.plane_normal = [0.0, 0.0, 1.0];
                self.plane_distance = 0.0;
            }
            if ui.button("XZ").clicked() {
                self.plane_normal = [0.0, 1.0, 0.0];
                self.plane_distance = 0.0;
            }
            if ui.button("YZ").clicked() {
                self.plane_normal = [1.0, 0.0, 0.0];
                self.plane_distance = 0.0;
            }
        });

        // Line-plane intersections
        if let Some(p) = self.line_plane_intersection(&line1) {
            if let Some((x, y, z)) = p.to_cartesian() {
                info_box(
                    ui,
                    &format!("L1-Plane intersection: ({:.3}, {:.3}, {:.3})", x, y, z),
                );
            }
        } else {
            info_box(ui, "L1-Plane: parallel (no intersection)");
        }

        if let Some(p) = self.line_plane_intersection(&line2) {
            if let Some((x, y, z)) = p.to_cartesian() {
                info_box(
                    ui,
                    &format!("L2-Plane intersection: ({:.3}, {:.3}, {:.3})", x, y, z),
                );
            }
        } else {
            info_box(ui, "L2-Plane: parallel (no intersection)");
        }

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_world_axes, "World Axes");
            ui.checkbox(&mut self.show_line1, "Line 1");
            ui.checkbox(&mut self.show_line2, "Line 2");
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_plane, "Plane");
            ui.checkbox(&mut self.show_intersection, "Intersections");
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
        ui.label("Plucker line coordinates in 3D PGA");
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(PROJECTIVE3_LINES_EDUCATION)
    }
}

impl VisualizationApp3D for Projective3LinesDemo {
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

        // World axes
        self.world_axes = Some(Axes::new(context, 0.02, 2.5));

        // Create cylinder mesh for lines
        let cylinder = CpuMesh::cylinder(16);

        // Line 1 (yellow)
        self.line1_mesh = Some(Gm::new(
            Mesh::new(context, &cylinder),
            ColorMaterial {
                color: Srgba::new(255, 220, 50, 255),
                ..Default::default()
            },
        ));

        // Line 2 (cyan)
        self.line2_mesh = Some(Gm::new(
            Mesh::new(context, &cylinder),
            ColorMaterial {
                color: Srgba::new(50, 220, 255, 255),
                ..Default::default()
            },
        ));

        // Point spheres
        let sphere = CpuMesh::sphere(16);

        // Points A, B (line 1 endpoints - orange)
        self.point_a_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(255, 150, 50, 255),
                ..Default::default()
            },
        ));
        self.point_b_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(255, 150, 50, 255),
                ..Default::default()
            },
        ));

        // Points C, D (line 2 endpoints - teal)
        self.point_c_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(50, 180, 200, 255),
                ..Default::default()
            },
        ));
        self.point_d_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(50, 180, 200, 255),
                ..Default::default()
            },
        ));

        // Reference plane (semi-transparent gray)
        let quad = CpuMesh::square();
        self.plane_mesh = Some(Gm::new(
            Mesh::new(context, &quad),
            ColorMaterial {
                color: Srgba::new(150, 150, 150, 100),
                is_transparent: true,
                ..Default::default()
            },
        ));

        // Intersection point (green)
        self.intersection_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(50, 255, 50, 255),
                ..Default::default()
            },
        ));

        // Lighting
        self.ambient_light = Some(AmbientLight::new(context, 0.4, Srgba::WHITE));
        self.key_light = Some(DirectionalLight::new(
            context,
            0.8,
            Srgba::WHITE,
            vec3(-1.0, -1.0, -1.0).normalize(),
        ));
    }

    fn render_3d(&mut self, frame: &mut FrameInput) {
        // Compute intersection FIRST before any mutable borrows
        let intersection_pos: Option<(f32, f32, f32)> = {
            let line1 = self.line1();
            self.line_plane_intersection(&line1)
                .and_then(|p| p.to_cartesian())
                .map(|(x, y, z)| (x as f32, y as f32, z as f32))
        };

        let camera = self.camera.as_mut().unwrap();
        let control = self.control.as_mut().unwrap();

        camera.set_viewport(frame.viewport);
        control.handle_events(camera, &mut frame.events);

        // Update line 1 mesh
        let p_a = vec3(self.point_a[0], self.point_a[1], self.point_a[2]);
        let p_b = vec3(self.point_b[0], self.point_b[1], self.point_b[2]);

        // Extend line beyond endpoints for visualization
        let dir1 = (p_b - p_a).normalize();
        let line1_start = p_a - dir1 * 3.0;
        let line1_end = p_b + dir1 * 3.0;

        if let Some(mesh) = &mut self.line1_mesh {
            mesh.set_transformation(Self::line_transform(line1_start, line1_end, 0.02));
        }

        // Update point meshes for line 1
        let point_scale = 0.08;
        if let Some(mesh) = &mut self.point_a_mesh {
            mesh.set_transformation(Mat4::from_translation(p_a) * Mat4::from_scale(point_scale));
        }
        if let Some(mesh) = &mut self.point_b_mesh {
            mesh.set_transformation(Mat4::from_translation(p_b) * Mat4::from_scale(point_scale));
        }

        // Update line 2 mesh
        let p_c = vec3(self.point_c[0], self.point_c[1], self.point_c[2]);
        let p_d = vec3(self.point_d[0], self.point_d[1], self.point_d[2]);

        let dir2 = (p_d - p_c).normalize();
        let line2_start = p_c - dir2 * 3.0;
        let line2_end = p_d + dir2 * 3.0;

        if let Some(mesh) = &mut self.line2_mesh {
            mesh.set_transformation(Self::line_transform(line2_start, line2_end, 0.02));
        }

        // Update point meshes for line 2
        if let Some(mesh) = &mut self.point_c_mesh {
            mesh.set_transformation(Mat4::from_translation(p_c) * Mat4::from_scale(point_scale));
        }
        if let Some(mesh) = &mut self.point_d_mesh {
            mesh.set_transformation(Mat4::from_translation(p_d) * Mat4::from_scale(point_scale));
        }

        // Update plane mesh
        if let Some(mesh) = &mut self.plane_mesh {
            let n = vec3(
                self.plane_normal[0],
                self.plane_normal[1],
                self.plane_normal[2],
            );
            let n_len = n.magnitude();
            if n_len > 1e-6 {
                let n_normalized = n / n_len;

                // Rotate from Z-up to normal direction
                let z_axis = vec3(0.0, 0.0, 1.0);
                let dot = z_axis.dot(n_normalized);

                let rotation = if dot > 0.9999 {
                    Mat4::identity()
                } else if dot < -0.9999 {
                    Mat4::from_angle_x(Rad(std::f32::consts::PI))
                } else {
                    let rot_axis = z_axis.cross(n_normalized);
                    let rot_axis_len = rot_axis.magnitude();
                    if rot_axis_len < 1e-6 {
                        Mat4::identity()
                    } else {
                        let rot_axis_normalized = rot_axis / rot_axis_len;
                        let angle = dot.clamp(-1.0, 1.0).acos();
                        Mat4::from_axis_angle(rot_axis_normalized, Rad(angle))
                    }
                };

                let translation = Mat4::from_translation(n_normalized * self.plane_distance);
                let scale = Mat4::from_scale(5.0);

                mesh.set_transformation(translation * rotation * scale);
            }
        }

        // Update intersection point (show line1-plane intersection)
        if let Some(mesh) = &mut self.intersection_mesh {
            if let Some((x, y, z)) = intersection_pos {
                mesh.set_transformation(
                    Mat4::from_translation(vec3(x, y, z)) * Mat4::from_scale(0.1),
                );
            } else {
                // Hide by scaling to zero
                mesh.set_transformation(Mat4::from_scale(0.0));
            }
        }

        // Collect objects to render
        let mut objects: Vec<&dyn Object> = Vec::new();

        if self.show_world_axes {
            if let Some(axes) = &self.world_axes {
                objects.push(axes);
            }
        }

        if self.show_line1 {
            if let Some(mesh) = &self.line1_mesh {
                objects.push(mesh);
            }
            if let Some(mesh) = &self.point_a_mesh {
                objects.push(mesh);
            }
            if let Some(mesh) = &self.point_b_mesh {
                objects.push(mesh);
            }
        }

        if self.show_line2 {
            if let Some(mesh) = &self.line2_mesh {
                objects.push(mesh);
            }
            if let Some(mesh) = &self.point_c_mesh {
                objects.push(mesh);
            }
            if let Some(mesh) = &self.point_d_mesh {
                objects.push(mesh);
            }
        }

        if self.show_plane {
            if let Some(mesh) = &self.plane_mesh {
                objects.push(mesh);
            }
        }

        if self.show_intersection {
            if let Some(mesh) = &self.intersection_mesh {
                objects.push(mesh);
            }
        }

        // Collect lights
        let lights: Vec<&dyn Light> = vec![
            self.ambient_light.as_ref().unwrap(),
            self.key_light.as_ref().unwrap(),
        ];

        // Render
        frame.screen().render(camera, objects, &lights);
    }
}

/// Educational content for the Plücker Lines visualization.
const PROJECTIVE3_LINES_EDUCATION: EducationalContent = EducationalContent {
    title: "Plucker Line Coordinates in 3D PGA",

    overview: "\
Plucker coordinates are a 6-component representation of 3D lines. In PGA, lines are \
grade-2 elements (bivectors) with the same 6 components, making the connection \
between classical Plucker coordinates and geometric algebra explicit.

A line is defined by joining two points: L = P1 ^ P2 (wedge product).",

    math_background: "\
Plucker coordinates consist of two 3-vectors:
  - Direction d: the line's direction vector
  - Moment m: encodes the line's position (m = p x d for point p on line)

Plucker constraint: d . m = 0 (valid lines satisfy this)

In PGA bivector representation:
  - Bulk part (e23, e31, e12): direction components
  - Weight part (e01, e02, e03): moment components

Line-line distance:
  For non-parallel lines L1, L2:
    dist = |d1 . m2 + d2 . m1| / |d1 x d2|

Line-plane intersection:
  L meet Pi = L v Pi (antiwedge product)
  Result is a point (finite or at infinity if parallel)",

    how_to_use: "\
- Drag points A, B to define Line 1 (yellow)
- Drag points C, D to define Line 2 (cyan)
- Observe Plucker coordinates update in real-time
- The d.m value shows Plucker constraint satisfaction
- Adjust the reference plane to see line-plane intersections
- Green sphere shows intersection point",

    key_concepts: "\
- Lines in 3D require 4 DOF (2 points define a line, minus 2 for scaling)
- Plucker coords have 6 components but 1 constraint: d.m = 0
- Join (wedge): Point ^ Point = Line
- Meet (antiwedge): Line v Plane = Point
- Skew lines have nonzero distance, coplanar lines can intersect
- Parallel lines: cross product of directions is zero",

    resources: &[
        (
            "Rigid GA Wiki - Lines",
            "https://rigidgeometricalgebra.org/wiki/index.php?title=Line",
        ),
        (
            "Plucker Coordinates (Wikipedia)",
            "https://en.wikipedia.org/wiki/Pl%C3%BCcker_coordinates",
        ),
        (
            "PGA for Computer Scientists",
            "https://bivector.net/PGA4CS.pdf",
        ),
    ],
};
