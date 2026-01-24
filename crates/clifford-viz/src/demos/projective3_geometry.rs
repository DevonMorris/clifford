//! 3D Point-Line-Plane Operations Visualization using native three-d rendering.
//!
//! This demo demonstrates join and meet operations in 3D PGA:
//! - Join (wedge): Point ^ Point = Line, Point ^ Line = Plane
//! - Meet (antiwedge): Plane v Plane = Line, Line v Plane = Point
//!
//! The user can interactively create points and see the resulting
//! geometric constructions.

use crate::common::prelude::*;
use clifford::ops::{Antiwedge, Wedge};
use clifford::specialized::projective::dim3::{Line, Plane, Point};
use three_d::*;

/// Geometric element that can be displayed.
#[derive(Clone, Copy, PartialEq, Default)]
enum GeometryMode {
    /// Show two points and their joining line.
    #[default]
    TwoPointsLine,
    /// Show point and line, and their joining plane.
    PointLinePlane,
    /// Show two planes and their meeting line.
    TwoPlanesLine,
    /// Show line and plane, and their meeting point.
    LinePlanePoint,
}

/// Demo for 3D Point-Line-Plane operations using native three-d rendering.
pub struct Projective3GeometryDemo {
    /// 3D camera for rendering.
    camera: Option<Camera>,
    /// Orbit control for camera interaction.
    control: Option<OrbitControl>,
    /// World coordinate axes.
    world_axes: Option<Axes>,

    /// Point 1 mesh (sphere).
    point1_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point 2 mesh (sphere).
    point2_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Point 3 mesh (sphere) - for results.
    point3_mesh: Option<Gm<Mesh, ColorMaterial>>,

    /// Line mesh (cylinder).
    line_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Result line mesh.
    result_line_mesh: Option<Gm<Mesh, ColorMaterial>>,

    /// Plane 1 mesh (quad).
    plane1_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Plane 2 mesh (quad).
    plane2_mesh: Option<Gm<Mesh, ColorMaterial>>,
    /// Result plane mesh.
    result_plane_mesh: Option<Gm<Mesh, ColorMaterial>>,

    /// Ambient lighting.
    ambient_light: Option<AmbientLight>,
    /// Key directional light.
    key_light: Option<DirectionalLight>,

    /// Current geometry mode.
    mode: GeometryMode,

    /// Point 1 coordinates.
    point1: [f32; 3],
    /// Point 2 coordinates.
    point2: [f32; 3],

    /// Line direction (for PointLinePlane and LinePlanePoint modes).
    line_direction: [f32; 3],
    /// Line point (for line definition).
    line_point: [f32; 3],

    /// Plane 1 normal.
    plane1_normal: [f32; 3],
    /// Plane 1 distance.
    plane1_distance: f32,

    /// Plane 2 normal.
    plane2_normal: [f32; 3],
    /// Plane 2 distance.
    plane2_distance: f32,

    /// Whether to show world axes.
    show_world_axes: bool,
}

impl Default for Projective3GeometryDemo {
    fn default() -> Self {
        Self {
            camera: None,
            control: None,
            world_axes: None,
            point1_mesh: None,
            point2_mesh: None,
            point3_mesh: None,
            line_mesh: None,
            result_line_mesh: None,
            plane1_mesh: None,
            plane2_mesh: None,
            result_plane_mesh: None,
            ambient_light: None,
            key_light: None,
            mode: GeometryMode::TwoPointsLine,
            point1: [1.0, 0.0, 0.0],
            point2: [0.0, 1.0, 1.0],
            line_direction: [1.0, 1.0, 0.0],
            line_point: [0.0, 0.0, 0.0],
            plane1_normal: [0.0, 0.0, 1.0],
            plane1_distance: 0.0,
            plane2_normal: [1.0, 0.0, 0.0],
            plane2_distance: 0.0,
            show_world_axes: true,
        }
    }
}

impl Projective3GeometryDemo {
    /// Create a PGA point from coordinates.
    fn make_point(&self, coords: &[f32; 3]) -> Point<f64> {
        Point::from_cartesian(
            f64::from(coords[0]),
            f64::from(coords[1]),
            f64::from(coords[2]),
        )
    }

    /// Create a PGA line from point and direction.
    fn make_line(&self) -> Line<f64> {
        let pt = self.make_point(&self.line_point);
        let dir = Point::ideal(
            f64::from(self.line_direction[0]),
            f64::from(self.line_direction[1]),
            f64::from(self.line_direction[2]),
        );
        pt.wedge(&dir)
    }

    /// Create a PGA plane from normal and distance.
    fn make_plane1(&self) -> Plane<f64> {
        Plane::from_normal_and_distance(
            f64::from(self.plane1_normal[0]),
            f64::from(self.plane1_normal[1]),
            f64::from(self.plane1_normal[2]),
            f64::from(self.plane1_distance),
        )
    }

    /// Create a PGA plane from normal and distance.
    fn make_plane2(&self) -> Plane<f64> {
        Plane::from_normal_and_distance(
            f64::from(self.plane2_normal[0]),
            f64::from(self.plane2_normal[1]),
            f64::from(self.plane2_normal[2]),
            f64::from(self.plane2_distance),
        )
    }

    /// Build cylinder transformation for a line.
    fn line_transform(p1: Vec3, p2: Vec3, radius: f32) -> Mat4 {
        let dir = p2 - p1;
        let length = dir.magnitude();
        if length < 1e-6 {
            return Mat4::identity();
        }
        let dir_normalized = dir / length;

        let scale = Mat4::from_nonuniform_scale(length, radius, radius);

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

        let translation = Mat4::from_translation(p1);
        translation * rotation * scale
    }

    /// Build plane transformation.
    fn plane_transform(normal: Vec3, distance: f32, size: f32) -> Mat4 {
        let n_len = normal.magnitude();
        if n_len < 1e-6 {
            return Mat4::identity();
        }
        let n_normalized = normal / n_len;

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

        let translation = Mat4::from_translation(n_normalized * distance);
        let scale = Mat4::from_scale(size);

        translation * rotation * scale
    }
}

impl VisualizationApp for Projective3GeometryDemo {
    fn name(&self) -> &'static str {
        "3D PGA Point-Line-Plane"
    }

    fn update(&mut self, _dt: f32) {
        // No animation
    }

    fn render(&mut self, _ui: &mut egui::Ui) {
        // 3D rendering is handled by render_3d()
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Mode Selection ===
        group_header(ui, "Operation Mode");

        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.mode, GeometryMode::TwoPointsLine, "P^P=L");
            ui.selectable_value(&mut self.mode, GeometryMode::PointLinePlane, "P^L=Pi");
        });
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.mode, GeometryMode::TwoPlanesLine, "Pi^Pi=L");
            ui.selectable_value(&mut self.mode, GeometryMode::LinePlanePoint, "L^Pi=P");
        });

        // === Mode-specific controls ===
        match self.mode {
            GeometryMode::TwoPointsLine => {
                section_separator(ui, Some("Point ^ Point = Line"));

                ui.horizontal(|ui| {
                    ui.label("Point 1:");
                    ui.add(
                        egui::DragValue::new(&mut self.point1[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.point1[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.point1[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });

                ui.horizontal(|ui| {
                    ui.label("Point 2:");
                    ui.add(
                        egui::DragValue::new(&mut self.point2[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.point2[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.point2[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });

                let p1 = self.make_point(&self.point1);
                let p2 = self.make_point(&self.point2);
                let line: Line<f64> = p1.wedge(&p2);
                let dir = line.direction();
                let mom = line.moment();

                info_box(
                    ui,
                    &format!(
                        "Result Line:\n  Dir: ({:.3}, {:.3}, {:.3})\n  Mom: ({:.3}, {:.3}, {:.3})",
                        dir.x(),
                        dir.y(),
                        dir.z(),
                        mom.x(),
                        mom.y(),
                        mom.z()
                    ),
                );
            }

            GeometryMode::PointLinePlane => {
                section_separator(ui, Some("Point ^ Line = Plane"));

                ui.horizontal(|ui| {
                    ui.label("Point:");
                    ui.add(
                        egui::DragValue::new(&mut self.point1[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.point1[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.point1[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });

                ui.horizontal(|ui| {
                    ui.label("Line pt:");
                    ui.add(
                        egui::DragValue::new(&mut self.line_point[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_point[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_point[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });

                ui.horizontal(|ui| {
                    ui.label("Line dir:");
                    ui.add(
                        egui::DragValue::new(&mut self.line_direction[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_direction[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_direction[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });

                let pt = self.make_point(&self.point1);
                let line = self.make_line();
                let plane: Plane<f64> = pt.wedge(&line);
                let normal = plane.normal();

                info_box(
                    ui,
                    &format!(
                        "Result Plane:\n  Normal: ({:.3}, {:.3}, {:.3})\n  Distance: {:.3}",
                        normal.x(),
                        normal.y(),
                        normal.z(),
                        plane.distance_from_origin()
                    ),
                );
            }

            GeometryMode::TwoPlanesLine => {
                section_separator(ui, Some("Plane v Plane = Line"));

                ui.horizontal(|ui| {
                    ui.label("Plane 1 n:");
                    ui.add(
                        egui::DragValue::new(&mut self.plane1_normal[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.plane1_normal[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.plane1_normal[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Plane 1 d:");
                    ui.add(egui::DragValue::new(&mut self.plane1_distance).speed(0.1));
                });

                ui.horizontal(|ui| {
                    ui.label("Plane 2 n:");
                    ui.add(
                        egui::DragValue::new(&mut self.plane2_normal[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.plane2_normal[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.plane2_normal[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Plane 2 d:");
                    ui.add(egui::DragValue::new(&mut self.plane2_distance).speed(0.1));
                });

                let plane1 = self.make_plane1();
                let plane2 = self.make_plane2();
                let line: Line<f64> = plane1.antiwedge(&plane2);
                let dir = line.direction();

                info_box(
                    ui,
                    &format!(
                        "Result Line:\n  Direction: ({:.3}, {:.3}, {:.3})",
                        dir.x(),
                        dir.y(),
                        dir.z()
                    ),
                );
            }

            GeometryMode::LinePlanePoint => {
                section_separator(ui, Some("Line v Plane = Point"));

                ui.horizontal(|ui| {
                    ui.label("Line pt:");
                    ui.add(
                        egui::DragValue::new(&mut self.line_point[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_point[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_point[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });

                ui.horizontal(|ui| {
                    ui.label("Line dir:");
                    ui.add(
                        egui::DragValue::new(&mut self.line_direction[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_direction[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.line_direction[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });

                ui.horizontal(|ui| {
                    ui.label("Plane n:");
                    ui.add(
                        egui::DragValue::new(&mut self.plane1_normal[0])
                            .speed(0.1)
                            .prefix("x: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.plane1_normal[1])
                            .speed(0.1)
                            .prefix("y: "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut self.plane1_normal[2])
                            .speed(0.1)
                            .prefix("z: "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Plane d:");
                    ui.add(egui::DragValue::new(&mut self.plane1_distance).speed(0.1));
                });

                let line = self.make_line();
                let plane = self.make_plane1();
                let point: Point<f64> = line.antiwedge(&plane);

                if let Some((x, y, z)) = point.to_cartesian() {
                    info_box(ui, &format!("Result Point: ({:.3}, {:.3}, {:.3})", x, y, z));
                } else {
                    info_box(ui, "Result: Point at infinity (parallel)");
                }
            }
        }

        // === Display Options ===
        section_separator(ui, Some("Display"));
        ui.checkbox(&mut self.show_world_axes, "World Axes");

        // === Camera Info ===
        section_separator(ui, Some("Camera"));
        ui.label(
            egui::RichText::new("Drag to orbit, scroll to zoom")
                .small()
                .weak(),
        );
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.label("Join (^) and Meet (v) operations in 3D PGA");
    }

    fn educational_content(&self) -> Option<EducationalContent> {
        Some(PROJECTIVE3_GEOMETRY_EDUCATION)
    }
}

impl VisualizationApp3D for Projective3GeometryDemo {
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

        let sphere = CpuMesh::sphere(16);
        let cylinder = CpuMesh::cylinder(16);
        let quad = CpuMesh::square();

        // Point meshes (orange, red, green for result)
        self.point1_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(255, 150, 50, 255),
                ..Default::default()
            },
        ));
        self.point2_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(255, 100, 50, 255),
                ..Default::default()
            },
        ));
        self.point3_mesh = Some(Gm::new(
            Mesh::new(context, &sphere),
            ColorMaterial {
                color: Srgba::new(50, 255, 100, 255),
                ..Default::default()
            },
        ));

        // Line meshes (yellow input, green result)
        self.line_mesh = Some(Gm::new(
            Mesh::new(context, &cylinder),
            ColorMaterial {
                color: Srgba::new(255, 220, 50, 255),
                ..Default::default()
            },
        ));
        self.result_line_mesh = Some(Gm::new(
            Mesh::new(context, &cylinder),
            ColorMaterial {
                color: Srgba::new(50, 255, 100, 255),
                ..Default::default()
            },
        ));

        // Plane meshes (blue/cyan transparent)
        self.plane1_mesh = Some(Gm::new(
            Mesh::new(context, &quad),
            ColorMaterial {
                color: Srgba::new(100, 150, 255, 100),
                is_transparent: true,
                ..Default::default()
            },
        ));
        self.plane2_mesh = Some(Gm::new(
            Mesh::new(context, &quad),
            ColorMaterial {
                color: Srgba::new(100, 255, 200, 100),
                is_transparent: true,
                ..Default::default()
            },
        ));
        self.result_plane_mesh = Some(Gm::new(
            Mesh::new(context, &quad),
            ColorMaterial {
                color: Srgba::new(50, 255, 100, 80),
                is_transparent: true,
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
        // Pre-compute all PGA results before any mutable borrows
        let mode = self.mode;
        let point1_coords = self.point1;
        let point2_coords = self.point2;
        let line_point_coords = self.line_point;
        let line_dir_coords = self.line_direction;
        let plane1_normal_coords = self.plane1_normal;
        let plane1_dist = self.plane1_distance;
        let plane2_normal_coords = self.plane2_normal;
        let plane2_dist = self.plane2_distance;

        // Pre-compute PGA operations
        let result_plane_normal_dist: Option<(Vec3, f32)> = match mode {
            GeometryMode::PointLinePlane => {
                let pga_pt = self.make_point(&point1_coords);
                let pga_line = self.make_line();
                let plane: Plane<f64> = pga_pt.wedge(&pga_line);
                let normal = plane.normal();
                Some((
                    vec3(normal.x() as f32, normal.y() as f32, normal.z() as f32),
                    plane.distance_from_origin() as f32,
                ))
            }
            _ => None,
        };

        let result_line_from_planes: Option<(Vec3, Vec3)> = match mode {
            GeometryMode::TwoPlanesLine => {
                let plane1 = self.make_plane1();
                let plane2 = self.make_plane2();
                let line: Line<f64> = plane1.antiwedge(&plane2);
                let dir = line.direction();
                if dir.norm() > 1e-6 {
                    let closest = line.closest_point(&Point::origin());
                    closest.to_cartesian().map(|(cx, cy, cz)| {
                        (
                            vec3(cx as f32, cy as f32, cz as f32),
                            vec3(dir.x() as f32, dir.y() as f32, dir.z() as f32),
                        )
                    })
                } else {
                    None
                }
            }
            _ => None,
        };

        let result_point_from_line_plane: Option<Vec3> = match mode {
            GeometryMode::LinePlanePoint => {
                let line = self.make_line();
                let plane = self.make_plane1();
                let point: Point<f64> = line.antiwedge(&plane);
                point
                    .to_cartesian()
                    .map(|(x, y, z)| vec3(x as f32, y as f32, z as f32))
            }
            _ => None,
        };

        let camera = self.camera.as_mut().unwrap();
        let control = self.control.as_mut().unwrap();

        camera.set_viewport(frame.viewport);
        control.handle_events(camera, &mut frame.events);

        let point_scale = 0.1;
        let line_radius = 0.03;
        let plane_size = 4.0;

        // Hide all meshes initially, then show based on mode
        let hide_transform = Mat4::from_scale(0.0);

        // Reset all to hidden
        if let Some(m) = &mut self.point1_mesh {
            m.set_transformation(hide_transform);
        }
        if let Some(m) = &mut self.point2_mesh {
            m.set_transformation(hide_transform);
        }
        if let Some(m) = &mut self.point3_mesh {
            m.set_transformation(hide_transform);
        }
        if let Some(m) = &mut self.line_mesh {
            m.set_transformation(hide_transform);
        }
        if let Some(m) = &mut self.result_line_mesh {
            m.set_transformation(hide_transform);
        }
        if let Some(m) = &mut self.plane1_mesh {
            m.set_transformation(hide_transform);
        }
        if let Some(m) = &mut self.plane2_mesh {
            m.set_transformation(hide_transform);
        }
        if let Some(m) = &mut self.result_plane_mesh {
            m.set_transformation(hide_transform);
        }

        match mode {
            GeometryMode::TwoPointsLine => {
                // Show two points and the resulting line
                let p1 = vec3(point1_coords[0], point1_coords[1], point1_coords[2]);
                let p2 = vec3(point2_coords[0], point2_coords[1], point2_coords[2]);

                if let Some(m) = &mut self.point1_mesh {
                    m.set_transformation(
                        Mat4::from_translation(p1) * Mat4::from_scale(point_scale),
                    );
                }
                if let Some(m) = &mut self.point2_mesh {
                    m.set_transformation(
                        Mat4::from_translation(p2) * Mat4::from_scale(point_scale),
                    );
                }

                // Extend line for visualization
                let dir = (p2 - p1).normalize();
                let line_start = p1 - dir * 3.0;
                let line_end = p2 + dir * 3.0;

                if let Some(m) = &mut self.result_line_mesh {
                    m.set_transformation(Self::line_transform(line_start, line_end, line_radius));
                }
            }

            GeometryMode::PointLinePlane => {
                // Show point, input line, and resulting plane
                let pt = vec3(point1_coords[0], point1_coords[1], point1_coords[2]);
                let lp = vec3(
                    line_point_coords[0],
                    line_point_coords[1],
                    line_point_coords[2],
                );
                let ld = vec3(line_dir_coords[0], line_dir_coords[1], line_dir_coords[2]);

                if let Some(m) = &mut self.point1_mesh {
                    m.set_transformation(
                        Mat4::from_translation(pt) * Mat4::from_scale(point_scale),
                    );
                }

                // Show input line
                let ld_norm = ld.normalize();
                let line_start = lp - ld_norm * 3.0;
                let line_end = lp + ld_norm * 3.0;
                if let Some(m) = &mut self.line_mesh {
                    m.set_transformation(Self::line_transform(line_start, line_end, line_radius));
                }

                // Show resulting plane (use precomputed value)
                if let Some((n, d)) = result_plane_normal_dist {
                    if let Some(m) = &mut self.result_plane_mesh {
                        m.set_transformation(Self::plane_transform(n, d, plane_size));
                    }
                }
            }

            GeometryMode::TwoPlanesLine => {
                // Show two planes and resulting line
                let n1 = vec3(
                    plane1_normal_coords[0],
                    plane1_normal_coords[1],
                    plane1_normal_coords[2],
                );
                let n2 = vec3(
                    plane2_normal_coords[0],
                    plane2_normal_coords[1],
                    plane2_normal_coords[2],
                );

                if let Some(m) = &mut self.plane1_mesh {
                    m.set_transformation(Self::plane_transform(n1, plane1_dist, plane_size));
                }
                if let Some(m) = &mut self.plane2_mesh {
                    m.set_transformation(Self::plane_transform(n2, plane2_dist, plane_size));
                }

                // Show intersection line (use precomputed value)
                if let Some((lp, ld)) = result_line_from_planes {
                    let ld_norm = ld.normalize();
                    let line_start = lp - ld_norm * 5.0;
                    let line_end = lp + ld_norm * 5.0;

                    if let Some(m) = &mut self.result_line_mesh {
                        m.set_transformation(Self::line_transform(
                            line_start,
                            line_end,
                            line_radius,
                        ));
                    }
                }
            }

            GeometryMode::LinePlanePoint => {
                // Show line, plane, and resulting intersection point
                let lp = vec3(
                    line_point_coords[0],
                    line_point_coords[1],
                    line_point_coords[2],
                );
                let ld = vec3(line_dir_coords[0], line_dir_coords[1], line_dir_coords[2]);
                let n = vec3(
                    plane1_normal_coords[0],
                    plane1_normal_coords[1],
                    plane1_normal_coords[2],
                );

                // Show input line
                let ld_norm = ld.normalize();
                let line_start = lp - ld_norm * 3.0;
                let line_end = lp + ld_norm * 3.0;
                if let Some(m) = &mut self.line_mesh {
                    m.set_transformation(Self::line_transform(line_start, line_end, line_radius));
                }

                // Show input plane
                if let Some(m) = &mut self.plane1_mesh {
                    m.set_transformation(Self::plane_transform(n, plane1_dist, plane_size));
                }

                // Show intersection point (use precomputed value)
                if let Some(p) = result_point_from_line_plane {
                    if let Some(m) = &mut self.point3_mesh {
                        m.set_transformation(
                            Mat4::from_translation(p) * Mat4::from_scale(point_scale * 1.5),
                        );
                    }
                }
            }
        }

        // Collect objects to render
        let mut objects: Vec<&dyn Object> = Vec::new();

        if self.show_world_axes {
            if let Some(axes) = &self.world_axes {
                objects.push(axes);
            }
        }

        // Add all meshes (hidden ones have scale 0)
        if let Some(m) = &self.point1_mesh {
            objects.push(m);
        }
        if let Some(m) = &self.point2_mesh {
            objects.push(m);
        }
        if let Some(m) = &self.point3_mesh {
            objects.push(m);
        }
        if let Some(m) = &self.line_mesh {
            objects.push(m);
        }
        if let Some(m) = &self.result_line_mesh {
            objects.push(m);
        }
        if let Some(m) = &self.plane1_mesh {
            objects.push(m);
        }
        if let Some(m) = &self.plane2_mesh {
            objects.push(m);
        }
        if let Some(m) = &self.result_plane_mesh {
            objects.push(m);
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

/// Educational content for the Point-Line-Plane operations visualization.
const PROJECTIVE3_GEOMETRY_EDUCATION: EducationalContent = EducationalContent {
    title: "Point-Line-Plane Operations in 3D PGA",

    overview: "\
Join (^) and Meet (v) are the fundamental operations in projective geometry. \
Join creates the smallest element containing both inputs, while meet finds their \
intersection.

These operations are dual to each other and form the basis of all geometric \
constructions in PGA.",

    math_background: "\
JOIN (Wedge Product ^):
  - Point ^ Point = Line (through both points)
  - Point ^ Line = Plane (containing point and line)
  - Line ^ Point = Plane

MEET (Antiwedge Product v):
  - Plane v Plane = Line (intersection)
  - Line v Plane = Point (intersection)
  - Plane v Line = Point

Grade relationships:
  - Join increases grade: 1 ^ 1 = 2, 1 ^ 2 = 3
  - Meet decreases grade: 3 v 3 = 2, 2 v 3 = 1

Degenerate cases:
  - Parallel elements meet at infinity
  - Coincident elements have zero result",

    how_to_use: "\
- Select operation mode (P^P, P^L, Pi^Pi, L^Pi)
- Adjust input element parameters
- Observe the resulting geometric element
- Green shows computed results
- Use presets for common configurations",

    key_concepts: "\
- Join and meet are coordinate-free operations
- Results automatically handle degenerate cases
- Points at infinity represent directions
- Planes at infinity don't exist in Euclidean space
- Operations are associative: (A ^ B) ^ C = A ^ (B ^ C)",

    resources: &[
        (
            "RGA Wiki - Join and Meet",
            "https://rigidgeometricalgebra.org/wiki/index.php?title=Join_and_meet",
        ),
        (
            "PGA for Computer Scientists",
            "https://bivector.net/PGA4CS.pdf",
        ),
    ],
};
