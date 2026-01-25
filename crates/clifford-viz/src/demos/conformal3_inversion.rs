//! Sphere Inversion - 3D Conformal GA Visualization
//!
//! This demo demonstrates sphere inversion (geometric inversion),
//! a fundamental conformal transformation implemented using 3D CGA.
//!
//! In sphere inversion through a sphere S with center O and radius r:
//! - A point P maps to P' such that |OP| * |OP'| = r^2 and P, O, P' are collinear
//! - Spheres map to spheres (or planes if they pass through the center)
//! - Planes map to spheres (or planes if they pass through the center)
//!
//! In CGA, inversion through sphere S is computed as:
//!     P' = S * P * S^(-1)
//! where * is the geometric product.

use crate::common::prelude::*;
use clifford::ops::Transform;
use clifford::specialized::conformal::dim3::{RoundPoint, Sphere};
use three_d::*;

/// Viewport bounds for 3D scene.
const SCENE_SIZE: f32 = 5.0;

/// Epsilon for numerical comparisons.
const EPSILON: f64 = 1e-10;

/// A point in the 3D scene.
#[derive(Clone)]
struct ScenePoint {
    /// The CGA round point.
    point: RoundPoint<f64>,
    /// Display name.
    name: String,
}

impl ScenePoint {
    /// Creates a new scene point from Euclidean coordinates.
    fn new(x: f64, y: f64, z: f64, name: &str) -> Self {
        Self {
            point: RoundPoint::from_euclidean(x, y, z),
            name: name.to_string(),
        }
    }

    /// Returns the Euclidean coordinates.
    fn euclidean(&self) -> Option<(f64, f64, f64)> {
        self.point.to_euclidean()
    }

    /// Sets the position from Euclidean coordinates.
    fn set_position(&mut self, x: f64, y: f64, z: f64) {
        self.point = RoundPoint::from_euclidean(x, y, z);
    }
}

/// A sphere in the 3D scene.
#[derive(Clone)]
struct SceneSphere {
    /// Center x.
    cx: f64,
    /// Center y.
    cy: f64,
    /// Center z.
    cz: f64,
    /// Radius.
    radius: f64,
    /// Display name.
    name: String,
}

impl SceneSphere {
    /// Creates a new scene sphere.
    fn new(cx: f64, cy: f64, cz: f64, radius: f64, name: &str) -> Self {
        Self {
            cx,
            cy,
            cz,
            radius,
            name: name.to_string(),
        }
    }

    /// Returns center and radius.
    fn center_radius(&self) -> (f64, f64, f64, f64) {
        (self.cx, self.cy, self.cz, self.radius)
    }
}

/// Result of inverting a sphere - can be a sphere or a plane.
#[derive(Clone)]
enum InvertedSphere {
    /// Sphere with center and radius.
    Sphere {
        /// Center x.
        cx: f64,
        /// Center y.
        cy: f64,
        /// Center z.
        cz: f64,
        /// Radius.
        radius: f64,
    },
    /// Plane (when original passes through inversion center).
    Plane {
        /// Normal x.
        _nx: f64,
        /// Normal y.
        _ny: f64,
        /// Normal z.
        _nz: f64,
        /// Distance from origin.
        _d: f64,
    },
}

/// Number of segments for wireframe sphere rings.
const SPHERE_RING_SEGMENTS: usize = 16;

/// Number of latitude/longitude rings for wireframe sphere.
const SPHERE_RINGS: usize = 4;

/// Demo state for 3D Sphere Inversion visualization.
pub struct Conformal3InversionDemo {
    // === Three-d 3D Resources ===
    /// Graphics context for mesh creation.
    context: Option<Context>,
    /// 3D camera for rendering.
    camera: Option<Camera>,
    /// Orbit control for camera interaction.
    control: Option<OrbitControl>,
    /// Ambient lighting for the scene.
    ambient_light: Option<AmbientLight>,
    /// Key light (main directional light).
    key_light: Option<DirectionalLight>,
    /// Coordinate axes visualization.
    world_axes: Option<Axes>,

    /// Inversion sphere wireframe segments (cylinders forming rings).
    inversion_sphere_segments: Vec<Gm<Mesh, ColorMaterial>>,
    /// Original point meshes (small spheres).
    point_meshes: Vec<Gm<Mesh, ColorMaterial>>,
    /// Inverted point meshes (small spheres).
    inv_point_meshes: Vec<Gm<Mesh, ColorMaterial>>,
    /// Original sphere wireframe segments (one Vec per sphere).
    sphere_segments: Vec<Vec<Gm<Mesh, ColorMaterial>>>,
    /// Inverted sphere wireframe segments (one Vec per sphere).
    inv_sphere_segments: Vec<Vec<Gm<Mesh, ColorMaterial>>>,
    /// Ray line meshes (cylinders connecting original to inverted).
    ray_meshes: Vec<Gm<Mesh, ColorMaterial>>,

    // === Inversion Sphere Parameters ===
    /// Inversion sphere center x.
    inversion_cx: f32,
    /// Inversion sphere center y.
    inversion_cy: f32,
    /// Inversion sphere center z.
    inversion_cz: f32,
    /// Inversion sphere radius.
    inversion_radius: f32,

    /// Points in the scene.
    points: Vec<ScenePoint>,
    /// Spheres in the scene.
    spheres: Vec<SceneSphere>,

    // === Display Options ===
    /// Whether to show coordinate axes.
    show_axes: bool,
    /// Whether to show the inversion sphere.
    show_inversion_sphere: bool,
    /// Whether to show the original objects.
    show_originals: bool,
    /// Whether to show the inverted objects.
    show_inverted: bool,
    /// Whether to show connecting rays.
    show_rays: bool,
}

impl Default for Conformal3InversionDemo {
    fn default() -> Self {
        // Initial configuration with some points and spheres
        let points = vec![
            ScenePoint::new(2.5, 1.0, 0.5, "P1"),
            ScenePoint::new(3.0, -1.5, 1.0, "P2"),
            ScenePoint::new(-2.0, 2.0, -0.5, "P3"),
        ];

        let spheres = vec![
            SceneSphere::new(-2.0, -1.0, 0.0, 0.8, "S1"),
            SceneSphere::new(3.0, 1.5, 0.5, 0.6, "S2"),
        ];

        Self {
            context: None,
            camera: None,
            control: None,
            ambient_light: None,
            key_light: None,
            world_axes: None,
            inversion_sphere_segments: Vec::new(),
            point_meshes: Vec::new(),
            inv_point_meshes: Vec::new(),
            sphere_segments: Vec::new(),
            inv_sphere_segments: Vec::new(),
            ray_meshes: Vec::new(),
            inversion_cx: 0.0,
            inversion_cy: 0.0,
            inversion_cz: 0.0,
            inversion_radius: 2.0,
            points,
            spheres,
            show_axes: true,
            show_inversion_sphere: true,
            show_originals: true,
            show_inverted: true,
            show_rays: true,
        }
    }
}

impl Conformal3InversionDemo {
    /// Get the CGA inversion sphere from current parameters.
    fn inversion_sphere(&self) -> Sphere<f64> {
        Sphere::from_center_radius(
            f64::from(self.inversion_cx),
            f64::from(self.inversion_cy),
            f64::from(self.inversion_cz),
            f64::from(self.inversion_radius),
        )
    }

    /// Invert a point through the inversion sphere using CGA.
    fn invert_point(&self, point: &RoundPoint<f64>) -> Option<(f64, f64, f64)> {
        let inv_sphere = self.inversion_sphere();
        let result = inv_sphere.transform(point);
        result.to_euclidean()
    }

    /// Invert a sphere through the inversion sphere.
    ///
    /// In CGA, sphere inversion can be computed as: S' = I * S * I^(-1)
    /// where I is the inversion sphere. For visualization, we use the
    /// classical geometric formula which is equivalent but easier to
    /// extract center/radius from.
    ///
    /// Classical inversion: For a sphere with center C and radius r,
    /// inverted through sphere with center O and radius R:
    /// - If |OC| = r (sphere passes through O), result is a plane
    /// - Otherwise, invert the near/far points and reconstruct
    fn invert_sphere(&self, sphere: &SceneSphere) -> Option<InvertedSphere> {
        let inv_cx = f64::from(self.inversion_cx);
        let inv_cy = f64::from(self.inversion_cy);
        let inv_cz = f64::from(self.inversion_cz);
        let inv_r_sq = f64::from(self.inversion_radius * self.inversion_radius);

        let (cx, cy, cz, r) = sphere.center_radius();

        // Distance from inversion center to sphere center
        let dist_to_center =
            ((cx - inv_cx).powi(2) + (cy - inv_cy).powi(2) + (cz - inv_cz).powi(2)).sqrt();

        if (dist_to_center - r).abs() < EPSILON {
            // Sphere passes through inversion center -> becomes a plane
            let dx = cx - inv_cx;
            let dy = cy - inv_cy;
            let dz = cz - inv_cz;
            let len = (dx * dx + dy * dy + dz * dz).sqrt();
            if len < EPSILON {
                return None;
            }
            return Some(InvertedSphere::Plane {
                _nx: dx / len,
                _ny: dy / len,
                _nz: dz / len,
                _d: inv_r_sq / (2.0 * dist_to_center),
            });
        }

        // General case: invert the near and far points of the sphere
        // (points along the line from inversion center through sphere center)
        let dir_x = (cx - inv_cx) / dist_to_center;
        let dir_y = (cy - inv_cy) / dist_to_center;
        let dir_z = (cz - inv_cz) / dist_to_center;

        // Near point (closest to inversion center)
        let near_x = cx - r * dir_x;
        let near_y = cy - r * dir_y;
        let near_z = cz - r * dir_z;

        // Far point (farthest from inversion center)
        let far_x = cx + r * dir_x;
        let far_y = cy + r * dir_y;
        let far_z = cz + r * dir_z;

        // Invert both points: P' = O + R^2 * (P - O) / |P - O|^2
        let near_dist_sq =
            (near_x - inv_cx).powi(2) + (near_y - inv_cy).powi(2) + (near_z - inv_cz).powi(2);
        let far_dist_sq =
            (far_x - inv_cx).powi(2) + (far_y - inv_cy).powi(2) + (far_z - inv_cz).powi(2);

        if near_dist_sq < EPSILON || far_dist_sq < EPSILON {
            return None;
        }

        let near_inv_x = inv_cx + inv_r_sq * (near_x - inv_cx) / near_dist_sq;
        let near_inv_y = inv_cy + inv_r_sq * (near_y - inv_cy) / near_dist_sq;
        let near_inv_z = inv_cz + inv_r_sq * (near_z - inv_cz) / near_dist_sq;

        let far_inv_x = inv_cx + inv_r_sq * (far_x - inv_cx) / far_dist_sq;
        let far_inv_y = inv_cy + inv_r_sq * (far_y - inv_cy) / far_dist_sq;
        let far_inv_z = inv_cz + inv_r_sq * (far_z - inv_cz) / far_dist_sq;

        // New sphere: center is midpoint, radius is half the distance
        Some(InvertedSphere::Sphere {
            cx: (near_inv_x + far_inv_x) / 2.0,
            cy: (near_inv_y + far_inv_y) / 2.0,
            cz: (near_inv_z + far_inv_z) / 2.0,
            radius: ((near_inv_x - far_inv_x).powi(2)
                + (near_inv_y - far_inv_y).powi(2)
                + (near_inv_z - far_inv_z).powi(2))
            .sqrt()
                / 2.0,
        })
    }

    /// Create a transformation matrix for a cylinder connecting two points.
    fn cylinder_transform(p1: Vector3<f32>, p2: Vector3<f32>, radius: f32) -> Mat4 {
        let diff = p2 - p1;
        let length = diff.magnitude();
        if length < 1e-6 {
            return Mat4::from_scale(0.0);
        }

        let dir = diff / length;
        let midpoint = (p1 + p2) * 0.5;

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

    /// Generate wireframe line segments for a sphere (latitude and longitude rings).
    /// Returns pairs of points representing line segment endpoints.
    fn sphere_wireframe_segments(
        center: Vector3<f32>,
        radius: f32,
        ring_segments: usize,
        num_rings: usize,
    ) -> Vec<(Vector3<f32>, Vector3<f32>)> {
        use std::f32::consts::PI;
        let mut segments = Vec::new();

        // Latitude rings (horizontal circles at different heights)
        for ring in 0..=num_rings {
            let phi = PI * (ring as f32) / (num_rings as f32);
            let r = radius * phi.sin();
            let z = radius * phi.cos();

            for i in 0..ring_segments {
                let theta1 = 2.0 * PI * (i as f32) / (ring_segments as f32);
                let theta2 = 2.0 * PI * ((i + 1) as f32) / (ring_segments as f32);

                let p1 = center + vec3(r * theta1.cos(), r * theta1.sin(), z);
                let p2 = center + vec3(r * theta2.cos(), r * theta2.sin(), z);
                segments.push((p1, p2));
            }
        }

        // Longitude rings (vertical great circles)
        for i in 0..num_rings {
            let theta = 2.0 * PI * (i as f32) / (num_rings as f32);

            for j in 0..ring_segments {
                let phi1 = PI * (j as f32) / (ring_segments as f32);
                let phi2 = PI * ((j + 1) as f32) / (ring_segments as f32);

                let p1 = center
                    + vec3(
                        radius * phi1.sin() * theta.cos(),
                        radius * phi1.sin() * theta.sin(),
                        radius * phi1.cos(),
                    );
                let p2 = center
                    + vec3(
                        radius * phi2.sin() * theta.cos(),
                        radius * phi2.sin() * theta.sin(),
                        radius * phi2.cos(),
                    );
                segments.push((p1, p2));
            }
        }

        segments
    }
}

impl VisualizationApp for Conformal3InversionDemo {
    fn name(&self) -> &'static str {
        "Conformal 3D - Sphere Inversion"
    }

    fn update(&mut self, _dt: f32) {
        // No animation in this demo
    }

    fn render(&mut self, _ui: &mut egui::Ui) {
        // 3D rendering is done in render_3d
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        // === Inversion Sphere ===
        group_header(ui, "Inversion Sphere");

        ui.horizontal(|ui| {
            ui.label("Center X:");
            ui.add(egui::DragValue::new(&mut self.inversion_cx).speed(0.1));
        });
        ui.horizontal(|ui| {
            ui.label("Center Y:");
            ui.add(egui::DragValue::new(&mut self.inversion_cy).speed(0.1));
        });
        ui.horizontal(|ui| {
            ui.label("Center Z:");
            ui.add(egui::DragValue::new(&mut self.inversion_cz).speed(0.1));
        });
        ui.horizontal(|ui| {
            ui.label("Radius:");
            ui.add(
                egui::DragValue::new(&mut self.inversion_radius)
                    .speed(0.1)
                    .range(0.1..=10.0),
            );
        });

        ui.add_space(12.0);

        // === Display Options ===
        group_header(ui, "Display");

        ui.checkbox(&mut self.show_axes, "Show Axes");
        ui.checkbox(&mut self.show_inversion_sphere, "Show Inversion Sphere");
        ui.checkbox(&mut self.show_originals, "Show Original Objects");
        ui.checkbox(&mut self.show_inverted, "Show Inverted Objects");
        ui.checkbox(&mut self.show_rays, "Show Inversion Rays");

        ui.add_space(12.0);

        // === Points ===
        group_header(ui, "Points");

        let num_points = self.points.len();
        for i in 0..num_points {
            let name = self.points[i].name.clone();
            let coords = self.points[i].euclidean();
            ui.horizontal(|ui| {
                ui.label(format!("{}:", name));
                if let Some((x, y, z)) = coords {
                    let mut px = x as f32;
                    let mut py = y as f32;
                    let mut pz = z as f32;
                    let changed_x = ui
                        .add(egui::DragValue::new(&mut px).speed(0.1).prefix("x:"))
                        .changed();
                    let changed_y = ui
                        .add(egui::DragValue::new(&mut py).speed(0.1).prefix("y:"))
                        .changed();
                    let changed_z = ui
                        .add(egui::DragValue::new(&mut pz).speed(0.1).prefix("z:"))
                        .changed();
                    if changed_x || changed_y || changed_z {
                        self.points[i].set_position(f64::from(px), f64::from(py), f64::from(pz));
                    }
                }
            });
        }

        ui.add_space(12.0);

        // === Spheres ===
        group_header(ui, "Spheres");

        let num_spheres = self.spheres.len();
        for i in 0..num_spheres {
            let name = self.spheres[i].name.clone();
            let mut sx = self.spheres[i].cx as f32;
            let mut sy = self.spheres[i].cy as f32;
            let mut sz = self.spheres[i].cz as f32;
            let mut sr = self.spheres[i].radius as f32;
            ui.horizontal(|ui| {
                ui.label(format!("{}:", name));
                let cx = ui
                    .add(egui::DragValue::new(&mut sx).speed(0.1).prefix("x:"))
                    .changed();
                let cy = ui
                    .add(egui::DragValue::new(&mut sy).speed(0.1).prefix("y:"))
                    .changed();
                let cz = ui
                    .add(egui::DragValue::new(&mut sz).speed(0.1).prefix("z:"))
                    .changed();
                let cr = ui
                    .add(
                        egui::DragValue::new(&mut sr)
                            .speed(0.1)
                            .prefix("r:")
                            .range(0.1..=5.0),
                    )
                    .changed();
                if cx || cy || cz || cr {
                    self.spheres[i].cx = f64::from(sx);
                    self.spheres[i].cy = f64::from(sy);
                    self.spheres[i].cz = f64::from(sz);
                    self.spheres[i].radius = f64::from(sr);
                }
            });
        }
    }

    fn info(&self, ui: &mut egui::Ui) {
        ui.label("Drag sliders to move points and spheres.");
        ui.label("Observe how they invert through the sphere.");
    }
}

impl VisualizationApp3D for Conformal3InversionDemo {
    fn init_3d(&mut self, context: &Context) {
        // Store context for dynamic mesh creation
        self.context = Some(context.clone());

        // Initialize camera
        self.camera = Some(Camera::new_perspective(
            Viewport::new_at_origo(1, 1),
            vec3(8.0, 6.0, 5.0),
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 0.0, 1.0),
            degrees(45.0),
            0.1,
            100.0,
        ));

        // Initialize orbit control
        self.control = Some(OrbitControl::new(
            vec3(0.0, 0.0, 0.0),
            1.0,
            SCENE_SIZE * 4.0,
        ));

        // Initialize world axes
        self.world_axes = Some(Axes::new(context, 0.02, 2.0));

        // Initialize lighting
        self.ambient_light = Some(AmbientLight::new(context, 0.4, Srgba::WHITE));

        // Key light (main directional light)
        self.key_light = Some(DirectionalLight::new(
            context,
            0.8,
            Srgba::WHITE,
            vec3(1.0, 0.5, 1.0).normalize(),
        ));

        // Create meshes
        let sphere_cpu = CpuMesh::sphere(16);
        let cylinder_cpu = CpuMesh::cylinder(8);

        // Inversion sphere wireframe segments (purple)
        // Calculate how many segments we need for the wireframe
        let num_wireframe_segments =
            (SPHERE_RINGS + 1) * SPHERE_RING_SEGMENTS + SPHERE_RINGS * SPHERE_RING_SEGMENTS;
        for _ in 0..num_wireframe_segments {
            self.inversion_sphere_segments.push(Gm::new(
                Mesh::new(context, &cylinder_cpu),
                ColorMaterial {
                    color: Srgba::new(180, 100, 220, 255),
                    ..Default::default()
                },
            ));
        }

        // Create point meshes (3 original + 3 inverted)
        for _ in 0..3 {
            // Original points (blue)
            self.point_meshes.push(Gm::new(
                Mesh::new(context, &sphere_cpu),
                ColorMaterial {
                    color: Srgba::new(100, 150, 255, 255),
                    ..Default::default()
                },
            ));
            // Inverted points (gold)
            self.inv_point_meshes.push(Gm::new(
                Mesh::new(context, &sphere_cpu),
                ColorMaterial {
                    color: Srgba::new(255, 200, 50, 255),
                    ..Default::default()
                },
            ));
            // Ray cylinders (gray)
            self.ray_meshes.push(Gm::new(
                Mesh::new(context, &cylinder_cpu),
                ColorMaterial {
                    color: Srgba::new(150, 150, 150, 150),
                    is_transparent: true,
                    ..Default::default()
                },
            ));
        }

        // Create sphere wireframe segments (2 original + 2 inverted)
        let num_sphere_segments =
            (SPHERE_RINGS + 1) * SPHERE_RING_SEGMENTS + SPHERE_RINGS * SPHERE_RING_SEGMENTS;
        for _ in 0..2 {
            // Original spheres (green wireframe)
            let mut orig_segs = Vec::with_capacity(num_sphere_segments);
            for _ in 0..num_sphere_segments {
                orig_segs.push(Gm::new(
                    Mesh::new(context, &cylinder_cpu),
                    ColorMaterial {
                        color: Srgba::new(100, 220, 100, 255),
                        ..Default::default()
                    },
                ));
            }
            self.sphere_segments.push(orig_segs);

            // Inverted spheres (orange wireframe)
            let mut inv_segs = Vec::with_capacity(num_sphere_segments);
            for _ in 0..num_sphere_segments {
                inv_segs.push(Gm::new(
                    Mesh::new(context, &cylinder_cpu),
                    ColorMaterial {
                        color: Srgba::new(255, 150, 50, 255),
                        ..Default::default()
                    },
                ));
            }
            self.inv_sphere_segments.push(inv_segs);
        }
    }

    fn render_3d(&mut self, frame: &mut FrameInput) {
        // Pre-compute all inversions BEFORE any mutable borrows
        let point_inversions: Vec<Option<(f64, f64, f64)>> = self
            .points
            .iter()
            .map(|pt| self.invert_point(&pt.point))
            .collect();

        let point_originals: Vec<Option<(f64, f64, f64)>> =
            self.points.iter().map(|pt| pt.euclidean()).collect();

        let sphere_originals: Vec<(f64, f64, f64, f64)> =
            self.spheres.iter().map(|s| s.center_radius()).collect();

        let sphere_inversions: Vec<Option<InvertedSphere>> =
            self.spheres.iter().map(|s| self.invert_sphere(s)).collect();

        // Display options (copy to avoid borrows)
        let show_inversion_sphere = self.show_inversion_sphere;
        let show_originals = self.show_originals;
        let show_inverted = self.show_inverted;
        let show_rays = self.show_rays;
        let show_axes = self.show_axes;
        let inv_cx = self.inversion_cx;
        let inv_cy = self.inversion_cy;
        let inv_cz = self.inversion_cz;
        let inv_radius = self.inversion_radius;

        // Update camera viewport and handle controls
        let camera = self.camera.as_mut().unwrap();
        let control = self.control.as_mut().unwrap();

        camera.set_viewport(frame.viewport);
        control.handle_events(camera, &mut frame.events);

        // Hide transform (scale to zero)
        let hide = Mat4::from_scale(0.0);

        // Update inversion sphere wireframe
        if show_inversion_sphere {
            let center = vec3(inv_cx, inv_cy, inv_cz);
            let segments = Self::sphere_wireframe_segments(
                center,
                inv_radius,
                SPHERE_RING_SEGMENTS,
                SPHERE_RINGS,
            );
            for (i, mesh) in self.inversion_sphere_segments.iter_mut().enumerate() {
                if let Some((p1, p2)) = segments.get(i) {
                    mesh.set_transformation(Self::cylinder_transform(*p1, *p2, 0.015));
                } else {
                    mesh.set_transformation(hide);
                }
            }
        } else {
            for mesh in &mut self.inversion_sphere_segments {
                mesh.set_transformation(hide);
            }
        }

        // Update point meshes
        let point_size = 0.12;
        for i in 0..point_originals.len() {
            // Original point
            if let Some(mesh) = self.point_meshes.get_mut(i) {
                if show_originals {
                    if let Some((x, y, z)) = point_originals[i] {
                        mesh.set_transformation(
                            Mat4::from_translation(vec3(x as f32, y as f32, z as f32))
                                * Mat4::from_scale(point_size),
                        );
                    } else {
                        mesh.set_transformation(hide);
                    }
                } else {
                    mesh.set_transformation(hide);
                }
            }

            // Inverted point
            if let Some(mesh) = self.inv_point_meshes.get_mut(i) {
                if show_inverted {
                    if let Some((inv_x, inv_y, inv_z)) = point_inversions[i] {
                        // Only show if within bounds
                        if inv_x.abs() < (SCENE_SIZE * 2.0) as f64
                            && inv_y.abs() < (SCENE_SIZE * 2.0) as f64
                            && inv_z.abs() < (SCENE_SIZE * 2.0) as f64
                        {
                            mesh.set_transformation(
                                Mat4::from_translation(vec3(
                                    inv_x as f32,
                                    inv_y as f32,
                                    inv_z as f32,
                                )) * Mat4::from_scale(point_size),
                            );
                        } else {
                            mesh.set_transformation(hide);
                        }
                    } else {
                        mesh.set_transformation(hide);
                    }
                } else {
                    mesh.set_transformation(hide);
                }
            }

            // Ray connecting original to inverted
            if let Some(mesh) = self.ray_meshes.get_mut(i) {
                if show_rays && show_originals && show_inverted {
                    if let (Some((ox, oy, oz)), Some((ix, iy, iz))) =
                        (point_originals[i], point_inversions[i])
                    {
                        if ix.abs() < (SCENE_SIZE * 2.0) as f64
                            && iy.abs() < (SCENE_SIZE * 2.0) as f64
                            && iz.abs() < (SCENE_SIZE * 2.0) as f64
                        {
                            mesh.set_transformation(Self::cylinder_transform(
                                vec3(ox as f32, oy as f32, oz as f32),
                                vec3(ix as f32, iy as f32, iz as f32),
                                0.02,
                            ));
                        } else {
                            mesh.set_transformation(hide);
                        }
                    } else {
                        mesh.set_transformation(hide);
                    }
                } else {
                    mesh.set_transformation(hide);
                }
            }
        }

        // Update sphere wireframes
        for i in 0..sphere_originals.len() {
            let (cx, cy, cz, r) = sphere_originals[i];
            let center = vec3(cx as f32, cy as f32, cz as f32);

            // Original sphere wireframe
            if let Some(segs) = self.sphere_segments.get_mut(i) {
                if show_originals {
                    let segments = Self::sphere_wireframe_segments(
                        center,
                        r as f32,
                        SPHERE_RING_SEGMENTS,
                        SPHERE_RINGS,
                    );
                    for (j, mesh) in segs.iter_mut().enumerate() {
                        if let Some((p1, p2)) = segments.get(j) {
                            mesh.set_transformation(Self::cylinder_transform(*p1, *p2, 0.012));
                        } else {
                            mesh.set_transformation(hide);
                        }
                    }
                } else {
                    for mesh in segs.iter_mut() {
                        mesh.set_transformation(hide);
                    }
                }
            }

            // Inverted sphere wireframe
            if let Some(segs) = self.inv_sphere_segments.get_mut(i) {
                if show_inverted {
                    if let Some(ref inverted) = sphere_inversions[i] {
                        match inverted {
                            InvertedSphere::Sphere {
                                cx: icx,
                                cy: icy,
                                cz: icz,
                                radius: ir,
                            } => {
                                if icx.abs() < (SCENE_SIZE * 3.0) as f64
                                    && icy.abs() < (SCENE_SIZE * 3.0) as f64
                                    && icz.abs() < (SCENE_SIZE * 3.0) as f64
                                    && *ir < (SCENE_SIZE * 3.0) as f64
                                {
                                    let inv_center = vec3(*icx as f32, *icy as f32, *icz as f32);
                                    let segments = Self::sphere_wireframe_segments(
                                        inv_center,
                                        *ir as f32,
                                        SPHERE_RING_SEGMENTS,
                                        SPHERE_RINGS,
                                    );
                                    for (j, mesh) in segs.iter_mut().enumerate() {
                                        if let Some((p1, p2)) = segments.get(j) {
                                            mesh.set_transformation(Self::cylinder_transform(
                                                *p1, *p2, 0.012,
                                            ));
                                        } else {
                                            mesh.set_transformation(hide);
                                        }
                                    }
                                } else {
                                    for mesh in segs.iter_mut() {
                                        mesh.set_transformation(hide);
                                    }
                                }
                            }
                            InvertedSphere::Plane { .. } => {
                                // For now, hide planes
                                for mesh in segs.iter_mut() {
                                    mesh.set_transformation(hide);
                                }
                            }
                        }
                    } else {
                        for mesh in segs.iter_mut() {
                            mesh.set_transformation(hide);
                        }
                    }
                } else {
                    for mesh in segs.iter_mut() {
                        mesh.set_transformation(hide);
                    }
                }
            }
        }

        // Collect objects to render
        let mut objects: Vec<&dyn Object> = Vec::new();

        if show_axes {
            if let Some(axes) = &self.world_axes {
                objects.push(axes);
            }
        }

        // Add inversion sphere wireframe
        for mesh in &self.inversion_sphere_segments {
            objects.push(mesh);
        }

        // Add point meshes
        for mesh in &self.point_meshes {
            objects.push(mesh);
        }
        for mesh in &self.inv_point_meshes {
            objects.push(mesh);
        }
        for mesh in &self.ray_meshes {
            objects.push(mesh);
        }

        // Add sphere wireframes
        for segs in &self.sphere_segments {
            for mesh in segs {
                objects.push(mesh);
            }
        }
        for segs in &self.inv_sphere_segments {
            for mesh in segs {
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
