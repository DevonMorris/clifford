//! 3D Euclidean Rotors Demo using three-d
//!
//! Run with: `cargo run -p clifford-viz --example euclidean3_three_d --features three-d --release`
//!
//! This demo shows GA rotors vs Euler angles for 3D rotation, with real 3D rendering.
//! All geometric computations use clifford types - three-d is only used for rendering.

use clifford::ops::{RightComplement, Transform};
use clifford::specialized::euclidean::dim3::{Bivector, Rotor, Vector};
use three_d::*;

/// Rotation mode for comparison
#[derive(Debug, Clone, Copy, PartialEq)]
enum RotationMode {
    Rotor,
    Euler,
}

fn main() {
    let window = Window::new(WindowSettings {
        title: "3D Euclidean Rotors (three-d)".to_string(),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();
    let mut gui = GUI::new(&context);

    // Camera
    let mut camera = Camera::new_perspective(
        window.viewport(),
        vec3(4.0, 3.0, 5.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        degrees(45.0),
        0.1,
        100.0,
    );
    let mut control = OrbitControl::new(vec3(0.0, 0.0, 0.0), 1.0, 15.0);

    // Coordinate axes
    let axes = Axes::new(&context, 0.015, 1.5);

    // Lighting - key light + fill light + ambient
    let ambient = AmbientLight::new(&context, 0.3, Srgba::WHITE);
    let key_light = DirectionalLight::new(&context, 2.5, Srgba::WHITE, vec3(-1.0, -1.0, -1.0));
    let fill_light = DirectionalLight::new(
        &context,
        1.0,
        Srgba::new_opaque(200, 200, 255),
        vec3(1.0, 0.5, 1.0),
    );

    // UI state
    let mut mode = RotationMode::Rotor;
    let mut rotor_angle: f32 = 0.0;
    // Rotation axis as clifford Vector
    let mut axis = Vector::new(0.0_f32, 0.0, 1.0);
    let mut euler_yaw: f32 = 0.0;
    let mut euler_pitch: f32 = 0.0;
    let mut euler_roll: f32 = 0.0;
    let mut auto_rotate = false;
    let mut show_plane = true;

    window.render_loop(move |mut frame_input| {
        let mut gui_consumed = false;

        // Axis components for UI (will update the Vector)
        let mut ax = axis.x();
        let mut ay = axis.y();
        let mut az = axis.z();

        // Update GUI
        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |ctx| {
                gui_consumed = ctx.wants_pointer_input() || ctx.wants_keyboard_input();

                egui::Window::new("3D Rotor Demo")
                    .default_pos([10.0, 10.0])
                    .default_width(300.0)
                    .show(ctx, |ui| {
                        ui.heading("Rotation Mode");
                        ui.horizontal(|ui| {
                            ui.selectable_value(&mut mode, RotationMode::Rotor, "Rotor (GA)");
                            ui.selectable_value(&mut mode, RotationMode::Euler, "Euler Angles");
                        });

                        ui.separator();

                        match mode {
                            RotationMode::Rotor => {
                                ui.heading("Rotor Parameters");

                                ui.label("Rotation Axis (Vector):");
                                ui.horizontal(|ui| {
                                    ui.colored_label(egui::Color32::from_rgb(255, 100, 100), "x:");
                                    ui.add(
                                        egui::Slider::new(&mut ax, -1.0..=1.0).show_value(false),
                                    );
                                    ui.label(format!("{:.2}", ax));
                                });
                                ui.horizontal(|ui| {
                                    ui.colored_label(egui::Color32::from_rgb(100, 255, 100), "y:");
                                    ui.add(
                                        egui::Slider::new(&mut ay, -1.0..=1.0).show_value(false),
                                    );
                                    ui.label(format!("{:.2}", ay));
                                });
                                ui.horizontal(|ui| {
                                    ui.colored_label(egui::Color32::from_rgb(100, 100, 255), "z:");
                                    ui.add(
                                        egui::Slider::new(&mut az, -1.0..=1.0).show_value(false),
                                    );
                                    ui.label(format!("{:.2}", az));
                                });

                                ui.horizontal(|ui| {
                                    if ui.small_button("e1").clicked() {
                                        ax = 1.0;
                                        ay = 0.0;
                                        az = 0.0;
                                    }
                                    if ui.small_button("e2").clicked() {
                                        ax = 0.0;
                                        ay = 1.0;
                                        az = 0.0;
                                    }
                                    if ui.small_button("e3").clicked() {
                                        ax = 0.0;
                                        ay = 0.0;
                                        az = 1.0;
                                    }
                                    if ui.small_button("e1+e2+e3").clicked() {
                                        ax = 1.0;
                                        ay = 1.0;
                                        az = 1.0;
                                    }
                                });

                                ui.separator();
                                ui.horizontal(|ui| {
                                    ui.label("Angle:");
                                    ui.add(
                                        egui::Slider::new(&mut rotor_angle, 0.0..=360.0)
                                            .suffix(" deg"),
                                    );
                                });

                                ui.checkbox(&mut auto_rotate, "Auto-rotate");
                                ui.checkbox(&mut show_plane, "Show rotation plane (bivector)");

                                // Compute and display GA values
                                ui.separator();
                                ui.heading("Geometric Algebra Values");

                                let axis_vec = Vector::new(ax, ay, az);
                                let axis_norm = axis_vec.norm();

                                if axis_norm > 1e-6 {
                                    let axis_unit = axis_vec.normalized();
                                    ui.label(format!(
                                        "Axis (unit): {:.2}e1 + {:.2}e2 + {:.2}e3",
                                        axis_unit.x(),
                                        axis_unit.y(),
                                        axis_unit.z()
                                    ));

                                    // Rotation plane is the Hodge dual of the axis
                                    let plane: Bivector<f32> = axis_unit.right_complement();
                                    ui.label(format!(
                                        "Plane (dual): {:.2}e23 + {:.2}e31 + {:.2}e12",
                                        plane.rx(),
                                        plane.ry(),
                                        plane.rz()
                                    ));

                                    let rotor =
                                        Rotor::from_angle_plane(rotor_angle.to_radians(), plane);
                                    ui.label(format!(
                                        "Rotor: {:.3} + {:.3}e23 + {:.3}e31 + {:.3}e12",
                                        rotor.s(),
                                        rotor.rx(),
                                        rotor.ry(),
                                        rotor.rz()
                                    ));

                                    ui.label(format!("|R| = {:.4}", rotor.norm()));
                                }
                            }
                            RotationMode::Euler => {
                                ui.heading("Euler Angles (YXZ)");

                                ui.horizontal(|ui| {
                                    ui.label("Yaw (Y):");
                                    ui.add(
                                        egui::Slider::new(&mut euler_yaw, -180.0..=180.0)
                                            .suffix(" deg"),
                                    );
                                });
                                ui.horizontal(|ui| {
                                    ui.label("Pitch (X):");
                                    ui.add(
                                        egui::Slider::new(&mut euler_pitch, -90.0..=90.0)
                                            .suffix(" deg"),
                                    );
                                });
                                ui.horizontal(|ui| {
                                    ui.label("Roll (Z):");
                                    ui.add(
                                        egui::Slider::new(&mut euler_roll, -180.0..=180.0)
                                            .suffix(" deg"),
                                    );
                                });

                                if (euler_pitch.abs() - 90.0).abs() < 5.0 {
                                    ui.separator();
                                    ui.colored_label(
                                        egui::Color32::from_rgb(255, 100, 100),
                                        "GIMBAL LOCK! Yaw and roll affect the same axis.",
                                    );
                                }

                                ui.separator();
                                if ui.button("Set pitch to 90 deg (gimbal lock)").clicked() {
                                    euler_pitch = 90.0;
                                }

                                // Show equivalent rotor
                                ui.separator();
                                ui.heading("Equivalent Rotor");
                                let r_yaw = Rotor::from_angle_plane(
                                    euler_yaw.to_radians(),
                                    Bivector::unit_ry(),
                                );
                                let r_pitch = Rotor::from_angle_plane(
                                    euler_pitch.to_radians(),
                                    Bivector::unit_rx(),
                                );
                                let r_roll = Rotor::from_angle_plane(
                                    euler_roll.to_radians(),
                                    Bivector::unit_rz(),
                                );
                                let combined = r_yaw * r_pitch * r_roll;
                                ui.label(format!(
                                    "R = {:.3} + {:.3}e23 + {:.3}e31 + {:.3}e12",
                                    combined.s(),
                                    combined.rx(),
                                    combined.ry(),
                                    combined.rz()
                                ));
                            }
                        }

                        ui.separator();
                        ui.label("Camera: drag to orbit, scroll to zoom");
                        if ui.button("Reset Camera").clicked() {
                            camera.set_view(
                                vec3(4.0, 3.0, 5.0),
                                vec3(0.0, 0.0, 0.0),
                                vec3(0.0, 1.0, 0.0),
                            );
                        }
                    });
            },
        );

        // Update axis from UI
        axis = Vector::new(ax, ay, az);

        // Camera control
        if !gui_consumed {
            control.handle_events(&mut camera, &mut frame_input.events);
        }

        // Auto-rotate
        if auto_rotate && mode == RotationMode::Rotor {
            rotor_angle += frame_input.elapsed_time as f32 * 0.03;
            if rotor_angle > 360.0 {
                rotor_angle -= 360.0;
            }
        }

        // Compute rotor using clifford
        let rotor = match mode {
            RotationMode::Rotor => {
                let norm = axis.norm();
                if norm > 1e-6 {
                    let plane = axis.normalized().right_complement();
                    Rotor::from_angle_plane(rotor_angle.to_radians(), plane)
                } else {
                    Rotor::new(1.0, 0.0, 0.0, 0.0) // identity
                }
            }
            RotationMode::Euler => {
                let r_yaw = Rotor::from_angle_plane(euler_yaw.to_radians(), Bivector::unit_ry());
                let r_pitch =
                    Rotor::from_angle_plane(euler_pitch.to_radians(), Bivector::unit_rx());
                let r_roll = Rotor::from_angle_plane(euler_roll.to_radians(), Bivector::unit_rz());
                r_yaw * r_pitch * r_roll
            }
        };

        // Transform cube vertices using clifford rotor
        let cube_verts = unit_cube_vertices();
        let transformed_verts: Vec<Vector<f32>> =
            cube_verts.iter().map(|v| rotor.transform(v)).collect();

        // Create cube mesh from transformed vertices
        let cube_color = if mode == RotationMode::Euler && (euler_pitch.abs() - 90.0).abs() < 5.0 {
            Srgba::new_opaque(255, 100, 100)
        } else {
            Srgba::new_opaque(100, 180, 255)
        };

        let cube = create_cube_mesh(&context, &transformed_verts, cube_color, &rotor);

        // Create local axes transformed by rotor
        let local_axes = create_transformed_axes(&context, &rotor);

        // Create rotation plane visualization (if enabled and in rotor mode)
        let plane_mesh = if show_plane && mode == RotationMode::Rotor && axis.norm() > 1e-6 {
            Some(create_plane_disk(&context, &axis.normalized()))
        } else {
            None
        };

        // Render
        let mut renderables: Vec<&dyn Object> =
            vec![&axes, &cube, &local_axes.0, &local_axes.1, &local_axes.2];

        if let Some(ref plane) = plane_mesh {
            renderables.push(plane);
        }

        let _ = frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.08, 0.08, 0.12, 1.0, 1.0))
            .render(&camera, renderables, &[&ambient, &key_light, &fill_light])
            .write(|| gui.render());

        FrameOutput::default()
    });
}

// =============================================================================
// Clifford-based geometry
// =============================================================================

/// Unit cube vertices as clifford Vectors (centered at origin, size 1)
fn unit_cube_vertices() -> [Vector<f32>; 8] {
    let s = 0.5;
    [
        Vector::new(-s, -s, -s),
        Vector::new(s, -s, -s),
        Vector::new(s, s, -s),
        Vector::new(-s, s, -s),
        Vector::new(-s, -s, s),
        Vector::new(s, -s, s),
        Vector::new(s, s, s),
        Vector::new(-s, s, s),
    ]
}

/// Create axes transformed by a rotor (all computation in clifford)
fn create_transformed_axes(
    context: &Context,
    rotor: &Rotor<f32>,
) -> (
    Gm<Mesh, ColorMaterial>,
    Gm<Mesh, ColorMaterial>,
    Gm<Mesh, ColorMaterial>,
) {
    let len = 0.8_f32;

    // Transform unit vectors with rotor
    let x_dir = rotor.transform(&(Vector::unit_x() * len));
    let y_dir = rotor.transform(&(Vector::unit_y() * len));
    let z_dir = rotor.transform(&(Vector::unit_z() * len));

    (
        create_arrow_from_vector(context, &x_dir, Srgba::new_opaque(255, 80, 80)),
        create_arrow_from_vector(context, &y_dir, Srgba::new_opaque(80, 255, 80)),
        create_arrow_from_vector(context, &z_dir, Srgba::new_opaque(80, 80, 255)),
    )
}

// =============================================================================
// Rendering helpers (convert clifford -> three-d at the last moment)
// =============================================================================

/// Convert clifford Vector to three-d Vec3
fn to_vec3(v: &Vector<f32>) -> Vec3 {
    vec3(v.x(), v.y(), v.z())
}

/// Create flat-shaded cube mesh from 8 clifford Vector vertices
/// Normals are transformed by the rotor to match the rotated faces
fn create_cube_mesh(
    context: &Context,
    verts: &[Vector<f32>],
    color: Srgba,
    rotor: &Rotor<f32>,
) -> Gm<Mesh, PhysicalMaterial> {
    // Face vertex indices and base normals (before rotation)
    let faces: [(usize, usize, usize, usize, Vector<f32>); 6] = [
        (4, 5, 6, 7, Vector::new(0.0, 0.0, 1.0)),  // Front +Z
        (1, 0, 3, 2, Vector::new(0.0, 0.0, -1.0)), // Back -Z
        (7, 6, 2, 3, Vector::new(0.0, 1.0, 0.0)),  // Top +Y
        (0, 1, 5, 4, Vector::new(0.0, -1.0, 0.0)), // Bottom -Y
        (5, 1, 2, 6, Vector::new(1.0, 0.0, 0.0)),  // Right +X
        (0, 4, 7, 3, Vector::new(-1.0, 0.0, 0.0)), // Left -X
    ];

    let mut positions = Vec::with_capacity(24);
    let mut normals = Vec::with_capacity(24);

    for (a, b, c, d, base_normal) in &faces {
        positions.push(to_vec3(&verts[*a]));
        positions.push(to_vec3(&verts[*b]));
        positions.push(to_vec3(&verts[*c]));
        positions.push(to_vec3(&verts[*d]));

        // Transform normal by rotor (same rotation as vertices)
        let rotated_normal = rotor.transform(base_normal);
        let normal = to_vec3(&rotated_normal);
        normals.push(normal);
        normals.push(normal);
        normals.push(normal);
        normals.push(normal);
    }

    let indices: Vec<u32> = (0..6)
        .flat_map(|face| {
            let base = face * 4;
            [base, base + 1, base + 2, base, base + 2, base + 3]
        })
        .collect();

    let mesh = CpuMesh {
        positions: Positions::F32(positions),
        indices: Indices::U32(indices),
        normals: Some(normals),
        ..Default::default()
    };

    Gm::new(
        Mesh::new(context, &mesh),
        PhysicalMaterial::new_opaque(
            context,
            &CpuMaterial {
                albedo: color,
                roughness: 0.4,
                metallic: 0.1,
                ..Default::default()
            },
        ),
    )
}

/// Create arrow from clifford Vector direction
fn create_arrow_from_vector(
    context: &Context,
    direction: &Vector<f32>,
    color: Srgba,
) -> Gm<Mesh, ColorMaterial> {
    let length = direction.norm();
    if length < 1e-6 {
        // Degenerate - return tiny mesh
        return Gm::new(
            Mesh::new(context, &CpuMesh::sphere(4)),
            ColorMaterial {
                color,
                ..Default::default()
            },
        );
    }

    let dir = to_vec3(&direction.normalized());
    let radius = 0.02;

    let mut mesh = CpuMesh::cylinder(8);
    mesh.transform(Mat4::from_nonuniform_scale(radius, length * 0.5, radius))
        .unwrap();

    // Rotate cylinder to point in direction
    let up = vec3(0.0, 1.0, 0.0);
    let rotation = if (up.dot(dir) - 1.0).abs() < 1e-6 {
        Mat4::identity()
    } else if (up.dot(dir) + 1.0).abs() < 1e-6 {
        Mat4::from_angle_z(radians(std::f32::consts::PI))
    } else {
        let axis = up.cross(dir).normalize();
        let angle = up.dot(dir).acos();
        Mat4::from_axis_angle(axis, radians(angle))
    };

    let translation = Mat4::from_translation(to_vec3(direction) * 0.5);
    mesh.transform(rotation).unwrap();
    mesh.transform(translation).unwrap();

    Gm::new(
        Mesh::new(context, &mesh),
        ColorMaterial {
            color,
            ..Default::default()
        },
    )
}

/// Create a disk representing the rotation plane (perpendicular to axis)
fn create_plane_disk(context: &Context, axis: &Vector<f32>) -> Gm<Mesh, ColorMaterial> {
    let segments = 32;
    let radius = 1.2;

    // Find two vectors perpendicular to axis using clifford cross product
    let perp1 = if axis.x().abs() < 0.9 {
        axis.cross(Vector::unit_x()).normalized()
    } else {
        axis.cross(Vector::unit_y()).normalized()
    };
    let perp2 = axis.cross(perp1).normalized();

    // Create circle vertices
    let mut positions = vec![vec3(0.0, 0.0, 0.0)]; // center
    for i in 0..segments {
        let angle = (i as f32 / segments as f32) * std::f32::consts::TAU;
        let (sin, cos) = angle.sin_cos();
        let point = perp1 * cos * radius + perp2 * sin * radius;
        positions.push(to_vec3(&point));
    }

    // Fan triangles
    let indices: Vec<u32> = (0..segments)
        .flat_map(|i| {
            let next = (i + 1) % segments;
            [0, i + 1, next + 1]
        })
        .collect();

    let mesh = CpuMesh {
        positions: Positions::F32(positions),
        indices: Indices::U32(indices),
        ..Default::default()
    };

    Gm::new(
        Mesh::new(context, &mesh),
        ColorMaterial {
            color: Srgba::new(180, 100, 255, 80), // Semi-transparent purple
            ..Default::default()
        },
    )
}
