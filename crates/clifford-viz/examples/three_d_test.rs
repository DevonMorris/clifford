//! Test example for three-d library integration with egui.
//!
//! Run with: `cargo run -p clifford-viz --example three_d_test --features three-d`
//!
//! This tests three-d for 3D rendering with egui UI overlay.

use three_d::*;

fn main() {
    // Create a window with an OpenGL context
    let window = Window::new(WindowSettings {
        title: "three-d + egui Test".to_string(),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();

    // Create egui GUI
    let mut gui = GUI::new(&context);

    // Create a perspective camera with orbit controls
    let mut camera = Camera::new_perspective(
        window.viewport(),
        vec3(3.0, 2.0, 4.0), // eye position
        vec3(0.0, 0.0, 0.0), // target
        vec3(0.0, 1.0, 0.0), // up
        degrees(45.0),       // fov
        0.1,                 // near
        100.0,               // far
    );
    let mut control = OrbitControl::new(vec3(0.0, 0.0, 0.0), 1.0, 10.0);

    // Create coordinate axes
    let axes = Axes::new(&context, 0.02, 1.0);

    // Lighting
    let ambient = AmbientLight::new(&context, 0.3, Srgba::WHITE);
    let light0 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, vec3(-1.0, -1.0, -1.0));
    let light1 = DirectionalLight::new(&context, 1.0, Srgba::WHITE, vec3(1.0, 1.0, 1.0));

    // UI state
    let mut rotation_angle: f32 = 0.0;
    let mut auto_rotate = false;

    // Render loop
    window.render_loop(move |mut frame_input| {
        // Handle camera controls (when not interacting with GUI)
        let mut gui_consumed = false;

        // Update GUI
        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |egui_ctx| {
                gui_consumed = egui_ctx.wants_pointer_input() || egui_ctx.wants_keyboard_input();

                egui::Window::new("3D Rotor Demo")
                    .default_pos([10.0, 10.0])
                    .show(egui_ctx, |ui| {
                        ui.heading("Rotation Controls");

                        ui.horizontal(|ui| {
                            ui.label("Angle:");
                            ui.add(
                                egui::Slider::new(&mut rotation_angle, 0.0..=360.0).suffix(" deg"),
                            );
                        });

                        ui.checkbox(&mut auto_rotate, "Auto-rotate");

                        ui.separator();

                        ui.label("Camera: drag to orbit, scroll to zoom");

                        if ui.button("Reset View").clicked() {
                            camera.set_view(
                                vec3(3.0, 2.0, 4.0),
                                vec3(0.0, 0.0, 0.0),
                                vec3(0.0, 1.0, 0.0),
                            );
                        }
                    });
            },
        );

        // Handle camera only if GUI didn't consume the input
        if !gui_consumed {
            control.handle_events(&mut camera, &mut frame_input.events);
        }

        // Auto-rotate
        if auto_rotate {
            rotation_angle += frame_input.elapsed_time as f32 * 0.05;
            if rotation_angle > 360.0 {
                rotation_angle -= 360.0;
            }
        }

        // Apply rotation to cube
        let angle_rad = rotation_angle.to_radians();
        let rotation = Mat4::from_angle_y(radians(angle_rad));

        // Create rotated cube for this frame
        let mut rotated_mesh = CpuMesh::cube();
        rotated_mesh.transform(Mat4::from_scale(0.5)).unwrap();
        rotated_mesh.transform(rotation).unwrap();

        let rotated_cube = Gm::new(
            Mesh::new(&context, &rotated_mesh),
            PhysicalMaterial::new_opaque(
                &context,
                &CpuMaterial {
                    albedo: Srgba::new_opaque(100, 200, 100),
                    ..Default::default()
                },
            ),
        );

        // Clear and render 3D scene
        let _ = frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.1, 0.1, 0.15, 1.0, 1.0))
            .render(
                &camera,
                axes.into_iter().chain(&rotated_cube),
                &[&ambient, &light0, &light1],
            )
            .write(|| gui.render());

        FrameOutput::default()
    });
}
