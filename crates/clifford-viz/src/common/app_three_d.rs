//! Three-d based application framework for visualization demos.
//!
//! This module provides two runners:
//! - [`run_three_d_app`] - for 2D demos using egui_plot
//! - [`run_three_d_app_3d`] - for native 3D demos using three-d rendering
//!
//! This module is only available with the `three-d` feature.

use super::app::{
    EducationalContent, VisualizationApp, configure_responsive_style, use_mobile_layout,
};
use three_d::*;

// =============================================================================
// 3D Visualization Trait
// =============================================================================

/// Extended trait for demos that use native three-d 3D rendering.
///
/// This trait extends [`VisualizationApp`] with methods for GPU-accelerated
/// 3D rendering. Implement this for demos that render actual 3D geometry
/// instead of 2D plots.
///
/// # Example
///
/// ```ignore
/// struct My3DDemo {
///     camera: Option<Camera>,
///     cube: Option<Gm<Mesh, ColorMaterial>>,
///     // ... other state
/// }
///
/// impl VisualizationApp for My3DDemo {
///     // ... standard methods for name, update, render (unused), controls
/// }
///
/// impl VisualizationApp3D for My3DDemo {
///     fn init_3d(&mut self, context: &Context) {
///         self.camera = Some(Camera::new_perspective(...));
///         self.cube = Some(Gm::new(Mesh::new(context, &CpuMesh::cube()), ...));
///     }
///
///     fn render_3d(&mut self, frame: &mut FrameInput) {
///         // Update camera, apply transformations, render objects
///     }
/// }
/// ```
pub trait VisualizationApp3D: VisualizationApp {
    /// Initialize 3D resources when the graphics context becomes available.
    ///
    /// This is called once before the first frame. Use it to create:
    /// - Camera and orbit controls
    /// - Meshes and materials
    /// - Lights
    ///
    /// # Arguments
    /// * `context` - The three-d graphics context for creating GPU resources
    fn init_3d(&mut self, context: &Context);

    /// Render 3D objects to the screen.
    ///
    /// This is called each frame after `update()` and before the egui overlay.
    /// The screen has already been cleared.
    ///
    /// Typical implementation:
    /// 1. Update camera viewport from `frame.viewport`
    /// 2. Handle camera events from `frame.events`
    /// 3. Apply transformations to objects
    /// 4. Call `frame.screen().render(&camera, objects, &lights)`
    ///
    /// # Arguments
    /// * `frame` - The current frame input with viewport, events, and screen
    fn render_3d(&mut self, frame: &mut FrameInput);
}

/// Run a 3D visualization app using native three-d rendering.
///
/// This function is similar to [`run_three_d_app`] but calls the
/// [`VisualizationApp3D`] methods for GPU-accelerated 3D rendering
/// instead of relying on egui_plot.
///
/// The render loop:
/// 1. Calculates delta time
/// 2. Calls `update(dt)`
/// 3. Clears the screen
/// 4. Calls `render_3d(frame)` for native 3D rendering
/// 5. Renders egui overlay (controls sidebar, info panel, learn popup)
///
/// # Type Parameters
/// * `T` - A type implementing [`VisualizationApp3D`] and [`Default`]
pub fn run_three_d_app_3d<T: VisualizationApp3D + Default + 'static>(title: &str) {
    #[cfg(not(target_arch = "wasm32"))]
    let max_size = Some((1920, 1080));
    #[cfg(target_arch = "wasm32")]
    let max_size = None;

    let window = Window::new(WindowSettings {
        title: title.to_string(),
        max_size,
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();
    let mut gui = GUI::new(&context);

    let mut demo = T::default();
    let mut learn_window_open = false;
    let mut sidebar_open = false;
    let mut last_time = 0.0;
    let mut initialized = false;

    window.render_loop(move |mut frame_input| {
        // Calculate delta time (accumulated_time is in milliseconds, convert to seconds)
        let current_time = frame_input.accumulated_time;
        let dt = if last_time > 0.0 {
            (((current_time - last_time) / 1000.0) as f32).min(0.1)
        } else {
            0.0
        };
        last_time = current_time;

        // Initialize 3D resources on first frame
        if !initialized {
            demo.init_3d(&context);
            initialized = true;
        }

        // Update demo logic
        demo.update(dt);

        // Clear screen with dark background
        frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.1, 0.1, 0.12, 1.0, 1.0));

        // Process egui FIRST to consume events over UI elements.
        // This prevents OrbitControl from capturing slider/button interactions.
        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |ctx| {
                render_app_ui_3d(ctx, &mut demo, &mut learn_window_open, &mut sidebar_open);
            },
        );

        // Render 3D objects (after egui consumes UI events)
        demo.render_3d(&mut frame_input);

        // Render egui to screen (on top of 3D)
        frame_input
            .screen()
            .write(|| gui.render())
            .expect("Failed to render GUI");

        FrameOutput::default()
    });
}

/// Render UI for 3D demos (no central panel - 3D renders behind).
fn render_app_ui_3d<T: VisualizationApp>(
    ctx: &egui::Context,
    demo: &mut T,
    learn_window_open: &mut bool,
    sidebar_open: &mut bool,
) {
    configure_responsive_style(ctx);
    let is_mobile = use_mobile_layout(ctx);

    // Mobile menu button
    if is_mobile && !*sidebar_open {
        egui::Area::new(egui::Id::new("menu_button_3d"))
            .fixed_pos(egui::pos2(8.0, 8.0))
            .show(ctx, |ui| {
                if ui
                    .add(egui::Button::new("\u{2630}").min_size(egui::vec2(36.0, 36.0)))
                    .on_hover_text("Open controls")
                    .clicked()
                {
                    *sidebar_open = true;
                }
            });
    }

    // Desktop: permanent side panel (semi-transparent for 3D)
    if !is_mobile {
        egui::SidePanel::left("controls_3d")
            .resizable(true)
            .default_width(250.0)
            .max_width(400.0)
            .frame(
                egui::Frame::side_top_panel(&ctx.style())
                    .fill(egui::Color32::from_rgba_unmultiplied(30, 30, 35, 230)),
            )
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Controls");
                    ui.separator();
                    demo.controls(ui);

                    if demo.educational_content().is_some() {
                        ui.add_space(16.0);
                        ui.separator();
                        if ui
                            .button("\u{1f4d6} Learn About This")
                            .on_hover_text("Open educational explanation")
                            .clicked()
                        {
                            *learn_window_open = true;
                        }
                    }
                });
            });
    }

    // Mobile: overlay panel
    if is_mobile && *sidebar_open {
        let screen_height = ctx.screen_rect().height();
        egui::Window::new("Controls")
            .id(egui::Id::new("mobile_controls_3d"))
            .fixed_pos(egui::pos2(0.0, 0.0))
            .fixed_size(egui::vec2(240.0, screen_height))
            .title_bar(false)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    if ui
                        .add(egui::Button::new("\u{2715}").min_size(egui::vec2(36.0, 36.0)))
                        .clicked()
                    {
                        *sidebar_open = false;
                    }
                    ui.heading("Controls");
                });
                ui.separator();

                egui::ScrollArea::vertical().show(ui, |ui| {
                    demo.controls(ui);

                    if demo.educational_content().is_some() {
                        ui.add_space(16.0);
                        ui.separator();
                        if ui
                            .button("\u{1f4d6} Learn About This")
                            .on_hover_text("Open educational explanation")
                            .clicked()
                        {
                            *learn_window_open = true;
                        }
                    }
                });
            });
    }

    // Educational window
    if let Some(content) = demo.educational_content() {
        egui::Window::new(format!("\u{1f4d6} {}", content.title))
            .open(learn_window_open)
            .default_width(if is_mobile { 300.0 } else { 500.0 })
            .default_height(400.0)
            .resizable(true)
            .collapsible(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    render_educational_content(ui, &content);
                });
            });
    }

    // Info panel at bottom (semi-transparent for 3D)
    if demo.show_info_panel() {
        egui::TopBottomPanel::bottom("info_3d")
            .resizable(true)
            .default_height(if is_mobile { 40.0 } else { 60.0 })
            .frame(
                egui::Frame::side_top_panel(&ctx.style())
                    .fill(egui::Color32::from_rgba_unmultiplied(30, 30, 35, 230)),
            )
            .show(ctx, |ui| {
                if is_mobile {
                    egui::ScrollArea::horizontal().show(ui, |ui| {
                        demo.info(ui);
                    });
                } else {
                    demo.info(ui);
                }
            });
    }

    // NOTE: No CentralPanel - 3D rendering fills the background
}

/// Run a visualization app using three-d.
///
/// This is the three-d equivalent of the old eframe-based `run_app`.
///
/// # Type Parameters
/// * `T` - A type implementing [`VisualizationApp`] and [`Default`]
pub fn run_three_d_app<T: VisualizationApp + Default + 'static>(title: &str) {
    // On native, cap window size. On WASM, use full browser window for responsiveness.
    #[cfg(not(target_arch = "wasm32"))]
    let max_size = Some((1920, 1080));
    #[cfg(target_arch = "wasm32")]
    let max_size = None;

    let window = Window::new(WindowSettings {
        title: title.to_string(),
        max_size,
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();
    let mut gui = GUI::new(&context);

    let mut demo = T::default();
    let mut learn_window_open = false;
    let mut sidebar_open = false;
    let mut last_time = 0.0;

    window.render_loop(move |mut frame_input| {
        // Calculate delta time (accumulated_time is in milliseconds, convert to seconds)
        let current_time = frame_input.accumulated_time;
        let dt = if last_time > 0.0 {
            (((current_time - last_time) / 1000.0) as f32).min(0.1)
        } else {
            0.0
        };
        last_time = current_time;

        // Update demo logic
        demo.update(dt);

        // Clear screen with dark background
        frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.1, 0.1, 0.12, 1.0, 1.0));

        // Render egui
        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |ctx| {
                render_app_ui(ctx, &mut demo, &mut learn_window_open, &mut sidebar_open);
            },
        );

        // Render egui to screen
        frame_input
            .screen()
            .write(|| gui.render())
            .expect("Failed to render GUI");

        FrameOutput::default()
    });
}

/// Render the standard app UI layout (sidebar, info panel, central area).
fn render_app_ui<T: VisualizationApp>(
    ctx: &egui::Context,
    demo: &mut T,
    learn_window_open: &mut bool,
    sidebar_open: &mut bool,
) {
    // Apply responsive styling
    configure_responsive_style(ctx);

    let is_mobile = use_mobile_layout(ctx);

    // Mobile menu button
    if is_mobile && !*sidebar_open {
        egui::Area::new(egui::Id::new("menu_button"))
            .fixed_pos(egui::pos2(8.0, 8.0))
            .show(ctx, |ui| {
                if ui
                    .add(egui::Button::new("\u{2630}").min_size(egui::vec2(36.0, 36.0)))
                    .on_hover_text("Open controls")
                    .clicked()
                {
                    *sidebar_open = true;
                }
            });
    }

    // Desktop: permanent side panel
    if !is_mobile {
        egui::SidePanel::left("controls")
            .resizable(true)
            .default_width(250.0)
            .max_width(400.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Controls");
                    ui.separator();
                    demo.controls(ui);

                    if demo.educational_content().is_some() {
                        ui.add_space(16.0);
                        ui.separator();
                        if ui
                            .button("\u{1f4d6} Learn About This")
                            .on_hover_text("Open educational explanation")
                            .clicked()
                        {
                            *learn_window_open = true;
                        }
                    }
                });
            });
    }

    // Mobile: overlay panel
    if is_mobile && *sidebar_open {
        let screen_height = ctx.screen_rect().height();
        egui::Window::new("Controls")
            .id(egui::Id::new("mobile_controls"))
            .fixed_pos(egui::pos2(0.0, 0.0))
            .fixed_size(egui::vec2(240.0, screen_height))
            .title_bar(false)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    if ui
                        .add(egui::Button::new("\u{2715}").min_size(egui::vec2(36.0, 36.0)))
                        .clicked()
                    {
                        *sidebar_open = false;
                    }
                    ui.heading("Controls");
                });
                ui.separator();

                egui::ScrollArea::vertical().show(ui, |ui| {
                    demo.controls(ui);

                    if demo.educational_content().is_some() {
                        ui.add_space(16.0);
                        ui.separator();
                        if ui
                            .button("\u{1f4d6} Learn About This")
                            .on_hover_text("Open educational explanation")
                            .clicked()
                        {
                            *learn_window_open = true;
                        }
                    }
                });
            });
    }

    // Educational window
    if let Some(content) = demo.educational_content() {
        egui::Window::new(format!("\u{1f4d6} {}", content.title))
            .open(learn_window_open)
            .default_width(if is_mobile { 300.0 } else { 500.0 })
            .default_height(400.0)
            .resizable(true)
            .collapsible(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    render_educational_content(ui, &content);
                });
            });
    }

    // Info panel at bottom
    if demo.show_info_panel() {
        egui::TopBottomPanel::bottom("info")
            .resizable(true)
            .default_height(if is_mobile { 40.0 } else { 60.0 })
            .show(ctx, |ui| {
                if is_mobile {
                    egui::ScrollArea::horizontal().show(ui, |ui| {
                        demo.info(ui);
                    });
                } else {
                    demo.info(ui);
                }
            });
    }

    // Central panel: visualization
    egui::CentralPanel::default().show(ctx, |ui| {
        demo.render(ui);
    });
}

/// Render educational content with nice formatting.
fn render_educational_content(ui: &mut egui::Ui, content: &EducationalContent) {
    // Overview section
    ui.heading("Overview");
    ui.label(content.overview);
    ui.add_space(12.0);

    // Mathematical Background section
    ui.heading("Mathematical Background");
    let code_bg = if ui.ctx().style().visuals.dark_mode {
        egui::Color32::from_rgb(40, 40, 50)
    } else {
        egui::Color32::from_rgb(240, 240, 245)
    };
    egui::Frame::none()
        .fill(code_bg)
        .inner_margin(12.0)
        .outer_margin(4.0)
        .rounding(4.0)
        .show(ui, |ui| {
            ui.monospace(content.math_background);
        });
    ui.add_space(12.0);

    // How to Use section
    ui.heading("How to Use");
    ui.label(content.how_to_use);
    ui.add_space(12.0);

    // Key Concepts section
    ui.heading("Key Concepts");
    ui.label(content.key_concepts);

    // Resources section (if any)
    if !content.resources.is_empty() {
        ui.add_space(12.0);
        ui.heading("Resources");
        for (name, url) in content.resources {
            ui.hyperlink_to(*name, *url);
        }
    }
}
