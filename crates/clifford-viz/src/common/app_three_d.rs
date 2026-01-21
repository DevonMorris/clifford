//! Three-d based application framework for visualization demos.
//!
//! This module provides [`run_three_d_app`] which runs a [`VisualizationApp`]
//! using the three-d renderer instead of eframe. This enables unified rendering
//! for both 2D (via egui_plot) and 3D demos.
//!
//! This module is only available with the `three-d` feature.

use super::app::{
    configure_responsive_style, use_mobile_layout, EducationalContent, VisualizationApp,
};
use three_d::*;

/// Run a visualization app using three-d.
///
/// This is the three-d equivalent of [`super::app::run_app`].
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
        // Calculate delta time
        let current_time = frame_input.accumulated_time;
        let dt = if last_time > 0.0 {
            ((current_time - last_time) as f32).min(0.1)
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
                render_app_ui(
                    ctx,
                    &mut demo,
                    &mut learn_window_open,
                    &mut sidebar_open,
                );
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
