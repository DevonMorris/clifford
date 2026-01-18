//! Landing page menu for selecting demos.
//!
//! This module provides the `DemoMenu` component which displays a list of
//! available demos and allows navigation between them in the web version.

use crate::common::colors;
use egui::ImageSource;

/// Demo entry with metadata for display.
struct DemoEntry {
    /// URL parameter name (e.g., "euclidean2").
    id: &'static str,
    /// Display name shown in the menu.
    name: &'static str,
    /// Brief description of what the demo shows.
    description: &'static str,
}

/// Euclidean geometry demos.
const EUCLIDEAN_DEMOS: &[DemoEntry] = &[DemoEntry {
    id: "euclidean2",
    name: "2D Rotors",
    description: "Rotation and dilation using rotors in 2D Euclidean space",
}];

/// Projective geometry demos.
const PROJECTIVE_DEMOS: &[DemoEntry] = &[
    DemoEntry {
        id: "projective2",
        name: "2D Point-Line Geometry",
        description: "Join/meet operations, motors, and transformations in 2D PGA",
    },
    DemoEntry {
        id: "projective2_robot",
        name: "Robot Arm (2D PGA)",
        description: "Forward kinematics with motor composition",
    },
];

/// Landing page showing available demos.
///
/// This implements `eframe::App` directly rather than `VisualizationApp`
/// since it's just a navigation menu, not a visualization demo.
pub struct DemoMenu;

impl DemoMenu {
    /// Create a new demo menu.
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self
    }
}

impl eframe::App for DemoMenu {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Install image loaders for egui (needed for include_image!)
        egui_extras::install_image_loaders(ctx);

        let is_mobile = ctx.screen_rect().width() < 600.0;

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add_space(if is_mobile { 16.0 } else { 32.0 });

                    // Logo
                    let logo_size = if is_mobile { 80.0 } else { 120.0 };
                    let logo: ImageSource<'_> = egui::include_image!("../../clifford.png");
                    ui.add(egui::Image::new(logo).max_size(egui::vec2(logo_size, logo_size)));

                    ui.add_space(if is_mobile { 8.0 } else { 12.0 });

                    // Title
                    ui.heading(egui::RichText::new("Clifford").size(if is_mobile { 28.0 } else { 36.0 }));
                    ui.label(egui::RichText::new("Geometric Algebra for Rust").size(if is_mobile { 14.0 } else { 16.0 }));

                    ui.add_space(if is_mobile { 12.0 } else { 16.0 });

                    // External links
                    ui.horizontal(|ui| {
                        ui.spacing_mut().item_spacing.x = if is_mobile { 8.0 } else { 16.0 };
                        ui.hyperlink_to("📦 Crates.io", "https://crates.io/crates/clifford");
                        ui.hyperlink_to("📚 Docs.rs", "https://docs.rs/clifford");
                        ui.hyperlink_to("🔗 GitHub", "https://github.com/DevonMorris/clifford");
                    });

                    ui.add_space(if is_mobile { 16.0 } else { 24.0 });
                });

                ui.separator();
                ui.add_space(if is_mobile { 12.0 } else { 16.0 });

                ui.vertical_centered(|ui| {
                    ui.label(egui::RichText::new("Interactive Demos").size(if is_mobile { 18.0 } else { 20.0 }).strong());
                    let desc_color = colors::text_secondary(ui.ctx());
                    ui.colored_label(desc_color, "Learn Geometric Algebra through visualization");
                });

                ui.add_space(if is_mobile { 12.0 } else { 16.0 });

                render_demo_category(ui, "Euclidean Geometry", EUCLIDEAN_DEMOS, is_mobile);
                ui.add_space(if is_mobile { 12.0 } else { 16.0 });
                render_demo_category(ui, "Projective Geometry (PGA)", PROJECTIVE_DEMOS, is_mobile);

                ui.add_space(if is_mobile { 20.0 } else { 32.0 });
                ui.separator();
                ui.add_space(8.0);

                ui.vertical_centered(|ui| {
                    let footer_color = colors::text_secondary(ui.ctx());
                    ui.colored_label(footer_color, "Built with Rust + egui • MIT License");
                });

                ui.add_space(16.0);
            });
        });
    }
}

/// Render a category of demos with a header and list of entries.
fn render_demo_category(ui: &mut egui::Ui, title: &str, demos: &[DemoEntry], is_mobile: bool) {
    ui.heading(title);
    ui.add_space(8.0);

    let description_color = colors::text_secondary(ui.ctx());

    for demo in demos {
        if is_mobile {
            // Stack vertically on mobile
            if ui.link(demo.name).clicked() {
                navigate_to_demo(demo.id);
            }
            ui.colored_label(description_color, demo.description);
            ui.add_space(8.0);
        } else {
            // Horizontal on desktop
            ui.horizontal(|ui| {
                if ui.link(demo.name).clicked() {
                    navigate_to_demo(demo.id);
                }
                ui.label(" - ");
                ui.colored_label(description_color, demo.description);
            });
            ui.add_space(4.0);
        }
    }
}

/// Navigate to a specific demo.
///
/// On web, this updates the URL query parameter to load the demo.
/// On native, this prints instructions for running the demo.
fn navigate_to_demo(demo_id: &str) {
    #[cfg(target_arch = "wasm32")]
    {
        if let Some(window) = web_sys::window() {
            let url = format!("?demo={}", demo_id);
            if let Err(e) = window.location().set_href(&url) {
                log::error!("Failed to navigate: {:?}", e);
            }
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        // On native, just print instructions (menu is mainly for web)
        eprintln!(
            "To run this demo natively:\n  cargo run -p clifford-viz --example {} --release",
            demo_id
        );
    }
}
