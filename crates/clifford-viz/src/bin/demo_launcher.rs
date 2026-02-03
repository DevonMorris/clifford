//! Three-d based demo launcher for web and native.
//!
//! This launcher uses three-d instead of eframe, enabling unified rendering
//! for both 2D and 3D demos.
//!
//! ## Running locally
//!
//! ```bash
//! cargo run -p clifford-viz --bin demo_launcher --features three-d --release
//! ```
//!
//! ## Building for web
//!
//! ```bash
//! cd crates/clifford-viz
//! trunk build --release --config Trunk.three-d.toml
//! ```

use clifford_viz::common::app::{
    EducationalContent, VisualizationApp, configure_responsive_style, use_mobile_layout,
};
use clifford_viz::common::app_three_d::run_three_d_app_3d;
use clifford_viz::demos::*;
use three_d::*;

fn main() {
    // Get demo name from URL query parameter (web) or command line (native)
    let demo_name = get_demo_name();

    match demo_name.as_str() {
        "euclidean2" => run_demo::<Euclidean2Demo>("Euclidean 2D - Rotor Animation"),
        "euclidean3" => run_three_d_app_3d::<Euclidean3Demo>("3D Euclidean Rotors"),
        "projective2" => run_demo::<Projective2Demo>("Projective 2D - Point-Line Geometry"),
        "projective2_robot" => run_demo::<RobotArmDemo>("Robot Arm - 2D PGA"),
        "projective3_motor" => {
            run_three_d_app_3d::<Projective3MotorDemo>("3D PGA Motors - Rigid Body Motion")
        }
        "projective3_lines" => run_three_d_app_3d::<Projective3LinesDemo>("3D PGA Plucker Lines"),
        "projective3_robot" => run_three_d_app_3d::<Projective3RobotDemo>("3D PGA Robot Arm"),
        "projective3_geometry" => {
            run_three_d_app_3d::<Projective3GeometryDemo>("3D PGA Point-Line-Plane")
        }
        "conformal2_circles" => run_demo::<Conformal2CirclesDemo>("Conformal 2D - Circles"),
        "conformal2_inversion" => {
            run_demo::<Conformal2InversionDemo>("Conformal 2D - Circle Inversion")
        }
        "conformal2_mobius" => {
            run_demo::<Conformal2MobiusDemo>("Conformal 2D - Mobius Transformations")
        }
        "conformal2_intersection" => {
            run_demo::<Conformal2IntersectionDemo>("Conformal 2D - Circle Intersection")
        }
        "conformal3_inversion" => {
            run_three_d_app_3d::<Conformal3InversionDemo>("Conformal 3D - Sphere Inversion")
        }
        "conformal3_circles" => {
            run_three_d_app_3d::<Conformal3CirclesDemo>("Conformal 3D - Circle Operations")
        }
        "complex_domain" => run_demo::<ComplexDomainDemo>("Complex Domain Coloring"),
        "complex_fractal" => run_demo::<ComplexFractalDemo>("Complex Fractals"),
        "dual_autodiff" => run_demo::<DualAutodiffDemo>("Dual Number Autodiff"),
        "minkowski2_diagram" => run_demo::<Minkowski2DiagramDemo>("Minkowski 1+1D Diagram"),
        "minkowski2_dilation" => run_demo::<Minkowski2DilationDemo>("Time Dilation"),
        "menu" | "" => run_menu(),
        _ => {
            eprintln!("Unknown demo: {}", demo_name);
            eprintln!("Available demos: euclidean2, euclidean3, projective2, projective2_robot,");
            eprintln!("  projective3_motor, projective3_lines, conformal2_circles,");
            eprintln!("  conformal2_inversion, conformal2_mobius, conformal2_intersection,");
            eprintln!(
                "  conformal3_inversion, conformal3_circles, complex_domain, complex_fractal,"
            );
            eprintln!("  dual_autodiff,");
            eprintln!("  minkowski2_diagram, minkowski2_dilation");
            run_menu();
        }
    }
}

/// Get demo name from URL (web) or default to menu.
fn get_demo_name() -> String {
    #[cfg(target_arch = "wasm32")]
    {
        if let Some(window) = web_sys::window() {
            if let Ok(search) = window.location().search() {
                if let Some(demo) = search.trim_start_matches('?').split('&').find_map(|param| {
                    let mut parts = param.splitn(2, '=');
                    if parts.next() == Some("demo") {
                        parts.next().map(String::from)
                    } else {
                        None
                    }
                }) {
                    return demo;
                }
            }
        }
        "menu".to_string()
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        std::env::args()
            .nth(1)
            .unwrap_or_else(|| "menu".to_string())
    }
}

/// Run a demo with the three-d framework.
fn run_demo<T: VisualizationApp + Default + 'static>(title: &str) {
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

        // Clear screen
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
                render_demo_ui(ctx, &mut demo, &mut learn_window_open, &mut sidebar_open);
            },
        );

        frame_input
            .screen()
            .write(|| gui.render())
            .expect("Failed to render GUI");

        FrameOutput::default()
    });
}

/// State for the menu screen.
struct MenuState {
    logo_texture: Option<egui::TextureHandle>,
    logo_load_attempted: bool,
}

impl MenuState {
    fn new() -> Self {
        Self {
            logo_texture: None,
            logo_load_attempted: false,
        }
    }

    fn get_logo(&mut self, ctx: &egui::Context) -> Option<&egui::TextureHandle> {
        if !self.logo_load_attempted {
            self.logo_load_attempted = true;
            let image_bytes = include_bytes!("../../clifford.png");
            if let Ok(image) = image::load_from_memory(image_bytes) {
                let image = image.to_rgba8();
                let size = [image.width() as usize, image.height() as usize];
                let pixels = image.into_raw();
                let color_image = egui::ColorImage::from_rgba_unmultiplied(size, &pixels);
                self.logo_texture = Some(ctx.load_texture(
                    "clifford-logo",
                    color_image,
                    egui::TextureOptions::LINEAR,
                ));
            }
        }
        self.logo_texture.as_ref()
    }
}

/// Run the demo menu.
fn run_menu() {
    // On native, cap window size. On WASM, use full browser window for responsiveness.
    #[cfg(not(target_arch = "wasm32"))]
    let max_size = Some((1920, 1080));
    #[cfg(target_arch = "wasm32")]
    let max_size = None;

    let window = Window::new(WindowSettings {
        title: "Clifford - Geometric Algebra Demos".to_string(),
        max_size,
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();
    let mut gui = GUI::new(&context);
    let mut menu_state = MenuState::new();

    window.render_loop(move |mut frame_input| {
        frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.1, 0.1, 0.12, 1.0, 1.0));

        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |ctx| {
                render_menu_ui(ctx, &mut menu_state);
            },
        );

        frame_input
            .screen()
            .write(|| gui.render())
            .expect("Failed to render GUI");

        FrameOutput::default()
    });
}

/// Render the demo UI with sidebar, info panel, and central area.
fn render_demo_ui<T: VisualizationApp>(
    ctx: &egui::Context,
    demo: &mut T,
    learn_window_open: &mut bool,
    sidebar_open: &mut bool,
) {
    configure_responsive_style(ctx);
    let is_mobile = use_mobile_layout(ctx);

    // Mobile menu button
    if is_mobile && !*sidebar_open {
        egui::Area::new(egui::Id::new("menu_button"))
            .fixed_pos(egui::pos2(8.0, 8.0))
            .show(ctx, |ui| {
                if ui
                    .add(egui::Button::new("\u{2630}").min_size(egui::vec2(36.0, 36.0)))
                    .clicked()
                {
                    *sidebar_open = true;
                }
            });
    }

    // Desktop sidebar
    if !is_mobile {
        egui::SidePanel::left("controls")
            .resizable(true)
            .default_width(250.0)
            .max_width(400.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    // Back to menu link
                    if ui.link("<- Back to Menu").clicked() {
                        navigate_to_demo("menu");
                    }
                    ui.separator();

                    ui.heading("Controls");
                    ui.separator();
                    demo.controls(ui);

                    if demo.educational_content().is_some() {
                        ui.add_space(16.0);
                        ui.separator();
                        if ui.button("\u{1f4d6} Learn About This").clicked() {
                            *learn_window_open = true;
                        }
                    }
                });
            });
    }

    // Mobile overlay sidebar
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
                    if ui.link("<- Back to Menu").clicked() {
                        navigate_to_demo("menu");
                    }
                    ui.separator();

                    demo.controls(ui);

                    if demo.educational_content().is_some() {
                        ui.add_space(16.0);
                        ui.separator();
                        if ui.button("\u{1f4d6} Learn About This").clicked() {
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

    // Info panel
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

    // Central panel
    egui::CentralPanel::default().show(ctx, |ui| {
        demo.render(ui);
    });
}

/// Render the menu UI.
fn render_menu_ui(ctx: &egui::Context, state: &mut MenuState) {
    configure_responsive_style(ctx);
    let is_mobile = use_mobile_layout(ctx);
    let sp = if is_mobile { 1.0 } else { 1.5 };

    // Get logo texture (may fail on some platforms)
    let logo = state.get_logo(ctx);
    let logo_size = if is_mobile {
        egui::vec2(80.0, 80.0)
    } else {
        egui::vec2(120.0, 120.0)
    };

    egui::CentralPanel::default().show(ctx, |ui| {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.vertical_centered(|ui| {
                ui.add_space(24.0 * sp);

                // Logo (only show if loaded successfully)
                if let Some(logo) = logo {
                    ui.add(egui::Image::new(logo).fit_to_exact_size(logo_size));
                    ui.add_space(8.0 * sp);
                }

                // Title
                ui.heading(egui::RichText::new("Clifford").strong());
                ui.label("Geometric Algebra for Rust");

                ui.add_space(12.0 * sp);

                // Links (use explicit URL opening for WASM compatibility)
                if is_mobile {
                    if ui.link("Crates.io").clicked() {
                        open_external_url("https://crates.io/crates/clifford");
                    }
                    if ui.link("Docs.rs").clicked() {
                        open_external_url("https://docs.rs/clifford");
                    }
                    if ui.link("GitHub").clicked() {
                        open_external_url("https://github.com/DevonMorris/clifford");
                    }
                } else {
                    ui.horizontal(|ui| {
                        if ui.link("Crates.io").clicked() {
                            open_external_url("https://crates.io/crates/clifford");
                        }
                        if ui.link("Docs.rs").clicked() {
                            open_external_url("https://docs.rs/clifford");
                        }
                        if ui.link("GitHub").clicked() {
                            open_external_url("https://github.com/DevonMorris/clifford");
                        }
                    });
                }

                ui.add_space(16.0 * sp);
            });

            ui.separator();
            ui.add_space(12.0 * sp);

            ui.vertical_centered(|ui| {
                ui.label(egui::RichText::new("Interactive Demos").strong());
                ui.label("Learn Geometric Algebra through visualization");
            });

            ui.add_space(12.0 * sp);

            render_demo_category(
                ui,
                "Euclidean Geometry",
                &[
                    ("euclidean2", "2D Rotors", "Rotation using rotors"),
                    ("euclidean3", "3D Rotors", "Native 3D rotation demo"),
                ],
            );

            ui.add_space(12.0 * sp);

            render_demo_category(
                ui,
                "Projective Geometry (PGA)",
                &[
                    (
                        "projective2",
                        "2D Point-Line",
                        "Join/meet operations in 2D PGA",
                    ),
                    (
                        "projective2_robot",
                        "Robot Arm",
                        "Forward kinematics with motors",
                    ),
                    (
                        "projective3_motor",
                        "3D Motors",
                        "Rigid body motion in 3D PGA",
                    ),
                    (
                        "projective3_lines",
                        "3D Plucker Lines",
                        "Line geometry with Plucker coords",
                    ),
                    (
                        "projective3_robot",
                        "3D Robot Arm",
                        "Forward kinematics with motors",
                    ),
                    (
                        "projective3_geometry",
                        "3D Point-Line-Plane",
                        "Join and meet operations",
                    ),
                ],
            );

            ui.add_space(12.0 * sp);

            render_demo_category(
                ui,
                "Conformal Geometry (CGA)",
                &[
                    (
                        "conformal2_circles",
                        "Circle from 3 Points",
                        "Create circles via wedge product",
                    ),
                    (
                        "conformal2_inversion",
                        "Circle Inversion",
                        "Conformal transformation",
                    ),
                    (
                        "conformal2_mobius",
                        "Mobius Transforms",
                        "Compose transformations",
                    ),
                    (
                        "conformal2_intersection",
                        "Circle Intersection",
                        "Find intersection points",
                    ),
                    (
                        "conformal3_inversion",
                        "3D Sphere Inversion",
                        "Sphere inversion in 3D CGA",
                    ),
                    (
                        "conformal3_circles",
                        "3D Circle from Points",
                        "Create circles in 3D CGA",
                    ),
                ],
            );

            ui.add_space(12.0 * sp);

            render_demo_category(
                ui,
                "Minkowski Spacetime",
                &[
                    (
                        "minkowski2_diagram",
                        "1+1D Diagram",
                        "Light cones and Lorentz boosts",
                    ),
                    (
                        "minkowski2_dilation",
                        "Time Dilation",
                        "The geometry of time dilation",
                    ),
                ],
            );

            ui.add_space(12.0 * sp);

            render_demo_category(
                ui,
                "Complex & Dual Numbers",
                &[
                    (
                        "complex_domain",
                        "Domain Coloring",
                        "Visualize complex functions",
                    ),
                    (
                        "complex_fractal",
                        "Mandelbrot & Julia",
                        "Explore complex fractals",
                    ),
                    (
                        "dual_autodiff",
                        "Dual Autodiff",
                        "Automatic differentiation",
                    ),
                ],
            );

            ui.add_space(24.0 * sp);
            ui.separator();
            ui.add_space(8.0);

            ui.vertical_centered(|ui| {
                ui.label("Built with Rust + three-d + egui");
            });

            ui.add_space(16.0);
        });
    });
}

/// Render a category of demos.
fn render_demo_category(ui: &mut egui::Ui, title: &str, demos: &[(&str, &str, &str)]) {
    ui.heading(title);
    ui.add_space(8.0);

    for (id, name, description) in demos {
        ui.horizontal(|ui| {
            if ui.link(*name).clicked() {
                navigate_to_demo(id);
            }
            ui.label(" - ");
            ui.label(*description);
        });
        ui.add_space(4.0);
    }
}

/// Open an external URL in a new tab.
///
/// In WASM, egui hyperlinks don't always work with three-d's canvas,
/// so we create an anchor element and click it. This works better on
/// mobile browsers which block window.open() when not triggered directly
/// by a user gesture.
fn open_external_url(url: &str) {
    #[cfg(target_arch = "wasm32")]
    {
        use wasm_bindgen::JsCast;
        if let Some(window) = web_sys::window() {
            if let Some(document) = window.document() {
                if let Ok(anchor) = document.create_element("a") {
                    let anchor: web_sys::HtmlAnchorElement = anchor.unchecked_into();
                    anchor.set_href(url);
                    anchor.set_target("_blank");
                    anchor.set_rel("noopener noreferrer");
                    anchor.click();
                }
            }
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        // On native, print the URL (or could use opener crate)
        eprintln!("Open in browser: {}", url);
    }
}

/// Navigate to a demo.
fn navigate_to_demo(demo_id: &str) {
    #[cfg(target_arch = "wasm32")]
    {
        if let Some(window) = web_sys::window() {
            let url = format!("?demo={}", demo_id);
            let _ = window.location().set_href(&url);
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        eprintln!(
            "To run: cargo run -p clifford-viz --bin demo_launcher --features three-d -- {}",
            demo_id
        );
    }
}

/// Render educational content.
fn render_educational_content(ui: &mut egui::Ui, content: &EducationalContent) {
    ui.heading("Overview");
    ui.label(content.overview);
    ui.add_space(12.0);

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

    ui.heading("How to Use");
    ui.label(content.how_to_use);
    ui.add_space(12.0);

    ui.heading("Key Concepts");
    ui.label(content.key_concepts);

    if !content.resources.is_empty() {
        ui.add_space(12.0);
        ui.heading("Resources");
        for (name, url) in content.resources {
            if ui.link(*name).clicked() {
                open_external_url(url);
            }
        }
    }
}
