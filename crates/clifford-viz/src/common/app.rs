//! Base application trait and runner for visualization demos.
//!
//! Provides a common structure for all visualization examples with:
//! - A main visualization panel
//! - A controls sidebar
//! - Optional info panel (status bar)
//! - Educational popup window with math explanations
//! - Consistent layout and styling

use egui;

/// Base trait for visualization demos.
///
/// Implementations provide the specific visualization content while
/// this framework handles the window setup and layout.
///
/// # Example
/// ```ignore
/// struct MyDemo {
///     angle: f32,
/// }
///
/// impl VisualizationApp for MyDemo {
///     fn name(&self) -> &'static str { "My Demo" }
///
///     fn update(&mut self, _dt: f32) {
///         // Update animation state
///     }
///
///     fn render(&mut self, ui: &mut egui::Ui) {
///         // Draw visualization and handle mouse interactions
///     }
///
///     fn controls(&mut self, ui: &mut egui::Ui) {
///         // Add sliders, buttons, etc.
///     }
/// }
/// ```
pub trait VisualizationApp {
    /// Display name shown in window title.
    fn name(&self) -> &'static str;

    /// Update logic called every frame.
    ///
    /// # Arguments
    /// * `dt` - Delta time since last frame in seconds
    fn update(&mut self, dt: f32);

    /// Render the main visualization area.
    ///
    /// This is typically where you'd draw to a plot or custom painter.
    /// Takes `&mut self` to allow handling mouse interactions.
    fn render(&mut self, ui: &mut egui::Ui);

    /// Render the control panel (sliders, buttons, etc.).
    fn controls(&mut self, ui: &mut egui::Ui);

    /// Optional: render an info panel with explanatory text.
    ///
    /// Default implementation shows nothing.
    fn info(&self, _ui: &mut egui::Ui) {}

    /// Optional: whether to show the info panel.
    ///
    /// Default is `true`.
    fn show_info_panel(&self) -> bool {
        true
    }

    /// Educational content for the "Learn" popup window.
    ///
    /// This should provide a comprehensive explanation of the mathematical
    /// concepts being visualized. Return `None` to disable the Learn button.
    ///
    /// # Sections
    ///
    /// The content is organized into sections:
    /// - **Overview**: What this visualization demonstrates
    /// - **Mathematical Background**: The underlying math (formulas, definitions)
    /// - **How to Use**: Instructions for interacting with the demo
    /// - **Key Concepts**: Important takeaways
    ///
    /// Use markdown-like formatting (will be rendered as rich text):
    /// - `**bold**` for emphasis
    /// - Formulas in plain text with Unicode subscripts (e₁₂, θ/2, etc.)
    fn educational_content(&self) -> Option<EducationalContent> {
        None
    }
}

/// Structured educational content for visualization demos.
///
/// This content appears in the "Learn" popup window to help users
/// understand the mathematical concepts being visualized.
#[derive(Debug, Clone, Default)]
pub struct EducationalContent {
    /// Brief title for the educational content.
    pub title: &'static str,
    /// Overview of what the visualization demonstrates.
    pub overview: &'static str,
    /// Mathematical background and formulas.
    pub math_background: &'static str,
    /// Instructions for using the visualization.
    pub how_to_use: &'static str,
    /// Key concepts to understand.
    pub key_concepts: &'static str,
    /// Optional: links to external resources.
    pub resources: &'static [(&'static str, &'static str)],
}

/// Wrapper that adapts a [`VisualizationApp`] to eframe's [`eframe::App`].
///
/// This is used internally by [`run_app`] and externally by the WASM
/// demo launcher to wrap demos in the standard visualization framework.
///
/// Only available with the `native` feature or on WASM.
#[cfg(any(feature = "native", target_arch = "wasm32"))]
pub struct AppWrapper<T: VisualizationApp> {
    /// The wrapped visualization app.
    app: T,
    /// Whether the educational window is open.
    learn_window_open: bool,
    /// Whether the sidebar is open (for mobile toggle).
    sidebar_open: bool,
}

#[cfg(any(feature = "native", target_arch = "wasm32"))]
impl<T: VisualizationApp + Default> Default for AppWrapper<T> {
    fn default() -> Self {
        Self {
            app: T::default(),
            learn_window_open: false,
            sidebar_open: false, // Closed by default (desktop ignores this)
        }
    }
}

#[cfg(any(feature = "native", target_arch = "wasm32"))]
impl<T: VisualizationApp> eframe::App for AppWrapper<T> {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Use egui's built-in timing (works on both native and WASM)
        // Cap at 0.1s to handle edge cases (first frame, sleeping, etc.)
        let dt = ctx.input(|i| i.stable_dt.min(0.1));

        // Update app logic
        self.app.update(dt);

        // Request continuous repaint for animations
        ctx.request_repaint();

        // Responsive layout based on screen size
        let screen_width = ctx.screen_rect().width();
        let is_mobile = screen_width < 600.0;

        // On mobile, show a floating menu button when sidebar is closed
        if is_mobile && !self.sidebar_open {
            egui::Area::new(egui::Id::new("menu_button"))
                .fixed_pos(egui::pos2(8.0, 8.0))
                .show(ctx, |ui| {
                    if ui
                        .add(egui::Button::new("\u{2630}").min_size(egui::vec2(36.0, 36.0)))
                        .on_hover_text("Open controls")
                        .clicked()
                    {
                        self.sidebar_open = true;
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
                        self.app.controls(ui);

                        if self.app.educational_content().is_some() {
                            ui.add_space(16.0);
                            ui.separator();
                            if ui
                                .button("\u{1f4d6} Learn About This")
                                .on_hover_text("Open educational explanation")
                                .clicked()
                            {
                                self.learn_window_open = true;
                            }
                        }
                    });
                });
        }

        // Mobile: overlay panel (doesn't push content)
        if is_mobile && self.sidebar_open {
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
                        // Close button on left (same position as hamburger menu)
                        if ui
                            .add(egui::Button::new("\u{2715}").min_size(egui::vec2(36.0, 36.0)))
                            .clicked()
                        {
                            self.sidebar_open = false;
                        }
                        ui.heading("Controls");
                    });
                    ui.separator();

                    egui::ScrollArea::vertical().show(ui, |ui| {
                        self.app.controls(ui);

                        if self.app.educational_content().is_some() {
                            ui.add_space(16.0);
                            ui.separator();
                            if ui
                                .button("\u{1f4d6} Learn About This")
                                .on_hover_text("Open educational explanation")
                                .clicked()
                            {
                                self.learn_window_open = true;
                            }
                        }
                    });
                });
        }

        // Educational window (popup)
        if let Some(content) = self.app.educational_content() {
            egui::Window::new(format!("\u{1f4d6} {}", content.title))
                .open(&mut self.learn_window_open)
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

        // Optional bottom panel: info (smaller on mobile)
        if self.app.show_info_panel() {
            egui::TopBottomPanel::bottom("info")
                .resizable(true)
                .default_height(if is_mobile { 40.0 } else { 60.0 })
                .show(ctx, |ui| {
                    self.app.info(ui);
                });
        }

        // Central panel: visualization
        egui::CentralPanel::default().show(ctx, |ui| {
            self.app.render(ui);
        });
    }
}

/// Render educational content with nice formatting.
#[cfg(any(feature = "native", target_arch = "wasm32"))]
fn render_educational_content(ui: &mut egui::Ui, content: &EducationalContent) {
    // Overview section
    ui.heading("Overview");
    ui.label(content.overview);
    ui.add_space(12.0);

    // Mathematical Background section
    ui.heading("Mathematical Background");
    // Render with monospace for formulas (theme-aware background)
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

/// Run a visualization app with default window options.
///
/// This function is only available on native platforms (not WASM) and requires
/// the `native` feature to be enabled.
///
/// # Type Parameters
/// * `T` - A type implementing [`VisualizationApp`] and [`Default`]
///
/// # Errors
/// Returns an error if the window creation fails.
///
/// # Example
/// ```ignore
/// fn main() -> eframe::Result<()> {
///     run_app::<MyDemo>()
/// }
/// ```
#[cfg(all(feature = "native", not(target_arch = "wasm32")))]
pub fn run_app<T: VisualizationApp + Default + 'static>() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1024.0, 768.0])
            .with_min_inner_size([640.0, 480.0]),
        ..Default::default()
    };

    eframe::run_native(
        T::default().name(),
        options,
        Box::new(|_cc| Ok(Box::new(AppWrapper::<T>::default()))),
    )
}

/// Run a visualization app with custom options.
///
/// This function is only available on native platforms (not WASM) and requires
/// the `native` feature to be enabled.
///
/// # Arguments
/// * `options` - Custom native options for the window
///
/// # Errors
/// Returns an error if the window creation fails.
#[cfg(all(feature = "native", not(target_arch = "wasm32")))]
pub fn run_app_with_options<T: VisualizationApp + Default + 'static>(
    options: eframe::NativeOptions,
) -> eframe::Result<()> {
    eframe::run_native(
        T::default().name(),
        options,
        Box::new(|_cc| Ok(Box::new(AppWrapper::<T>::default()))),
    )
}

/// Configuration for window creation.
#[derive(Debug, Clone)]
pub struct WindowConfig {
    /// Window title
    pub title: String,
    /// Initial window size (width, height)
    pub size: (f32, f32),
    /// Minimum window size
    pub min_size: (f32, f32),
    /// Whether the window is resizable
    pub resizable: bool,
}

impl Default for WindowConfig {
    fn default() -> Self {
        Self {
            title: "Clifford Visualization".to_string(),
            size: (1024.0, 768.0),
            min_size: (640.0, 480.0),
            resizable: true,
        }
    }
}

impl WindowConfig {
    /// Create native options from this config.
    ///
    /// This method is only available on native platforms (not WASM) and requires
    /// the `native` feature to be enabled.
    #[cfg(all(feature = "native", not(target_arch = "wasm32")))]
    #[must_use]
    pub fn to_native_options(&self) -> eframe::NativeOptions {
        eframe::NativeOptions {
            viewport: egui::ViewportBuilder::default()
                .with_inner_size([self.size.0, self.size.1])
                .with_min_inner_size([self.min_size.0, self.min_size.1])
                .with_resizable(self.resizable),
            ..Default::default()
        }
    }
}
