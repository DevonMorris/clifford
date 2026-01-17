//! Base application trait and runner for visualization demos.
//!
//! Provides a common structure for all visualization examples with:
//! - A main visualization panel
//! - A controls sidebar
//! - Optional info panel (status bar)
//! - Educational popup window with math explanations
//! - Consistent layout and styling

use eframe::egui;

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
///     fn render(&self, ui: &mut egui::Ui) {
///         // Draw visualization
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
    fn render(&self, ui: &mut egui::Ui);

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

/// Internal wrapper that adapts a [`VisualizationApp`] to eframe's [`eframe::App`].
struct AppWrapper<T: VisualizationApp> {
    /// The wrapped visualization app.
    app: T,
    /// Timestamp of the last frame, used to compute delta time.
    last_time: Option<std::time::Instant>,
    /// Whether the educational window is open.
    learn_window_open: bool,
}

impl<T: VisualizationApp + Default> Default for AppWrapper<T> {
    fn default() -> Self {
        Self {
            app: T::default(),
            last_time: None,
            learn_window_open: false,
        }
    }
}

impl<T: VisualizationApp> eframe::App for AppWrapper<T> {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Calculate delta time
        let now = std::time::Instant::now();
        let dt = self
            .last_time
            .map(|t| now.duration_since(t).as_secs_f32())
            .unwrap_or(0.0);
        self.last_time = Some(now);

        // Update app logic
        self.app.update(dt);

        // Request continuous repaint for animations
        ctx.request_repaint();

        // Left panel: controls
        egui::SidePanel::left("controls")
            .resizable(true)
            .default_width(250.0)
            .show(ctx, |ui| {
                ui.heading("Controls");
                ui.separator();
                self.app.controls(ui);

                // Learn button at bottom of controls (if educational content exists)
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

        // Educational window (popup)
        if let Some(content) = self.app.educational_content() {
            egui::Window::new(format!("\u{1f4d6} {}", content.title))
                .open(&mut self.learn_window_open)
                .default_width(500.0)
                .default_height(400.0)
                .resizable(true)
                .collapsible(true)
                .show(ctx, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        render_educational_content(ui, &content);
                    });
                });
        }

        // Optional bottom panel: info
        if self.app.show_info_panel() {
            egui::TopBottomPanel::bottom("info")
                .resizable(true)
                .default_height(60.0)
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
fn render_educational_content(ui: &mut egui::Ui, content: &EducationalContent) {
    // Overview section
    ui.heading("Overview");
    ui.label(content.overview);
    ui.add_space(12.0);

    // Mathematical Background section
    ui.heading("Mathematical Background");
    // Render with monospace for formulas
    egui::Frame::none()
        .fill(egui::Color32::from_rgba_unmultiplied(40, 40, 50, 255))
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
/// # Arguments
/// * `options` - Custom native options for the window
///
/// # Errors
/// Returns an error if the window creation fails.
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
