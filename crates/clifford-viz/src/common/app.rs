//! Base application trait and runner for visualization demos.
//!
//! Provides a common structure for all visualization examples with:
//! - A main visualization panel
//! - A controls sidebar
//! - Optional info panel
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
}

/// Internal wrapper that adapts a [`VisualizationApp`] to eframe's [`eframe::App`].
struct AppWrapper<T: VisualizationApp> {
    /// The wrapped visualization app.
    app: T,
    /// Timestamp of the last frame, used to compute delta time.
    last_time: Option<std::time::Instant>,
}

impl<T: VisualizationApp + Default> Default for AppWrapper<T> {
    fn default() -> Self {
        Self {
            app: T::default(),
            last_time: None,
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
            });

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
