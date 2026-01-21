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

// =============================================================================
// Responsive Styling
// =============================================================================

/// Screen size category for responsive design.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScreenSize {
    /// Mobile: < 600px
    Mobile,
    /// Tablet: 600-1024px
    Tablet,
    /// Desktop: > 1024px
    Desktop,
}

impl ScreenSize {
    /// Determine screen size category from width.
    #[must_use]
    pub fn from_width(width: f32) -> Self {
        if width < 600.0 {
            Self::Mobile
        } else if width < 1024.0 {
            Self::Tablet
        } else {
            Self::Desktop
        }
    }

    /// Check if this is a mobile screen.
    #[must_use]
    pub fn is_mobile(self) -> bool {
        self == Self::Mobile
    }
}

/// Check if mobile layout should be used.
///
/// Returns true if:
/// - Width < 600px (narrow screen), OR
/// - Height > Width (portrait orientation)
///
/// This ensures the hamburger menu appears on narrow screens
/// and when the window is taller than it is wide.
#[must_use]
pub fn use_mobile_layout(ctx: &egui::Context) -> bool {
    let rect = ctx.screen_rect();
    rect.width() < 600.0 || rect.height() > rect.width()
}

/// Configure responsive styles based on screen size.
///
/// Call this once per frame to set up font sizes and spacing that adapt
/// to the current screen width. This eliminates the need for manual
/// `if is_mobile` checks throughout the UI code.
///
/// # Example
/// ```ignore
/// fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
///     configure_responsive_style(ctx);
///     // ... rest of UI code uses normal egui widgets
/// }
/// ```
pub fn configure_responsive_style(ctx: &egui::Context) {
    use egui::{FontFamily, FontId, TextStyle};

    let screen_size = ScreenSize::from_width(ctx.screen_rect().width());

    // Base font size scales with screen size
    let base_size = match screen_size {
        ScreenSize::Mobile => 14.0,
        ScreenSize::Tablet => 15.0,
        ScreenSize::Desktop => 16.0,
    };

    let heading_size = match screen_size {
        ScreenSize::Mobile => 20.0,
        ScreenSize::Tablet => 24.0,
        ScreenSize::Desktop => 28.0,
    };

    // Configure text styles
    let text_styles: std::collections::BTreeMap<TextStyle, FontId> = [
        (
            TextStyle::Small,
            FontId::new(base_size * 0.75, FontFamily::Proportional),
        ),
        (
            TextStyle::Body,
            FontId::new(base_size, FontFamily::Proportional),
        ),
        (
            TextStyle::Button,
            FontId::new(base_size, FontFamily::Proportional),
        ),
        (
            TextStyle::Heading,
            FontId::new(heading_size, FontFamily::Proportional),
        ),
        (
            TextStyle::Monospace,
            FontId::new(base_size * 0.9, FontFamily::Monospace),
        ),
    ]
    .into();

    // Apply styles
    ctx.style_mut(|style| {
        style.text_styles = text_styles;

        // Adjust spacing for mobile
        if screen_size.is_mobile() {
            style.spacing.item_spacing = egui::vec2(6.0, 4.0);
            style.spacing.button_padding = egui::vec2(6.0, 3.0);
        }
    });
}

/// Get the current screen size category.
#[must_use]
pub fn screen_size(ctx: &egui::Context) -> ScreenSize {
    ScreenSize::from_width(ctx.screen_rect().width())
}
