//! Color palette for consistent styling across all visualization demos.
//!
//! Colors are organized by their semantic meaning:
//! - Geometric objects (points, lines, planes, etc.)
//! - Transformations (rotors, motors)
//! - Coordinate axes
//! - UI states

use egui::Color32;

/// Primary color palette for geometric objects and UI elements.
pub mod palette {
    use super::*;

    // Primary geometric objects
    /// Color for points in visualizations.
    pub const POINT: Color32 = Color32::from_rgb(66, 133, 244); // Blue
    /// Color for lines in visualizations.
    pub const LINE: Color32 = Color32::from_rgb(234, 67, 53); // Red
    /// Color for planes in visualizations.
    pub const PLANE: Color32 = Color32::from_rgb(52, 168, 83); // Green
    /// Color for circles in visualizations.
    pub const CIRCLE: Color32 = Color32::from_rgb(251, 188, 5); // Yellow
    /// Color for spheres in visualizations.
    pub const SPHERE: Color32 = Color32::from_rgb(154, 160, 166); // Gray

    // Transformations
    /// Color for rotors in visualizations.
    pub const ROTOR: Color32 = Color32::from_rgb(171, 71, 188); // Purple
    /// Color for motors in visualizations.
    pub const MOTOR: Color32 = Color32::from_rgb(0, 172, 193); // Cyan

    // Axes
    /// Color for the X axis.
    pub const X_AXIS: Color32 = Color32::from_rgb(244, 67, 54); // Red
    /// Color for the Y axis.
    pub const Y_AXIS: Color32 = Color32::from_rgb(76, 175, 80); // Green
    /// Color for the Z axis.
    pub const Z_AXIS: Color32 = Color32::from_rgb(33, 150, 243); // Blue
    /// Color for the time axis (Minkowski).
    pub const T_AXIS: Color32 = Color32::from_rgb(255, 193, 7); // Amber

    // UI states
    /// Color for selected objects.
    pub const SELECTED: Color32 = Color32::from_rgb(255, 235, 59); // Highlight yellow
    /// Color for hovered objects.
    pub const HOVERED: Color32 = Color32::from_rgb(255, 255, 255); // White
    /// Color for grid lines.
    pub const GRID: Color32 = Color32::from_rgb(66, 66, 66); // Dark gray
    /// Background color.
    pub const BACKGROUND: Color32 = Color32::from_rgb(30, 30, 30); // Near black

    // Secondary geometric objects (for multiple objects of same type)
    /// Secondary point color.
    pub const POINT_SECONDARY: Color32 = Color32::from_rgb(100, 181, 246); // Light blue
    /// Secondary line color.
    pub const LINE_SECONDARY: Color32 = Color32::from_rgb(239, 154, 154); // Light red
    /// Secondary plane color.
    pub const PLANE_SECONDARY: Color32 = Color32::from_rgb(129, 199, 132); // Light green
}

/// Convert HSL (Hue, Saturation, Lightness) to RGB Color32.
///
/// # Arguments
/// * `hue` - Hue value from 0.0 to 1.0 (representing 0° to 360°)
/// * `saturation` - Saturation from 0.0 to 1.0
/// * `lightness` - Lightness from 0.0 to 1.0
///
/// This is useful for domain coloring in complex number visualizations.
#[must_use]
pub fn hsl_to_color(hue: f32, saturation: f32, lightness: f32) -> Color32 {
    let hue = hue.rem_euclid(1.0);
    let saturation = saturation.clamp(0.0, 1.0);
    let lightness = lightness.clamp(0.0, 1.0);

    let c = (1.0 - (2.0 * lightness - 1.0).abs()) * saturation;
    let x = c * (1.0 - ((hue * 6.0) % 2.0 - 1.0).abs());
    let m = lightness - c / 2.0;

    let (r, g, b) = match (hue * 6.0) as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };

    Color32::from_rgb(
        ((r + m) * 255.0) as u8,
        ((g + m) * 255.0) as u8,
        ((b + m) * 255.0) as u8,
    )
}

/// Create a color from a normalized value (0.0 to 1.0) using a rainbow gradient.
///
/// Useful for visualizing scalar fields or parametric curves.
#[must_use]
pub fn rainbow(t: f32) -> Color32 {
    hsl_to_color(t, 0.8, 0.5)
}

/// Create a color with alpha transparency.
#[must_use]
pub fn with_alpha(color: Color32, alpha: u8) -> Color32 {
    Color32::from_rgba_unmultiplied(color.r(), color.g(), color.b(), alpha)
}

/// Linearly interpolate between two colors.
#[must_use]
pub fn lerp_color(a: Color32, b: Color32, t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let inv_t = 1.0 - t;
    Color32::from_rgb(
        (f32::from(a.r()) * inv_t + f32::from(b.r()) * t) as u8,
        (f32::from(a.g()) * inv_t + f32::from(b.g()) * t) as u8,
        (f32::from(a.b()) * inv_t + f32::from(b.b()) * t) as u8,
    )
}
