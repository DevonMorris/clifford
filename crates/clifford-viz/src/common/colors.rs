//! Color palette for consistent styling across all visualization demos.
//!
//! This palette uses the "Soft Technical" design language:
//! - Desaturated colors that don't fatigue the eyes
//! - Warm neutral backgrounds
//! - Clear visual hierarchy through saturation and opacity
//!
//! Colors are organized by their semantic meaning:
//! - Backgrounds and surfaces
//! - Grid elements (subtle, don't compete with content)
//! - Coordinate axes (desaturated, harmonious)
//! - Geometric objects (soft, pleasant)
//! - Transformations
//! - Interactive states (can be more saturated)

use egui::Color32;

/// Primary color palette for geometric objects and UI elements.
///
/// # Design Principles
///
/// 1. **No pure saturated primaries** - They clash and tire the eyes
/// 2. **Supporting elements fade back** - Grid, labels use low-contrast colors
/// 3. **Focus elements can pop** - But use saturation sparingly
/// 4. **Warm neutrals** - Background has slight warmth, not pure gray
pub mod palette {
    use super::*;

    // === Backgrounds ===
    /// Main background color (warm dark gray).
    pub const BACKGROUND: Color32 = Color32::from_rgb(30, 30, 34); // #1e1e22
    /// Elevated surface color (for panels, cards).
    pub const SURFACE: Color32 = Color32::from_rgb(37, 37, 42); // #25252a

    // === Grid (very subtle, shouldn't compete with content) ===
    /// Default grid line color.
    pub const GRID: Color32 = Color32::from_rgb(50, 50, 56); // #323238
    /// Minor grid lines (most subtle).
    pub const GRID_MINOR: Color32 = Color32::from_rgb(42, 42, 48); // #2a2a30
    /// Major grid lines (slightly more visible).
    pub const GRID_MAJOR: Color32 = Color32::from_rgb(58, 58, 64); // #3a3a40

    // === Axes (desaturated, harmonious - no Christmas tree effect) ===
    /// Color for the X axis (muted brick red).
    pub const X_AXIS: Color32 = Color32::from_rgb(184, 84, 80); // #b85450
    /// Color for the Y axis (muted sage green).
    pub const Y_AXIS: Color32 = Color32::from_rgb(90, 154, 94); // #5a9a5e
    /// Color for the Z axis (muted steel blue).
    pub const Z_AXIS: Color32 = Color32::from_rgb(80, 133, 176); // #5085b0
    /// Color for the time axis in Minkowski space (muted amber).
    pub const T_AXIS: Color32 = Color32::from_rgb(200, 160, 80); // #c8a050

    // === Geometric Objects (soft, pleasant) ===
    /// Color for points in visualizations (soft blue).
    pub const POINT: Color32 = Color32::from_rgb(106, 155, 212); // #6a9bd4
    /// Color for lines in visualizations (soft coral).
    pub const LINE: Color32 = Color32::from_rgb(196, 112, 104); // #c47068
    /// Color for planes in visualizations (soft green).
    pub const PLANE: Color32 = Color32::from_rgb(104, 168, 120); // #68a878
    /// Color for circles in visualizations (soft gold).
    pub const CIRCLE: Color32 = Color32::from_rgb(212, 180, 100); // #d4b464
    /// Color for spheres in visualizations (soft gray).
    pub const SPHERE: Color32 = Color32::from_rgb(148, 152, 160); // #9498a0

    // === Transformations ===
    /// Color for rotors in visualizations (soft purple).
    pub const ROTOR: Color32 = Color32::from_rgb(152, 120, 168); // #9878a8
    /// Color for motors in visualizations (soft cyan).
    pub const MOTOR: Color32 = Color32::from_rgb(100, 160, 170); // #64a0aa

    // === Secondary (lighter variants for multiple objects of same type) ===
    /// Secondary point color (lighter blue).
    pub const POINT_SECONDARY: Color32 = Color32::from_rgb(140, 175, 212); // #8cafd4
    /// Secondary line color (lighter coral).
    pub const LINE_SECONDARY: Color32 = Color32::from_rgb(212, 148, 140); // #d4948c
    /// Secondary plane color (lighter green).
    pub const PLANE_SECONDARY: Color32 = Color32::from_rgb(140, 188, 152); // #8cbc98

    // === Interactive States (can be more saturated for visibility) ===
    /// Color for selected objects (highlight gold).
    pub const SELECTED: Color32 = Color32::from_rgb(255, 213, 79); // #ffd54f
    /// Color for hovered objects (soft highlight).
    pub const HOVERED: Color32 = Color32::from_rgb(200, 220, 255); // #c8dcff
    /// Color for active/animating objects.
    pub const ACTIVE: Color32 = Color32::from_rgb(174, 213, 129); // #aed581

    // === Text ===
    /// Primary text color (off-white, less harsh than pure white).
    pub const TEXT_PRIMARY: Color32 = Color32::from_rgb(224, 224, 224); // #e0e0e0
    /// Secondary text color (for labels, hints).
    pub const TEXT_SECONDARY: Color32 = Color32::from_rgb(158, 158, 158); // #9e9e9e
}

/// Line weight constants for consistent visual hierarchy.
///
/// Use semantically appropriate weights:
/// - `HAIRLINE` for subtle supporting elements
/// - `THIN` for secondary elements
/// - `NORMAL` for standard shapes
/// - `THICK` for primary/focus elements
/// - `HEAVY` for strong emphasis
pub mod line_weights {
    /// Hairline weight for grid minor lines (0.5px).
    pub const HAIRLINE: f32 = 0.5;
    /// Thin weight for grid major lines, supporting elements (1.0px).
    pub const THIN: f32 = 1.0;
    /// Normal weight for standard shapes, axes (1.5px).
    pub const NORMAL: f32 = 1.5;
    /// Thick weight for primary/focus elements (2.0px).
    pub const THICK: f32 = 2.0;
    /// Heavy weight for strong emphasis, selected items (3.0px).
    pub const HEAVY: f32 = 3.0;
}

/// Spacing constants for consistent layout.
pub mod spacing {
    /// Extra small spacing (4px).
    pub const XS: f32 = 4.0;
    /// Small spacing (8px).
    pub const SM: f32 = 8.0;
    /// Medium spacing (12px).
    pub const MD: f32 = 12.0;
    /// Large spacing (16px).
    pub const LG: f32 = 16.0;
    /// Extra large spacing (24px).
    pub const XL: f32 = 24.0;
}

/// Transition duration constants for animations.
pub mod transitions {
    /// Fast transition for UI feedback (150ms).
    pub const FAST: f32 = 0.15;
    /// Normal transition for UI state changes (250ms).
    pub const NORMAL: f32 = 0.25;
    /// Slow transition for educational animations (500ms).
    pub const SLOW: f32 = 0.5;
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
/// Uses slightly reduced saturation for a softer look.
#[must_use]
pub fn rainbow(t: f32) -> Color32 {
    hsl_to_color(t, 0.65, 0.55) // Reduced saturation, slightly lighter
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

/// Desaturate a color by a given amount (0.0 = no change, 1.0 = grayscale).
#[must_use]
pub fn desaturate(color: Color32, amount: f32) -> Color32 {
    let amount = amount.clamp(0.0, 1.0);
    let gray = (f32::from(color.r()) * 0.299
        + f32::from(color.g()) * 0.587
        + f32::from(color.b()) * 0.114) as u8;
    lerp_color(color, Color32::from_rgb(gray, gray, gray), amount)
}

/// Lighten a color for hover states.
#[must_use]
pub fn highlight(color: Color32, amount: f32) -> Color32 {
    let amount = amount.clamp(0.0, 1.0);
    lerp_color(color, Color32::WHITE, amount * 0.3)
}

/// Darken a color for pressed states.
#[must_use]
pub fn darken(color: Color32, amount: f32) -> Color32 {
    let amount = amount.clamp(0.0, 1.0);
    lerp_color(color, Color32::BLACK, amount * 0.3)
}
