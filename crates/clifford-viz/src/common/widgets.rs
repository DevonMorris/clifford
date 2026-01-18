//! Reusable UI widgets for visualization demos.
//!
//! These widgets provide consistent UI patterns across all demos:
//! - Angle sliders (displaying degrees while storing radians)
//! - Vector input fields
//! - GA value displays with proper notation
//! - Info boxes and tooltips
//!
//! # Design Principles
//!
//! - Use spacing constants from [`spacing`] for consistent layout
//! - Use color palette from [`palette`] for visual consistency
//! - Group related controls visually
//! - Provide clear visual hierarchy

use egui::Color32;

use super::colors::{palette, spacing};

/// A slider for angle input that displays in degrees but stores radians.
///
/// # Arguments
/// * `ui` - The egui UI context
/// * `label` - Label displayed next to the slider
/// * `radians` - Mutable reference to the angle value (in radians)
pub fn angle_slider(ui: &mut egui::Ui, label: &str, radians: &mut f32) {
    let mut degrees = radians.to_degrees();
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(
            egui::Slider::new(&mut degrees, -180.0..=180.0)
                .suffix("\u{00b0}")
                .fixed_decimals(1),
        );
    });
    *radians = degrees.to_radians();
}

/// A slider for angle input with a custom range.
pub fn angle_slider_range(
    ui: &mut egui::Ui,
    label: &str,
    radians: &mut f32,
    min_deg: f32,
    max_deg: f32,
) {
    let mut degrees = radians.to_degrees();
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(
            egui::Slider::new(&mut degrees, min_deg..=max_deg)
                .suffix("\u{00b0}")
                .fixed_decimals(1),
        );
    });
    *radians = degrees.to_radians();
}

/// Input fields for a 2D vector.
///
/// # Arguments
/// * `ui` - The egui UI context
/// * `label` - Label for the vector
/// * `x`, `y` - Mutable references to the vector components
pub fn vector2_input(ui: &mut egui::Ui, label: &str, x: &mut f32, y: &mut f32) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(x).prefix("x: ").speed(0.1));
        ui.add(egui::DragValue::new(y).prefix("y: ").speed(0.1));
    });
}

/// Input fields for a 2D vector with custom range.
pub fn vector2_input_range(
    ui: &mut egui::Ui,
    label: &str,
    x: &mut f32,
    y: &mut f32,
    range: std::ops::RangeInclusive<f32>,
) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(
            egui::DragValue::new(x)
                .prefix("x: ")
                .speed(0.1)
                .range(range.clone()),
        );
        ui.add(
            egui::DragValue::new(y)
                .prefix("y: ")
                .speed(0.1)
                .range(range),
        );
    });
}

/// Input fields for a 3D vector.
pub fn vector3_input(ui: &mut egui::Ui, label: &str, x: &mut f32, y: &mut f32, z: &mut f32) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(x).prefix("x: ").speed(0.1));
        ui.add(egui::DragValue::new(y).prefix("y: ").speed(0.1));
        ui.add(egui::DragValue::new(z).prefix("z: ").speed(0.1));
    });
}

/// Display a GA (Geometric Algebra) value with proper basis notation.
///
/// # Arguments
/// * `ui` - The egui UI context
/// * `label` - Name of the value (e.g., "Rotor")
/// * `components` - Slice of (basis_label, coefficient) pairs
///
/// # Example
/// ```ignore
/// ga_value_display(ui, "R", &[("1", 0.707), ("e\u{2081}\u{2082}", 0.707)]);
/// // Displays: R = 0.707·1 + 0.707·e₁₂
/// ```
pub fn ga_value_display(ui: &mut egui::Ui, label: &str, components: &[(&str, f32)]) {
    ui.horizontal(|ui| {
        ui.label(format!("{} =", label));

        let text = components
            .iter()
            .enumerate()
            .filter(|(_, (_, val))| val.abs() > 1e-6) // Skip near-zero terms
            .map(|(i, (basis, val))| {
                let sign = if i == 0 || *val < 0.0 { "" } else { "+ " };
                if *basis == "1" {
                    format!("{}{:.3}", sign, val)
                } else {
                    format!("{}{:.3}{}", sign, val, basis)
                }
            })
            .collect::<Vec<_>>()
            .join(" ");

        ui.monospace(if text.is_empty() { "0" } else { &text });
    });
}

/// Display a formatted value with a label.
pub fn value_display(ui: &mut egui::Ui, label: &str, value: f32, precision: usize) {
    ui.horizontal(|ui| {
        ui.label(format!("{}:", label));
        ui.monospace(format!("{:.prec$}", value, prec = precision));
    });
}

/// Display a 2D point/vector value.
pub fn point2_display(ui: &mut egui::Ui, label: &str, x: f32, y: f32) {
    ui.horizontal(|ui| {
        ui.label(format!("{}:", label));
        ui.monospace(format!("({:.3}, {:.3})", x, y));
    });
}

/// Display a 3D point/vector value.
pub fn point3_display(ui: &mut egui::Ui, label: &str, x: f32, y: f32, z: f32) {
    ui.horizontal(|ui| {
        ui.label(format!("{}:", label));
        ui.monospace(format!("({:.3}, {:.3}, {:.3})", x, y, z));
    });
}

/// An info box with styled background for explanatory text.
///
/// Uses the design system's surface color with warm tint for better readability.
pub fn info_box(ui: &mut egui::Ui, text: &str) {
    egui::Frame::none()
        .fill(Color32::from_rgba_unmultiplied(25, 25, 30, 220)) // Warm dark with slight blue tint
        .inner_margin(spacing::SM)
        .outer_margin(spacing::XS)
        .rounding(4.0)
        .show(ui, |ui| {
            ui.label(
                egui::RichText::new(text)
                    .color(palette::TEXT_PRIMARY)
                    .size(12.0),
            );
        });
}

/// A collapsible section with a header.
pub fn collapsible_section<R>(
    ui: &mut egui::Ui,
    title: &str,
    default_open: bool,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::CollapsingResponse<R> {
    egui::CollapsingHeader::new(title)
        .default_open(default_open)
        .show(ui, add_contents)
}

/// A section separator with an optional title.
///
/// Creates visual separation between control groups with proper spacing.
pub fn section_separator(ui: &mut egui::Ui, title: Option<&str>) {
    ui.add_space(spacing::SM);
    if let Some(title) = title {
        ui.horizontal(|ui| {
            ui.separator();
            ui.label(
                egui::RichText::new(title)
                    .small()
                    .color(palette::TEXT_SECONDARY),
            );
            ui.separator();
        });
    } else {
        ui.separator();
    }
    ui.add_space(spacing::XS);
}

/// A help tooltip that appears on hover.
///
/// Uses a styled question mark that's more universally readable.
pub fn help_marker(ui: &mut egui::Ui, text: &str) {
    ui.label(
        egui::RichText::new("?")
            .small()
            .color(palette::TEXT_SECONDARY),
    )
    .on_hover_text(text);
}

/// A labeled help marker next to a control.
pub fn with_help<R>(
    ui: &mut egui::Ui,
    help_text: &str,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> R {
    ui.horizontal(|ui| {
        let r = add_contents(ui);
        help_marker(ui, help_text);
        r
    })
    .inner
}

/// Radio button group for mode selection.
pub fn mode_selector<T: PartialEq + Copy>(
    ui: &mut egui::Ui,
    current: &mut T,
    options: &[(T, &str)],
) {
    ui.horizontal(|ui| {
        for (value, label) in options {
            ui.selectable_value(current, *value, *label);
        }
    });
}

/// A toggle button (checkbox styled as button).
pub fn toggle_button(ui: &mut egui::Ui, label: &str, value: &mut bool) -> egui::Response {
    let response = ui.selectable_label(*value, label);
    if response.clicked() {
        *value = !*value;
    }
    response
}

/// A compact row of toggle checkboxes.
///
/// Use this instead of vertical checkbox lists for display options.
pub fn toggle_row(ui: &mut egui::Ui, toggles: &mut [(&str, &mut bool)]) {
    ui.horizontal(|ui| {
        for (label, value) in toggles.iter_mut() {
            ui.checkbox(value, *label);
        }
    });
}

/// A read-only value display with a subtle background.
///
/// Use for computed values that the user can't edit.
pub fn readonly_value(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.horizontal(|ui| {
        ui.label(format!("{}:", label));
        egui::Frame::none()
            .fill(Color32::from_rgba_unmultiplied(40, 40, 46, 180))
            .inner_margin(egui::vec2(spacing::XS, 2.0))
            .rounding(2.0)
            .show(ui, |ui| {
                ui.monospace(value);
            });
    });
}

/// A styled header for a control group.
pub fn group_header(ui: &mut egui::Ui, title: &str) {
    ui.add_space(spacing::SM);
    ui.label(
        egui::RichText::new(title)
            .strong()
            .color(palette::TEXT_PRIMARY),
    );
    ui.add_space(spacing::XS);
}
