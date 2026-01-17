//! Common utilities shared across all visualization demos.
//!
//! This module provides the building blocks for creating interactive
//! geometric algebra visualizations:
//!
//! - [`app`]: Base application trait and window runner
//! - [`animation`]: Animation timing and playback controls
//! - [`colors`]: Consistent color palette
//! - [`grid`]: 2D grid and axis rendering
//! - [`shapes`]: 2D shape primitives (points, arrows, circles, etc.)
//! - [`widgets`]: Reusable UI components (sliders, displays, etc.)
//!
//! # Quick Start
//!
//! ```ignore
//! use clifford_viz::common::prelude::*;
//!
//! struct MyDemo {
//!     angle: f32,
//!     animation: Animation,
//! }
//!
//! impl Default for MyDemo {
//!     fn default() -> Self {
//!         Self {
//!             angle: 0.0,
//!             animation: Animation::default(),
//!         }
//!     }
//! }
//!
//! impl VisualizationApp for MyDemo {
//!     fn name(&self) -> &'static str { "My Demo" }
//!
//!     fn update(&mut self, dt: f32) {
//!         self.animation.update(dt);
//!         if self.animation.playing {
//!             self.angle = self.animation.angle();
//!         }
//!     }
//!
//!     fn render(&self, ui: &mut egui::Ui) {
//!         // Use egui_plot to draw
//!     }
//!
//!     fn controls(&mut self, ui: &mut egui::Ui) {
//!         angle_slider(ui, "Rotation", &mut self.angle);
//!         animation_controls(ui, &mut self.animation);
//!     }
//! }
//!
//! fn main() -> eframe::Result<()> {
//!     run_app::<MyDemo>()
//! }
//! ```

pub mod animation;
pub mod app;
pub mod colors;
pub mod grid;
pub mod shapes;
pub mod widgets;

/// Prelude for convenient imports.
///
/// Imports the most commonly used types and functions.
pub mod prelude {
    pub use super::animation::{Animation, animation_controls, easing, progress_slider};
    pub use super::app::{
        EducationalContent, VisualizationApp, WindowConfig, run_app, run_app_with_options,
    };
    pub use super::colors::{hsl_to_color, lerp_color, palette, rainbow, with_alpha};
    pub use super::grid::{axes_2d, grid_2d, grid_2d_with_minor, polar_grid};
    pub use super::shapes::{
        arc_2d, arc_with_arrow, arrow_2d, bivector_2d, circle_2d, infinite_line_2d, labeled_arrow,
        labeled_line_segment, labeled_point, line_from_homogeneous, line_segment, point_marker,
        polygon, regular_polygon,
    };
    pub use super::widgets::{
        angle_slider, angle_slider_range, collapsible_section, ga_value_display, help_marker,
        info_box, mode_selector, point2_display, point3_display, section_separator, toggle_button,
        value_display, vector2_input, vector2_input_range, vector3_input, with_help,
    };
}
