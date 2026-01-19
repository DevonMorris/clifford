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
//! - [`camera3d`]: 3D camera with perspective projection
//! - [`shapes3d`]: 3D shape primitives (wireframes, axes, etc.)
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
pub mod camera3d;
pub mod colors;
pub mod grid;
pub mod shapes;
pub mod shapes3d;
pub mod widgets;

/// Prelude for convenient imports.
///
/// Imports the most commonly used types and functions.
pub mod prelude {
    pub use super::animation::{Animation, animation_controls, easing, progress_slider};
    #[cfg(any(feature = "native", target_arch = "wasm32"))]
    pub use super::app::AppWrapper;
    pub use super::app::{
        EducationalContent, ScreenSize, VisualizationApp, WindowConfig, configure_responsive_style,
        screen_size,
    };
    #[cfg(all(feature = "native", not(target_arch = "wasm32")))]
    pub use super::app::{run_app, run_app_with_options};
    // 3D camera and shapes
    pub use super::camera3d::{Camera3D, camera_controls, camera_response};
    pub use super::shapes3d::{
        arrow_3d, circle_3d, coordinate_axes, coordinate_axes_labeled, line_3d, plane_3d,
        point_3d, unit_cube_vertices, wireframe_box, wireframe_box_vertices, wireframe_sphere,
    };
    pub use super::colors::{
        // Theme-aware color functions
        active,
        background,
        circle,
        // Modules
        dark,
        // Color utilities
        darken,
        desaturate,
        grid,
        grid_major,
        grid_minor,
        highlight,
        hovered,
        hsl_to_color,
        is_dark_mode,
        lerp_color,
        light,
        line,
        line_secondary,
        line_weights,
        motor,
        palette,
        plane,
        plane_secondary,
        point,
        point_secondary,
        rainbow,
        rotor,
        selected,
        spacing,
        sphere,
        surface,
        t_axis,
        text_primary,
        text_secondary,
        transitions,
        with_alpha,
        x_axis,
        y_axis,
        z_axis,
    };
    pub use super::grid::{axes_2d, grid_2d, grid_2d_with_minor, origin_marker, polar_grid};
    pub use super::shapes::{
        arc_2d, arc_with_arrow, arrow_2d, bivector_2d, circle_2d, infinite_line_2d, labeled_arrow,
        labeled_line_segment, labeled_point, line_from_homogeneous, line_segment, point_marker,
        polygon, regular_polygon,
    };
    pub use super::widgets::{
        angle_slider, angle_slider_range, collapsible_section, ga_value_display, group_header,
        help_marker, info_box, mode_selector, point2_display, point3_display, readonly_value,
        section_separator, toggle_button, toggle_row, value_display, vector2_input,
        vector2_input_range, vector3_input, with_help,
    };
}
