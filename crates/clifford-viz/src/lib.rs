//! Interactive visualizations for Clifford Geometric Algebra.
//!
//! This crate provides tools for creating interactive visualizations
//! that demonstrate geometric algebra concepts. It is built on top of
//! [egui](https://github.com/emilk/egui) for immediate-mode GUI and
//! [egui_plot](https://docs.rs/egui_plot) for 2D plotting.
//!
//! # Features
//!
//! - **Common Framework**: Base application trait, animation utilities,
//!   color palette, and reusable widgets
//! - **2D Shapes**: Points, arrows, circles, arcs, lines, and polygons
//! - **Grid System**: Cartesian and polar grids with customizable styling
//! - **Animation**: Play/pause controls, speed adjustment, easing functions
//!
//! # Getting Started
//!
//! Create a demo by implementing the [`common::app::VisualizationApp`] trait:
//!
//! ```ignore
//! use clifford_viz::common::prelude::*;
//! use egui_plot::{Plot, PlotPoints, Line};
//!
//! #[derive(Default)]
//! struct RotorDemo {
//!     angle: f32,
//!     animation: Animation,
//! }
//!
//! impl VisualizationApp for RotorDemo {
//!     fn name(&self) -> &'static str {
//!         "2D Rotor Demo"
//!     }
//!
//!     fn update(&mut self, dt: f32) {
//!         self.animation.update(dt);
//!         if self.animation.playing {
//!             self.angle = self.animation.angle();
//!         }
//!     }
//!
//!     fn render(&self, ui: &mut egui::Ui) {
//!         Plot::new("rotor_plot")
//!             .data_aspect(1.0)
//!             .show(ui, |plot_ui| {
//!                 // Draw grid
//!                 for line in grid_2d(5.0, 1.0) {
//!                     plot_ui.line(line);
//!                 }
//!
//!                 // Draw rotated vector
//!                 let x = self.angle.cos() as f64;
//!                 let y = self.angle.sin() as f64;
//!                 for arrow in arrow_2d(0.0, 0.0, x * 2.0, y * 2.0, palette::POINT) {
//!                     plot_ui.line(arrow);
//!                 }
//!             });
//!     }
//!
//!     fn controls(&mut self, ui: &mut egui::Ui) {
//!         angle_slider(ui, "Angle", &mut self.angle);
//!         animation_controls(ui, &mut self.animation);
//!     }
//! }
//!
//! fn main() -> eframe::Result<()> {
//!     run_app::<RotorDemo>()
//! }
//! ```
//!
//! # Module Organization
//!
//! - [`common`]: Shared utilities for all visualizations
//!   - [`common::app`]: Application framework
//!   - [`common::animation`]: Animation timing
//!   - [`common::colors`]: Color palette
//!   - [`common::grid`]: Grid rendering
//!   - [`common::shapes`]: 2D shape primitives
//!   - [`common::widgets`]: UI components
//! - `testing` (feature-gated): Testing utilities for verifying visualization correctness
//!   - `testing::invariants`: Geometric invariant helpers
//!   - `testing::scene`: Scene graph for assertions
//!
//! # Examples
//!
//! Run the built-in examples with:
//!
//! ```bash
//! cargo run -p clifford-viz --example euclidean2 --release
//! ```

#![warn(missing_docs)]

pub mod common;

/// Demo implementations (available with `three-d` feature or on WASM).
#[cfg(any(feature = "three-d", target_arch = "wasm32"))]
pub mod demos;

#[cfg(feature = "testing")]
pub mod testing;

/// Re-export egui for convenience.
pub use egui;

/// Re-export egui_plot for convenience.
pub use egui_plot;


/// Prelude module for convenient imports.
///
/// This re-exports the most commonly used items from the crate.
pub mod prelude {
    pub use crate::common::prelude::*;
    pub use egui;
    pub use egui_plot;
}
