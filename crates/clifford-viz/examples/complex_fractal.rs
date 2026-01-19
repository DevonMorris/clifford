//! Mandelbrot and Julia Set Visualization
//!
//! Run with: `cargo run -p clifford-viz --example complex_fractal --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::ComplexFractalDemo;

fn main() -> eframe::Result<()> {
    run_app::<ComplexFractalDemo>()
}
