//! Mandelbrot and Julia Set Visualization
//!
//! Run with: `cargo run -p clifford-viz --example complex_fractal --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::ComplexFractalDemo;

fn main() {
    run_three_d_app::<ComplexFractalDemo>("Mandelbrot and Julia Set");
}
