//! Time Dilation and Twin Paradox Visualization
//!
//! Run with: `cargo run -p clifford-viz --example minkowski2_dilation --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Minkowski2DilationDemo;

fn main() -> eframe::Result<()> {
    run_app::<Minkowski2DilationDemo>()
}
