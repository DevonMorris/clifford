//! Time Dilation and Twin Paradox Visualization
//!
//! Run with: `cargo run -p clifford-viz --example minkowski2_dilation --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::Minkowski2DilationDemo;

fn main() {
    run_three_d_app::<Minkowski2DilationDemo>("Time Dilation");
}
