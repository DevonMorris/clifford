//! Dual Number Automatic Differentiation Visualization
//!
//! Run with: `cargo run -p clifford-viz --example dual_autodiff --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::DualAutodiffDemo;

fn main() {
    run_three_d_app::<DualAutodiffDemo>("Dual Number Automatic Differentiation");
}
