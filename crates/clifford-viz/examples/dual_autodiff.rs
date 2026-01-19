//! Dual Number Automatic Differentiation Visualization
//!
//! Run with: `cargo run -p clifford-viz --example dual_autodiff --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::DualAutodiffDemo;

fn main() -> eframe::Result<()> {
    run_app::<DualAutodiffDemo>()
}
