//! Minkowski 1+1D Spacetime Diagram Visualization
//!
//! Run with: `cargo run -p clifford-viz --example minkowski2_diagram --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Minkowski2DiagramDemo;

fn main() -> eframe::Result<()> {
    run_app::<Minkowski2DiagramDemo>()
}
