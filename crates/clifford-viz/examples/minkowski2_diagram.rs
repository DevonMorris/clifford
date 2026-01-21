//! Minkowski 1+1D Spacetime Diagram Visualization
//!
//! Run with: `cargo run -p clifford-viz --example minkowski2_diagram --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::Minkowski2DiagramDemo;

fn main() {
    run_three_d_app::<Minkowski2DiagramDemo>("Minkowski 1+1D Spacetime Diagram");
}
