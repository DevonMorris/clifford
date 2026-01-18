//! 2D Euclidean Rotor Visualization
//!
//! Run with: `cargo run -p clifford-viz --example euclidean2 --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Euclidean2Demo;

fn main() -> eframe::Result<()> {
    run_app::<Euclidean2Demo>()
}
