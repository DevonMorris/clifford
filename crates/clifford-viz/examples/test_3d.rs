//! 3D Visualization Test
//!
//! Run with: `cargo run -p clifford-viz --example test_3d --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Test3DDemo;

fn main() -> eframe::Result<()> {
    run_app::<Test3DDemo>()
}
