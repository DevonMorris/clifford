//! 2D Projective Geometry (PGA) Visualization
//!
//! Run with: `cargo run -p clifford-viz --example projective2 --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Projective2Demo;

fn main() -> eframe::Result<()> {
    run_app::<Projective2Demo>()
}
