//! Circle-Circle Intersection - 2D Conformal GA Visualization
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_intersection --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Conformal2IntersectionDemo;

fn main() -> eframe::Result<()> {
    run_app::<Conformal2IntersectionDemo>()
}
