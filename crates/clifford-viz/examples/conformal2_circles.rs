//! Circle from Three Points - 2D Conformal GA Visualization
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_circles --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Conformal2CirclesDemo;

fn main() -> eframe::Result<()> {
    run_app::<Conformal2CirclesDemo>()
}
