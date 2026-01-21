//! Circle from Three Points - 2D Conformal GA Visualization
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_circles --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::Conformal2CirclesDemo;

fn main() {
    run_three_d_app::<Conformal2CirclesDemo>("Conformal 2D - Circles");
}
