//! Circle-Circle Intersection - 2D Conformal GA Visualization
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_intersection --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::Conformal2IntersectionDemo;

fn main() {
    run_three_d_app::<Conformal2IntersectionDemo>("Circle-Circle Intersection - 2D Conformal GA");
}
