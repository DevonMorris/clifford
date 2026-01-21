//! 2D Euclidean Rotor Visualization
//!
//! Run with: `cargo run -p clifford-viz --example euclidean2 --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::Euclidean2Demo;

fn main() {
    run_three_d_app::<Euclidean2Demo>("2D Euclidean Rotor Visualization");
}
