//! 2D Projective Geometry (PGA) Visualization
//!
//! Run with: `cargo run -p clifford-viz --example projective2 --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::Projective2Demo;

fn main() {
    run_three_d_app::<Projective2Demo>("2D Projective Geometry");
}
