//! Mobius Transformations Visualization
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_mobius --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::Conformal2MobiusDemo;

fn main() {
    run_three_d_app::<Conformal2MobiusDemo>("Conformal 2D - Mobius Transformations");
}
