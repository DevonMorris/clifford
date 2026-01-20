//! 3D Euclidean Rotors Demo
//!
//! Run with: `cargo run -p clifford-viz --example euclidean3 --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Euclidean3Demo;

fn main() -> eframe::Result<()> {
    run_app::<Euclidean3Demo>()
}
