//! Möbius Transformations Visualization
//!
//! Run with: `cargo run -p clifford-viz --example conformal2_mobius --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::Conformal2MobiusDemo;

fn main() -> eframe::Result<()> {
    run_app::<Conformal2MobiusDemo>()
}
