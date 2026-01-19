//! Complex Domain Coloring Visualization
//!
//! Run with: `cargo run -p clifford-viz --example complex_domain --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::ComplexDomainDemo;

fn main() -> eframe::Result<()> {
    run_app::<ComplexDomainDemo>()
}
