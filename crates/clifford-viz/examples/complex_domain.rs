//! Complex Domain Coloring Visualization
//!
//! Run with: `cargo run -p clifford-viz --example complex_domain --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::ComplexDomainDemo;

fn main() {
    run_three_d_app::<ComplexDomainDemo>("Complex Domain Coloring");
}
