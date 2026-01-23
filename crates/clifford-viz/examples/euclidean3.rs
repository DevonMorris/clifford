//! 3D Euclidean rotor visualization using native three-d rendering.
//!
//! Run with:
//! ```bash
//! cargo run -p clifford-viz --example euclidean3 --features three-d --release
//! ```

use clifford_viz::common::app_three_d::run_three_d_app_3d;
use clifford_viz::demos::Euclidean3Demo;

fn main() {
    run_three_d_app_3d::<Euclidean3Demo>("3D Euclidean Rotors");
}
