//! 3D Projective Motor visualization using native three-d rendering.
//!
//! Run with:
//! ```bash
//! cargo run -p clifford-viz --example projective3_motor --features three-d --release
//! ```

use clifford_viz::common::app_three_d::run_three_d_app_3d;
use clifford_viz::demos::Projective3MotorDemo;

fn main() {
    run_three_d_app_3d::<Projective3MotorDemo>("3D PGA Motors - Rigid Body Motion");
}
