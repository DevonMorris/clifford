//! 3D Robot Arm visualization using native three-d rendering.
//!
//! Run with:
//! ```bash
//! cargo run -p clifford-viz --example projective3_robot --features three-d --release
//! ```

use clifford_viz::common::app_three_d::run_three_d_app_3d;
use clifford_viz::demos::Projective3RobotDemo;

fn main() {
    run_three_d_app_3d::<Projective3RobotDemo>("3D PGA Robot Arm");
}
