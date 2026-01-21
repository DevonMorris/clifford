//! 2D Robot Arm Visualization with PGA Motors
//!
//! Run with: `cargo run -p clifford-viz --example projective2_robot --features three-d --release`

use clifford_viz::common::app_three_d::run_three_d_app;
use clifford_viz::demos::RobotArmDemo;

fn main() {
    run_three_d_app::<RobotArmDemo>("2D Robot Arm Visualization with PGA Motors");
}
