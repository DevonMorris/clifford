//! 2D Robot Arm Visualization with PGA Motors
//!
//! Run with: `cargo run -p clifford-viz --example projective2_robot --release`

use clifford_viz::common::app::run_app;
use clifford_viz::demos::RobotArmDemo;

fn main() -> eframe::Result<()> {
    run_app::<RobotArmDemo>()
}
